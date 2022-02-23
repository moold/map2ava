use clap::{AppSettings, Arg, Command};
use hashbrown::HashMap;
use rust_htslib::bam::{
    ext::BamRecordExtensions, record::Cigar, record::CigarStringView, Read, Reader, Record,
};
use std::{cmp::min, collections::VecDeque, rc::Rc, string::String};

//version number
const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));
const BIN: u32 = 6; // 2^^6 == 64

//qs, qe, rev, tid, ts, te
type Paf = (i32, i32, bool, u32, i32, i32);

// a continuous aligned block without RefSkip (`N`) cigar
#[derive(Debug)]
struct Block {
    qs: i32,
    qe: i32,
    ts: i64,
    te: i64,
    cs: usize, //cigar start index
    ce: usize, //cigar end index
}

// an alignment block
static REVERSE: u8 = 0b0000_0001; // reverse flag
#[derive(Debug)]
struct Aln {
    flag: u8,
    tid: i32, //ref id
    ts: i64,  //ref. start
    te: i64,  //ref. end
    qid: u32,
    cigar: CigarStringView, // cigars
    blocks: Vec<Block>,
}

// convert qname:String => qname:usize
struct Name {
    name: Rc<String>,
    len: usize,
}

struct Q {
    names: Vec<Name>,
    idxs: HashMap<Rc<String>, u32>, //(name, index in names)
}

impl Q {
    fn new() -> Q {
        Q {
            names: Vec::new(),
            idxs: HashMap::new(),
        }
    }

    fn get_idx(&mut self, name: impl Into<String>, len: usize) -> u32 {
        let name = Rc::new(name.into());
        if self.idxs.contains_key(&name) {
            *self.idxs.get(&name).unwrap()
        } else {
            let idx = self.names.len();
            self.names.push(Name {
                name: name.clone(),
                len,
            });
            self.idxs.insert(name, idx.try_into().unwrap());
            idx as u32
        }
    }
}

fn split_to_blocks(p: i64, cigar: &CigarStringView) -> Vec<Block> {
    let mut blocks: Vec<Block> = Vec::new();
    let mut qs = 0;
    let mut qe = 0;
    let mut ts = p;
    let mut te = p;
    let mut cs = 0;
    for (i, c) in cigar.iter().enumerate() {
        match c {
            Cigar::SoftClip(l) | Cigar::HardClip(l) => {
                if i == 0 {
                    qe += l;
                    qs += l;
                    cs += 1;
                } else {
                    blocks.push(Block {
                        qs: qs as i32,
                        qe: qe as i32 - 1,
                        ts,
                        te: te - 1,
                        cs,
                        ce: i,
                    });
                    cs = i + 1;
                }
            }
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                qe += l;
                te += *l as i64;
            }
            Cigar::Ins(l) => {
                qe += l;
            }
            Cigar::Del(l) => {
                te += *l as i64;
            }
            Cigar::RefSkip(l) => {
                // split alignment to Block for full-length iso-seq data
                blocks.push(Block {
                    qs: qs as i32,
                    qe: qe as i32 - 1,
                    ts,
                    te: te - 1,
                    cs,
                    ce: i,
                });
                te += *l as i64;
                ts = te;
                qs = qe;
                cs = i + 1;
            }
            _ => (),
        }
    }
    let i = cigar.len();
    if cs < i {
        blocks.push(Block {
            qs: qs as i32,
            qe: qe as i32 - 1,
            ts,
            te: te - 1,
            cs,
            ce: i,
        });
    }
    blocks
}

fn get_read_matched_pos(tp: i64, qs: i32, ts: i64, cigar: &[Cigar], is_start: bool) -> i32 {
    let mut qe = qs - 1;
    let mut te = ts - 1;
    for c in cigar {
        match c {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                qe += *l as i32;
                let l = *l as i64;
                te += l;
                if te - l <= tp && tp <= te {
                    return qe - (te - tp) as i32; //return offset
                }
            }
            Cigar::Ins(l) => {
                let l = *l as i32;
                qe += l;
            }
            Cigar::Del(l) => {
                let l = *l as i64;
                te += l;
                if te - l <= tp && tp <= te {
                    return if is_start {
                        (tp - te - 1) as i32
                    } else {
                        (te - l - tp) as i32
                    }; //return offset
                }
            }
            _ => unreachable!(),
        };
    }
    -1
}

fn get_ava(
    bk1: &[Block],
    cg1: &CigarStringView,
    bk2: &[Block],
    cg2: &CigarStringView,
) -> Vec<(i32, i32, i32, i32)> {
    let mut ovls: Vec<(i32, i32, i32, i32)> = Vec::new();
    for c1 in bk1 {
        let mut ts = 0;
        let mut te = 0;
        let mut qs1 = 0;
        let mut qe1 = 0;
        let mut qs2 = 0;
        let mut qe2 = 0;

        for c2 in bk2 {
            if c2.te <= c1.ts {
                continue;
            } else if c2.ts <= c1.ts && c1.ts <= c2.te && c2.te <= c1.te {
                ts = c1.ts;
                te = c2.te;
                qs1 = c1.qs;
                qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                qe2 = c2.qe;
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else if c2.ts <= c1.ts && c2.te >= c1.te {
                ts = c1.ts;
                te = c1.te;
                qs1 = c1.qs;
                qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = c1.qe;
                qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else if c1.ts <= c2.ts && c2.te <= c1.te {
                ts = c2.ts;
                te = c2.te;
                qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                qs2 = c2.qs;
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                qe2 = c2.qe;
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else if c1.ts <= c2.ts && c2.ts <= c1.te && c2.te >= c1.te {
                ts = c2.ts;
                te = c1.te;
                qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                qs2 = c2.qs;
                // println!("ts: {} qs1:{:?} qs2:{}",ts, qs1, qs2 );
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = c1.qe;
                qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else {
                break;
            }
        }

        if ts != te && qs1 >= 0 && qe1 >= 0 && qs2 >= 0 && qe2 >= 0 {
            ovls.push((qs1, qe1, qs2, qe2));
        }
    }

    // for i in &ovls{
    //     println!("{:?}", i);
    // }

    let mut res: Vec<(i32, i32, i32, i32)> = Vec::new();
    if !ovls.is_empty() {
        let (mut qs1, mut qe1, mut qs2, mut qe2) = ovls[0];
        for (qs1_, qe1_, qs2_, qe2_) in ovls.into_iter().skip(1) {
            if ((qs1_ - qe1) - (qs2_ - qe2)).abs() < 30 {
                qe1 = qe1_;
                qe2 = qe2_;
            } else {
                res.push((qs1, qe1, qs2, qe2));
                qs1 = qs1_;
                qe1 = qe1_;
                qs2 = qs2_;
                qe2 = qe2_;
            }
        }
        res.push((qs1, qe1, qs2, qe2));
    }
    res
}

fn get_min_depth(s: usize, e: usize, depth: &mut [u32]) -> u32 {
    let mut d = u32::MAX;
    // for i in s..e+1 {
    for v in depth.iter_mut().take(e + 1).skip(s) {
        *v += 1;
        if *v < d {
            d = *v;
        }
    }
    d
}

fn out_ava(
    buf: &mut VecDeque<Aln>,
    avas: &mut HashMap<u32, Vec<Paf>>,
    q: &mut Q,
    flen: i32,
    ffra: f32,
    fcount: u32,
    fdepth: u32,
) {
    let f = buf.pop_front().unwrap();
    let fid = f.qid;
    let fname = &q.names[fid as usize];
    for b in buf {
        let bid = b.qid;
        let bname = &q.names[bid as usize];
        if b.tid != f.tid || b.ts > f.te {
            break;
        } else if fid == bid {
            continue;
        }

        for (qs1, qe1, qs2, qe2) in get_ava(&f.blocks, &f.cigar, &b.blocks, &b.cigar) {
            if flen > 0 && min(qe1 - qs1 + 1, qe2 - qs2 + 1) < flen {
                continue;
            }

            if ffra > 0.
                && (qe1 - qs1 < (ffra * (fname.len as f32)) as i32
                    || qe2 - qs2 < (ffra * (bname.len as f32)) as i32)
            {
                continue;
            }

            let (qs1, qe1) = if f.flag & REVERSE == 1 {
                (fname.len as i32 - qe1 - 1, fname.len as i32 - qs1 - 1)
            } else {
                (qs1, qe1)
            };

            let (qs2, qe2) = if b.flag & REVERSE == 1 {
                (bname.len as i32 - qe2 - 1, bname.len as i32 - qs2 - 1)
            } else {
                (qs2, qe2)
            };
            let std = f.flag & REVERSE == b.flag & REVERSE;
            let ava = (qs1, qe1 + 1, std, bid, qs2, qe2 + 1);
            let rev_ava = (qs2, qe2 + 1, std, fid, qs1, qe1 + 1);
            avas.entry(fid).or_insert(Vec::new()).push(ava);
            avas.entry(bid).or_insert(Vec::new()).push(rev_ava);
        }
    }

    if let Some(mut v) = avas.remove(&fid) {
        v.sort_unstable_by(|a, b| (b.1 - b.0).cmp(&(a.1 - a.0))); //reversed sort by overlapping length
        let mut c = 0;
        let mut depth: Vec<u32> = if fdepth > 0 {
            vec![0; (q.names[fid as usize].len >> BIN) + 1]
        } else {
            Vec::new()
        };
        for i in v {
            if fdepth > 0
                && get_min_depth(i.0 as usize >> BIN, i.1 as usize >> BIN, &mut depth) > fdepth
            {
                continue;
            }
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                q.names[fid as usize].name,
                q.names[fid as usize].len,
                i.0,
                i.1,
                if i.2 { '-' } else { '+' },
                q.names[i.3 as usize].name,
                q.names[i.3 as usize].len,
                i.4,
                i.5,
            );
            c += 1;
            if fcount > 0 && c >= fcount {
                break;
            }
        }
    }
}

fn main() {
    let args = Command::new("map2ava")
        .version(VERSION)
        .about("A tool to convert read-to-ref mapping to read-vs-read overlapping")
        .arg_required_else_help(true)
        .global_setting(AppSettings::DeriveDisplayOrder)
        .arg(Arg::new("input").required(true).help("input file"))
        .arg(
            Arg::new("length")
                .short('l')
                .long("length")
                .value_name("INT.FLOAT")
                .help("discard a record if overlap length < min(INT, FLOAT * read_length)")
                .takes_value(true),
        )
        .arg(
            Arg::new("count")
                .short('c')
                .long("count")
                .value_name("INT")
                .help("only output the longest INT records for each query")
                .takes_value(true),
        )
        .arg(
            Arg::new("depth")
                .short('d')
                .long("depth")
                .value_name("INT")
                .help("only output the longest INT fold records for each query")
                .takes_value(true),
        )
        .get_matches();

    let infile = args.value_of("input").expect("Missing input file!");
    let _len: f32 = args.value_of_t("length").unwrap_or(0.);
    let (flen, ffra) = (_len as i32, _len.fract());
    let fcount: u32 = args.value_of_t("count").unwrap_or(0);
    let fdepth: u32 = args.value_of_t("depth").unwrap_or(0);

    let mut bam = Reader::from_path(infile).unwrap();
    let mut r = Record::new();
    let mut buf: VecDeque<Aln> = VecDeque::new();
    let mut q = Q::new();
    let (mut pre_tid, mut pre_pos) = (0, 0);

    let mut avas: HashMap<u32, Vec<Paf>> = HashMap::with_capacity(5000);
    while let Some(ret) = bam.read(&mut r) {
        ret.expect("BAM/SAM parsing failed!");
        if r.is_unmapped() {
            continue;
        }
        assert!(
            r.tid() > pre_tid || r.reference_start() >= pre_pos,
            "Unsorted input file!"
        );
        if !r.is_secondary() {
            //only consider the primary alignment
            let cigar = r.cigar();
            let blocks = split_to_blocks(r.reference_start(), &cigar);
            let qid = q.get_idx(
                String::from_utf8(r.qname().to_vec()).unwrap(),
                r.seq_len_from_cigar(true),
            );
            let mut b = Aln {
                flag: 0,
                tid: r.tid(),
                ts: r.reference_start(),
                te: r.reference_end(),
                qid,
                cigar,
                blocks,
            };
            if r.is_reverse() {
                b.flag |= REVERSE;
            }
            buf.push_back(b);

            let f = buf.front().unwrap();
            let b = buf.back().unwrap();
            if b.tid != f.tid || b.ts > f.te {
                out_ava(&mut buf, &mut avas, &mut q, flen, ffra, fcount, fdepth);
            }
        }
        pre_tid = r.tid();
        pre_pos = r.reference_start();
    }

    while !buf.is_empty() {
        out_ava(&mut buf, &mut avas, &mut q, flen, ffra, fcount, fdepth);
    }
}
