use std::{
    env::args, cmp::min, collections::VecDeque,
};
use rust_htslib::bam::{
    ext::BamRecordExtensions, record::Cigar, record::CigarStringView, Read, Reader, Record,
};

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
    ts: i64, //ref. start
    te: i64, //ref. end
    qlen: usize, //read length
    qname: String, //read name
    cigar: CigarStringView, // cigars
    blocks: Vec<Block>,
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
            Cigar::RefSkip(l) => {// split alignment to Block for full-length iso-seq data
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
                    return if is_start { (tp - te - 1) as i32 } else { (te - l - tp) as i32 }; //return offset
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

fn out_ava(buf: &mut VecDeque<Aln>) {
    let f = buf.pop_front().unwrap();
    for b in buf {
        if b.tid != f.tid || b.ts > f.te {
            break;
        } else if f.qname == b.qname {
            continue;
        }

        for (qs1, qe1, qs2, qe2) in get_ava(&f.blocks, &f.cigar, &b.blocks, &b.cigar) {
            let (qs1, qe1) = if f.flag & REVERSE == 1 {
                (f.qlen as i32 - qe1 - 1, f.qlen as i32 - qs1 - 1)
            } else {
                (qs1, qe1)
            };

            let (qs2, qe2) = if b.flag & REVERSE == 1 {
                (b.qlen as i32 - qe2 - 1, b.qlen as i32 - qs2 - 1)
            } else {
                (qs2, qe2)
            };

            let std = if f.flag & REVERSE == b.flag & REVERSE {
                '+'
            } else {
                '-'
            };

            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                f.qname,
                f.qlen,
                qs1,
                qe1 + 1,
                std,
                b.qname,
                b.qlen,
                qs2,
                qe2 + 1
            );
        }
    }
}

fn main() {
    let infile = args().nth(1).expect("missing input file");
    let mut bam = Reader::from_path(infile).unwrap();
    let mut r = Record::new();
    let mut buf: VecDeque<Aln> = VecDeque::new();
    let (mut pre_tid, mut pre_pos) = (0, 0);
    while let Some(ret) = bam.read(&mut r) {
        ret.expect("BAM/SAM parsing failed!");
        if r.is_unmapped() {
            continue;
        }
        assert!(r.tid() > pre_tid || r.reference_start() >= pre_pos, "Unsorted input file!");
        if !r.is_secondary() {//only consider the primary alignment
            let cigar = r.cigar();
            let blocks = split_to_blocks(r.reference_start(), &cigar);
            let mut b = Aln {
                flag: 0,
                tid: r.tid(),
                ts: r.reference_start(),
                te: r.reference_end(),
                qlen: r.seq_len_from_cigar(true),
                qname: String::from_utf8(r.qname().to_vec()).unwrap(),
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
                out_ava(&mut buf);
            }
        }
        pre_tid = r.tid();
        pre_pos = r.reference_start();
    }

    while !buf.is_empty() {
        out_ava(&mut buf);
    }
}
