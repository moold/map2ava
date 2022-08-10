use clap::{AppSettings, Arg, Command};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::{atomic::AtomicCell, thread};
use hashbrown::HashMap;
use rust_htslib::bam::{
    ext::BamRecordExtensions, record::Cigar, record::CigarStringView, Read, Reader, Record,
};
use std::{
    cmp::min,
    io::{self, Write},
    string::String,
    sync::Arc,
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

//version number
const VERSION: &str = include_str!(concat!(env!("OUT_DIR"), "/VERSION"));
const BIN: u32 = 6; // 2^^6 == 64
                    //number of records for each batch reading
const BATCH: usize = 2_000_000;
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
const REVERSE: u8 = 0b0000_0001; // reverse flag
const SPLIT: u8 = 0b0000_0010; // split mapping flag
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
    name: Arc<String>,
    len: usize,
}

struct Q {
    names: Vec<Name>,
    idxs: HashMap<Arc<String>, u32>, //(name, index in names)
}

impl Q {
    fn new() -> Q {
        Q {
            names: Vec::new(),
            idxs: HashMap::new(),
        }
    }

    fn get_idx(&mut self, name: impl Into<String>, len: usize) -> u32 {
        let name = Arc::new(name.into());
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

// split an alignment block to continuous aligned blocks without RefSkip (`N`) cigar
// for iso-seq data
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

//get read matched pos using ref pos, skip gap
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
        for c2 in bk2 {
            let (mut ts, mut te, mut qs1, mut qe1, mut qs2, mut qe2);
            if c2.te <= c1.ts {
                continue;
            } else if c2.ts <= c1.ts && c1.ts <= c2.te && c2.te <= c1.te {
                ts = c1.ts;
                te = c2.te;
                qs1 = match cg1.0[c1.cs] {
                    Cigar::Match(_) => c1.qs,
                    _ => get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true),
                };
                qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                qe2 = match cg2.0[c2.ce - 1] {
                    Cigar::Match(_) => c2.qe,
                    _ => get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false),
                };
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else if c2.ts <= c1.ts && c2.te >= c1.te {
                ts = c1.ts;
                te = c1.te;
                qs1 = match cg1.0[c1.cs] {
                    Cigar::Match(_) => c1.qs,
                    _ => get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true),
                };
                qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = match cg1.0[c1.ce - 1] {
                    Cigar::Match(_) => c1.qe,
                    _ => get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false),
                };
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
                qs2 = match cg2.0[c2.cs] {
                    Cigar::Match(_) => c2.qs,
                    _ => get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true),
                };
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                qe2 = match cg2.0[c2.ce - 1] {
                    Cigar::Match(_) => c2.qe,
                    _ => get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false),
                };
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else if c1.ts <= c2.ts && c2.ts <= c1.te && c2.te >= c1.te {
                ts = c2.ts;
                te = c1.te;
                qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                qs2 = match cg2.0[c2.cs] {
                    Cigar::Match(_) => c2.qs,
                    _ => get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true),
                };
                // println!("ts: {} qs1:{:?} qs2:{}",ts, qs1, qs2 );
                while (qs1 < 0 && qs2 >= 0) || (qs1 >= 0 && qs2 < 0) {
                    ts -= min(qs1, qs2) as i64;
                    qs1 = get_read_matched_pos(ts, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], true);
                    qs2 = get_read_matched_pos(ts, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], true);
                }

                qe1 = match cg1.0[c1.ce - 1] {
                    Cigar::Match(_) => c1.qe,
                    _ => get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false),
                };
                qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                while (qe1 < 0 && qe2 >= 0) || (qe1 >= 0 && qe2 < 0) {
                    te += min(qe1, qe2) as i64;
                    qe1 = get_read_matched_pos(te, c1.qs, c1.ts, &cg1.0[c1.cs..c1.ce], false);
                    qe2 = get_read_matched_pos(te, c2.qs, c2.ts, &cg2.0[c2.cs..c2.ce], false);
                }
            } else {
                break;
            }

            if ts != te && qs1 >= 0 && qe1 >= qs1 && qs2 >= 0 && qe2 >= qs2 {
                ovls.push((qs1, qe1, qs2, qe2));
            }
        }
    }

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

fn accumulate_depth(s: usize, e: usize, depth: &mut [u32]) {
    for v in depth.iter_mut().take(e + 1).skip(s) {
        *v += 1;
    }
}

fn get_min_depth(s: usize, e: usize, depth: &[u32]) -> u32 {
    *depth[s..e + 1].iter().min().unwrap_or(&0)
}

fn out_sorted_paf(
    mut pafs: Vec<Paf>,
    q: &Q,
    fid: u32,
    fcount: u32,
    fdepth: u32,
    depth: &mut Vec<u32>,
    s: Option<&Sender<(u32, Vec<Paf>)>>,
) {
    pafs.sort_unstable_by_key(|k| k.0 - k.1); //reversed sort by overlapping length
    let mut c = 0;
    if fdepth > 0 {
        depth.resize((q.names[fid as usize].len >> BIN) + 1, 0);
        depth.fill(0);
    }

    let mut ret = (fid, Vec::new()); //TODO change to map
    for i in pafs {
        if fdepth > 0 {
            let (s, e) = (i.0 as usize >> BIN, i.1 as usize >> BIN);
            if get_min_depth(s, e, depth) < fdepth {
                accumulate_depth(s, e, depth);
            } else {
                continue;
            }
        }
        if s.is_some() {
            ret.1.push(i);
        } else {
            let f = &q.names[fid as usize];
            let r = &q.names[i.3 as usize];
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                f.name,
                f.len,
                i.0,
                i.1,
                if i.2 { '+' } else { '-' },
                r.name,
                r.len,
                i.4,
                i.5,
            );
        }
        c += 1;
        if fcount > 0 && c >= fcount {
            break;
        }
    }
    if !ret.1.is_empty() {
        s.unwrap().send(ret).unwrap();
    }
}

fn filt_pafs_by_count(pafs: &mut Vec<Paf>, fcount: u32) {
    pafs.sort_unstable_by_key(|k| k.0 - k.1); //reversed sort by overlapping length
    (*pafs).truncate(fcount as usize);
}

fn filt_pafs_by_depth(pafs: &mut Vec<Paf>, len: usize, fdepth: u32, depth: &mut Vec<u32>) {
    pafs.sort_unstable_by_key(|k| k.0 - k.1); //reversed sort by overlapping length
    depth.resize((len >> BIN) + 1, 0);
    depth.fill(0);
    pafs.retain(|&i| {
        let (s, e) = (i.0 as usize >> BIN, i.1 as usize >> BIN);
        if get_min_depth(s, e, depth) < fdepth {
            accumulate_depth(s, e, depth);
            true
        } else {
            false
        }
    });
}

fn caculate_ava(
    f: &Aln,
    b: &Aln,
    fname: &Name,
    bname: &Name,
    pafs: &mut HashMap<u32, Vec<Paf>>,
    depth_buf: &mut Vec<u32>, // buffer to caculate overlapping depth for each query
    flen: i32,
    ffra: f32,
    fcount: u32,
    fdepth: u32,
    s: &Sender<(u32, Vec<Paf>)>,
) {
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

        let (qs1, qe1) = if f.flag & REVERSE > 0 {
            (fname.len as i32 - qe1 - 1, fname.len as i32 - qs1 - 1)
        } else {
            (qs1, qe1)
        };

        let (qs2, qe2) = if b.flag & REVERSE > 0 {
            (bname.len as i32 - qe2 - 1, bname.len as i32 - qs2 - 1)
        } else {
            (qs2, qe2)
        };
        let std = f.flag & REVERSE == b.flag & REVERSE;
        if fcount == 0 && fdepth == 0 {
            s.send((f.qid, vec![(qs1, qe1 + 1, std, b.qid, qs2, qe2 + 1)]))
                .unwrap();
        } else {
            let paf = (qs1, qe1 + 1, std, b.qid, qs2, qe2 + 1);
            let pafs_v = pafs.entry(f.qid).or_insert(Vec::new());
            pafs_v.push(paf);
            if fcount > 0 && pafs_v.len() as u32 > fcount * 5 {
                // filter to reduce RAM
                filt_pafs_by_count(pafs_v, fcount);
            }
            if fdepth > 0 && pafs_v.len() as u32 > fdepth * 10 {
                // filter to reduce RAM
                filt_pafs_by_depth(pafs_v, fname.len, fdepth, depth_buf);
            }
        }
    }
}

// find the first Aln index overlapping with f
fn find_starti(buf: &[Aln], f: &Aln, fsi: usize) -> usize {
    for (i, b) in buf.iter().enumerate().skip(fsi) {
        if b.tid == f.tid && b.te >= f.ts {
            return i;
        }
    }
    fsi
}

fn out_ava_thread(
    buf: &[Aln],
    g_pafs: &mut HashMap<u32, Vec<Paf>>,
    q: &Q,
    flen: i32,
    ffra: f32,
    fcount: u32,
    fdepth: u32,
    thread: usize,
    buf_i: usize,
    is_last: bool,
) -> usize {
    let no_overlap = fcount == 0 && fdepth == 0;
    let valid_len = if is_last {
        buf.len()
    } else {
        find_starti(buf, buf.last().unwrap(), 0)
    };

    thread::scope(|work| {
        let buf_i = Arc::new(AtomicCell::new(buf_i));
        let (ou_s, ou_r): (Sender<(u32, Vec<Paf>)>, Receiver<(u32, Vec<Paf>)>) = bounded(1024);
        work.spawn(move |_| {
            // output thread
            let stdout = io::stdout();
            let lock = stdout.lock();
            let mut w = io::BufWriter::with_capacity(1024000, lock);
            while let Ok((fid, pafs)) = ou_r.recv() {
                let f = &q.names[fid as usize];
                for paf in pafs {
                    let r = &q.names[paf.3 as usize];
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        f.name,
                        f.len,
                        paf.0,
                        paf.1,
                        if paf.2 { '+' } else { '-' },
                        r.name,
                        r.len,
                        paf.4,
                        paf.5,
                    )
                    .unwrap();
                }
            }
        });

        work.spawn(move |_| {
            //work thread
            thread::scope(|scoped| {
                let mut handles = Vec::with_capacity(thread);
                for _i in 0..thread {
                    let buf_i = buf_i.clone();
                    let ou_s = ou_s.clone();
                    let handle = scoped.spawn(move |_| {
                        let mut fsi = 0;
                        let mut depth_buf: Vec<u32> = Vec::with_capacity(500);
                        // reads with split mappings
                        let mut pafs: HashMap<u32, Vec<Paf>> = HashMap::with_capacity(500);
                        // qids with pafs in g_pafs from previous batch
                        let mut idx = buf_i.fetch_add(1);
                        while idx < valid_len {
                            let f = &buf[idx];
                            if no_overlap {
                                fsi = idx + 1;
                            }
                            let mut fsi_updated = false;
                            for (j, b) in buf.iter().enumerate().skip(fsi) {
                                //find the leftmost boundary for the next f
                                if (!fsi_updated) && (b.tid == f.tid && b.te >= f.ts) {
                                    fsi = j;
                                    fsi_updated = true;
                                }

                                if b.tid > f.tid || (b.tid == f.tid && b.ts > f.te) {
                                    break;
                                } else if f.qid == b.qid
                                    || b.tid < f.tid
                                    || (b.tid == f.tid && b.te < f.ts)
                                {
                                    continue;
                                }

                                caculate_ava(
                                    f,
                                    b,
                                    &q.names[f.qid as usize],
                                    &q.names[b.qid as usize],
                                    &mut pafs,
                                    &mut depth_buf,
                                    flen,
                                    ffra,
                                    fcount,
                                    fdepth,
                                    &ou_s,
                                );
                            }
                            if f.flag & SPLIT == 0 {
                                if let Some(v) = pafs.remove(&f.qid) {
                                    out_sorted_paf(
                                        v,
                                        q,
                                        f.qid,
                                        fcount,
                                        fdepth,
                                        &mut depth_buf,
                                        Some(&ou_s),
                                    );
                                }
                            }
                            idx = buf_i.fetch_add(1);
                        }
                        pafs
                    });
                    handles.push(handle);
                }

                // wait all threads finish and get results
                let res: Vec<HashMap<u32, Vec<Paf>>> =
                    handles.into_iter().map(|h| h.join().unwrap()).collect();
                for pafs in res {
                    // save records without output
                    for (k, v) in pafs {
                        g_pafs.entry(k).or_insert(Vec::new()).extend(v);
                    }
                }
            })
            .unwrap();
        });
    })
    .unwrap();

    valid_len
}

fn main() {
    let args = Command::new("map2ava")
        .version(VERSION)
        .about("A tool to convert read-to-ref mapping to read-vs-read overlapping")
        .arg_required_else_help(true)
        .global_setting(AppSettings::DeriveDisplayOrder)
        .arg(
            Arg::new("input")
                .required(true)
                .help("input sorted bam/sam file"),
        )
        .arg(
            Arg::new("thread")
                .short('t')
                .long("thread")
                .value_name("INT")
                .default_value("3")
                .help("number of threads.")
                .takes_value(true),
        )
        .arg(
            Arg::new("mapq")
                .short('q')
                .long("mapq")
                .value_name("INT")
                .default_value("1")
                .help(
                    "minimum mapping quality, alignments with a lower mapping quality are ignored.",
                )
                .takes_value(true),
        )
        .arg(
            Arg::new("len")
                .short('l')
                .long("len")
                .value_name("INT")
                .default_value("0")
                .help(
                    "minimum read length, reads with length < INT are ignored.",
                )
                .takes_value(true),
        )
        .arg(
            Arg::new("maplen")
                .short('L')
                .long("maplen")
                .value_name("INT.FLOAT")
                // .default_value("10.")
                .help("discard a record if overlap length < min(INT, FLOAT * read_length)")
                .takes_value(true),
        )
        .arg(
            Arg::new("count")
                .short('c')
                .long("count")
                .value_name("INT")
                .help("only output the longest INT records for each query")
                .hide(true)
                .takes_value(true),
        )
        .arg(
            Arg::new("depth")
                .short('d')
                .long("depth")
                .hide(true)
                .value_name("INT")
                .help("only output the longest INT fold records for each query")
                .takes_value(true),
        )
        .get_matches();

    let infile = args.value_of("input").expect("Missing input file!");
    let thread: usize = args.value_of_t("thread").unwrap();
    let fqual: u8 = args.value_of_t("mapq").unwrap();
    let frlen: usize = args.value_of_t("len").unwrap();
    let _len: f32 = args.value_of_t("maplen").unwrap_or(0.);
    let (flen, ffra) = (_len as i32, _len.fract());
    let fcount: u32 = args.value_of_t("count").unwrap_or(0);
    let fdepth: u32 = args.value_of_t("depth").unwrap_or(0);

    let mut bam = Reader::from_path(infile).unwrap();
    let mut r = Record::new();
    let mut buf: Vec<Aln> = Vec::with_capacity(BATCH);
    let mut buf_i = 0;
    let mut q = Q::new();
    let (mut pre_tid, mut pre_pos) = (0, 0);

    let mut pafs: HashMap<u32, Vec<Paf>> = HashMap::with_capacity(5000);
    while let Some(ret) = bam.read(&mut r) {
        ret.expect("BAM/SAM parsing failed!");
        let rlen = r.seq_len_from_cigar(true);
        if r.is_unmapped() || r.mapq() < fqual || rlen < frlen {
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
                rlen,
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
            if r.aux(b"SA").is_ok() {
                b.flag |= SPLIT;
            }
            buf.push(b);
            if buf.len() > BATCH {
                buf_i = out_ava_thread(
                    &buf, &mut pafs, &q, flen, ffra, fcount, fdepth, thread, buf_i, false,
                );
                let skip_len = find_starti(&buf, &buf[buf_i], 0);
                buf.drain(..skip_len);
                buf_i -= skip_len;
            }
            pre_tid = r.tid();
            pre_pos = r.reference_start();
        }
    }
    if !buf.is_empty() {
        out_ava_thread(
            &buf, &mut pafs, &q, flen, ffra, fcount, fdepth, thread, buf_i, true,
        );
    }
    // for records with split mapping
    let mut depth_buf: Vec<u32> = Vec::with_capacity(500);
    for (fid, pafs) in pafs.into_iter() {
        out_sorted_paf(pafs, &q, fid, fcount, fdepth, &mut depth_buf, None);
    }
}
