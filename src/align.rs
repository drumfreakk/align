
use std::string::String;
use std::cmp::{max, min, Ordering};
use std::fmt;

use crate::score;
use crate::util::known_i;

pub struct Alignment {
    score: Option<i8>,
    seq: [String; 2],
    alignments: Vec<Position>,
}

#[derive(PartialEq)]
pub struct Position {
    length: usize,
    start: [usize; 2],
    strictness: AAMatchType,
    aligned: bool,
}

#[derive(PartialEq)]
pub enum AAMatchType {
    Exact,
    Class,
}

struct AACompare {
    equal: bool,
    strictness: AAMatchType,
}

#[derive(Debug)]
struct SegmentExtremes {
    start: usize,
    end: usize,
}

impl Alignment {
    fn lcs(&mut self, classes: Option<&[String]>) {
        let n = self.seq[0].len();
        let m = self.seq[1].len();
    
        let mut dp_in = Vec::new();
        for _i in 0..m+1 { dp_in.push(0); }
        let mut dp = Vec::new();
        for _i in 0..n+1 { dp.push(dp_in.clone()); }

        let mut out: Vec<Position> = Vec::new();

        for i in 0..n {
            'outer: for j in 0..m {
                let cmp = matching_aa(known_i(&self.seq[0], i), known_i(&self.seq[1], j), classes);
                if cmp.equal {
                    dp[i+1][j+1] =  dp[i][j] + 1;
                    if dp[i+1][j+1] > 0 {
                        let length = dp[i+1][j+1];
                        let start = [i+1-length, j+1-length];
    
                        for i in 0..out.len() {
                            if out[i].start == start {
                                if out[i].length < length {
                                    out[i].length = length;
                                    if out[i].strictness == AAMatchType::Exact {
                                        out[i].strictness = cmp.strictness;
                                    }
                                }
                                continue 'outer;
                            }
                        }
                        out.push(Position{length, start, strictness: cmp.strictness, aligned: start[0] == start[1]});
                    }
                } else {
                    dp[i+1][j+1] = 0;
                }
            }
        }
        if out.len() == 0 {
            self.alignments = Vec::new();
        } else{
            out.sort_unstable_by(|a,b| a.partial_cmp(b).unwrap());
            out.reverse();
            self.alignments = out;
        }
    }

    fn min_pos<I>(&self, mut index: I) -> SegmentExtremes
        where I: Iterator<Item = usize>, {
        let mut extr_i = index.next().unwrap();
        let mut min_s = self.alignments[extr_i].start[0];
        for i in index {
            if min_s > self.alignments[i].start[0] {
                min_s = self.alignments[i].start[0];
                extr_i = i;
            }
        }
        SegmentExtremes {
            start: min_s,
            end: min_s + self.alignments[extr_i].length,
        }
    }

    fn max_pos<I>(&self, index: I) -> SegmentExtremes
        where I: Iterator<Item = usize>, {
        let mut max_s = 0;
        let mut extr_l = 0;
        for i in index {
            //TODO is index 0 actually always the right one?
            if max_s < self.alignments[i].start[0] || max_s == 0 {
                max_s = self.alignments[i].start[0];
                extr_l = self.alignments[i].length;
            }
        }
        SegmentExtremes {
            start: max_s,
            end: max_s + extr_l,
        }
    }
       
    fn match_lengths(&mut self, from_start: bool) {
        while self.seq[0].len() != self.seq[1].len() {
            if self.seq[0].len() < self.seq[1].len() {
                if from_start {
                    self.seq[0].insert(0, '-');
                } else {
                    self.seq[0].push('-');
                }
            } else{
                if from_start {
                    self.seq[1].insert(0,'-');
                } else {
                    self.seq[1].push('-');
                }
            } 
        }
        self.trim_gaps();
    }

    fn trim_gaps(&mut self) {
        let mut to_pop = Vec::new();
        for i in 0..min(self.seq[0].len(), self.seq[1].len()) {
            if known_i(&self.seq[0], i) == known_i(&self.seq[1], i) && known_i(&self.seq[0], i) == '-' {
                to_pop.push(i);
            }
        }
        to_pop.reverse();
        for p in to_pop {
            self.seq[0].remove(p);
            self.seq[1].remove(p);
        }
    }

    fn score(&mut self, scores: &score::Scores) -> Result<(), String> {
        match score::score_subs(scores, &self.seq[0], &self.seq[1]) {
            Ok(i)   => {
                self.score = Some(i);
                return Ok(());
            },
            Err(e)  => return Err(e),
        }
    }

    fn insert_gap(&mut self, on_left: bool, align_on_match: usize, pos: usize) {
        let insert_on = match on_left {
            true    => 0,
            false   => 1,
        };
        let substr = &self.alignments[align_on_match];
        if substr.start[0] > substr.start[1] {
            self.seq[insert_on].insert(pos, '-');
        } else if substr.start[0] < substr.start[1] {
            self.seq[(insert_on+1)%2].insert(pos, '-');
        }
    }
}

impl Position {
    fn max(&self) -> usize { max(self.start[0], self.start[1]) }
    fn min(&self) -> usize { min(self.start[0], self.start[1]) }

    fn to_segment_max(&self) -> SegmentExtremes {
        SegmentExtremes {
            start: self.max(),
            end: self.max() + self.length,
        }
    }

    fn to_segment_min(&self) -> SegmentExtremes {
        SegmentExtremes {
            start: self.min(),
            end: self.min() + self.length,
        }
    }
}

impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.length > other.length {
            return Some(Ordering::Greater);
        } else if self.length < other.length {
            return Some(Ordering::Less);
        } else {
            let dself = self.max() - self.min();
            let dother = other.max() - other.min();

            if dself < dother {
                return Some(Ordering::Greater);
            } else if dself > dother {
                return Some(Ordering::Less);
            } else {
                if self.strictness == other.strictness {
                    return Some(Ordering::Equal);
                } else if self.strictness == AAMatchType::Exact {
                    return Some(Ordering::Greater);
                } else {
                    return Some(Ordering::Less);
                }
            }
        }
    }
}

fn matching_aa(a: char, b: char, classes: Option<&[String]>) -> AACompare {
    if a == b {
        return AACompare{
            equal: true,
            strictness: AAMatchType::Exact,
        };
    } else {
        match classes {
            Some(c) => {
                for class in c {
                    if class.contains(a) && class.contains(b) {
                        return AACompare{
                            equal: true,
                            strictness: AAMatchType::Class,
                        };
                    }
                }
            },
            None    => ()
        }
    }
    AACompare{
        equal: false,
        strictness: AAMatchType::Exact,
    }
}

pub fn align_seqs(scores: &score::Scores, seq0: &String, seq1: &String, classes: &[String]) -> Alignment {
  
    for _i in 0..50 {
        print!("-");
    }
    print!("\n");

    let mut to_align = Alignment{score: None, seq: [seq0.to_string(), seq1.to_string()], alignments: Vec::new()};
    to_align.lcs(Some(classes));
    display_pos_list(&to_align.alignments, Some(5));
    println!("{}", to_align);

    //TODO put this in the middle part
    // Align based on top sequence
    // Only creates gaps at the start & end
  /*  let substr = &to_align.alignments[0];
    let offset = substr.max() - substr.min();
    for _i in 0..offset {
        if substr.start[0] < substr.start[1] {
            to_align.seq[0].insert(0, '-');
        } else {
            to_align.seq[1].insert(0, '-');
        }
    }
    to_align.match_lengths(false);
    to_align.lcs(Some(classes));

    // Align based on next best sequence
    // Can introduce gaps
    display_pos_list(&to_align.alignments, Some(5));
    to_align.score(scores).unwrap();
*/    

    for i in 0..1 {
        println!("WORKING ON {}", i);
        let substr = &to_align.alignments[i];

        let fix_min;
        let fix_max;
        if i == 0 {
            let len = max(to_align.seq[0].len(), to_align.seq[1].len());
            fix_min = SegmentExtremes{start: 0, end: 0};
            fix_max = SegmentExtremes{start: len, end: len};
        } else {
            fix_min = to_align.min_pos(0..i);
            fix_max = to_align.max_pos(0..i);
        }

        let align_min = substr.to_segment_min();
        let align_max = substr.to_segment_max();

        //TODO this right?
        if !(align_min.end < fix_min.start && align_max.start >= fix_max.end) && !substr.aligned {
            println!("no overlap, adjusting");
        
            // TODO <= & > correct?
            if fix_max.end <= align_min.start {
                // insert gap on right
                println!("right");
                for j in align_min.start..align_max.start {
                    to_align.insert_gap(false, i, j);
                }
                to_align.match_lengths(false);
            } else if fix_min.start > align_min.end {
                // Insert gap on the left of the reference align
                println!("left");
                for _j in align_min.end..align_max.end {
                    to_align.insert_gap(true, i, align_max.end);
                }
                to_align.match_lengths(true);
            } else {
                // overlap or middle?, continue
                println!("middle");
                let i_move = match substr.start[0] > substr.start[1] {
                    true    => 1,
                    false   => 0,
                };

                let i_stay = (i_move + 1) % 2;
                for _i in substr.start[i_move]..substr.start[i_stay] {
                    to_align.seq[i_move].insert(substr.start[i_move], '-');
                    to_align.seq[i_stay].insert(align_max.end, '-');
                }
                to_align.match_lengths(true);
            }
            // TODO check if its better
        } else {
            println!("overlap, skipping");
        }

        to_align.lcs(Some(classes));
        to_align.score(scores).unwrap();
        display_pos_list(&to_align.alignments, Some(5));
        println!("{}", to_align);
    }

    to_align
}


/* Display */

fn display_pos_list(pos: &Vec<Position>, amount: Option<usize>) {
    println!("{} matches:\nIndex\tLength\tStart0\tStart1\tAligned\tStrict", pos.len());
    let amount = match amount {
        Some(i) => i,
        None    => pos.len(),
    };
    for i in 0..amount {
        println!("{}\t{}\t{}\t{}\t{}\t{}", 
            i,
            pos[i].length,
            pos[i].start[0],
            pos[i].start[1],
            pos[i].aligned,
            pos[i].strictness
            );
    }
}

impl fmt::Display for Alignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut res = Vec::new();
        res.push(write!(f, "Alignment"));
        
        match self.score {
            Some(i) => res.push(write!(f, " (Score: {}):\n", i)),
            None    => res.push(write!(f, ":\n")),
        };

        for i in 0..max(self.seq[0].len(), self.seq[1].len()) {
            if i % 5 == 0 {
                res.push(write!(f, "{}", i));
            } else {
                if i % 5 != 1 || i < 10 {
                    res.push(write!(f, " "));
                }
            }
        }
        res.push(write!(f, "\n"));

        for i in 0..self.seq.len() {
            for j in 0..self.seq[i].len() {
                for substr in &self.alignments {
                    if substr.start[0] == substr.start[1] {
                        if substr.start[i] == j {
                            res.push(write!(f, "\u{1b}[1m"));
                        } 
                        if substr.start[i] + substr.length == j {
                            res.push(write!(f, "\u{1b}[22m"));
                        }
                    }
                }
                res.push(write!(f, "{}", known_i(&self.seq[i], j)));
            }
            res.push(write!(f, "\u{1b}[22m\n"));
        }

        for r in res {
            match r {
                Err(e)  => return Err(e),
                _       => (),
            };
        }
        return Ok(())
    }
}

impl fmt::Display for AAMatchType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AAMatchType::Exact => write!(f, "Exact"),
            AAMatchType::Class => write!(f, "Class")
        }
    }
}
