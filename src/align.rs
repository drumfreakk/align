
use std::string::String;
use std::cmp::{max, min, Ordering};

use crate::score;
use crate::util::get_known_index;


#[derive(Debug)]
pub struct Alignment {
    score: Option<i8>,
    seq: [String; 2],
}

#[derive(Debug, PartialEq)]
struct Position{
    length: usize,
    start: [usize; 2],
    strictness: AAMatchType,
}

#[derive(Debug, PartialEq)]
enum AAMatchType {
    Exact,
    Class,
}

struct AACompare {
    equal: bool,
    strictness: AAMatchType,
}

impl Alignment {
    pub fn lcs(&self, classes: &[String]) -> Vec<Position> {
        let n = self.seq[0].len();
        let m = self.seq[1].len();
    
        let mut dp_in = Vec::new();
        for _i in 0..m+1 { dp_in.push(0); }
        let mut dp = Vec::new();
        for _i in 0..n+1 { dp.push(dp_in.clone()); }

        let mut out: Vec<Position> = Vec::new();

        for i in 0..n {
            'outer: for j in 0..m {
                let cmp = matching_aa(get_known_index(&self.seq[0], i), get_known_index(&self.seq[1], j), classes);
                if cmp.equal {
                    dp[i+1][j+1] =  dp[i][j] + 1;
                    if dp[i+1][j+1] > 0 {
                        let length = dp[i+1][j+1];
                        let start = [i+1-length, j+1-length];
    
                        for i in 0..out.len() {
                            if out[i].start == start {
                                if out[i].length < length {
                                    out[i].length = length;
                                }
                                continue 'outer;
                            }
                        }
                        out.push(Position{length, start, strictness: cmp.strictness});
                    }
                } else {
                    dp[i+1][j+1] = 0;
                }
            }
        }
        out.sort_unstable_by(|a,b| a.partial_cmp(b).unwrap());
        out.reverse();
        out
    }

    pub fn match_lengths(&mut self) {
       while self.seq[0].len() != self.seq[1].len() {
           if self.seq[0].len() < self.seq[1].len() {
               self.seq[0].push('-');
           } else{
               self.seq[1].push('-');
           }
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
            let dself: usize;
            let dother: usize;
            if self.start[0] < self.start[1]{
                dself = self.start[1] - self.start[0];
            } else {
                dself = self.start[0] - self.start[1];
            }
            if other.start[0] < other.start[1]{
                dother = other.start[1] - other.start[0];
            } else {
                dother = other.start[0] - other.start[1];
            }

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

fn matching_aa(a: char, b: char, classes: &[String]) -> AACompare {
    if a == b {
        return AACompare{
            equal: true,
            strictness: AAMatchType::Exact,
        };
    } else {
        for class in classes {
            if class.contains(a) && class.contains(b) {
                return AACompare{
                    equal: true,
                    strictness: AAMatchType::Class,
                };
            }
        }
    }
    AACompare{
        equal: false,
        strictness: AAMatchType::Class,
    }
}

pub fn align_seqs(scores: &score::Scores, seq1: &String, seq2: &String, classes: &[String]) -> Alignment {
   
    /* Algorithm:
     *  Find longest matching sequences
     *   Alt: longest high-ranking/positive sequences
     *  Align based on longest sequence
     *  Add gaps to align other sequences
     *  Repeat using other sequences to get alternatives?
     */
    let mut to_align = Alignment{score: None, seq: [seq1.to_string(), seq2.to_string()]};
    let substrs = to_align.lcs(classes);
    println!("{:#?}", substrs); 

    let substr = &substrs[0];


    let offset = max(substr.start[0], substr.start[1]) - min(substr.start[0], substr.start[1]);
    for _i in 0..offset {
        if substr.start[0] < substr.start[1] {
            to_align.seq[0].insert(0, '-');
        } else {
            to_align.seq[1].insert(0, '-');
        }
    }
    to_align.match_lengths();

    to_align.score = Some(score::score_subs(scores, &to_align.seq[0], &to_align.seq[1]).unwrap());
    to_align
}

