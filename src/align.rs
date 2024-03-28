
use std::string::String;
use std::cmp::{max, min, Ordering};
use std::fmt;

use crate::score;
use crate::util::get_known_index;

pub struct Alignment {
    score: Option<i8>,
    seq: [String; 2],
}

#[derive(PartialEq)]
struct Position{
    length: usize,
    start: [usize; 2],
    strictness: AAMatchType,
}

#[derive(PartialEq)]
enum AAMatchType {
    Exact,
    Class,
}

struct AACompare {
    equal: bool,
    strictness: AAMatchType,
}

impl Alignment {
    fn lcs(&self, classes: Option<&[String]>) -> Vec<Position> {
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
                                    if out[i].strictness == AAMatchType::Exact {
                                        out[i].strictness = cmp.strictness;
                                    }
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

    fn match_lengths(&mut self) {
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
            let dself = max(self.start[0], self.start[1]) - min(self.start[0], self.start[1]);
            let dother = max(other.start[0], other.start[1]) - min(other.start[0], other.start[1]);

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
   
    /* Algorithm:
     *  Find & rank matching sequences
     *      Based on length, closeness in original sequences & exactness of match
     *  Align based on best match
     *  Add gaps to align other sequences
     *      Check to make sure adding a gap actually improves the score
     *  Repeat using other sequences to get alternatives?
     */

    let mut to_align = Alignment{score: None, seq: [seq0.to_string(), seq1.to_string()]};
    let substrs = to_align.lcs(Some(classes));
    println!("{}", to_align);
    display_pos_list(&substrs);

    // Align based on top sequence
    // Only creates gaps at the start & end
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

    // Align based on next best sequence
    // Can introduce gaps
    // TODO:

    to_align.score = Some(score::score_subs(scores, &to_align.seq[0], &to_align.seq[1]).unwrap());
    to_align
}


/* Display */

fn display_pos_list(pos: &Vec<Position>) {
    println!("{} matches:\nIndex\tLength\tStart0\tStart1\tStrict", pos.len());
    for i in 0..pos.len(){
        println!("{}\t{}\t{}\t{}\t{}", 
            i, 
            pos[i].length, 
            pos[i].start[0], 
            pos[i].start[1], 
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
        
        //TODO: matching class
        let mut matches = Vec::new();
        for substr in self.lcs(None) {
            if substr.start[0] == substr.start[1] {
                matches.push(substr);
            }
        }

        for i in 0..self.seq.len(){
            for j in 0..self.seq[i].len(){
                for k in 0..matches.len(){
                    if matches[k].start[i] == j {
                        res.push(write!(f, "\u{1b}[1m"));
                    } 
                    if matches[k].start[i] + matches[k].length == j {
                        res.push(write!(f, "\u{1b}[22m"));
                    }
                }
                res.push(write!(f, "{}", get_known_index(&self.seq[i], j)));
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

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Match: (Length: {}, Strictness: {}) Starts: [{}, {}]", 
            self.length, 
            self.strictness,
            self.start[0],
            self.start[1])
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
