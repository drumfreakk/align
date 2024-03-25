
use std::string::String;
use crate::score::Scores;
use crate::util::get_known_index;
use std::cmp::{max, min};


#[derive(Debug)]
pub struct Alignment {
    score: Option<i8>,
    seq: [String; 2],
}

#[derive(Debug)]
pub struct Position<'a>{
    from: &'a mut Alignment,
    length: usize,
    start: [usize; 2],
}

impl Alignment {
    pub fn score(&self) -> Option<i8> { self.score }
    pub fn seq(&self, i: usize) -> &String { &self.seq[i] }

    //TODO: return all matching sequences
    pub fn lcs(&mut self) -> Position {
        let n = self.seq[0].len();
        let m = self.seq[1].len();
    
        let mut dp_in = Vec::new();
        for _i in 0..m+1 { dp_in.push(0); }
        let mut dp = Vec::new();
        for _i in 0..n+1 { dp.push(dp_in.clone()); }

        let mut length: usize = 0;
        let mut start: [usize; 2] = [0,0];
    
        for i in 0..n {
            for j in 0..m {
                if get_known_index(&self.seq[0], i) == get_known_index(&self.seq[1], j) {
                    dp[i+1][j+1] =  dp[i][j] + 1;
                    if dp[i+1][j+1] > length {
                        length = dp[i+1][j+1];
                        start = [i+1-length, j+1-length];
                    }
                } else {
                    dp[i+1][j+1] = 0;
                }
            }
        }
        Position{
            from: self,
            length,
            start
        }
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

pub fn align_seqs(scores: &Scores, seq1: &String, seq2: &String) -> Alignment {
   
    /* Algorithm:
     *  Find longest matching sequences
     *   Alt: longest high-ranking/positive sequences
     *  Align based on longest sequence
     *  Add gaps to align other sequences
     *  Repeat using other sequences to get alternatives?
     */
    let mut to_align = Alignment{score: None, seq: [seq1.to_string(), seq2.to_string()]};
    let substr = to_align.lcs();

    let offset = max(substr.start[0], substr.start[1]) - min(substr.start[0], substr.start[1]);
    for _i in 0..offset {
        if substr.start[0] < substr.start[1] {
            substr.from.seq[0].insert(0, '-');
        } else {
            substr.from.seq[1].insert(0, '-');
        }
    }
    substr.from.match_lengths();

    println!("{:#?}", substr); 
    
    to_align.score = Some(crate::score::score_subs(scores, &to_align.seq[0], &to_align.seq[1]).unwrap());
    to_align
}

