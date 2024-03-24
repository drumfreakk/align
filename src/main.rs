
mod scoring;
use crate::scoring as score;
use crate::scoring::Scores;

#[derive(Debug)]
struct Alignment {
    score: i8,
    seq1: String,
    seq2: String,
}

fn align(scores: &Scores, seq1: &String, seq2: &String) -> Alignment {
    let seq1 = String::from("TRASMKASELEVHGH-");
    let seq2 = String::from("-----HGS-IQVKGHL");
    
    Alignment{score: score::score_subs(scores, &seq1, &seq2).expect("hi"), seq1: seq1, seq2: seq2}
}

fn main() {
    let seq1 = String::from("HGSIQVKGHL");
    let seq2 = String::from("TRASMKASELEVHGH");

    let classes = [String::from("CSTPAG"), String::from("NDEQ"), String::from("HRK"), String::from("MILV"), String::from("FYW")];
    let scores_array = [[1,0,-1,0,-1], [0,1,0,-2,-2], [-1,0,1,-2,-2], [0,-2,-2,2,0], [-1,-2,-2,0,2]];
    let scores = score::init_scores(&classes, &scores_array, 4, 6, -4, -2);
    
    let align = align(&scores, &seq1, &seq2);
    println!("{:#?}", &align);
}
