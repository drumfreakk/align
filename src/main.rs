
mod util;
mod score;
mod align;

fn main() {
    let classes = [String::from("CSTPAG"), String::from("NDEQ"), String::from("HRK"), String::from("MILV"), String::from("FYW")];
    let scores_array = [[1,0,-1,0,-1], [0,1,0,-2,-2], [-1,0,1,-2,-2], [0,-2,-2,2,0], [-1,-2,-2,0,2]];
    let scores = crate::score::init_scores(&classes, &scores_array, 4, 6, -4, -2);
    
    let seq1 = String::from("HTRANN");
    let seq2 = String::from("TRANNY");
    
    let align = crate::align::align_seqs(&scores, &seq1, &seq2, &classes);
    println!("{:#?}", &align);


    let seq1 = String::from("HGSIQVKGHL");
    let seq2 = String::from("TRASMKASELEVHGH");
    
    let align = crate::align::align_seqs(&scores, &seq1, &seq2, &classes);
    println!("{:#?}", &align);


}
