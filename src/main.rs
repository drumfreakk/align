
mod util;
mod score;
mod align;

fn main() {
    let classes = [String::from("CSTPAG"), String::from("NDEQ"), String::from("HRK"), String::from("MILV"), String::from("FYW")];
    let scores_array = [[1,0,-1,0,-1], [0,1,0,-2,-2], [-1,0,1,-2,-2], [0,-2,-2,2,0], [-1,-2,-2,0,2]];
    let scores = crate::score::init_scores(&classes, &scores_array, 4, 6, -4, -2);
    
/*    let seq0 = String::from("HTRANN");
    let seq1 = String::from("TRANNY");
    
    let align = crate::align::align_seqs(&scores, &seq0, &seq1, &classes);
    println!("{}", &align);
*/

    let seq0 = String::from("TRASMKASELEVHGH");
    let seq1 = String::from("HGSIQVKGHL");
    
    let align = crate::align::align_seqs(&scores, &seq0, &seq1, &classes);
    println!("{}", &align);


}
