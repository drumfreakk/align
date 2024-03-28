use std::string::String;
use std::collections::HashMap;
use crate::util::known_i;

pub type Scores = HashMap<char, HashMap<char, i8>>;

pub fn init_scores(classes: &[String], scores_array: &[[i8; 5]], score_same: i8, score_cys: i8, score_gap_start: i8, score_gap_extend: i8) -> Scores {
    let mut scores = HashMap::new();
    for i in 0..classes.len(){
        let mut scores_inner = HashMap::new();
        for j in 0..scores_array[i].len(){
            for k in 0..classes[j].len() {
                scores_inner.insert(known_i(&classes[j], k), scores_array[j][i]);
            }
            scores_inner.insert('-', score_gap_extend);
            scores_inner.insert('_', score_gap_start);
        }
        for j in 0..classes[i].len(){
            let mut temp_inner = scores_inner.clone();
            let val = known_i(&classes[i], j);
            if val == 'C'{
                temp_inner.insert(val, score_cys);
            } else {
                temp_inner.insert(val, score_same);
            }
            scores.insert(known_i(&classes[i], j), temp_inner);
        }
    }
    let mut scores_g_s = HashMap::new();
    let mut scores_g_e = HashMap::new();
    for j in 0..scores_array[0].len(){
        for k in 0..classes[j].len() {
            scores_g_s.insert(known_i(&classes[j], k), score_gap_start);
            scores_g_e.insert(known_i(&classes[j], k), score_gap_extend);
        }
    }
    scores.insert('_', scores_g_s);
    scores.insert('-', scores_g_e);
    scores
}

fn check_gap(seq: &str, index: usize, gap: &mut bool) -> char {
    let mut aa = known_i(seq, index);
    if !*gap && aa == '-' {
        aa = '_';
        *gap = true;
    } else if aa != '-' {
        *gap = false;
    }
    aa
}

pub fn score_subs(scores: &Scores, seq1: &str, seq2: &str) -> Result<i8, String> {
    let mut score = 0;
    if seq1.len() != seq2.len(){
        return Err(String::from("Sequences have unequal lengths"));
    }
    
    let mut gap1 = false;
    let mut gap2 = false;
    for i in 0..seq1.len(){
        score += scores[&check_gap(seq1, i, &mut gap1)][&check_gap(seq2, i, &mut gap2)];
    }
    Ok(score)
}
