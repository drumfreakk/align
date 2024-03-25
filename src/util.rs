
pub fn get_known_index(s: &str, i: usize) -> char {
    s.chars().nth(i).expect("Invalid index")
}

