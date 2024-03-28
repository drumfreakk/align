
pub fn known_i(s: &str, i: usize) -> char {
    s.chars().nth(i).expect("Invalid index")
}

