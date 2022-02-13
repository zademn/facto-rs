use algorithms::L;
use dixon::{dixon, generate_factor_base};
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("3233180591", 10).unwrap();
    let q = Integer::from_str_radix("2452395553", 10).unwrap();
    let n = p.clone() * &q;
    // let n = Integer::from_str_radix("2760531481 ", 10).unwrap();

    let bound = (L(n.clone()) as f64).sqrt().round() as u64 + 1;

    println!("bound = {}", bound);
    let fb = generate_factor_base(4 * bound);
    // println!("factor base: {:?}", fb);
    println!("factor base len: {}", fb.len());

    dixon(n, &fb);
}
