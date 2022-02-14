use algorithms::big_l;
use dixon::{dixon, generate_factor_base};
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("4163498621", 10).unwrap();
    let q = Integer::from_str_radix("3359099369", 10).unwrap();
    let n = p.clone() * &q;

    let bound = (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;
    //let bound = 7;

    println!("bound = {}", bound);
    let fb = generate_factor_base(bound);
    // println!("factor base: {:?}", fb);
    println!("factor base len: {}", fb.len());

    dixon(n, &fb);
}
