use algorithms::L;
use quadratic_sieve::{generate_factor_base_qs, quadratic_sieve};
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("3233180591", 10).unwrap();
    let q = Integer::from_str_radix("2452395553", 10).unwrap();
    //let n = p.clone() * &q;
    let n = Integer::from_str_radix("2760531481 ", 10).unwrap();

    let bound = (L(n.clone()) as f64).sqrt().round() as u64 + 1;

    println!("bound = {}", bound);
    let fb = generate_factor_base_qs(&n, bound);
    // println!("factor base: {:?}", fb);
    println!("factor base len: {}", fb.len());
    // 2760531481
    // 1000000000
    quadratic_sieve(n, 1000000, &fb);
}
