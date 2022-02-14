use algorithms::big_l;
use quadratic_sieve::{generate_factor_base_qs, quadratic_sieve};
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("17771068455963924973", 10).unwrap();
    let q = Integer::from_str_radix("14911841925380075107", 10).unwrap();
    let n = p * &q; // 264999363660442213444840706597252947111

    let bound = (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;

    println!("bound = {}", bound);
    let fb = generate_factor_base_qs(&n, bound);
    // println!("factor base: {:?}", fb);
    println!("factor base len: {}", fb.len());
    // 2760531481
    // 1000000000
    quadratic_sieve(n, 100000, &fb);
}
