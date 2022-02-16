use facto_rs::Factorizer;
use facto_rs::QuadraticSieve;
use rug::Integer;

fn main() {
    let n = Integer::from(15347u32);
    let factor_base = vec![2, 17, 23, 29];
    let extra_relations = 2;
    let sieve_size = 1000;
    let bound = 30;
    let qs = QuadraticSieve::new(
        n,
        bound,
        sieve_size,
        factor_base,
        extra_relations,
        true, // verbose
        None, // no multithreading
    );
    let res = qs.factor();
    assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
}
