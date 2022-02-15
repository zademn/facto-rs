use facto_rs::Factorizer;
use facto_rs::QuadraticSieveBuilder;
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("3507360361", 10).unwrap();
    let q = Integer::from_str_radix("3916272539", 10).unwrap();
    let n = p.clone() * &q; // 13735779066161426579

    let qs_builder = QuadraticSieveBuilder::new(n);
    let qs = qs_builder.build();
    let res = qs.factor();
    assert_eq!(res, Some((p, q)));
}
