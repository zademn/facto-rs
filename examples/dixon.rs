use facto_rs::DixonBuilder;
use facto_rs::Factorizer;
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("42509", 10).unwrap();
    let q = Integer::from_str_radix("63299", 10).unwrap();
    let n = p.clone() * &q;

    let dixon_builder = DixonBuilder::new(n);
    let dixon_struct = dixon_builder.build();
    let res = dixon_struct.factor();
    assert_eq!(res, Some((p, q)));
}
