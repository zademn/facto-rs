use facto_rs::Factorizer;
use facto_rs::QuadraticSieveBuilder;
use rug::Integer;

fn main() {
    let p = Integer::from_str_radix("9270271495854497069", 10).unwrap();
    let q = Integer::from_str_radix("15658485116535926689", 10).unwrap();
    let n = p.clone() * &q; // 13735779066161426579

    // let p = Integer::from_str_radix("629449496776611089395111", 10).unwrap();
    // let q = Integer::from_str_radix("1157146192290239025720557", 10).unwrap();
    // let n =
    //     Integer::from_str_radix("728365088434062605443262763506401413148047996827", 10).unwrap();
    let qs_builder = QuadraticSieveBuilder::new(n);
    let qs = qs_builder.build();
    let res = qs.factor();
    assert_eq!(res, Some((p, q)));
}
