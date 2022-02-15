use facto_rs::fermat::fermat;
use facto_rs::Factorizer;
use facto_rs::Fermat;
use rug::Integer;

fn main() {
    let n = Integer::from_str_radix("5959", 10).unwrap();
    let res = fermat(&n);
    assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
    let f = Fermat::new(n.clone());
    let res = f.factor();
    assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
}
