use facto_rs::traits::FullFactorizer;
use facto_rs::trial_division::trial_division;
use facto_rs::TrialDivision;
use rug::Integer;

fn main() {
    //let n = Integer::from(100u32);
    let n = Integer::from_str_radix("100", 10).unwrap();
    let fs = trial_division(&n);
    let mut n_ = Integer::from(1u32);
    for f in fs {
        n_ *= f.n();
    }
    assert_eq!(n, n_);

    let td = TrialDivision::new(n.clone());
    let fs = td.factor_full();
    let mut n_ = Integer::from(1u32);
    for f in fs {
        n_ *= f.n();
    }
    assert_eq!(n, n_);
}
