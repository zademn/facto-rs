use rug::Integer;
use crate::traits::{Factor, Factorizer, FullFactorizer};

pub struct TrialDivision {
    n: Integer,
}
impl TrialDivision {
    fn new(n: Integer) -> Self {
        Self { n }
    }
}
impl Factorizer for TrialDivision {
    fn factor(&self) -> Option<(Integer, Integer)> {
        let n = self.n.clone();
        if n == 1u32 || n == 2u32 || n == 3u32 {
            return Some((n, Integer::from(1u32)));
        }
        let bound = n.clone().sqrt();
        let mut i = Integer::from(5u32);
        while i <= bound {
            if n.clone() % &i == 0u32 {
                return Some((i.clone(), n / i));
            }
            i += 2;
        }
        None
    }
}
impl FullFactorizer for TrialDivision {
    fn factor_full(&self) -> Vec<Factor> {
        trial_division(&self.n)
    }
}

/// Trial division function.
fn trial_division(n: &Integer) -> Vec<Factor> {
    let n = n.clone();
    let mut res: Vec<Factor> = vec![];
    // Remove factors of 2 if they exist
    let (mut n, c) = n.remove_factor(&2u32.into());
    if c != 0 {
        res.push(Factor::new(2u32.into(), c));
    }
    let bound = n.clone().sqrt();
    let mut i = Integer::from(3u32);
    while i <= bound {
        let (n_, c) = n.remove_factor(&i);
        if c != 0 {
            res.push(Factor::new(i.clone(), c));
        }
        n = n_;
        i += 2;
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;
    use std::str::FromStr;
    #[test]
    fn test_trial_division() {
        //let n = Integer::from(100u32);
        let n = Integer::from_str("100").unwrap();
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
}
