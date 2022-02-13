use rug::Integer;
use types::Factor;

pub fn trial_division(n: Integer) -> Vec<Factor> {
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
        i = i + 2;
    }

    res
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;
    use rug::Integer;
    #[test]
    fn test_trial_division() {
        //let n = Integer::from(100u32);
        let n = Integer::from_str("100").unwrap();

        let fs = trial_division(n.clone());
        dbg!(&fs);
        let mut n_ = Integer::from(1);
        for f in fs {
            n_ *= f.n();
        }
        assert_eq!(n, n_)
    }
}
