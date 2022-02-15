//! Module that provides Fermat's factorizzation method and additional functionalities.

use crate::traits::Factorizer;
use rug::{ops::Pow, Integer};

/// Structure to handle Fermat's factorization.
/// After creating a [Fermat] structure call `.factor()` to start factorizing the number.
/// When factoring, it calls the [fermat] function internally.
/// # Example
/// ```rust
/// # use facto_rs::fermat::Fermat;
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(5959u32);
/// let fermat_struct = Fermat::new(n);
/// let res = fermat_struct.factor();
/// assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
/// ```
pub struct Fermat {
    n: Integer,
}
impl Fermat {
    /// Creates a new [Fermat] struct.
    pub fn new(n: Integer) -> Self {
        Self { n }
    }
}
impl Factorizer for Fermat {
    fn factor(&self) -> Option<(Integer, Integer)> {
        fermat(&self.n)
    }
}

/// Fermat's factorization method.
/// # Example
/// ```no_run
/// # use facto_rs::fermat::fermat;
/// # use rug::Integer;
/// let n = Integer::from(5959u32);
/// let res = fermat(&n);
/// assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
/// ```
pub fn fermat(n: &Integer) -> Option<(Integer, Integer)> {
    let n = n.clone();
    let mut a: Integer = n.clone().sqrt() + 1;
    let mut b2 = a.clone().pow(2) - &n;
    while !b2.is_perfect_square() {
        b2 = b2 + 2u32 * &a + 1u32;
        a += 1u32;
    }
    let b = b2.sqrt();
    Some((a.clone() - &b, a + &b))
}

#[allow(unused)]
fn fermat_big(n: Integer) -> (Integer, Integer) {
    let a: Integer = n.clone().sqrt() + 1u32;
    let mut i: Integer = Integer::from(0);
    loop {
        let x: Integer = a.clone() + &i;
        let y2 = x.clone().pow(2) - &n;
        let y = y2.sqrt();
        if (x.clone() - &y) != 1 {
            return (x.clone() - &y, x + y);
        }
        i += 1;
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;
    use std::str::FromStr;
    #[test]
    fn test_fermat() {
        let n = Integer::from_str("5959").unwrap();
        let res = fermat(&n);
        assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
        let f = Fermat::new(n.clone());
        let res = f.factor();
        assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
    }
}
