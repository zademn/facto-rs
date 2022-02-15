//! Module that provides traits and structs used in the crate.

use std::collections::HashMap;

use lazy_static::lazy_static;
use rug::{ops::Pow, Integer};

/// A Factor struct for keeping the prime and exponent.
#[derive(Debug, Clone)]
pub struct Factor {
    prime: Integer,
    exp: u32,
}
impl Factor {
    /// Create a new [Factor];
    pub fn new(prime: Integer, exp: u32) -> Self {
        Self { prime, exp }
    }
    /// Compute `p^e`;
    pub fn n(&self) -> Integer {
        let p = self.prime.clone();
        p.pow(self.exp)
    }
    /// Compute `p^e mod n` for a given `n`;
    pub fn n_mod(&self, n: &Integer) -> Integer {
        let p = self.prime.clone();
        p.pow_mod(&Integer::from(self.exp), n).unwrap()
    }
}

/// Types that implement this trait should factor a number into 2 integers.
pub trait Factorizer {
    /// Return 2 factors of a number
    fn factor(&self) -> Option<(Integer, Integer)>;
}
/// Types that implement this trait should fully factor a number
pub trait FullFactorizer {
    /// Return the full factorization of a number.
    fn factor_full(&self) -> Vec<Factor>;
}

/// Number of precomputed primes and their log2's.
pub static N_PRECOMPUTED_PRIMES: usize = 200000;
lazy_static! {
    #[derive(Debug)]
    pub static ref PRIMES: Vec<u64> = {
        let mut t = Vec::with_capacity(N_PRECOMPUTED_PRIMES);
        let mut i = Integer::from(2u32);
        for _ in 0..N_PRECOMPUTED_PRIMES {
            t.push(i.to_u64().unwrap());
            i = i.next_prime()
        }
        t
    };
    /// First primes and their logarithms
    #[derive(Debug)]
    pub static ref LOG_PRIMES: HashMap<u64, u8> = {
        let mut t: HashMap<u64, u8> = HashMap::new();
        for &p in PRIMES.iter(){
            t.insert(p, (p as f64).log2().round() as u8);
        }
        t
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;
    #[test]
    fn test_new_factor() {
        let f = Factor::new(Integer::from(3), 4);
        assert_eq!(f.n(), 81);
    }
}
