use std::collections::HashMap;

use lazy_static::lazy_static;
use rug::{ops::Pow, Integer};
#[derive(Debug, Clone)]
pub struct Factor {
    prime: Integer,
    exp: u32,
}
impl Factor {
    pub fn new(prime: Integer, exp: u32) -> Self {
        Self { prime, exp }
    }
    pub fn n(&self) -> Integer {
        let p = self.prime.clone();
        p.pow(self.exp)
    }
}

pub static N_PRECOMPUTED_PRIMES: usize = 10000;
/// First primes and their logarithms
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

    #[test]
    fn test_primes() {
        // for p in LOG_PRIMES.iter() {
        //     dbg!(p);
        // }

        assert_eq!(PRIMES.len(), 1000);
    }
}
