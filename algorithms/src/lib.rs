use std::collections::HashMap;

use nalgebra::{DMatrix, DVector};
use rug::{Complete, Integer};

/// xor 2 arrays of u8 against each other.
fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    a.iter().zip(b).map(|(aa, bb)| aa ^ bb).collect()
}

/// From here
/// https://stackoverflow.com/questions/45918104/how-to-do-arithmetic-modulo-another-number-without-overflow
pub fn mul_mod(mut x: u64, mut y: u64, n: u64) -> u64 {
    let msb = 0x8000_0000_0000_0000;
    let mut d = 0;
    let mp2 = n >> 1;
    x %= n;
    y %= n;

    if n & msb == 0 {
        for _ in 0..64 {
            d = if d > mp2 { (d << 1) - n } else { d << 1 };
            if x & msb != 0 {
                d += y;
            }
            if d >= n {
                d -= n;
            }
            x <<= 1;
        }
        d
    } else {
        for _ in 0..64 {
            d = if d > mp2 {
                d.wrapping_shl(1).wrapping_sub(n)
            } else {
                // the case d == m && x == 0 is taken care of
                // after the end of the loop
                d << 1
            };
            if x & msb != 0 {
                let (mut d1, overflow) = d.overflowing_add(y);
                if overflow {
                    d1 = d1.wrapping_sub(n);
                }
                d = if d1 >= n { d1 - n } else { d1 };
            }
            x <<= 1;
        }
        if d >= n {
            d - n
        } else {
            d
        }
    }
}

/// Simple left to right exponentiation
pub fn pow_mod(a: u64, e: u64, n: u64) -> u64 {
    if e == 0 {
        return 1;
    }
    if e == 1 {
        return a % n;
    }
    let l = 64 - e.leading_zeros();
    let mut x = a;
    for i in (0..l - 1).rev() {
        x = mul_mod(x, x, n);
        if (e >> i) & 1 == 1 {
            x = mul_mod(a, x, n);
        }
    }
    x
}

/// Removes the `factor` from a number `n` by repeated division
/// Returns the remainder and the count.
pub fn remove_factor(mut n: u64, factor: u64) -> (u64, u64) {
    let mut count = 0;
    while n % factor == 0 {
        n = n / factor;
        count += 1;
    }
    return (n, count);
}

/// Computes the prime factorization of an integer `n` given a factor base by repeated division
/// If at the end `n!=1` then the factor base is toos small and the function will return None
pub fn fb_factorization(mut n: Integer, fb: &[u64]) -> Option<HashMap<u64, u32>> {
    let mut res = HashMap::new();
    for &p in fb {
        let (r, c) = n.remove_factor(&Integer::from(p));
        res.insert(p, c);
        n = r;
    }
    if n == 1u32 {
        return Some(res);
    }
    None
}
/// Given a positive odd integer m and integer a the algorithms returns jacobi symbol (a / m)
/// Algorithm 2.3.5 from here:
/// http://thales.doa.fmph.uniba.sk/macaj/skola/teoriapoli/primes.pdf
pub fn jacobi(mut a: u64, mut m: u64) -> i64 {
    // 1. Reduction loops
    a = a % m;
    let mut t = 1;
    while a != 0 {
        while (a & 1) == 0 {
            // even
            a = a / 2;
            let r = m % 8;
            if r == 3 || r == 5 {
                t = -t;
            }
        }
        std::mem::swap(&mut a, &mut m);
        if a % 4 == 3 && m % 4 == 3 {
            t = -t;
        }
        a = a % m;
    }
    // 2. Termination
    if m == 1 {
        return t;
    }
    return 0;
}

/// Given an odd prime p and an integer a with jacobi(a, p) = 1
/// the algorithm  a solution x to x^2 ≡ a (mod p)
/// Resources:
/// Algorithm 2.3.8 from http://thales.doa.fmph.uniba.sk/macaj/skola/teoriapoli/primes.pdf
/// Wikipedia page: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
pub fn tonelli_shanks(a: Integer, p: u64) -> u64 {
    let a = (a % &p).to_u64().unwrap();
    assert_eq!(jacobi(a, p), 1);
    // Check easy cases
    let r = p & 8;
    if r == 3 || r == 7 {
        let x = pow_mod(a, (p + 1) / 4, p);
        return x;
    }
    if r == 5 {
        let mut x = pow_mod(a, (p + 3) / 8, p);
        let c = mul_mod(x, x, p);
        if c != a {
            x = mul_mod(x, 1 << ((p - 1) / 4), p)
        }
        return x;
    }

    // Find randominteger d in [2, p-1] with jacobi(d, p) = -1
    let mut d = 2;
    while jacobi(d, p) != -1 {
        d += 1;
    }

    // Write p - 1 = 2^s * t with t odd.
    let (t, s) = remove_factor(p - 1, 2);

    let big_a = pow_mod(a, t, p);
    let big_d = pow_mod(d, t, p);
    let mut m = 0;
    for i in 0..s {
        if pow_mod(mul_mod(big_a, pow_mod(big_d, m, p), p), 1 << (s - 1 - i), p) + 1 == 0 % p {
            m = m + (1 << i);
        }
    }
    let x = mul_mod(pow_mod(a, (t + 1) / 2, p), pow_mod(big_d, m / 2, p), p);
    x
}

/// Given an odd prime p and an integer a with jacobi(a, p) = 1
/// the algorithm  a solution x to x^2 ≡ a (mod p)
/// Resources:
/// Algorithm 2.3.8 from http://thales.doa.fmph.uniba.sk/macaj/skola/teoriapoli/primes.pdf
/// Wikipedia page: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
fn tonelli_shanks_big(a: Integer, p: Integer) -> Integer {
    assert_eq!(a.jacobi(&p), 1);

    let one = Integer::from(1u32);
    let two = Integer::from(2u32);

    // Check easy cases
    let r = p.mod_u(8);
    if r == 3 || r == 7 {
        let e = (&p + 1u32).complete() / 4u32;
        let x = a.pow_mod(&e, &p).unwrap();
        return x;
    }
    if r == 5 {
        let e = (p.clone() + 3u32) / 8u32;
        let mut x = a.clone().pow_mod(&e, &p).unwrap();
        let c = x.clone().pow_mod(&two, &p).unwrap();
        if c != a {
            let e = (&p - 1u32).complete() / 4u32;
            x = x * two.pow_mod(&e, &p).unwrap() % &p;
        }
        return x;
    }

    // Find randominteger d in [2, p-1] with jacobi(d, p) = -1
    let mut d = Integer::from(2u32);
    while d.jacobi(&p) != -1 {
        d += 1u32;
    }

    // Write p - 1 = 2^s * t with t odd.
    let (t, s) = (p.clone() - 1u32).remove_factor(&two);

    let big_a = a.clone().pow_mod(&t, &p).unwrap();
    let big_d = d.clone().pow_mod(&t, &p).unwrap();
    let mut m = Integer::ZERO;
    for i in 0..s {
        let e = Integer::from(1u32) << (s - 1 - i); // 2 ^ (s - 1 - i)
        let ad = &big_a * big_d.clone().pow_mod(&m, &p).unwrap();
        if ((ad).pow_mod(&e, &p).unwrap() + 1u32) % &p == 0u32 {
            m = &m + (one.clone() << i);
        }
    }
    let e = (t + 1u32) / 2u32;
    let x = a.pow_mod(&e, &p).unwrap() * d.pow_mod(&(m / 2), &p).unwrap() % &p;
    x
}

/// Fast gaussian elimination in GF(2) implemented after
/// https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
/// Takes ownership of the given matrix.
/// Returns a tuple composed of the modified matrix and a Vec<bool> which marks the rows with pivots
pub fn gaussian_elimination_gf2(mut m: DMatrix<u8>) -> (DMatrix<u8>, Vec<bool>) {
    let mut marks = vec![false; m.nrows()];
    for j in 0..m.ncols() {
        // println!("{}", m);
        // println!("{:?}", marks);

        // Iterate through columns
        for i in 0..m.nrows() {
            // Search for pivot
            if m[(i, j)] == 1 {
                // println!("Pivot found: {:?}", (i + 1, j + 1));
                marks[i] = true; // if we found a pivot found mark row i
                for k in 0..m.ncols() {
                    if k != j {
                        if m[(i, k)] == 1 {
                            // println!("Other row found: {:?}", (i + 1, k + 1));
                            let t = DVector::from_vec(xor(
                                m.column(k).as_slice(),
                                m.column(j).as_slice(),
                            ));
                            m.set_column(k, &t);
                            // println!("{}", m);
                            // println!("{:?}", marks);
                        }
                    }
                }
                break; // no reason to search for another pivot so break
            }
        }
    }
    (m, marks)
}

///2^sqrt(log2(N) * log2(log2(N)))
pub fn L(n: Integer) -> u64 {
    let ln_n = n.significant_bits();
    let ln_ln_n = Integer::from(ln_n).significant_bits();
    let e = (Integer::from(ln_n) * ln_ln_n).sqrt().to_u32().unwrap();
    return 1 << e; // 2^e
}
#[cfg(test)]
mod tests {
    use std::num::IntErrorKind;

    use nalgebra::{DMatrix, RowDVector};

    use super::*;

    #[test]
    fn test_xor() {
        let bv1 = vec![0, 0, 1, 1];
        let bv2 = &vec![0, 1, 0, 1];

        let xor_res = xor(&bv1, &bv2);
        dbg!(&xor_res);
        assert_eq!(xor_res, vec![0, 1, 1, 0]);
    }
    #[test]
    fn test_pow_mod() {
        let a = 13;
        let e = 1234567890;
        let n = 7831289736129812936;
        let res = pow_mod(a, e, n);
        dbg!(res);
        assert_eq!(res, 5138850779699736345);
    }

    #[test]
    fn test_jacobi() {
        let ks = [1, 1, 1, 2, 2, 2, 3, 3, 3, 23, 23, 23, 30, 30, 30];
        let ns = [1, 3, 5, 1, 3, 5, 9, 11, 13, 1, 17, 31, 1, 3, 5];
        let rs = [1, 1, 1, 1, -1, -1, 0, 1, 1, 1, -1, -1, 1, 0, 0];
        for ((&k, n), r) in ks.iter().zip(ns).zip(rs) {
            assert_eq!(jacobi(k, n), r);
        }
    }
    #[test]
    fn test_remove_factor() {
        let v = 2u64.pow(3) * 3u64.pow(4) * 5u64.pow(10);
        let (res, c1) = remove_factor(v, 2);
        let (res, c2) = remove_factor(res, 3);
        let (res, c3) = remove_factor(res, 5);
        let (res, c4) = remove_factor(res, 2);

        assert_eq!(c1, 3);
        assert_eq!(c2, 4);
        assert_eq!(c3, 10);
        assert_eq!(c4, 0);
        assert_eq!(res, 1);
    }
    #[test]
    fn test_mul_mod() {
        let half = 1 << 16;
        let max = std::u64::MAX;

        assert_eq!(mul_mod(0, 0, 2), 0);
        assert_eq!(mul_mod(1, 0, 2), 0);
        assert_eq!(mul_mod(0, 1, 2), 0);
        assert_eq!(mul_mod(1, 1, 2), 1);
        assert_eq!(mul_mod(42, 1, 2), 0);
        assert_eq!(mul_mod(1, 42, 2), 0);
        assert_eq!(mul_mod(42, 42, 2), 0);
        assert_eq!(mul_mod(42, 42, 42), 0);
        assert_eq!(mul_mod(42, 42, 41), 1);
        assert_eq!(mul_mod(1239876, 2948635, 234897), 163320);

        assert_eq!(mul_mod(1239876, 2948635, half), 18476);
        assert_eq!(mul_mod(half, half, half), 0);
        assert_eq!(mul_mod(half + 1, half + 1, half), 1);

        assert_eq!(mul_mod(max, max, max), 0);
        assert_eq!(mul_mod(1239876, 2948635, max), 3655941769260);
        assert_eq!(mul_mod(1239876, max, max), 0);
        assert_eq!(mul_mod(1239876, max - 1, max), max - 1239876);
        assert_eq!(mul_mod(max, 2948635, max), 0);
        assert_eq!(mul_mod(max - 1, 2948635, max), max - 2948635);
        assert_eq!(mul_mod(max - 1, max - 1, max), 1);
        assert_eq!(mul_mod(2, max / 2, max - 1), 0);
    }

    #[test]
    fn test_tonelli_shanks() {
        let a = Integer::from(44u32);
        let p = 83;
        assert_eq!(tonelli_shanks(a, p), 25);
    }

    #[test]
    fn test_fb_factorization() {
        let n = 2u64.pow(3) * 3u64.pow(4) * 5u64.pow(10);
        let mut good = HashMap::new();
        good.insert(2, 3);
        good.insert(3, 4);
        good.insert(5, 10);
        let n = Integer::from(n);
        let res = fb_factorization(n.clone(), &[2, 3, 5]);

        assert_eq!(res, Some(good));
        let res = fb_factorization(n.clone(), &[2, 3, 5]);

        assert_eq!(res, None);
    }

    #[test]
    fn test_gaussian_elimination() {
        let m = DMatrix::from_rows(&[
            RowDVector::from_row_slice(&[1u8, 1, 0, 0]),
            RowDVector::from_row_slice(&[1, 1, 0, 1]),
            RowDVector::from_row_slice(&[0, 1, 1, 1]),
            RowDVector::from_row_slice(&[0, 0, 1, 0]),
            RowDVector::from_row_slice(&[0, 0, 0, 1]),
        ]);

        let m2 = DMatrix::from_rows(&[
            RowDVector::from_row_slice(&[1u8, 0, 0, 0]),
            RowDVector::from_row_slice(&[0, 0, 0, 1]),
            RowDVector::from_row_slice(&[0, 1, 0, 0]),
            RowDVector::from_row_slice(&[0, 0, 1, 0]),
            RowDVector::from_row_slice(&[1, 0, 0, 1]),
        ]);
        let (res, marks) = gaussian_elimination_gf2(m);
        assert_eq!(marks, [true, true, true, true, false]);
        assert_eq!(m2, res);
    }
}
