use algorithms::gaussian_elimination_gf2;
use algorithms::{fb_factorization, tonelli_shanks};
use indicatif::ProgressBar;
use nalgebra::DMatrix;
use rug::{ops::Pow, Complete, Integer};
use std::{collections::HashMap, vec};
use types::LOG_PRIMES;

pub fn generate_factor_base_qs(n: &Integer, bound: u64) -> Vec<u64> {
    let mut p = Integer::from(2u32);
    let mut fb = vec![];
    while p < bound {
        if n.legendre(&p) == 1 {
            fb.push(p.to_u64().unwrap());
        }
        p = p.next_prime();
    }
    fb
}

/// Find roots +-a_i such that  a_i^2 ≡ n (mod p_i)
/// for a given `n` and a given factor_base
fn find_roots(n: &Integer, factor_base: &[u64]) -> HashMap<u64, (u64, u64)> {
    let mut roots = HashMap::new();
    for &pi in factor_base.iter().skip(1) {
        let root = tonelli_shanks(n.clone(), pi);
        // use algorithms::pow_mod;
        // assert_eq!(
        //     pow_mod(root, 2, pi),
        //     (n.clone() % pi).to_u64().unwrap(),
        //     "{root}, {pi}, {n}"
        // );
        roots.insert(pi, (root, pi - root));
    }
    roots
}

/// n =  integer we want to factor
/// s = size of the sieve interval
/// fb = factor base
pub fn quadratic_sieve(n: Integer, s: usize, factor_base: &[u64]) -> Option<(Integer, Integer)> {
    //let k = 2 * factor_base.len(); // minimum of relations we want
    let k = factor_base.len() + 10; // minimum of relations we want

    //let mut relations = Hash
    let max_factor = *factor_base.iter().max().unwrap();

    // 1. Initialization
    // Find roots x_p^2 ≡ n (mod p)
    println!("Searching roots...");
    let mut roots: HashMap<u64, (u64, u64)> = find_roots(&n, factor_base);
    roots.insert(2, (1, 1)); // add 2 to factor base
    println!("Roots computed");

    // 2. Sieving
    // let mut relations = HashMap::new();
    // let mut exponents = HashMap::new();
    let mut relations_x = Vec::new();
    let mut relations_qx = Vec::new();
    let mut exponents = Vec::new();

    // interval count

    println!("Starting the sieving process");
    let mut i = 0;
    let bar = ProgressBar::new(k as u64);
    loop {
        // Set interval bounds. We move in steps of `s`
        let start = n.clone().sqrt() + 1u32 + i * s as u32;
        let end = start.clone() + s as u32;
        if end > n {
            panic!("not enough relations found")
        }

        // Init sieve vector
        let mut sieve = vec![0u32; s]; // sieve array
                                       // let mut sieve = HashMap::new();

        // Threshold.
        let t = (start.clone().pow(2) - &n).significant_bits();
        //let epsilon = if t < 22 { 2 } else { t - 20 };
        let epsilon = t - (max_factor as f64).log2().round() as u32; // ugh.

        // Log implementation
        // Iterate through the primes and roots
        for (&p, &(xp, xp_)) in roots.iter() {
            let xp = Integer::from(xp);
            let xp_ = Integer::from(xp_);

            let mut j = (((xp.clone() - &start) % p + p) % p).to_usize().unwrap();
            while j < s {
                sieve[j] += *LOG_PRIMES.get(&p).unwrap() as u32;
                j += p as usize;
            }
            if xp != 1u32 {
                let mut j = (((xp_.clone() - &start) % p + p) % p).to_usize().unwrap();
                while j < s {
                    sieve[j] += *LOG_PRIMES.get(&p).unwrap() as u32;
                    j += p as usize;
                }
            }
        }

        // Check for B-smooth numbers in the interval
        for (i, &log_sum) in sieve.iter().enumerate() {
            // Check if the log sum is above the threshold
            if log_sum as u32 >= epsilon {
                let xi = start.clone() + i as u64;
                let qxi = (xi.clone() * &xi - &n) % &n;
                // Check for B-smooth and add to relations
                if let Some(factorization) = fb_factorization(qxi.clone(), factor_base) {
                    relations_x.push(xi);
                    relations_qx.push(qxi % &n);
                    exponents.push(factorization);
                }
            }
        }
        bar.set_position(relations_x.len() as u64);
        bar.set_message(format!("Sieving interval: [{}, {}]", start, end));

        if relations_x.len() > k {
            bar.finish();
            break;
        }
        i += 1;
    }
    println!("Number of relations found: {}", relations_x.len());
    println!("Building relations and factorization...");

    // Transform exponends into gf2 matrix
    let nrows = exponents.len();
    let ncols = exponents[0].len();
    // Gather all them into a vec and %2 them
    let all_exponents = exponents
        .iter()
        .map(|v| v.values().map(|e| *e as u8 % 2).collect::<Vec<u8>>())
        .flatten()
        .collect::<Vec<u8>>();
    // Create the matrix
    let matrix = DMatrix::from_row_slice(nrows, ncols, &all_exponents);
    // Compute the gaussian elimination.
    // We only care about the marked columns since they contain the linearly independent vectors
    let (matrix, marked) = gaussian_elimination_gf2(matrix);

    // println!("{:?}", marked);
    // Get the index of the l.i. and l.d. vectors
    let independent_rows = marked
        .iter()
        .enumerate()
        .filter(|(_, &v)| v)
        .map(|(i, _)| i)
        .collect::<Vec<usize>>();
    let dependent_rows = marked
        .iter()
        .enumerate()
        .filter(|(_, &v)| !v)
        .map(|(i, _)| i)
        .collect::<Vec<usize>>();
    // println!("{:?}", independent_rows);
    // println!("{:?}", dependent_rows);

    //3. Factorziation
    // For each dependent row try to form relations and factor.
    println!("Starting to search factorization through dependent rows:");
    for dependent_row in dependent_rows {
        // a. Get the column index for each one in our dependent row.
        let ones_idx = matrix
            .row(dependent_row)
            .iter()
            .enumerate()
            .filter(|(_, &v)| v == 1)
            .map(|(i, _)| i)
            .collect::<Vec<usize>>();
        // b. For each column that has a one in our dependent row
        // find the first row in the matrix that has a one on that position.
        let mut rows_idx = Vec::with_capacity(ones_idx.len() + 1);
        for idx in ones_idx {
            for &i in &independent_rows {
                if matrix[(i, idx)] == 1 {
                    rows_idx.push(i);
                }
            }
        }
        // add the dependent row
        rows_idx.push(dependent_row);

        // Check dependency
        // let mut v = DVector::from_iterator(ncols, matrix.row(rows_idx[0]).iter().cloned());
        // for &row in rows_idx.iter().skip(1) {
        //     let v2 = DVector::from_iterator(ncols, matrix.row(row).iter().cloned());
        //     v = v + v2;
        // }
        // assert_eq!(v.iter().all(|e| e % 2 == 0), true, "rows are not linearly dependent")

        // c. multiply all relations and compute the gcd
        let x: Integer =
            Integer::product(rows_idx.iter().map(|&idx| &relations_x[idx])).complete() % &n;
        let y = {
            // Fast check square root
            // Collect primes and add exponents to create factorization
            let mut factorization = HashMap::new();
            for &idx in rows_idx.iter() {
                for (&p, &e) in exponents[idx].iter() {
                    *factorization.entry(p).or_insert(0) += e;
                }
            }

            let mut res = Integer::from(1u32);
            for (&p, &e) in factorization.iter() {
                let p = Integer::from(p);
                let e = Integer::from(e);
                assert_eq!(e.clone() % 2u32, 0u32, "Not a square");
                res = res * p.pow_mod(&(e / 2u32), &n).unwrap() % &n;
            }
            res
        };
        let p = if x > y {
            (x - y).gcd(&n)
        } else {
            (y - x).gcd(&n)
        };

        // check if we found something
        if p != 1u32 && p != n {
            let q = n.clone() / &p;
            println!("p = {}", p);
            println!("q = {}", q);
            println!("n % p {}", n % &p);
            return Some((p, q));
        }
    }
    None
}

#[cfg(test)]
mod tests {

    use super::*;
    use algorithms::big_l;
    #[test]
    fn test_quadratic_sieve() {
        let p = Integer::from_str_radix("3507360361", 10).unwrap();
        let q = Integer::from_str_radix("3916272539", 10).unwrap();
        let n = p.clone() * &q; // 13735779066161426579

        let bound = 2 * (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;
        let fb = generate_factor_base_qs(&n, bound);
        let s = 100000;
        let res = quadratic_sieve(n, s, &fb);
        assert_eq!(res, Some((p, q)));
    }
}
