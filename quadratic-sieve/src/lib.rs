use algorithms::gaussian_elimination_gf2;
use algorithms::{fb_factorization, tonelli_shanks};
use nalgebra::{DMatrix, DVector};
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
        let xp = tonelli_shanks(n.clone(), pi);
        roots.insert(pi, (xp, pi - xp));
    }
    roots
}

/// n =  integer we want to factor
/// s = size of the sieve interval
/// fb = factor base
pub fn quadratic_sieve(n: Integer, s: usize, factor_base: &[u64]) {
    //let k = 2 * factor_base.len(); // minimum of relations we want
    let k = factor_base.len() + 5; // minimum of relations we want

    //let mut relations = Hash
    let max_factor = *factor_base.iter().max().unwrap();

    // 1. Initialization
    // Find roots x_p^2 ≡ n (mod p)
    println!("Searching roots...");
    let mut roots: HashMap<u64, (u64, u64)> = find_roots(&n, &factor_base);
    roots.insert(2, (1, 1));
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
    loop {
        // Set interval bounds. We move in steps of `s`
        let start = n.clone().sqrt() + 1u32 + i * s as u32;
        let end = start.clone() + s as u32;
        println!("Sieving interval: [{}, {}]", start, end);
        if end > n {
            panic!("not enough relations found")
        }

        // Init sieve vector
        let mut sieve = vec![0u8; s]; // sieve array

        // Threshold.
        let t = (start.clone().pow(2) - &n).significant_bits();
        let epsilon = if t < 22 { 2 } else { t - 20 }; //- (max_factor as f64).log2().round() as u32; // ugh.
        let epsilon = 2;
        println!("epsilon = {}", epsilon);

        // Iterate through the primes and roots
        for (&p, &(xp, xp_)) in roots.iter() {
            let mut xp = Integer::from(xp);
            let mut xp_ = Integer::from(xp_);

            // Section 5.2 from
            // http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf
            // let mut p_multiple = ((start.clone() + 1u32) / p) * p;
            // let t_sub = p_multiple.clone() - &xp;
            // let t_add = p_multiple.clone() + &xp;
            // while t_sub <= end {
            //     if start < t_sub && t_sub < end {
            //         sieve[(p_multiple.clone() - &start).to_usize().unwrap()] +=
            //             LOG_PRIMES.get(&p).unwrap();
            //     }

            //     if start < t_add && t_add < end {
            //         sieve[(p_multiple.clone() + &start).to_usize().unwrap()] +=
            //             LOG_PRIMES.get(&p).unwrap();
            //     }
            //     p_multiple += p;
            // }
            // // Mark in the array

            while xp < end {
                sieve[(((xp.clone() - &start) % p + p) % p).to_usize().unwrap()] +=
                    LOG_PRIMES.get(&p).unwrap();
                xp += p;
            }
            if xp_ != 1 {
                while xp_ < end {
                    sieve[(((xp_.clone() - &start) % p + p) % p).to_usize().unwrap()] +=
                        LOG_PRIMES.get(&p).unwrap();
                    xp_ += p;
                }
            }
        }

        // There are values.
        println!(
            "{:?}",
            sieve
                .iter()
                .enumerate()
                .filter(|(_, &v)| v > 0)
                .collect::<Vec<(usize, &u8)>>()
        );

        for (i, &log_sum) in sieve.iter().enumerate() {
            // Check if the log sum is above the threshold
            if log_sum as u32 >= epsilon {
                let xi = start.clone() + i as u64;
                let qxi = xi.clone() * &xi - &n;
                // Check for B-smooth and add to relations
                if let Some(factorization) = fb_factorization(qxi.clone(), &factor_base) {
                    relations_x.push(xi);
                    relations_qx.push(qxi % &n);
                    exponents.push(factorization);
                    println!("wanted relaions: {}", k);
                    println!("relations found: {}", relations_x.len());
                }
            }
        }

        if relations_x.len() > k {
            break;
        }
        i += 1;
    }

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
    let (matrix, marked) = gaussian_elimination_gf2(matrix.clone());

    println!("{:?}", marked);
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
    println!("{:?}", independent_rows);
    println!("{:?}", dependent_rows);

    //3. Factorziation
    // For each dependent row try to form relations and factor.
    for dependent_row in dependent_rows {
        // a. Get the column index for each one in our dependent row.
        let ones_idx = matrix
            .row(dependent_row)
            .iter()
            .enumerate()
            .filter(|(_, &v)| v == 1)
            .map(|(i, _)| i)
            .collect::<Vec<usize>>();
        println!("{:?}", ones_idx);
        // b. For each one_idx find the first column in the matrix
        // that has a one on that position.
        let mut cols_idx = Vec::with_capacity(ones_idx.len() + 1);
        for idx in ones_idx {
            for &i in &independent_rows {
                if matrix[(i, idx)] == 1 {
                    cols_idx.push(i);
                }
            }
        }
        // add the dependent row
        cols_idx.push(dependent_row);

        // c. multiply all relations and compute the gcd
        let x: Integer = Integer::product(cols_idx.iter().map(|&idx| &relations_x[idx])).complete();
        let y2: Integer =
            Integer::product(cols_idx.iter().map(|&idx| &relations_qx[idx])).complete();
        let y = y2.sqrt() % &n;
        let p = (x - y).gcd(&n);

        // check if we found something
        if p != 1u32 && p != n {
            println!("p = {}", p);
            println!("q = {}", n.clone() / &p);
            println!("n % p {}", n.clone() % &p);
            break;
        } else {
            println!("p = {}", p);
            println!("n = {}", n);
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
