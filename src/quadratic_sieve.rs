//! Module that provides the quadratic sieve factorizzation method and additional functionalities.

use crate::algorithms::{big_l, fb_factorization, gaussian_elimination_gf2, tonelli_shanks};
use crate::traits::{Factorizer, LOG_PRIMES};
use crossbeam;
use indicatif::ProgressBar;
use nalgebra::DMatrix;
use rayon::prelude::*;
use rug::{ops::Pow, Complete, Integer};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::{collections::HashMap, vec};

/// A builder used to configure a [QuadraticSieve] struct
/// # Example
/// To use with computed defaults provide a Rug Integer in the constructor and call `.build()`:
/// ```rust
/// # use facto_rs::quadratic_sieve::{QuadraticSieveBuilder, QuadraticSieve};
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(15347u32);
/// let builder = QuadraticSieveBuilder::new(n).bound(30);
/// let qs: QuadraticSieve = builder.build();
/// let res = qs.factor();
/// assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
/// ```
///
/// Otherwise you can configure each parameter by calling a function with the name of the parameter:
/// ```rust
/// # use facto_rs::quadratic_sieve::{QuadraticSieveBuilder, QuadraticSieve};
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(15347u32);
/// let builder = QuadraticSieveBuilder::new(n)
///     .bound(30)
///     .factor_base(vec![2, 17, 23, 29])
///     .extra_relations(1)
///     .sieve_size(1000)
///     .verbose(true);
/// let qs: QuadraticSieve = builder.build();
/// let res = qs.factor();
/// assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
/// ```
pub struct QuadraticSieveBuilder {
    n: Integer,
    bound: Option<u64>,
    sieve_size: usize,
    factor_base: Option<Vec<u64>>,
    extra_relations: usize,
    verbose: bool,
    num_cores: Option<usize>,
}

impl QuadraticSieveBuilder {
    /// Creates a new [QuadraticSieveBuilder]. You must provide an integer because it's used in `build` to compute defaults.

    pub fn new(n: Integer) -> Self {
        Self {
            n,
            bound: None,
            factor_base: None,
            sieve_size: 10000,
            extra_relations: 5,
            verbose: true,
            num_cores: None,
        }
    }
    /// Manually set the bound.
    pub fn bound(mut self, bound: u64) -> Self {
        self.bound = Some(bound);
        self
    }
    /// Manually set the sieve size
    pub fn sieve_size(mut self, sieve_size: usize) -> Self {
        self.sieve_size = sieve_size;
        self
    }
    /// Manually set the factor_base.
    pub fn factor_base(mut self, factor_base: Vec<u64>) -> Self {
        self.factor_base = Some(factor_base);
        self
    }
    /// Manually set the number of extra_relations to search.
    pub fn extra_relations(mut self, extra_relations: usize) -> Self {
        self.extra_relations = extra_relations;
        self
    }
    /// Manually set verbosity.
    pub fn verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
    }
    /// Manually set verbosity.
    pub fn num_cores(mut self, num_cores: usize) -> Self {
        self.num_cores = Some(num_cores);
        self
    }
    /// Buld a [QuadraticSieve] using the provided configuration.
    /// The default configs:
    /// - computes the bound as sqrt(L(n)) where L(n) = exp(sqrt(ln(n) * ln(ln(n))))
    /// - creates the factor base using all primes up to `bound` where n has a square root mod p
    /// - sieve size: 10000 if 10000 > n else n
    /// - extra_relations = 5
    pub fn build(self) -> QuadraticSieve {
        let bound = self
            .bound
            .unwrap_or(2 * (big_l(self.n.clone()) as f64).sqrt().round() as u64 + 1);
        let factor_base = self
            .factor_base
            .unwrap_or(generate_factor_base_qs(&self.n, bound));
        let sieve_size = if self.n < 10000u32 {
            self.n.to_usize().unwrap() + 1
        } else {
            self.sieve_size
        };
        QuadraticSieve::new(
            self.n,
            bound,
            sieve_size,
            factor_base,
            self.extra_relations,
            self.verbose,
            self.num_cores,
        )
    }
}

/// A builder used to configure a [QuadraticSieve] struct. After creating a QuadraticSieve structure call `.factor()` to start factorizing the number.
/// When factoring, it calls the [quadratic_sieve] function internally.
/// # Example
/// To use with computed defaults provide a Rug Integer in the constructor and call `.build()`:
/// ```rust
/// # use facto_rs::quadratic_sieve::QuadraticSieve;
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(15347u32);
/// let factor_base = vec![2, 17, 23, 29];
/// let extra_relations = 2;
/// let sieve_size = 1000;
/// let bound = 30;
/// let qs = QuadraticSieve::new(n, bound, sieve_size, factor_base, extra_relations, true, None);
/// let res = qs.factor();
/// assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
/// ```
pub struct QuadraticSieve {
    n: Integer,
    bound: u64,
    factor_base: Vec<u64>,
    sieve_size: usize,
    extra_relations: usize,
    verbose: bool,
    num_cores: Option<usize>,
}

impl QuadraticSieve {
    /// Creates a new [QuadraticSieve] struct.
    pub fn new(
        n: Integer,
        bound: u64,
        sieve_size: usize,
        factor_base: Vec<u64>,
        extra_relations: usize,
        verbose: bool,
        num_cores: Option<usize>,
    ) -> Self {
        Self {
            n,
            bound,
            sieve_size,
            factor_base,
            extra_relations,
            verbose,
            num_cores,
        }
    }
}

impl Factorizer for QuadraticSieve {
    fn factor(&self) -> Option<(Integer, Integer)> {
        quadratic_sieve(
            &self.n,
            self.sieve_size,
            &self.factor_base,
            self.extra_relations,
            self.verbose,
            self.num_cores,
        )
    }
}

/// Generates all primes `p` up to a given `bound` where `n` has a square root mod `p`
/// ```rust
/// # use facto_rs::quadratic_sieve::generate_factor_base_qs;
/// # use rug::Integer;
/// let n = Integer::from(15347u32);
/// let factor_base = generate_factor_base_qs(&n, 30);
/// assert_eq!(factor_base, vec![2, 17, 23, 29]);
/// ```
pub fn generate_factor_base_qs(n: &Integer, bound: u64) -> Vec<u64> {
    let mut fb = vec![2];
    let mut p = Integer::from(3u32);
    while p < bound {
        if n.jacobi(&p) == 1 {
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

/// Sieves for possible relations. Returns:
/// - A vector of x
/// - A vector of qx
/// - A vector of the factorization of each qx
fn sieve(
    n: &Integer,
    sieve_size: usize,
    start: &Integer,
    factor_base: &[u64],
    roots: &HashMap<u64, (u64, u64)>,
    epsilon: u32,
) -> (Vec<Integer>, Vec<Integer>, Vec<BTreeMap<u64, u32>>) {
    let mut sieve_vec = vec![0u32; sieve_size + 1]; // sieve array
                                                    // Init thread vectors to keep relations
    let mut relations_x = Vec::new();
    let mut relations_qx = Vec::new();
    let mut exponents = Vec::new();
    for (&p, &(xp, xp_)) in roots {
        let xp = Integer::from(xp);
        let xp_ = Integer::from(xp_);
        let mut j = (((xp.clone() - start) % p + p) % p).to_usize().unwrap();
        while j < sieve_size {
            sieve_vec[j] += *LOG_PRIMES.get(&p).unwrap() as u32;
            j += p as usize;
        }
        if xp != 1u32 {
            let mut j = (((xp_.clone() - start) % p + p) % p).to_usize().unwrap();
            while j < sieve_size {
                sieve_vec[j] += *LOG_PRIMES.get(&p).unwrap() as u32;
                j += p as usize;
            }
        }
    }
    // Check for B-smooth numbers in the interval
    for (i, &log_sum) in sieve_vec.iter().enumerate() {
        // Check if the log sum is above the threshold
        if log_sum as u32 >= epsilon {
            let xi = start.clone() + i as u64;
            let qxi = (xi.clone() * &xi - n) % n;
            // Check for B-smooth and add to relations
            if let Some(factorization) = fb_factorization(qxi.clone(), factor_base) {
                relations_x.push(xi);
                relations_qx.push(qxi % n);
                exponents.push(factorization);
            }
        }
    }
    (relations_x, relations_qx, exponents)
}

/// Quadratic sieve factorization algorithm.
/// ```no_run
/// # use facto_rs::quadratic_sieve::quadratic_sieve;
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(15347u32);
/// let factor_base = vec![2, 17, 23, 29];
/// let extra_relations = 2;
/// let sieve_size = 1000;
/// let res = quadratic_sieve(&n, sieve_size, &factor_base, extra_relations, true, None);
/// assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
/// ```
pub fn quadratic_sieve(
    n: &Integer,
    sieve_size: usize,
    factor_base: &[u64],
    extra_relations: usize,
    verbose: bool,
    num_cores: Option<usize>,
) -> Option<(Integer, Integer)> {
    let n = n.clone();
    let k = factor_base.len() + extra_relations; // minimum of relations we want

    //let mut relations = Hash
    let max_factor = *factor_base.iter().max().unwrap();

    // 1. Initialization
    // Find roots x_p^2 ≡ n (mod p)
    if verbose {
        println!("Searching roots...");
    }
    let mut roots: HashMap<u64, (u64, u64)> = find_roots(&n, factor_base);
    roots.insert(2, (1, 1)); // add 2 to factor base
    if verbose {
        println!("Roots computed");
    }

    // 2. Sieving
    // let mut relations = HashMap::new();
    // let mut exponents = HashMap::new();
    let mut relations_x = Vec::new();
    let mut relations_qx = Vec::new();
    let mut exponents = Vec::new();

    // interval count

    let mut bar = None;
    if verbose {
        println!("Starting the sieving process");
        bar = Some(ProgressBar::new(k as u64));
        bar.as_mut().unwrap().set_position(0);
    }
    let mut i = 0;
    loop {
        // Set interval bounds. We move in steps of `s`
        let start = n.clone().sqrt() + 1u32 + i * sieve_size as u32;
        let end = start.clone() + sieve_size as u32;
        if end > n {
            panic!("not enough relations found")
        }
        // Threshold.
        let t = (start.clone().pow(2) - &n).significant_bits();
        //let epsilon = if t < 22 { 2 } else { t - 20 };
        let epsilon = t - (max_factor as f64).log2().round() as u32; // ugh.

        // Check for multithreading.
        if let Some(num_cores) = num_cores {
            let interval_length = (end.clone() - &start).to_usize().unwrap() / num_cores;
            let res: Vec<(Vec<Integer>, Vec<Integer>, Vec<BTreeMap<u64, u32>>)> = (0..num_cores)
                .into_par_iter()
                .map(|core| {
                    let sieve_size_interval = sieve_size / num_cores;
                    let start_interval = start.clone() + interval_length as u64 * core as u64;
                    // Init thread vectors to keep relations
                    sieve(
                        &n,
                        sieve_size_interval,
                        &start_interval,
                        factor_base,
                        &roots,
                        epsilon,
                    )
                })
                .collect();
            for (rel_x, rel_qx, exp) in res {
                relations_x.extend(rel_x);
                relations_qx.extend(rel_qx);
                exponents.extend(exp);
            }
        } else {
            // No multithreading
            let (rel_x, rel_qx, exp) = sieve(&n, sieve_size, &start, factor_base, &roots, epsilon);
            relations_x.extend(rel_x);
            relations_qx.extend(rel_qx);
            exponents.extend(exp);
        }

        if let Some(bar) = bar.as_mut() {
            bar.set_position(relations_x.len() as u64);
            bar.set_message(format!("Sieving interval: [{}, {}]", start, end));
        }

        // Check if we found enough relations
        if relations_x.len() > k {
            if let Some(bar) = bar {
                bar.finish();
            }
            break;
        }
        i += 1;
    }
    if verbose {
        println!("Number of relations found: {}", relations_x.len());
        println!("Building relations...");
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
    if verbose {
        println!("Starting gaussian elimination...");
    }
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
    if verbose {
        println!("Starting to search factorization through dependent rows:");
    }
    let bar = ProgressBar::new(dependent_rows.len() as u64);
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
            if verbose {
                println!("p = {}", p);
                println!("q = {}", q);
                println!("n % p {}", n % &p);
            }
            if p < q {
                return Some((p, q));
            } else {
                return Some((q, p));
            }
        }
        bar.tick();
    }
    None
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::algorithms::big_l;
    #[test]
    fn test_quadratic_sieve() {
        let p = Integer::from_str_radix("3507360361", 10).unwrap();
        let q = Integer::from_str_radix("3916272539", 10).unwrap();
        let n = p.clone() * &q; // 13735779066161426579

        let bound = 2 * (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;
        let fb = generate_factor_base_qs(&n, bound);
        let s = 100000;
        let res = quadratic_sieve(&n, s, &fb, 5, false, None);
        assert_eq!(res, Some((p, q)));
    }

    #[test]
    fn test_quadratic_sieve_structs() {
        let p = Integer::from_str_radix("3507360361", 10).unwrap();
        let q = Integer::from_str_radix("3916272539", 10).unwrap();
        let n = p.clone() * &q; // 13735779066161426579

        let qs_builder = QuadraticSieveBuilder::new(n);
        let qs = qs_builder.build();
        let res = qs.factor();
        assert_eq!(res, Some((p, q)));
    }
}
