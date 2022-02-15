//! Module that provides Dixon's factorizzation method and additional functionalities

use crate::algorithms::{big_l, fb_factorization, gaussian_elimination_gf2};
use crate::traits::Factorizer;
use indicatif::ProgressBar;
use nalgebra::DMatrix;
use nalgebra::RowDVector;
use rug::Complete;
use rug::Integer;
use std::collections::HashMap;

/// A builder used to configure a Dixon struct
/// # Example
/// To use with computed defaults provide a Rug Integer in the constructor and call `.build()`:
/// ```rust
/// # use facto_rs::dixon::{DixonBuilder, Dixon};
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(84923u32);
/// let builder = DixonBuilder::new(n);
/// let dixon_struct: Dixon = builder.build();
/// let res = dixon_struct.factor();
/// assert_eq!(res, Some((Integer::from(163u32), Integer::from(521u32))));
/// ```
///
/// Otherwise you can configure each parameter by calling a function with the name of the parameter:
/// ```rust
/// # use facto_rs::dixon::{DixonBuilder, Dixon};
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(84923u32);
/// let builder = DixonBuilder::new(n)
///     .bound(8)
///     .factor_base(vec![2, 3, 5, 7])
///     .extra_relations(3);
/// let dixon_struct: Dixon = builder.build();
/// let res = dixon_struct.factor();
/// assert_eq!(res, Some((Integer::from(163u32), Integer::from(521u32))));
/// ```
pub struct DixonBuilder {
    n: Integer,
    bound: Option<u64>,
    factor_base: Option<Vec<u64>>,
    extra_relations: Option<usize>,
}

impl DixonBuilder {
    /// Creates a new [DixonBuilder]. You must provide an integer because it's used in `build` to compute defaults.
    pub fn new(n: Integer) -> Self {
        Self {
            n,
            bound: None,
            factor_base: None,
            extra_relations: None,
        }
    }
    /// Manually set the bound.
    pub fn bound(mut self, bound: u64) -> DixonBuilder {
        self.bound = Some(bound);
        self
    }
    /// Manually set the factor_base.
    pub fn factor_base(mut self, factor_base: Vec<u64>) -> DixonBuilder {
        self.factor_base = Some(factor_base);
        self
    }
    /// Manually set the number of extra_relations to search.
    pub fn extra_relations(mut self, extra_relations: usize) -> DixonBuilder {
        self.extra_relations = Some(extra_relations);
        self
    }
    /// Buld a [Dixon] using the provided configuration.
    /// The default configs:
    /// - computes the bound as sqrt(L(n)) where L(n) = exp(sqrt(ln(n) * ln(ln(n))))
    /// - creates the factor base using all primes up to `bound`
    /// - extra_relations = 2
    pub fn build(self) -> Dixon {
        let bound = self
            .bound
            .unwrap_or((big_l(self.n.clone()) as f64).sqrt().round() as u64 + 1);
        let factor_base = self.factor_base.unwrap_or(generate_factor_base(bound));
        let extra_relations = self.extra_relations.unwrap_or(2);
        Dixon::new(self.n, bound, factor_base, extra_relations)
    }
}

/// Structure to handle Dixon's factorization. 
/// After creating a Dixon structure call `.factor()` to start factorizing the number.
/// When factoring, it calls the [dixon] function internally.
/// # Example
/// To use with computed defaults provide a Rug Integer in the constructor and call `.build()`:
/// ```rust
/// # use facto_rs::dixon::Dixon;
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(84923u32);
/// let factor_base = vec![2, 3, 5, 7];
/// let extra_relations = 2;
/// let bound = 8;
/// let dixon_struct = Dixon::new(n, bound, factor_base, extra_relations);
/// let res = dixon_struct.factor();
/// assert_eq!(res, Some((Integer::from(163u32), Integer::from(521u32))));
/// ```
pub struct Dixon {
    n: Integer,
    bound: u64,
    factor_base: Vec<u64>,
    extra_relations: usize,
}
impl Dixon {
    /// Creates a new [Dixon] struct.
    pub fn new(n: Integer, bound: u64, factor_base: Vec<u64>, extra_relations: usize) -> Self {
        Self {
            n,
            bound,
            factor_base,
            extra_relations,
        }
    }
}

impl Factorizer for Dixon {
    fn factor(&self) -> Option<(Integer, Integer)> {
        dixon(&self.n, &self.factor_base, self.extra_relations)
    }
}

/// Generates all primes up to a given bound.
/// ```rust
/// # use facto_rs::dixon::generate_factor_base;
/// let factor_base = generate_factor_base(8);
/// assert_eq!(factor_base, vec![2, 3, 5, 7]);
/// ```
pub fn generate_factor_base(bound: u64) -> Vec<u64> {
    let mut p = Integer::from(2u32);
    let mut fb = vec![];
    while p < bound {
        fb.push(p.to_u64().unwrap());
        p = p.next_prime();
    }
    fb
}

/// Dixon factorization algorithm.
/// ```no_run
/// # use facto_rs::dixon::dixon;
/// # use facto_rs::traits::Factorizer;
/// # use rug::Integer;
/// let n = Integer::from(84923u32);
/// let factor_base = vec![2, 3, 5, 7];
/// let extra_relations = 2;
/// let res = dixon(&n, &factor_base, extra_relations);
/// assert_eq!(res, Some((Integer::from(163u32), Integer::from(521u32))));
/// ```
// # Algorithm from here:
// # https://dspace.cvut.cz/bitstream/handle/10467/94585/F8-DP-2021-Vladyka-Ondrej-DP_Vladyka_Ondrej_2021.pdf?sequence=-1&isAllowed=y
pub fn dixon(
    n: &Integer,
    factor_base: &[u64],
    extra_relations: usize,
) -> Option<(Integer, Integer)> {
    let n = n.clone();
    // number of relations we want
    let k = factor_base.len() + extra_relations;
    // Relation and exponents arrays
    let mut relations_x = Vec::with_capacity(k);
    let mut relations_y2 = Vec::with_capacity(k);
    let mut exponents = Vec::with_capacity(k);

    let n_sqrt = n.clone().sqrt();
    // search for possible relations
    let bar = ProgressBar::new(k as u64);
    let mut j: u32 = 1;
    loop {
        let x_j = n_sqrt.clone() + j;
        let y2_j = x_j.clone() * &x_j % &n;

        // Check if y^2 is b smooth. If it is, add the relation.
        if let Some(factorization) = fb_factorization(y2_j.clone(), factor_base) {
            relations_x.push(x_j);
            relations_y2.push(y2_j);
            exponents.push(factorization);
            // println!("wanted relaions: {}", k);
            // println!("relations found: {}", relations_x.len());
            bar.inc(1);
        };
        j += 1;
        if relations_x.len() >= k {
            bar.finish();
            break;
        }
    }
    println!("Number of relations found: {}", relations_x.len());

    // Transform exponends into gf2 matrix
    // let nrows = exponents.len();
    // let ncols = exponents[0].len();
    // Gather all them into a vec and %2 them
    // let all_exponents = exponents
    //     .iter()
    //     .map(|v| v.values().map(|e| (*e % 2) as u8).collect::<Vec<u8>>())
    //     .flatten()
    //     .collect::<Vec<u8>>();
    // // Create the matrix
    // let matrix = DMatrix::from_row_slice(nrows, ncols, &all_exponents);

    let matrix = DMatrix::from_rows(
        &exponents
            .iter()
            .map(|v| {
                RowDVector::from_row_slice(&v.values().map(|e| (*e % 2) as u8).collect::<Vec<u8>>())
            })
            .collect::<Vec<RowDVector<u8>>>(),
    );

    // Compute the gaussian elimination.
    // We only care about the marked columns since they contain the linearly independent vectors
    let (matrix, marked) = gaussian_elimination_gf2(matrix);

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
            if p < q {
                return Some((p, q));
            } else {
                return Some((q, p));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::algorithms::big_l;
    #[test]
    fn test_dixon() {
        let p = Integer::from_str_radix("42509", 10).unwrap();
        let q = Integer::from_str_radix("63299", 10).unwrap();
        let n = p.clone() * &q;

        let bound = 2 * (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;
        let fb = generate_factor_base(bound);
        let res = dixon(&n, &fb, 10);
        assert_eq!(res, Some((p, q)));
    }
    #[test]
    fn test_dixon_structs() {
        let p = Integer::from_str_radix("42509", 10).unwrap();
        let q = Integer::from_str_radix("63299", 10).unwrap();
        let n = p.clone() * &q;

        let dixon_builder = DixonBuilder::new(n);
        let dixon_struct = dixon_builder.build();
        let res = dixon_struct.factor();
        assert_eq!(res, Some((p, q)));
    }
}
