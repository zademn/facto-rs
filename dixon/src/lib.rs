use algorithms::{fb_factorization, gaussian_elimination_gf2};
use indicatif::ProgressBar;
use nalgebra::DMatrix;
use nalgebra::RowDVector;
use rug::Complete;
use rug::Integer;
use std::collections::HashMap;

pub fn generate_factor_base(bound: u64) -> Vec<u64> {
    let mut p = Integer::from(2u32);
    let mut fb = vec![];
    while p < bound {
        fb.push(p.to_u64().unwrap());
        p = p.next_prime();
    }
    fb
}

/// Algorithm from here:
///  https://dspace.cvut.cz/bitstream/handle/10467/94585/F8-DP-2021-Vladyka-Ondrej-DP_Vladyka_Ondrej_2021.pdf?sequence=-1&isAllowed=y
pub fn dixon(n: Integer, factor_base: &[u64]) -> Option<(Integer, Integer)> {
    // number of relations we want
    let k = 2 * factor_base.len();
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
    fn test_dixon() {
        let p = Integer::from_str_radix("42509", 10).unwrap();
        let q = Integer::from_str_radix("63299", 10).unwrap();
        let n = p.clone() * &q;

        let bound = 2 * (big_l(n.clone()) as f64).sqrt().round() as u64 + 1;
        let fb = generate_factor_base(bound);
        let res = dixon(n, &fb);
        assert_eq!(res, Some((p, q)));
    }
}
