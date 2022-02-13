use algorithms::{fb_factorization, gaussian_elimination_gf2, L};
use nalgebra::DMatrix;
use nalgebra::{DVector, RowDVector};
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
pub fn dixon(n: Integer, factor_base: &[u64]) /*-> (Integer, Integer)*/
{
    // number of relations we want
    let k = 2 * factor_base.len();
    // Relation and exponents arrays
    let mut relations_x = Vec::with_capacity(k);
    let mut relations_y2 = Vec::with_capacity(k);
    let mut exponents = Vec::with_capacity(k);

    let n_sqrt = n.clone().sqrt();
    // search for possible relations
    let mut j: u32 = 1;
    loop {
        let x_j = n_sqrt.clone() + j;
        let y2_j = x_j.clone() * &x_j % &n;

        // Check if y^2 is b smooth. If it is add the relation
        if let Some(factorization) = fb_factorization(y2_j.clone(), &factor_base) {
            // relations.insert(xj.clone(), yj);
            // exponents.insert(xj, factorization);
            relations_x.push(x_j);
            relations_y2.push(y2_j);
            exponents.push(factorization);
            println!("wanted relaions: {}", k);
            println!("relations found: {}", relations_x.len());
        };
        j += 1;
        if relations_x.len() >= k {
            break;
        }
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
            Integer::product(cols_idx.iter().map(|&idx| &relations_y2[idx])).complete();
        let y = y2.sqrt() % &n;
        let p = (x - y).gcd(&n);

        // check if we found something
        if p != 1 && p != n {
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

    #[test]
    fn test_dixon() {
        let p = Integer::from_str_radix("17355072078961481401", 10).unwrap();
        let q = Integer::from_str_radix("17355072078961481401", 10).unwrap();
        let n = p.clone() * &q;

        let bound = (L(n.clone()) as f64).sqrt().round() as u64 + 1;
        let fb = generate_factor_base(10000000);
        dixon(n, &fb);
        assert!(false);
    }
}
