# FactoRS

**FactoRS** crate that provides different methods of factorization. 
[Rug](https://crates.io/crates/rug) is used to handle big integer operations.  
Some algorithms such as [Dixon] and [QuadraticSieve] have a linear algebra step. For these algorithms [nalgebra](https://docs.rs/nalgebra/latest/nalgebra/) is used to handle the operations.

## How to use
The crate provides the following structs: `TrialDivision`, `Fermat`, `Dixon`, `QuadraticSieve`. Each one implements the `Factorizer` trait:
```rust
pub trait Factorizer {
    /// Return 2 factors of a number
    fn factor(&self) -> Option<(Integer, Integer)>;
}
```
To factor a number you must pass the number into the constructor (`new::(n: Integer, ...)`) alongside extra configurations.

For `QuadraticSieve` and `Dixon` I provided 2 builder classes `QuadraticSieveBuilder` and `DixonBuilder` that should help with the configuration of the algorithm.

### Examples
Factorization using Fermat's algorithm
```rust
// Using the function
let n = Integer::from_str("5959").unwrap();
let res = fermat(&n);
assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
// Using the struct
let f = Fermat::new(n.clone());
let res = f.factor();
assert_eq!(res, Some((Integer::from(59u32), Integer::from(101u32))));
```

Factorization using Quadratic sieve algorithm
```rust
let p = Integer::from_str_radix("3507360361", 10).unwrap();
let q = Integer::from_str_radix("3916272539", 10).unwrap();
let n = p.clone() * &q; // 13735779066161426579

let qs_builder = QuadraticSieveBuilder::new(n);
let qs = qs_builder.build();
let res = qs.factor();
assert_eq!(res, Some((p, q)));
```

An example using the builder API:
```rust
let n = Integer::from(15347u32);
let builder = QuadraticSieveBuilder::new(n)
    .bound(30) 
    .factor_base(vec![2, 17, 23, 29])
    .extra_relations(1)
    .sieve_size(1000);
let qs: QuadraticSieve = builder.build();
let res = qs.factor();
assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
```

The builders provide defaults. They are computed when we call the `.build()` method.  
For example here we set only the bound parameter:
```rust
let n = Integer::from(15347u32);
let builder = QuadraticSieveBuilder::new(n)
    .bound(30) 
let qs: QuadraticSieve = builder.build();
let res = qs.factor();
assert_eq!(res, Some((Integer::from(103u32), Integer::from(149u32))));
```

For more examples check out the `examples` directory.

### Remarks
For bigger numbers (> 50 bits) use `--release` mode. Gaussian elimination takes too long for a reason unknown to me yet. 