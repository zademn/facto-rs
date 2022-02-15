#![warn(missing_docs)]

// # facto_rs
//! **facto_rs** crate that provides different methods of factorization.
//! [Rug](https://crates.io/crates/rug) is used to handle big integer operations.
//! Some algorithms such as [Dixon] and [QuadraticSieve] have a linear algebra step. For that  algorithms [nalgebra](https://docs.rs/nalgebra/latest/nalgebra/) is used to handle the operations.
//! 
pub mod algorithms;
pub mod dixon;
pub mod fermat;
pub mod quadratic_sieve;
pub mod traits;
pub mod trial_division;

pub use dixon::{Dixon, DixonBuilder};
pub use fermat::Fermat;
pub use quadratic_sieve::{QuadraticSieve, QuadraticSieveBuilder};
pub use traits::Factorizer;
pub use trial_division::TrialDivision;
