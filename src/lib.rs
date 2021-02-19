//! A library for modeling Bolin & Rodhes reservoirs.

#![warn(missing_docs)]
pub mod reservoir;
pub mod utils;
pub mod plot;

pub mod prelude {
    // pub use crate::utils;
    pub use crate::reservoir::{Reservoir};
    // pub use crate::plot;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
