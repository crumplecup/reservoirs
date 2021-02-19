/*!
# Reservoirs - A library for modeling Bolin & Rodhe reservoirs.
Bolin & Rodhe (1973) describe methods for characterizing the turnover rate of mass accumulating in a reservoir.
The distribution of ages of particles in the reservoir constrain the possible input and output rates that could
produce the observed record.  The functions in this crate allow the user to compare synthetic accumulation records to
observed records using the K-S and Kuiper statistics, to determine the best-fitting input/output pair for an observed record.

In my research at Oregon State University, I estimate the transit times of stream sediments moving through headwater valleys
of the Coast Range by fitting reservoir models to a record of charcoal ages sampled from stream bank deposits.  Inherited age
refers to the age of charcoal when it enters a stream deposit.  If we do not account for inherited age in the model, then transit times
become artificially inflated.  I have added an inherited age capacity to reservoirs, and while this is not a traditional feature
of Bolin & Rodhe reservoirs, it is useful for dealing with charcoal ages.

This library includes the full code base used to estimate transit times for my ongoing dissertation, published here in the interest
of academic transparency.

 - Please direct questions, comments or insults to the [github repository](https://github.com/crumplecup/reservoirs).
 - View the crate documentation on [docs.rs](https://docs.rs/reservoirs/).

 ## Quick Start

To use reservoirs, add it to your `Cargo.toml`
```toml
[dependencies]
reservoirs = "^0.1.1"
```

Create reservoirs using a builder pattern.  First make a 'blank' reservoir using new, then assign it
features using the input, output and inherit methods.

```rust
use reservoirs::prelude::*;

// build step by step
let mut res = Reservoir::new();
res = res.input(&0.78)?;
res = res.output(&0.78)?;
res = res.inherit(&vec![10.0, 20.0, 27.0, 100.3, 7000.0, 10000.0]);

// or inline, same result
let res_b = Reservoir::new()
    .input(&0.78)?
    .output(&0.78)?
    .inherit(&vec![10.0, 20.0, 27.0, 100.3, 7000.0, 10000.0]);

assert_eq!(res, res_b);

```
*/

#![warn(missing_docs)]
pub mod reservoir;
pub mod utils;
pub mod plot;

pub mod prelude {
    // pub use crate::utils;
    pub use crate::reservoir::{Reservoir, Sample};
    pub use crate::plot;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
