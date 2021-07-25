use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let debris_flows = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.
    let index = 0..20000i32;

    let mut cdf = utils::cdf_rng(&debris_flows, &index);
    utils::record(&mut cdf, "/home/erik/output/debris_flow_deposits_cdf.csv").unwrap();
}
