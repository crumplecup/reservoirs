use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FG")
        .map(|x| x.age)
        .collect(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000)
        .obs(&fg) // Observations to fit.
        .period(40000.0) // Time period of individual simulations in years.
        .range(777) // Seed for rng for reproducibility.
        .runs(100); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/debris_flow_transits_ad.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.12)
        .storage_rate_gravels(0.12)
        .turnover(&318.0) // Set turnover period from the Anderson-Darling test.
        .manager(&model); // Load model parameters.

    let mut rec = fluvial.transit_times();
    utils::record(&mut rec, "/home/erik/output/transits_cdf_gr_ad.csv").unwrap();
}
