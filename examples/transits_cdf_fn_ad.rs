use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Anderson-Darling test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let ff: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FF")
        .map(|x| x.age + 50.0)
        .collect(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&ff) // Observations to fit.
        .obs_len(&ff) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(10000); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/gravels_cdf_ad.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.175)
        .storage_rate_gravels(0.108)
        .turnover(&318.0)
        .manager(&model); // Load model parameters.

    let mut rec = fluvial.transit_times();
    utils::record(&mut rec, "/home/erik/output/fines_cdf_ad_1000.csv").unwrap();
}
