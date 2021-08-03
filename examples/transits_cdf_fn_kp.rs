use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kuiper test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let ff: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FF")
        .map(|x| x.age + 50.0)
        .collect(); // Mean gravel deposit ages.
    let df = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age + 50.0)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&ff) // Observations to fit.
        .obs_len(&df) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(777) // Seed for rng for reproducibility.
        .runs(10000); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/gravels_transits_ks.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.301)
        .storage_rate_gravels(0.04)
        .turnover(&191.0) // Set turnover period from the Kuiper test.
        .manager(&model); // Load model parameters.

    let mut rec = fluvial.transit_times();
    utils::record(&mut rec, "/home/erik/output/transits_cdf_fn_kp.csv").unwrap();
}
