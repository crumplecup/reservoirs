use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kuiper test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FG")
        .map(|x| x.age + 50.0)
        .collect(); // Mean gravel deposit ages.


    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&fg) // Observations to fit.
        .obs_len(&fg) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(10000); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/debris_flow_transits_kp.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.10)
        .storage_rate_gravels(0.011)
        .turnover(&191.0) // Set turnover period from the Kuiper test.
        .manager(&model); // Load model parameters.

    let mut rec = fluvial.transit_times();
    utils::record(&mut rec, "/home/erik/output/gravels_cdf_kp.csv").unwrap();
}
