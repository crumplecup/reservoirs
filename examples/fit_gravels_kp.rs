use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Chi-squared test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FG")
        .map(|x| x.age + 50.0)
        .collect(); // Mean gravel deposit ages.
    let df = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age + 50.0)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .capture_gravels(0.0..1.0) // Range of gravel capture rates to model.
        .duration(10000) // Duration of timed() searches in hours.
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&fg) // Observations to fit.
        .obs_len(&df) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(1000) // Number of times to run the model per sampling point.
        .storage_gravels(0.0..1.0); // Range of gravel storage rates to model.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/debris_flow_transits_kp.csv")
        .unwrap() // Set source as debris-flow deposits.
        .turnover(&191.0) // Set turnover period from the Chi-squared test.
        .manager(&model); // Load model parameters.

    // Fit model to observed deposit ages for specified duration.
    // Change directory path for user, panics on invalid path

    fluvial
        .fit_rates_timed(&fg, "/home/erik/output/gravels_kp_1kx_1000.csv")
        .unwrap();
}
