use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("/home/erik/data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep.iter().filter(|x| x.facies == "FG").map(|x| x.age).collect(); // Mean gravel deposit ages.
    let iat = Sample::read("/home/erik/data/iat.csv").unwrap(); // Mean inherited ages of charcoal in deposits.
    let ia: Vec<f64> = iat.iter().map(|x| x.age).collect(); // Vector of inherited ages from all classes of deposits.

    // Set model parameters.
    let model = ModelManager::new()
        .capture_gravels(0.0..0.4) // Range of gravel capture rates to model.
        .duration(10000) // Duration of timed() searches in hours.
        .obs(&fg) // Observations to fit.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(40) // Number of times to run the model per sampling point.
        .source_runs(1) // Number of times to run the model per gravel source.
        .storage_gravels(0.0..0.4); // Range of gravel storage rates to model.

    // Source deposits for gravels.
    let debris_flows = Reservoir::new()
        .input(&0.89).unwrap() // Input rate for debris-flow deposits from the Kolmogorov-Smirnov test.
        .output(&0.89).unwrap() // Output rate for fluvial removal of deposits from the Kolmogorov-Smirnov test.
        .inherit(&ia) // Inherited ages of charcoal in debris-flow deposits.
        .model(&model); // Load model parameters.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source(&[debris_flows.clone()]) // Set source as debris-flow deposits.
        .turnover(&208.0) // Set turnover period from the Kolmogorov-Smirnov test.
        .manager(&model); // Load model parameters.

    // Fit model to observed deposit ages for specified duration.
    // Change directory path for user, panics on invalid path
    fluvial.fit_rates_timed(&fg, "/home/erik/output/gravels_ks_40x_1000.csv").unwrap();
}