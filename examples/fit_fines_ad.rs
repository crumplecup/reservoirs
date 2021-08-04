use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Anderson-Darling test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("/home/erik/data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let ff: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FF")
        .map(|x| x.age)
        .collect(); // Mean fines deposit ages.
    let iat = Sample::read("/home/erik/data/iat.csv").unwrap(); // Mean inherited ages of charcoal in deposits.
    let ia: Vec<f64> = iat.iter().map(|x| x.age).collect(); // Vector of inherited ages from all classes of deposits.

    // Set model parameters.
    let model = ModelManager::new()
        .capture_fines(0.0..0.5) // Range of fines capture rates to model.
        .capture_gravels(0.29..0.30) // Range of gravel capture rates to model, set to AD test.
        .duration(10000) // Duration of timed() searches in hours.
        .fines(true) // Fit for fines deposits.
        .obs(&ff) // Observations to fit.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(200) // Number of times to run the model per sampling point.
        .source_runs(100) // Number of times to run the model per gravel source.
        .storage_fines(0.0..0.5) // Range of fines storage rates to model.
        .storage_gravels(0.28..0.29); // Range of gravel storage rates to model, set to AD test.

    // Source deposits for gravels.
    let _debris_flows = Reservoir::new()
        .input(&0.35)
        .unwrap() // Input rate for debris-flow deposits from the Anderson-Darling test.
        .output(&0.35)
        .unwrap() // Output rate for fluvial removal of deposits from the Anderson-Darling test.
        .inherit(&ia) // Inherited ages of charcoal in debris-flow deposits.
        .model(&model); // Load model parameters.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        // .source(&debris_flows) // Set source as debris-flow deposits.
        .turnover(&318.0) // Set turnover period from the Anderson-Darling test.
        .manager(&model); // Load model parameters.

    // Fit model to observed deposit ages for specified duration.
    // Change directory path for user, panics on invalid path
    fluvial
        .fit_rates_timed(&ff, "/home/erik/output/fines_ad_20x_1000.csv")
        .unwrap();
}
