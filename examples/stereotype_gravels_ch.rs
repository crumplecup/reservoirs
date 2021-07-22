use reservoirs::prelude::*;
/// Produces a csv file with a vector of fluvial gravel ages fit to observed deposit ages using the Chi-squared test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("/home/erik/data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep.iter().filter(|x| x.facies == "FG").map(|x| x.age).collect(); // Mean gravel deposit ages.
    let iat = Sample::read("/home/erik/data/iat.csv").unwrap(); // Inherited ages of charcoal in deposits.
    let ia: Vec<f64> = iat.iter().map(|x| x.age).collect(); // Vector of inherited ages from all classes of deposits.

    // Set model parameters.
    let model = ModelManager::new()
        .obs(&fg)
        .period(40000.0) // Time period of individual simulations in years.
        .range(777) // Seed for rng for reproducibility.
        .runs(20) // Number of times to run the model per sampling point.
        .source_runs(1); // Number of times to run the model per gravel source.

    // Source deposits for gravels.
    let debris_flows = Reservoir::new()
        .input(&0.46).unwrap() // Input rate for debris-flow deposits from the Chi-squared test.
        .output(&0.46).unwrap() // Output rate for fluvial removal of deposits from the Chi-squared test.
        .inherit(&ia) // Inherited ages of charcoal in debris-flow deposits.
        .model(&model); // Load model parameters.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source(&[debris_flows.clone()]) // Set source as debris-flow deposits.
        .turnover(&293.0) // Set turnover period from the Chi-squared test.
        .manager(&model); // Load model parameters.

    // Stereotype gravel deposit age.
    let mut stereo = fluvial
        .capture_rate_gravels(0.177) // Set capture rate for gravels in flux from the Chi-squared test.
        .storage_rate_gravels(0.114) // Set storage rate for gravels in storage from the Chi-squared test.
        .stereotype_rate(); // Return run in model.runs most characteristic of the distribution using the K-S test.
    // Change directory path for user, panics on invalid path
    utils::record(&mut stereo, "/home/erik/output/gravels_stereo_ch.csv").unwrap();
}