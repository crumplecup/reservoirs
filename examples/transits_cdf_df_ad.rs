use reservoirs::prelude::*;
/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let iat = Sample::read("data/iat.csv").unwrap(); // Mean inherited ages of charcoal in deposits.
    let _ia: Vec<f64> = iat.iter().map(|x| x.age).collect(); // Vector of inherited ages from all classes of deposits.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000)
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(10000); // Number of times to run the model per sampling point.

    // Source deposits for gravels.
    let debris_flows = Reservoir::new()
        .input(&0.35)
        .unwrap() // Input rate for debris-flow deposits from the Kolmogorov-Smirnov test.
        .output(&0.35)
        .unwrap() // Output rate for fluvial removal of deposits from the Kolmogorov-Smirnov test.
        //.inherit(&ia) // Inherited ages of charcoal in debris-flow deposits.
        .model(&model); // Load model parameters.

    let mut rec = debris_flows.transit_times();
    utils::record(&mut rec, "/home/erik/output/transits_cdf_df_ad.csv").unwrap();
}
