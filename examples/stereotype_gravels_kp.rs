use reservoirs::prelude::*;

/// Produces a csv file of model fit to observed deposit ages.
fn main() {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv").unwrap(); // Mean charcoal ages of deposits.
    let fg: Vec<f64> = dep
        .iter()
        .filter(|x| x.facies == "FG")
        .map(|x| x.age + 50.0)
        .collect(); // Mean gravel deposit ages.
    let _df = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age + 50.0)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&fg) // Observations to fit.
        .obs_len(&fg) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1017) // Seed for rng for reproducibility.
        .runs(100000); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        // .source(df)
        .source_from_csv("data/debris_flow_transits_kp.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.08333428)
        .storage_rate_gravels(0.01669015)
        .turnover(&191.0)
        .manager(&model.clone()); // Load model parameters.

    let mut rec = fluvial.cherry_pick();
    utils::record(&mut rec, "/home/erik/output/stereotype_gravels_kp_1000.csv").unwrap();
}
