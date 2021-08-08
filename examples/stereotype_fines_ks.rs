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
        .obs_len(&ff) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1012) // Seed for rng for reproducibility.
        .runs(70000); // Number of times to run the model per sampling point.

    let mut source_flux = df.to_vec();
    let mut flux = fg.to_vec();
    source_flux.append(&mut flux);
    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        // .source(source_flux)
        .source_from_csv("data/gravels_cdf_ks.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.025)
        .storage_rate_gravels(0.425)
        .turnover(&208.00)
        .manager(&model); // Load model parameters.

    let mut rec = fluvial.cherry_pick();
    utils::record(&mut rec, "/home/erik/output/stereotype_fines_ks_1000.csv").unwrap();
}
