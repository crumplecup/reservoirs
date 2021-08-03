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
    let df = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age + 50.0)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&fg) // Observations to fit.
        .obs_len(&df) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(777) // Seed for rng for reproducibility.
        .runs(1); // Number of times to run the model per sampling point.


    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source_from_csv("data/debris_flow_transits_ch.csv")
        .unwrap() // Set source as debris-flow deposits.
        .capture_rate_gravels(0.3421237)
        .storage_rate_gravels(0.02205183)
        .turnover(&306.3833)
        .manager(&model.clone()); // Load model parameters.

    let mut max = rand_distr::num_traits::Float::max_value();
    let mut best = Vec::new();
    for i in 0..10000 {
        let fit = fluvial.clone().manager(&model.clone().range(i)).transit_times();
        let gof = utils::gof(&fit, &fg);
        if gof[4] < max {
            best = fit;
            max = gof[4];
        }
    }
    utils::record(&mut best, "/home/erik/output/stereotype_gravels_ch.csv").unwrap();
}
