use reservoirs::prelude::*;


/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
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

    let mut source_flux = df.to_vec();
    let mut flux = fg.to_vec();
    source_flux.append(&mut flux);

    // Set model parameters.
    let model = ModelManager::new()
        .duration(10000)
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&ff) // Observations to fit.
        .obs_len(&source_flux) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .thresholds(5.5, 14800.0, 0.33, 0.29)
        .runs(200); // Number of times to run the model per sampling point.



    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source(df)
        // .unwrap() // Set source as debris-flow deposits.
        // .capture_rate_gravels(0.1867216)
        .capture_rate_gravels(0.3285686)
        // .storage_rate_gravels(0.1223394)
        .storage_rate_gravels(0.957134)
        .turnover(&309.6175) // Set turnover period from the Kolmogorov-Smirnov test.
        .manager(&model.clone()); // Load model parameters.

    fluvial.hit_rates_timed("/home/erik/output/fines_hits_200_1000.csv").unwrap();
}
