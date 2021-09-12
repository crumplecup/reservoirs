use reservoirs::prelude::*;

/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() -> Result<(), ResError> {
    // Load charcoal age data.
    // Change directory path for user, panics on invalid path
    let dep = Sample::read("data/dep.csv")?; // Mean charcoal ages of deposits.
    let iat = Sample::read("data/iat.csv")?;
    let ia = iat.iter().map(|x| x.age).collect::<Vec<f64>>();

    let df = dep
        .iter()
        .filter(|x| x.facies == "DF")
        .map(|x| x.age + 50.0)
        .collect::<Vec<f64>>(); // Mean gravel deposit ages.

    // Set model parameters.
    let model = ModelManager::new()
        .batch(1) // Number of samples between saves.
        .duration(10000)
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&df) // Observations to fit.
        .obs_len(&df) // Number of samples to collect from source.
        .period(1000.0) // Time period of individual simulations in years.
        .periods(1000.0..50000.0) // Range of periods to fit.
        .range(1000) // Seed for rng for reproducibility.
        .rates(0.01..2.0) // Range of rates to fit.
        .runs(3); // Number of times to run the model per sampling point.

    let reservoir = Reservoir::new()
        .input(&0.35)
        .unwrap()
        .output(&0.35)
        .unwrap()
        .model(&model)
        .inherit(&ia)
        .sim();

    let fit = reservoir.fit_reservoir();

    println!("{:?}", fit);
    Ok(())
}
