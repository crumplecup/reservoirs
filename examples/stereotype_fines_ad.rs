use reservoirs::prelude::*;
/// Produces a csv file with a vector of fluvial gravel ages fit to observed deposit ages using the Chi-squared test.
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
    let mut more = fg.to_vec();
    source_flux.append(&mut flux);
    source_flux.append(&mut more);


    // Set model parameters.
    let model = ModelManager::new()
        .index(0..20000) // Range of years to fit transit time probabilities.
        .obs(&ff) // Observations to fit.
        .obs_len(&ff) // Number of samples to collect from source.
        .period(40000.0) // Time period of individual simulations in years.
        .range(1000) // Seed for rng for reproducibility.
        .runs(100000); // Number of times to run the model per sampling point.

    // Reservoir for gravel deposits.
    let fluvial = Fluvial::new()
        .source(source_flux)
        // .source(&debris_flows) // Set source as debris-flow deposits.
        .turnover(&204.0) // Set turnover period from the Anderson-Darling test.
        .manager(&model); // Load model parameters.

    // Stereotype gravel deposit age.
    let mut stereo = fluvial
        .capture_rate_fines(0.1) // Set capture rate for fines in flux from the Anderson-Darling test.
        .capture_rate_gravels(0.295) // Set capture rate for gravels in flux from the Anderson-Darling test.
        .storage_rate_fines(0.1) // Set storage rate for fines in storage from the Anderson-Darling test.
        .storage_rate_gravels(0.065) // Set storage rate for gravels in storage from the Anderson-Darling test.
        .cherry_pick(); // Return run in model.runs most characteristic of the distribution using the K-S test.
                            // Change directory path for user, panics on invalid path
    utils::record(&mut stereo, "/home/erik/output/fines_stereo_ad_1000.csv").unwrap();
}
