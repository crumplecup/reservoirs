use reservoirs::prelude::*;

/// Produces a csv file of model fit to observed deposit ages using the Kolmogorov-Smirnov test.
fn main() {
    let gravels = utils::read_f64("data/gravels_cdf_ad.csv").unwrap();
    let debris_flows = utils::read_f64("data/debris_flow_transits_ad.csv").unwrap();
    let mut source = utils::convo(&gravels, &debris_flows, 20000);
    utils::record(&mut source, "home/erik/output/fines_source_ad.csv").unwrap();
}
