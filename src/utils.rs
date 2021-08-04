use crate::errors;
use rand::seq::IteratorRandom;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex;
use serde::Serialize;

/// Anderson-Darling Two-Sample Test
pub fn ad_dual(sample: &[f64], other: &[f64]) -> f64 {
    // join the two vectors and sort
    let mut x = sample.to_vec();
    let mut y = other.to_vec();
    x.append(&mut y);
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lnx = sample.len();
    let lny = other.len();
    let k = lnx + lny;
    let mut ad_i = Vec::new();
    for i in x {
        ad_i.push(sample.iter().filter(|z| **z <= i).count() as f64);
    }
    let mut adi = Vec::new();

    for (i, val) in ad_i.iter().take(k - 1).enumerate() {
        let step = f64::powi((k as f64 * val) - (lnx * (i + 1)) as f64, 2)
            / ((i + 1) * (k - (i + 1))) as f64;
        adi.push(step);
    }
    let mut ad = adi.iter().sum::<f64>();
    ad /= (lnx * lny) as f64;
    ad
}

/// Bootstrap a synthetic dataset from observed samples.
pub fn bootstrap(obs: &[f64]) -> Vec<f64> {
    let ln = obs.len();
    let mut boot: Vec<f64> = Vec::with_capacity(ln);
    for _ in 0..ln {
        boot.push(obs.iter().cloned().choose(&mut rand::thread_rng()).unwrap());
    }
    boot
}

/// Estimate the cumulative distribution function (CDF) of a vector of observations.
///  - `x` is a reference to a slice of f64 values
///  - returns a vector of tuples(val, cdf)
pub fn cdf(x: &[f64]) -> Vec<(f64, f64)> {
    let ln = x.len() as f64;
    let mut cx = x.to_vec();
    cx.sort_by(|a, b| a.partial_cmp(b).unwrap());
    cx.dedup();
    let mut cdf = Vec::new();
    for i in cx {
        let num = x.iter().filter(|x| **x <= i).count() as f64 / ln;
        cdf.push((i, num));
    }
    cdf
}

/// Estimate the cumulative distribution function (CDF) of a vector of observations `x`,
/// binned into `bins` number of obs to preserve memory on large vectors.
///  - `x` is a reference to a slice of f64 values.
///  - `bins` is the number of observations to subsample from `x`.
///  - Returns a vector of the cdf.
pub fn cdf_bin(obs: &[f64], bins: usize) -> Vec<f64> {
    let cdf: Vec<(f64, f64)> = cdf(obs);
    let mut cdf_bin = Vec::new();
    for i in 0..bins {
        let res: Vec<(f64, f64)> = cdf
            .iter()
            .cloned()
            .filter(|z| z.1 <= (i + 1) as f64 / bins as f64)
            .collect();
        let mut res_max = 0.0;
        if !res.is_empty() {
            res_max = res[res.len() - 1].0;
        }
        cdf_bin.push(res_max);
    }
    cdf_bin
}

/// Appends one vector to another and generates the CDF of observations.
pub fn cdf_dual(obs: &[f64], other: &[f64]) -> Vec<(f64, f64)> {
    // join the two vectors, sort and deduplicate
    let mut x = obs.to_vec();
    let mut y = other.to_vec();
    let lnx = x.len() as f64; // note the original lengths
    let lny = y.len() as f64;
    let xo = x.clone(); // clone the originals for later use
    let yo = y.clone();
    x.append(&mut y);
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x.dedup();
    // construct cdf of each
    let mut cdf = Vec::new();
    for i in x {
        let numx = xo.iter().filter(|z| **z <= i).count() as f64 / lnx;
        let numy = yo.iter().filter(|z| **z <= i).count() as f64 / lny;
        cdf.push((numx, numy));
    }
    cdf
}

/// Returns the cdf given a pmf, implementing a cumulative sum over f64 values.
///
/// # Examples
///
/// ```{rust}
/// use::reservoirs::prelude::*;
/// fn main () -> Result<(), ResError> {
///     let two = vec![1.0, 2.0];
///     let three = vec![3.0, 4.0];
///     let four = vec![1.0, 2.0, 3.0, 4.0];
///
///     // convert to cdf
///     let cdf_two = utils::cdf(&two).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_four = utils::cdf(&four).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_two_of_four = utils::cdf_rng(&two, &(0..4));
///     let cdf_three_of_four = utils::cdf_rng(&three, &(0..4));
///
///     // convert to pmf
///     let pmf_two = utils::pmf_from_cdf(&cdf_two);
///     let pmf_four = utils::pmf_from_cdf(&cdf_four);
///     let pmf_two_of_four = utils::pmf_from_cdf(&cdf_two_of_four);
///     let pmf_three_of_four = utils::pmf_from_cdf(&cdf_three_of_four);
///
///     // convert back to cdf
///     let two_cdf = utils::cdf_from_pmf(&pmf_two);
///     let four_cdf = utils::cdf_from_pmf(&pmf_four);
///     let two_of_four_cdf = utils::cdf_from_pmf(&pmf_two_of_four);
///     let three_of_four_cdf = utils::cdf_from_pmf(&pmf_three_of_four);
///
///     let thresh = 0.0001;
///     for i in 0..pmf_two.len() {
///         assert_eq!((two_cdf[i] - cdf_two[i]) < thresh, true);
///     }
///     for i in 0..pmf_four.len() {
///         assert_eq!((four_cdf[i] - cdf_four[i]) < thresh, true);
///         assert_eq!((two_of_four_cdf[i] - cdf_two_of_four[i]) < thresh, true);
///         assert_eq!((three_of_four_cdf[i] - cdf_three_of_four[i]) < thresh, true);
///     }
///     Ok(())
/// }
/// ```
pub fn cdf_from_pmf(pmf: &[f64]) -> Vec<f64> {
    let mut acc = 0.0;
    let res = pmf
        .iter()
        .map(|x| {
            acc += *x;
            acc
        })
        .collect::<Vec<f64>>();
    res
}

/// take the cdf of observations over a range of discrete values
///  - `obs` is a reference to a slice of f64 values.
///  - `range` is range of discrete f64 values.
///  - Returns a vector of the cdf along range.
///
/// # Examples
///
/// ```{rust}
/// use reservoirs::prelude::*;
/// fn main() -> Result<(), ResError> {
///     let two = vec![1.0, 2.0];
///     let whole = vec![0.5, 1.0];
///     let four = vec![1.0, 2.0, 3.0, 4.0];
///     let half = vec![0.5, 1.0, 1.0, 1.0];
///     let quarter = vec![0.25, 0.5, 0.75, 1.0];
///
///     let three = vec![3.0, 4.0];
///     let third = vec![0.0, 0.0, 0.5, 1.0];
///
///     let cdf_two = utils::cdf(&two).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_four = utils::cdf(&four).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_two_of_four = utils::cdf_rng(&two, &(1..4));
///     let cdf_three_of_four = utils::cdf_rng(&three, &(1..4));
///
///     let thresh = 0.0001;
///     for i in 0..cdf_two.len() {
///         assert_eq!((whole[i] - cdf_two[i]) < thresh, true);
///     }
///     for i in 0..cdf_four.len() {
///         assert_eq!((quarter[i] - cdf_four[i]) < thresh, true);
///         assert_eq!((half[i] - cdf_two_of_four[i]) < thresh, true);
///         assert_eq!((third[i] - cdf_three_of_four[i]) < thresh, true);
///     }
///     Ok(())
/// }
/// ```
pub fn cdf_rng(obs: &[f64], range: &std::ops::Range<i32>) -> Vec<f64> {
    let cdf: Vec<(f64, f64)> = cdf(obs);
    let mut cdf_rng = Vec::new();
    for r in range.start..=range.end {
        let res: Vec<(f64, f64)> = cdf
            .iter()
            .cloned()
            .filter(|(a, _)| *a <= r as f64)
            .collect();
        let mut res_max = 0.0;
        if !res.is_empty() {
            res_max = res[res.len() - 1].1;
        }
        cdf_rng.push(res_max);
    }
    cdf_rng
}

/// Chi-squared goodness-of-fit test.
pub fn chi_squared(obs: &[f64], other: &[f64]) -> f64 {
    // join the two vectors, sort and deduplicate
    let mut x = obs.to_vec();
    let mut y = other.to_vec();
    let xo = x.clone(); // clone the originals for later use
    let yo = y.clone();
    x.append(&mut y);
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x.dedup();
    // construct cdf of each
    let mut difs = Vec::new();
    for i in x {
        let numx = xo.iter().filter(|z| **z <= i).count() as f64;
        let numy = yo.iter().filter(|z| **z <= i).count() as f64;
        if numy > 0.0 {
            difs.push(f64::powi(numx - numy, 2) / numy);
        }
    }
    difs.iter().sum::<f64>()
}


/// Convolve frequency distributions along an index of values.
///
/// # Examples
///
/// ```{rust}
/// use reservoirs::prelude::*;
/// fn main () -> Result<(), ResError> {
///     let two = vec![1.0, 2.0];
///     let whole = vec![0.5, 0.5];
///     let four = vec![1.0, 2.0, 3.0, 4.0];
///     let half = vec![0.0, 0.5, 0.5, 0.0, 0.0];
///     let quarter = vec![0.0, 0.25, 0.25, 0.25, 0.25];
///
///     let three = vec![3.0, 4.0];
///     let third = vec![0.0, 0.0, 0.0, 0.5, 0.5];
///
///     let cdf_two = utils::cdf(&two).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_four = utils::cdf(&four).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_two_of_four = utils::cdf_rng(&two, &(0..4));
///     let cdf_three_of_four = utils::cdf_rng(&three, &(0..4));
///
///     let pmf_two = utils::pmf_from_cdf(&cdf_two);
///     let pmf_four = utils::pmf_from_cdf(&cdf_four);
///     let pmf_two_of_four = utils::pmf_from_cdf(&cdf_two_of_four);
///     let pmf_three_of_four = utils::pmf_from_cdf(&cdf_three_of_four);
///
///     let thresh = 0.0001;
///     for i in 0..pmf_two.len() {
///         assert_eq!((whole[i] - pmf_two[i]) < thresh, true);
///     }
///     for i in 0..pmf_four.len() {
///         assert_eq!((quarter[i] - pmf_four[i]) < thresh, true);
///         assert_eq!((half[i] - pmf_two_of_four[i]) < thresh, true);
///         assert_eq!((third[i] - pmf_three_of_four[i]) < thresh, true);
///     }
///     let convo_two_three = utils::convo(&two, &three, 4);
///     let out_two_three = vec![0.0, 0.25, 0.25, 0.25, 0.25];
///     println!("convo of two and three is {:?}", convo_two_three);
///     for i in 0..convo_two_three.len() {
///         assert_eq!((out_two_three[i] - convo_two_three[i]) < thresh, true);
///     }
///
/// Ok(())
/// }
/// ```
pub fn convo(x: &[f64], y: &[f64], rng: i32) -> Vec<f64> {
    let len = (rng + 1) as usize;
    let range = 0..rng;

    let cdf_x = cdf_rng(x, &range);
    let mut pmf_x = pmf_from_cdf(&cdf_x);
    let cdf_y = cdf_rng(y, &range);
    let mut pmf_y = pmf_from_cdf(&cdf_y);

    let mut real_planner = RealFftPlanner::<f64>::new();
    // create a FFT
    let r2c = real_planner.plan_fft_forward(len);
    // make input and output vectors
    // let mut in_data = r2c.make_input_vec();
    let mut spectrum_x = r2c.make_output_vec();
    let mut spectrum_y = r2c.make_output_vec();

    // Forward transform the input data
    r2c.process(&mut pmf_x, &mut spectrum_x).unwrap();
    r2c.process(&mut pmf_y, &mut spectrum_y).unwrap();
    // Add the frequency distributions together
    let mut spectrum = spectrum_x
        .iter()
        .zip(spectrum_y.iter())
        .map(|(a, b)| a + b)
        .collect::<Vec<Complex<f64>>>();

    // create an iFFT and an output vector
    let c2r = real_planner.plan_fft_inverse(len);
    let mut out_data = c2r.make_output_vec();

    c2r.process(&mut spectrum, &mut out_data).unwrap();
    out_data = out_data
        .iter()
        .map(|a| a / (rng + 1) as f64)
        .collect::<Vec<f64>>();
    out_data
}

/// Generates expected values from a Poisson distribution for a single event envelope.
pub fn fish(rate: f64, t: f64) -> f64 {
    let mut val = -rate * t;
    val = val.exp() * rate * t;
    val
}

/// Goodness of fit test suite.
pub fn gof(obs: &[f64], other: &[f64]) -> Vec<f64> {
    let cdf = cdf_dual(obs, other);
    let ad = ad_dual(obs, other);
    // anderson-darling test
    let lnx = obs.len();
    let lny = other.len();
    let k = lnx + lny;
    let ada = (2.492 - 1.0) * (1.0 - (1.55 / k as f64)) + 1.0;
    // chi-squared pearsons test
    let ch = chi_squared(obs, other);
    // let ch = cdf
    //     .iter()
    //     .filter(|x| x.0 > 0.0)
    //     .map(|x| f64::powi(x.1 - x.0, 2) / x.0)
    //     .sum::<f64>();
    // kuiper test
    let kp1 = cdf.iter().map(|x| x.0 - x.1).fold(0.0, f64::max);
    let kp2 = cdf.iter().map(|x| x.1 - x.0).fold(0.0, f64::max);
    let kp = kp1 + kp2;
    // kolmogorov-smirnov test
    let ks = cdf
        .iter()
        .map(|x| rand_distr::num_traits::abs(x.0 - x.1))
        .fold(0.0, f64::max);
    let ksa = 1.35810 * f64::sqrt(k as f64 / (lnx * lny) as f64); // k-s crit value at 0.05
    vec![ad, ada, ch, kp, ks, ksa]
}

/// Produce integer-like index of f64 values from generic range.
///
/// # Examples
///
/// ```{rust}
/// use::reservoirs::prelude::*;
/// fn main() -> Result<(), ResError> {
///     let five = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
///     let index = utils::index_f64(0.0, 5.0, 1.0);
///
///     let thresh = 0.00001;
///     for i in 0..=5 {
///         assert_eq!((five[i] - index[i]) < thresh, true)
///     }
///     Ok(())
/// }
/// ```
pub fn index_f64(start: f64, end: f64, step: f64) -> Vec<f64> {
    let mut k = start;
    let mut index = Vec::new();
    while k <= end {
        index.push(k);
        k += step;
    }
    index
}

/// Calculates the low point along `y` and returns the value of `x` at the low point.
pub fn low_point(x: &[f64], y: &[f64]) -> f64 {
    let mut y_slope: Vec<f64> = Vec::new();
    for i in 0..(y.len() - 1) {
        y_slope.push(y[i + 1] - y[i]);
    }
    let mut y_momentum: Vec<f64> = vec![y_slope.iter().sum::<f64>().abs()];
    for i in 0..(y_slope.len() - 1) {
        let before: f64 = y_slope[0..(i + 1)].iter().sum::<f64>().abs();
        let after: f64 = y_slope[(i + 1)..y_slope.len()].iter().sum::<f64>().abs();
        y_momentum.push(before + after);
    }
    y_momentum.push(y_slope.iter().sum::<f64>().abs());
    let mut low: Vec<f64> = Vec::new();
    let y_max: f64 = y_momentum.iter().cloned().fold(f64::NAN, |m, v| v.max(m));
    for (i, val) in y_momentum.iter().enumerate() {
        if (y_max - *val) < 0.00001 {
            low.push(x[i])
        }
    }
    low[low.len() - 1]
}

/// Calculate the mean of a slice of f64 values.
///  - `numbers` is a reference to a slice of f64 values.
///  - Returns the mean of `numbers`.
///
/// # Examples
///
/// ```rust
/// let numbers = vec![1.0, 1.5, 2.0, 2.5, 3.0];
/// let mn = reservoirs::utils::mean(&numbers);
/// assert_eq!(2.0, mn);
/// ```
pub fn mean(numbers: &[f64]) -> f64 {
    let sum: f64 = numbers.iter().sum();

    sum as f64 / numbers.len() as f64
}

/// Calculate the median of a slice of f64 values.
///  - `numbers` is a reference to a slice of f64 values.
///  - Returns the median of `numbers`.
///
/// # Examples
///
/// ```rust
/// let numbers = vec![1.0, 3.0, 7.0, 10.0];
/// let med = reservoirs::utils::median(&numbers);
/// assert_eq!(5.0, med);
/// ```
pub fn median(numbers: &[f64]) -> f64 {
    let len = numbers.len();
    let mut x: Vec<f64> = numbers.to_vec();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = len / 2;
    if len % 2 == 0 {
        mean(&x[(mid - 1)..(mid + 1)].to_vec())
    } else {
        x[mid]
    }
}

/// pmf from observations
///
/// calculates the pmf from the cdf
///
/// # Examples
///
/// ```{rust}
/// use reservoirs::prelude::*;
/// fn main () -> Result<(), ResError> {
///     let two = vec![1.0, 2.0];
///     let whole = vec![0.5, 0.5];
///     let four = vec![1.0, 2.0, 3.0, 4.0];
///     let half = vec![0.5, 0.5, 0.0, 0.0];
///     let quarter = vec![0.25, 0.25, 0.25, 0.25];
///
///     let three = vec![3.0, 4.0];
///     let third = vec![0.0, 0.0, 0.5, 0.5];
///
///     let cdf_two = utils::cdf(&two).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_four = utils::cdf(&four).iter().map(|(a, b)| *b).collect::<Vec<f64>>();
///     let cdf_two_of_four = utils::cdf_rng(&two, &(1..4));
///     let cdf_three_of_four = utils::cdf_rng(&three, &(1..4));
///
///     let pmf_two = utils::pmf_from_cdf(&cdf_two);
///     let pmf_four = utils::pmf_from_cdf(&cdf_four);
///     let pmf_two_of_four = utils::pmf_from_cdf(&cdf_two_of_four);
///     let pmf_three_of_four = utils::pmf_from_cdf(&cdf_three_of_four);
///
///     let thresh = 0.0001;
///     for i in 0..pmf_two.len() {
///         assert_eq!((whole[i] - pmf_two[i]) < thresh, true);
///     }
///     for i in 0..pmf_four.len() {
///         assert_eq!((quarter[i] - pmf_four[i]) < thresh, true);
///         assert_eq!((half[i] - pmf_two_of_four[i]) < thresh, true);
///         assert_eq!((third[i] - pmf_three_of_four[i]) < thresh, true);
///     }
/// Ok(())
/// }
/// ```
pub fn pmf_from_cdf(cdf: &[f64]) -> Vec<f64> {
    let mut f = Vec::new();
    let mut k = 0;
    while k < cdf.len() {
        if k == 0 {
            f.push(cdf[k]);
        }
        if k > 0 {
            f.push(cdf[k] - cdf[k - 1]);
        }
        k += 1;
    }
    f
}

/// Calculate the value of the CDF of `obs` at a given threshold `thresh`.
/// Called by [quantiles](#method.quantiles).
pub fn quantile(obs: &[f64], thresh: &f64) -> f64 {
    let cdf: Vec<(f64, f64)> = cdf(obs);
    let sub: Vec<(f64, f64)> = cdf.iter().cloned().filter(|x| x.1 <= *thresh).collect();
    let mut res = 0.0;
    if !sub.is_empty() {
        res = sub[sub.len() - 1].0;
    }
    res
}

/// Calculate the quantiles for a box-and-whisker plot (2.5%, 25%, 50%, 75%, 97.5%).
///  - `obs` is a reference to a slice of f64 observed values.
///  - Returns the observed value at each quantile of the CDF.
pub fn quantiles(obs: &[f64]) -> Vec<f64> {
    let quants = vec![
        quantile(obs, &0.025),
        quantile(obs, &0.25),
        quantile(obs, &0.5),
        quantile(obs, &0.75),
        quantile(obs, &0.975),
    ];
    quants
}

/// Read data from csv file.
pub fn read_f64(path: &str) -> Result<Vec<f64>, errors::ResError> {
    let mut dat = Vec::new();
    let var = std::fs::File::open(path)?;
    let mut rdr = csv::Reader::from_reader(var);
    for result in rdr.records() {
        let row = result?;
        let row = row.deserialize(None)?;
        dat.push(row);
    }
    Ok(dat)
}

/// Write statistical results to csv file.
pub fn record<T: Serialize>(rec: &mut Vec<T>, path: &str) -> Result<(), errors::ResError> {
    let mut wtr = csv::Writer::from_path(path)?;
    for i in rec {
        wtr.serialize(i)?;
    }
    wtr.flush()?;
    Ok(())
}
