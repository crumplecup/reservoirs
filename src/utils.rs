use crate::errors;
use rand::seq::IteratorRandom;
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
        let step: Vec<f64> = sample.iter().filter(|z| **z <= i).map(|z| *z).collect();
        ad_i.push(step.len() as f64);
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
        let num: Vec<&f64> = x.iter().filter(|x| **x <= i).collect();
        let res = num.len() as f64 / ln;
        cdf.push((i, res));
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
        let numx: Vec<&f64> = xo.iter().filter(|z| **z <= i).collect();
        let numy: Vec<&f64> = yo.iter().filter(|z| **z <= i).collect();
        let resx = numx.len() as f64 / lnx;
        let resy = numy.len() as f64 / lny;
        cdf.push((resx, resy));
    }
    cdf
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
