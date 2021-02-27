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
///  - Returns a vector of tuples(val, cdf).
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

/// Calculates the low point along `y` and returns the value of `x` at the low point.
pub fn low_point(x: Vec<f64>, y: Vec<f64>) -> f64 {
    let mut x_slope: Vec<f64> = Vec::new();
    println!("x is len {}", x.len());
    for i in 0..(x.len()-1) {
        x_slope.push(x[i+1] - x[i]);
    }
    println!("x_slope is len {}", x_slope.len());
    let mut x_momentum: Vec<f64> = vec![x_slope.iter().sum::<f64>().abs()];
    println!("x_momentum is {:?}", x_momentum);
    for i in 0..(x_slope.len()-1) {
        println!("i is {}", i);
        let before: f64 = x_slope[0..(i+1)].iter().sum::<f64>().abs();
        let after: f64 = x_slope[(i+1)..x_slope.len()].iter().sum::<f64>().abs();
        x_momentum.push(before + after);
    }
    x_momentum.push(x_slope.iter().sum::<f64>().abs());
    println!("x_momemtum is len {}", x_momentum.len());
    let mut low: Vec<f64> = Vec::new();
    let x_max: f64 = x_momentum.iter().fold(0.0, |acc, z| f64::max(acc, *z));
    for (i, val) in x_momentum.iter().enumerate() {
        if (*val - x_max) < 0.0001 {
            low.push(y[i])
        }
    }
    y[0]
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
    let mid = len / 2;
    if len % 2 == 0 {
        mean(&numbers[(mid - 1)..(mid + 1)].to_vec())
    } else {
        numbers[mid]
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
