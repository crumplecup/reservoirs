use rand::distributions::Uniform;
use rand::thread_rng;
use rand_distr::num_traits::abs;
use rand_distr::{Distribution, Exp};
use std::ops::Range;

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

pub fn mean(numbers: &[f64]) -> f64 {
    let sum: f64 = numbers.iter().sum();

    sum as f64 / numbers.len() as f64
}

pub fn median(numbers: &[f64]) -> f64 {
    let len = numbers.len();
    let mid = len / 2;
    if len % 2 == 0 {
        mean(&numbers[(mid - 1)..(mid + 1)].to_vec())
    } else {
        numbers[mid]
    }
}
