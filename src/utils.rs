
use std::ops::Range;
use rand_distr::{Exp, Distribution};
use rand::thread_rng;
use rand::distributions::Uniform;
use rand_distr::num_traits::abs;


pub fn cdf(x: Vec<f64>) -> Vec<(f64, f64)> {
    let ln = x.len() as f64;
    let mut cx = x.clone();
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

pub fn cdf_bin(obs: &Vec<f64>, bins: usize) -> Vec<f64> {
    let cdf: Vec<(f64, f64)>= cdf(obs.clone());
    let mut cdf_bin = Vec::new();
    for i in 0..bins {
        let res: Vec<(f64, f64)> = cdf.iter().cloned().filter(|z| z.1 <= (i+1) as f64 / bins as f64).collect();
        let mut res_max = 0.0;
        if res.len() > 0 {
            res_max = res[res.len() - 1].0;
        }
        cdf_bin.push(res_max);
    }
    cdf_bin
}


pub fn gof(a: &Vec<f64>, b: &Vec<f64>) -> (f64, f64) {
    let mut x = a.clone();
    let mut y = b.clone();
    let lnx = x.len() as f64;
    let lny = y.len() as f64;
    let xo = x.clone();
    let yo = y.clone();
    x.append(&mut y);
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x.dedup();
    let mut cdf = Vec::new();
    for i in x {
        let numx: Vec<&f64> = xo.iter().filter(|z| **z <= i).collect();
        let numy: Vec<&f64> = yo.iter().filter(|z| **z <= i).collect();
        let resx = numx.len() as f64 / lnx;
        let resy = numy.len() as f64 / lny;
        cdf.push((resx, resy));
    }
    let ks = cdf.iter().map(|x| abs(x.0 - x.1)).fold(0.0, f64::max);
    let kp1 = cdf.iter().map(|x| x.0 - x.1).fold(0.0, f64::max);
    let kp2 = cdf.iter().map(|x| x.1 - x.0).fold(0.0, f64::max);
    let kp = kp1 + kp2;
    (ks, kp)
}

pub fn mean(numbers: &Vec<f64>) -> f64 {
    let sum: f64 = numbers.iter().sum();

    sum as f64 / numbers.len() as f64
}

pub fn median(numbers: &Vec<f64>) -> f64 {
    let len = numbers.len();
    let mid = len / 2;
    if len % 2 == 0 {
        mean(&numbers[(mid - 1)..(mid + 1)].to_vec())
    } else {
        f64::from(numbers[mid])
    }
}

