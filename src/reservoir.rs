//! Structs and methods for Bolin & Rodhe reservoir models.
use crate::errors;
use crate::plot;
use crate::utils;
use log::*;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Exp};
use rayon::prelude::*;
use realfft::RealFftPlanner;
use rustfft::num_complex::Complex;
use serde::{Deserialize, Serialize};

/// Holder struct for passing test values between functions.
/// Makes functions calls and returns less confusing.
#[derive(Clone, Debug)]
pub struct Fit {
    ad: f64,
    ada: f64,
    chs: f64,
    kp: f64,
    ks: f64,
    ksa: f64,
}

/// Summary statistics for several Fit objects.
#[derive(Clone, Debug)]
pub struct Fits {
    fit: Fit,
    n: f64,
}

/// Holder struct for goodness-of-fit statistics.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Gof {
    input: f64,
    output: f64,
    ad: f64,
    ada: f64,
    chs: f64,
    kp: f64,
    ks: f64,
    ksa: f64,
    n: f64,
}

impl Gof {
    /// Convert csv record to Gof struct.
    pub fn read(path: &str) -> Result<Vec<Gof>, errors::ResError> {
        let mut gof = Vec::new();
        let var = std::fs::File::open(path)?;
        let mut rdr = csv::Reader::from_reader(var);
        for result in rdr.records() {
            let row = result?;
            let row: Gof = row.deserialize(None)?;
            gof.push(row);
        }
        Ok(gof)
    }

    /// Write statistical results to csv file.
    pub fn record(rec: &mut Vec<Gof>, title: &str) -> Result<(), errors::ResError> {
        let mut wtr = csv::Writer::from_path(title)?;
        for i in rec {
            wtr.serialize(i)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Get the values associated with fields of the struct.
    pub fn values(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        (
            self.input,
            self.output,
            self.ad,
            self.ada,
            self.chs,
            self.kp,
            self.ks,
            self.ksa,
            self.n,
        )
    }
}

/// Holder struct for goodness-of-fit statistics.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ReservoirFit {
    rate: f64,
    period: f64,
    ad: f64,
    ch: f64,
    kp: f64,
    ks: f64,
    n: f64,
}

impl ReservoirFit {
    /// Convert csv record to Gof struct.
    pub fn read(path: &str) -> Result<Vec<ReservoirFit>, errors::ResError> {
        let mut fit = Vec::new();
        let var = std::fs::File::open(path)?;
        let mut rdr = csv::Reader::from_reader(var);
        for result in rdr.records() {
            let row = result?;
            let row: ReservoirFit = row.deserialize(None)?;
            fit.push(row);
        }
        Ok(fit)
    }

    /// Write statistical results to csv file.
    pub fn record(rec: &mut Vec<ReservoirFit>, title: &str) -> Result<(), errors::ResError> {
        let mut wtr = csv::Writer::from_path(title)?;
        for i in rec {
            wtr.serialize(i)?;
        }
        wtr.flush()?;
        Ok(())
    }

    /// Get the values associated with fields of the struct.
    pub fn values(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        (
            self.rate,
            self.period,
            self.ad,
            self.ch,
            self.kp,
            self.ks,
            self.n,
        )
    }
}

/// Holds parameters for bootstrapping confidence intervals of reservoir models.
#[derive(Clone, Debug)]
pub struct Bootstrap {
    model: Model,
    bins: usize,
    samples: usize,
}

impl Bootstrap {
    /// Sets the model bin number for comparing ranges in [search](#method.search).
    pub fn bins(mut self, size: usize) -> Self {
        self.bins = size;
        self
    }

    /// Creates a Bootstrap struct from a given Model.
    pub fn new(model: Model) -> Self {
        Bootstrap {
            model,
            bins: 100,
            samples: 3,
        }
    }

    /// Sets the number of samples desired for calculating the mean in a range,
    /// called by [search](#method.search).
    pub fn samples(mut self, count: usize) -> Self {
        self.samples = count;
        self
    }

    /// Fits range of rates to an observed distribution using a steady state reservoir.
    pub fn search(
        &mut self,
        rate: std::ops::Range<f64>,
        obs: &[f64],
        path: &str,
    ) -> Result<(), errors::ResError> {
        // initialize time variables
        let dur = std::time::Duration::new(60 * 60 * self.model.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        // load existing record if any
        let exists = std::path::Path::new(path).exists();
        match exists {
            true => {
                let mut fits: Vec<f64> = utils::read_f64(path)?;
                rec.append(&mut fits);
            }
            false => {}
        }
        // perform bootstrap
        let mut count: usize = 0;
        let rec_path = format!("{}{}", path, "rec.csv");
        while std::time::SystemTime::now() < now + dur {
            let mut gof = Vec::new();
            let boot = utils::bootstrap(obs);
            while gof.len() < (self.bins * self.samples) {
                let mut new = self.model.steady(rate.clone(), &boot);
                gof.append(&mut new);
                println!(
                    "{}% complete.",
                    (gof.len() as f64 / (self.bins * self.samples) as f64 * 100.0).round()
                );
            }
            let gof_path = format!("{}{}{}{}", path, "gof_", count, ".csv");
            let png_path = format!("{}{}{}{}", path, "gof_", count, ".png");
            utils::record(&mut gof, &gof_path)?;
            count += 1;
            let (x, y) = Record::bin_ave(&gof, self.bins);
            plot::xy(&x, &y, &png_path)?;
            rec.push(utils::low_point(&x, &y));
            utils::record(&mut rec, &rec_path)?;
        }

        Ok(())
    }
}

/// Holds model characteristics associated with a reservoir.
#[derive(Clone, Debug)]
pub struct Model {
    reservoir: Reservoir,
    period: f64,
    runs: usize,
    batch: usize,
    duration: u64,
    range: i32,
}

impl Model {
    /// Sets the model batch size to desired number of samples.
    pub fn batch(mut self, samples: usize) -> Self {
        self.batch = samples;
        self
    }

    /// Sets the model period to desired time in years.
    pub fn duration(mut self, hours: &u64) -> Self {
        self.duration = hours.to_owned();
        self
    }

    /// Fits a range of input and output rates to an observed accumulation record.
    /// A wrapper that calls [fit_rng](struct.Reservoir.html#method.fit_rng) runs for a specified period of time `dur`,
    /// and writes the results to a csv file.
    ///  - `self` - set inherited age before running.
    ///  - `input` is the range of possible input rates.
    ///  - `output` is the range of possible output rates.
    ///  - `obs` is the vector of observed values in years.
    ///  - `title` labels the csv file, and must be a valid path name ending in ".csv"
    ///  - Returns a csv written to the path specified by `title` containing goodness-of-fit statistics.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use reservoirs::prelude::*;
    ///
    /// fn main() -> Result<(), ResError> {
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter()
    ///     .filter(|x| x.facies == "DF")
    ///     .map(|x| x.age)
    ///     .collect();
    ///
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter()
    ///     .map(|x| x.age)
    ///     .collect();
    ///
    /// let mut debris_flows = Reservoir::new()
    ///     .input(&0.687)?
    ///     .output(&0.687)?
    ///     .inherit(&ia);
    ///
    /// // model parameters
    /// let batch = 10;  // fit 10 input/output pairs at a time using rayon
    /// let duration = 1; // fit model for 1 hour
    /// let period = 30000.0; // run simulations for 30000 years
    /// let runs = 1000; // run 1000 simulated accumulations per candidate pair for goodness-of-fit
    ///
    /// // create reservoir model using builder pattern
    /// let mut model = Model::new(debris_flows)
    ///     .batch(batch)
    ///     .duration(&duration)
    ///     .period(&period)
    ///     .runs(runs);
    /// // for one hour, fit batches of 10 randomly selected rate pairs (from range 0.01 to 1.0)
    /// // to observed debris flows
    /// // by running 1000 simulations for 30000 years for each pair
    /// model.fit_range(0.01..1.0, 0.01..1.0, &df, "examples/df_fit_1k.csv");
    /// Ok(())
    /// }
    /// ```
    pub fn fit_range(
        &mut self,
        input: std::ops::Range<f64>,
        output: std::ops::Range<f64>,
        obs: &[f64],
        title: &str,
    ) {
        let dur = std::time::Duration::new(60 * 60 * self.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(title).exists();
        match exists {
            true => {
                let mut gof: Vec<Gof> = Gof::read(title).unwrap();
                rec.append(&mut gof);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let mut new = self.fit_rng(input.clone(), output.clone(), obs);
            {
                rec.append(&mut new);
            }
            Gof::record(&mut rec, title).unwrap();
        }
    }

    /// Randomly selects rate pairs from ranges `input` and `output`, and simulates accumulation records
    /// in batches using [fit_rate](#method.fit_rate).  Returns the selected input/output pair and the mean
    /// goodness-of-fit statistics for each pair.  Called by [fit_range](method.fit_range).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let ages = vec![10.0, 42.0, 77.7, 12.0, 99.9, 10000.0, 777.7];
    ///
    ///     let mut debris_flows = Reservoir::new()
    ///         .input(&0.862)?
    ///         .output(&0.862)?;
    ///
    ///     // model parameters
    ///     let batch = 10;  // fit 10 input/output pairs at a time using rayon
    ///     let period = 3000.0; // run simulations for 3000 years
    ///     let runs = 20; // run 20 simulated accumulations per candidate pair for goodness-of-fit
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(debris_flows)
    ///         .batch(batch)
    ///         .period(&period)
    ///         .runs(runs);
    ///     // fit batches of 10 randomly selected rate pairs (from range 0.01 to 1.0)
    ///     // to observed debris flows
    ///     // by running 20 simulations for 3000 years for each pair
    ///     let gofs = model.fit_rng(0.01..1.0, 0.01..1.0, &ages);
    ///     Ok(())
    /// }
    /// ```
    pub fn fit_rng(
        &mut self,
        input: std::ops::Range<f64>,
        output: std::ops::Range<f64>,
        obs: &[f64],
    ) -> Vec<Gof> {
        let mut inputs = Vec::with_capacity(self.batch);
        let mut outputs = Vec::with_capacity(self.batch);
        let mut fits = Vec::with_capacity(self.batch);
        for i in 0..self.batch {
            inputs.push(
                rand::distributions::Uniform::from(input.clone()).sample(&mut self.reservoir.range),
            );
            outputs.push(
                rand::distributions::Uniform::from(
                    output.clone().start.max(inputs[i] * 0.975)
                        ..output.clone().end.min(inputs[i] * 1.0125),
                )
                .sample(&mut self.reservoir.range),
            );
            fits.push(
                self.clone().reservoir(
                    self.reservoir
                        .clone()
                        .input(&inputs[i])
                        .unwrap()
                        .output(&outputs[i])
                        .unwrap(),
                ),
            );
        }
        let gof: Vec<Fits> = fits.par_iter().map(|x| x.clone().fit_rate(obs)).collect();
        let mut gofs = Vec::with_capacity(self.batch);
        for i in 0..self.batch {
            gofs.push(Gof {
                input: inputs[i],
                output: outputs[i],
                ad: gof[i].fit.ad,
                ada: gof[i].fit.ada,
                chs: gof[i].fit.chs,
                kp: gof[i].fit.kp,
                ks: gof[i].fit.ks,
                ksa: gof[i].fit.ksa,
                n: gof[i].n,
            })
        }
        gofs
    }

    /// Runs simulations on a reservoir,
    /// returns the mean goodness-of-fit statistics compared to accumulation record `other`.
    /// Called by [fit_rng](#method.fit_rng) and [steady](#method.steady).  To use,
    /// set characteristics of the reservoir and model before running.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let ages = vec![10.0, 42.0, 77.7, 12.0, 99.9, 10000.0, 777.7];
    ///
    ///     let mut debris_flows = Reservoir::new()
    ///         .input(&0.862)?
    ///         .output(&0.862)?;
    ///
    ///     // model parameters
    ///     let period = 3000.0; // run simulations for 3000 years
    ///     let runs = 20; // run 20 simulated accumulations per candidate pair for goodness-of-fit
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(debris_flows)
    ///         .period(&period)
    ///         .runs(runs);
    ///     // fit selected rate pairs
    ///     // to observed debris flows
    ///     // by running 20 simulations for 3000 years for each pair
    ///     let fits = model.fit_rate(&ages);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn fit_rate(&mut self, other: &[f64]) -> Fits {
        let mut res: Vec<Reservoir> = Vec::with_capacity(self.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds: Vec<u64> = seeder
            .sample_iter(&mut self.reservoir.range)
            .take(self.runs)
            .collect();
        for seed in seeds {
            res.push(self.reservoir.clone().range(seed));
        }
        res = res.par_iter().cloned().map(|x| x.sim()).collect();
        let fits: Vec<Fit> = res.par_iter().cloned().map(|x| x.gof1(other)).collect();
        let ads: Vec<f64> = fits.iter().map(|x| x.ad).collect();
        let ad = utils::median(&ads);
        let adas: Vec<f64> = fits.iter().map(|x| x.ada).collect();
        let ada = utils::median(&adas);
        let chss: Vec<f64> = fits.iter().map(|x| x.chs).collect();
        let chs = utils::median(&chss);
        let kps: Vec<f64> = fits.iter().map(|x| x.kp).collect();
        let kp = utils::median(&kps);
        let kss: Vec<f64> = fits.iter().map(|x| x.ks).collect();
        let ks = utils::median(&kss);
        let ksas: Vec<f64> = fits.iter().map(|x| x.ksa).collect();
        let ksa = utils::median(&ksas);
        let ns: Vec<f64> = res.iter().map(|x| x.mass.len() as f64).collect();
        let n = utils::median(&ns);

        Fits {
            fit: Fit {
                ad,
                ada,
                chs,
                kp,
                ks,
                ksa,
            },
            n,
        }
    }

    /// Fits a range of steady state reservoirs to an observed accumulation record.
    /// A wrapper that calls [steady](struct.Reservoir.html#method.steady) runs for a specified period of time `dur`,
    /// and writes the results to a csv file.
    ///  - `self` - set inherited age before running.
    ///  - `rate` is the range of possible rates.
    ///  - `obs` is the vector of observed values in years.
    ///  - `title` labels the csv file, and must be a valid path name ending in ".csv"
    ///  - Returns a csv written to the path specified by `title` containing goodness-of-fit statistics.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter()
    ///     .filter(|x| x.facies == "DF")
    ///     .map(|x| x.age)
    ///     .collect();
    ///
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter()
    ///     .map(|x| x.age)
    ///     .collect();
    ///
    /// let mut debris_flows = Reservoir::new()
    ///     .input(&0.687)?
    ///     .output(&0.687)?
    ///     .inherit(&ia);
    ///
    /// // model parameters
    /// let batch = 10;  // fit 10 input/output pairs at a time using rayon
    /// let duration = 1; // fit model for 1 hour
    /// let period = 30000.0; // run simulations for 30000 years
    /// let runs = 1000; // run 1000 simulated accumulations per candidate pair for goodness-of-fit
    ///
    /// // create reservoir model using builder pattern
    /// let mut model = Model::new(debris_flows)
    ///     .batch(batch)
    ///     .duration(&duration)
    ///     .period(&period)
    ///     .runs(runs);
    /// // for one hour, fit batches of 10 randomly selected rate pairs (from range 0.01 to 1.0)
    /// // to observed debris flows
    /// // by running 1000 simulations for 30000 years for each pair
    /// model.fit_steady(0.01..1.0, &df, "examples/df_fit_1k.csv");
    ///
    /// Ok(())
    /// }
    /// ```
    pub fn fit_steady(&mut self, rate: std::ops::Range<f64>, obs: &[f64], title: &str) {
        let dur = std::time::Duration::new(60 * 60 * self.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(title).exists();
        match exists {
            true => {
                let mut gof: Vec<Gof> = Gof::read(title).unwrap();
                rec.append(&mut gof);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let mut new = self.steady(rate.clone(), obs);
            {
                rec.append(&mut new);
            }
            Gof::record(&mut rec, title).unwrap();
        }
    }

    /// Return quantile statistics for transit times at randomly selected rates for a given duration.
    /// Calls [transits](#method.transits).
    pub fn fit_transits(
        &mut self,
        rate: std::ops::Range<f64>,
        title: &str,
    ) -> Result<(), errors::ResError> {
        let dur = std::time::Duration::new(60 * 60 * self.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(title).exists();
        match exists {
            true => {
                let mut stats: Vec<Transits> = Transits::read(title)?;
                rec.append(&mut stats);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let mut new = self.transits(rate.clone());
            {
                rec.append(&mut new);
            }
            utils::record(&mut rec, title)?;
        }
        Ok(())
    }

    /// Return weighted mean transit time of particles in the reservoir.
    ///
    /// # Examples
    ///
    /// ```{rust}
    /// use::reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     // model parameters
    ///     let batch = 10;  // fit 10 input/output pairs at a time using rayon
    ///     let period = 3000.0; // run simulations for 3000 years
    ///     let runs = 10; // run 10 simulated accumulations per candidate pair for goodness-of-fit
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(Reservoir::new().input(&0.73)?.output(&0.73)?)
    ///         .batch(batch)
    ///         .period(&period)
    ///         .runs(runs)
    ///         .range(period as i32);
    ///
    ///     let transits = model.fit_transit_time();
    ///     Ok(())
    /// }
    /// ```
    pub fn fit_transit_time(&mut self) -> f64 {
        let pmf = self.clone().transit_times();
        let index = utils::index_f64(0.0, pmf.len() as f64, 1.0);
        let wgt_mn = index
            .iter()
            .zip(pmf.iter())
            .map(|(a, b)| a * b)
            .collect::<Vec<f64>>();
        let mean = wgt_mn.iter().fold(0.0, |acc, x| acc + x);
        mean
    }

    /// Returns quantile statistics on transit times for rates in a given range.
    /// Calls [fit_transit_rate](#method.fit_transit_rate).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     // model parameters
    ///     let batch = 10;  // fit 10 input/output pairs at a time using rayon
    ///     let period = 3000.0; // run simulations for 3000 years
    ///     let runs = 10; // run 10 simulated accumulations per candidate pair for goodness-of-fit
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(Reservoir::new().input(&0.73)?.output(&0.73)?)
    ///         .batch(batch)
    ///         .period(&period)
    ///         .runs(runs)
    ///         .range(period as i32);
    ///
    ///     let transits = model.transits((0.01..1.5));
    ///     Ok(())
    /// }
    /// ```
    pub fn transits(&mut self, rate: std::ops::Range<f64>) -> Vec<Transits> {
        info!("Randomly select rates from the given rates to fit.");
        let mut rates = Vec::with_capacity(self.batch);
        info!("Clone the parent model and give each copy a different reservoir rate.");
        let mut fits = Vec::with_capacity(self.batch);
        for i in 0..self.batch {
            rates.push(
                rand::distributions::Uniform::from(rate.clone()).sample(&mut self.reservoir.range),
            );
            fits.push(
                self.clone().reservoir(
                    self.reservoir
                        .clone()
                        .input(&rates[i])
                        .unwrap()
                        .output(&rates[i])
                        .unwrap(),
                ),
            );
        }
        let stats: Vec<f64> = fits.iter().map(|x| x.clone().fit_transit_time()).collect();
        let transits: Vec<Transits> = stats
            .iter()
            .zip(rates.iter())
            .map(|(a, b)| Transits { rate: *a, mean: *b })
            .collect();
        transits
    }

    /// Creates a new model from a given `reservoir` with default `period`, `runs` and `batch` values.
    pub fn new(reservoir: Reservoir) -> Self {
        Model {
            reservoir,
            period: 30000.0,
            runs: 100,
            batch: 10,
            duration: 1,
            range: 10000,
        }
    }

    /// Sets the model period to desired time in years.
    pub fn period(mut self, years: &f64) -> Self {
        self.period = years.to_owned();
        self
    }

    /// Sets the number of runs per sample to estimate goodness-of-fit.
    pub fn range(mut self, limit: i32) -> Self {
        self.range = limit;
        self
    }

    /// Sets the model reservoir.
    pub fn reservoir(mut self, res: Reservoir) -> Self {
        self.reservoir = res;
        self
    }

    /// Sets the number of runs per sample to estimate goodness-of-fit.
    pub fn runs(mut self, times: usize) -> Self {
        self.runs = times;
        self
    }

    /// Randomly selects a rate from ranges `rate` for a steady state reservoir,
    /// and simulates accumulation records
    /// in batches using [fit_rate](#method.fit_rate).  Returns the selected input/output pair and the mean
    /// goodness-of-fit statistics compared to `obs` for each pair.
    /// Called by [fit_steady](#method.fit_steady).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let ages = vec![10.0, 42.0, 77.7, 12.0, 99.9, 10000.0, 777.7];
    ///     let mut debris_flows = Reservoir::new()
    ///         .input(&0.862)?
    ///         .output(&0.862)?;
    ///
    ///     // model parameters
    ///     let batch = 10;  // fit 10 input/output pairs at a time using rayon
    ///     let period = 3000.0; // run simulations for 30000 years
    ///     let runs = 20; // run 1000 simulated accumulations per candidate pair for goodness-of-fit
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(debris_flows)
    ///         .batch(batch)
    ///         .period(&period)
    ///         .runs(runs);
    ///     // fit a batch of 10 randomly selected rate pairs (from range 0.01 to 1.0)
    ///     // to observed debris flows
    ///     // by running 20 simulations for 3000 years for each pair
    ///     let gofs = model.steady(0.01..1.0, &ages);
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn steady(&mut self, rate: std::ops::Range<f64>, obs: &[f64]) -> Vec<Gof> {
        info!("Randomly select rates from the given rates to fit.");
        let mut rates = Vec::with_capacity(self.batch);
        info!("Clone the parent model and give each copy a different reservoir rate.");
        let mut fits = Vec::with_capacity(self.batch);
        for i in 0..self.batch {
            rates.push(
                rand::distributions::Uniform::from(rate.clone()).sample(&mut self.reservoir.range),
            );
            fits.push(
                self.clone().reservoir(
                    self.reservoir
                        .clone()
                        .input(&rates[i])
                        .unwrap()
                        .output(&rates[i])
                        .unwrap(),
                ),
            );
        }
        let gof: Vec<Fits> = fits.par_iter().map(|x| x.clone().fit_rate(obs)).collect();
        let mut gofs = Vec::with_capacity(self.batch);
        for i in 0..self.batch {
            gofs.push(Gof {
                input: rates[i],
                output: rates[i],
                ad: gof[i].fit.ad,
                ada: gof[i].fit.ada,
                chs: gof[i].fit.chs,
                kp: gof[i].fit.kp,
                ks: gof[i].fit.ks,
                ksa: gof[i].fit.ksa,
                n: gof[i].n,
            })
        }
        gofs
    }

    /// Selects from among a number of simulated accumulation records, the record most
    /// characteristic of the average number and age distribution of deposits.
    ///  - `self` is a model with reservoir characteristics set.
    ///  - `bins` is the number of bins with which to construct a cdf.
    ///  - Returns the record with the lowest loss statistic.
    ///
    /// For large numbers of runs, the cdf becomes a large vector.  Selecting a smaller
    /// number of `bins` reduces memory consumption and increases speed, but can reduce accuracy.
    ///
    /// The loss statistic is the K-S stat plus a normalized 'distance from average number
    /// of deposits'.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///
    ///     let mut debris_flows = Reservoir::new()
    ///         .input(&0.862)?
    ///         .output(&0.862)?;
    ///
    ///     // model parameters
    ///     let period = 3000.0; // run simulations for 30000 years
    ///     let runs = 20; // run 1000 simulated accumulations per candidate pair for goodness-of-fit
    ///     let bins = 500; // split observation into bins for deriving CDF
    ///
    ///     // create reservoir model using builder pattern
    ///     let mut model = Model::new(debris_flows)
    ///         .period(&period)
    ///         .runs(runs);
    ///
    ///     // sample a stereotypical record from 1000 runs of 30000 years
    ///     let eg = model.stereotype(bins);
    ///     let egx = model.stereotype(bins);
    ///     // compare the CDF of the synthetic example to the observed debris-flow deposit record
    ///     plot::comp_cdf(&eg, &egx, "examples/transit_stereotype.png");
    ///
    ///     Ok(())
    /// }
    ///```
    pub fn stereotype(&mut self, bins: usize) -> Vec<f64> {
        let mut res: Vec<Model> = Vec::with_capacity(self.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds: Vec<u64> = seeder
            .sample_iter(&mut self.reservoir.range)
            .take(self.runs)
            .collect();

        for seed in seeds {
            // make boot number copies of reservoir
            res.push(self.clone().reservoir(self.reservoir.clone().range(seed)));
        }
        let res: Vec<Reservoir> = res.par_iter().cloned().map(|x| x.reservoir.sim()).collect(); // simulate accumulation record for each copy
        let mut ns: Vec<f64> = res
            .par_iter()
            .cloned()
            .map(|x| x.mass.len() as f64)
            .collect(); // number of deposits in reservoir
        let mid_n = utils::median(&ns); // median number of deposits
        ns = ns
            .iter()
            .map(|x| rand_distr::num_traits::abs((x / mid_n) - 1.0))
            .collect(); // distance from median length
                        // collect reservoir masses into single vector and calculate the cdf
        let mut rec = Vec::new(); // vector of mass
        for r in res.clone() {
            rec.append(&mut r.mass.clone()); // add each run to make one long vector
        }
        let cdf = utils::cdf_bin(&rec, bins); // subsample vector to length bins

        // TODO:  parallelize
        let gof: Vec<Fit> = res.par_iter().cloned().map(|x| x.gof1(&cdf)).collect(); // ks and kp values
        let ks: Vec<f64> = gof.par_iter().cloned().map(|x| x.ks).collect(); // clip to just ks values
        let mut least = 1.0; // test for lowest fit (set to high value)
        let mut low = Reservoir::new(); // initialize variable to hold lowest fit
        for (i, val) in ns.iter().enumerate() {
            let loss = ks[i] + val; // loss function
            if loss < least {
                // if lowest value
                low = res[i].clone(); // copy to low
                least = loss; // set least to new low value
            }
        }

        low.mass
    }

    /// Produces a raw vector of observations of simulated transit times from multiple runs.
    /// Used to produce statistics on the mean, median and statistical moments.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///
    ///     let period: f64 = 2000.0;  // length of time to simulate accumulation in the reservoir in years
    ///     let runs: usize = 20; // number of simulations to run for estimation
    ///     let df_rt: f64 = 0.72; // rate for steady state debris-flow accumulation
    ///     let until = 4095; // range limit for tracking transit times
    ///
    ///     // to estimate transit times, omit inherit age from charcoal
    ///     let debris_flows = Reservoir::new().input(&df_rt)?.output(&df_rt)?;
    ///     let mut df_mod = Model::new(debris_flows).period(&period).runs(runs).range(until);
    ///
    ///     // vector of transit times
    ///     let df_t = df_mod.transit_times();
    ///     println!("Debris-flow transit time quantiles are {:?}", utils::quantiles(&df_t));
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn transit_times(&mut self) -> Vec<f64> {
        let mut res = Vec::with_capacity(self.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds: Vec<u64> = seeder
            .sample_iter(&mut self.reservoir.range)
            .take(self.runs)
            .collect();
        for seed in seeds {
            res.push(self.reservoir.clone().range(seed));
        }
        res = res.par_iter().cloned().map(|x| x.sim()).collect();

        let index = 0..self.range;
        let mut out = vec![Complex::new(0.0, 0.0); index.len()];
        let mut real_planner = RealFftPlanner::<f64>::new();
        info!("Constructing FftPlanner for fourier transforms.");
        let r2c = real_planner.plan_fft_forward(index.len() + 1);
        for r in res {
            let cdf = utils::cdf_rng(&r.mass, &index);
            let mut pmf = utils::pmf_from_cdf(&cdf);
            let mut spectrum = r2c.make_output_vec();
            info!("Forward transforming the pmf using fft.");
            r2c.process(&mut pmf, &mut spectrum).unwrap();
            info!("Summing transformed pmfs.");
            out = out
                .iter()
                .zip(spectrum.iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<Complex<f64>>>();
        }
        info!("Constructing FftPlanner for inverse fourier transforms.");
        let c2r = real_planner.plan_fft_inverse(index.len() + 1);
        info!("Inverse fourier transform using fft.");
        let mut out_data = c2r.make_output_vec();
        c2r.process(&mut out, &mut out_data).unwrap();
        info!("Normalize output by dividing by length.");
        out_data = out_data
            .par_iter()
            .map(|a| a / index.len() as f64)
            .collect::<Vec<f64>>();
        info!("Normalize output by dividing by sum of pmfs (having added several pmfs together).");
        let sum_out = out_data.iter().fold(0.0, |acc, x| acc + *x);
        out_data = out_data
            .par_iter()
            .map(|a| a / sum_out)
            .collect::<Vec<f64>>();

        out_data
    }
}

impl Default for Model {
    fn default() -> Self {
        Model {
            reservoir: Reservoir::new(),
            period: 30000.0,
            runs: 100,
            batch: 10,
            duration: 1,
            range: 10000,
        }
    }
}

/// Model manager struct for control flow of modeling operations.
#[derive(Clone, Debug)]
pub struct ModelManager {
    batch: usize,
    capture_gravels: std::ops::Range<f64>,
    capture_fines: std::ops::Range<f64>,
    duration: u64,
    fines: bool,
    index: std::ops::Range<i32>,
    obs: Vec<f64>,
    obs_len: usize,
    period: f64,
    periods: std::ops::Range<f64>,
    range: rand::rngs::StdRng,
    rates: std::ops::Range<f64>,
    runs: usize,
    source_runs: usize,
    storage_gravels: std::ops::Range<f64>,
    storage_fines: std::ops::Range<f64>,
    thresholds: Thresholds,
    turnover: std::ops::Range<f64>,
}

/// Methods for modeling.
impl ModelManager {
    /// Reference a FluvialModel type to create a new ModelManager instance.
    pub fn new() -> ModelManager {
        ModelManager {
            batch: 1,
            capture_fines: 0.0..1.0,
            capture_gravels: 0.0..1.0,
            duration: 1,
            fines: false,
            index: 0..1000,
            obs: Vec::new(),
            obs_len: 0,
            period: 100.0,
            periods: 100.0..1000.0,
            range: rand::SeedableRng::seed_from_u64(777),
            rates: 0.01..2.0,
            runs: 10,
            source_runs: 10,
            storage_fines: 0.0..1.0,
            storage_gravels: 0.0..1.0,
            thresholds: Thresholds::new(0.0, 0.0, 0.0, 0.0),
            turnover: 200.0..300.0,
        }
    }

    /// Set the batch size for runs.
    pub fn batch(mut self, batch: usize) -> Self {
        self.batch = batch;
        self
    }

    /// Set capture rate of fines.
    pub fn capture_fines(mut self, capture_fines: std::ops::Range<f64>) -> Self {
        self.capture_fines = capture_fines;
        self
    }

    /// Set capture rate of gravels.
    pub fn capture_gravels(mut self, capture_gravels: std::ops::Range<f64>) -> Self {
        self.capture_gravels = capture_gravels;
        self
    }

    /// Set duration of timed runs.
    pub fn duration(mut self, duration: u64) -> Self {
        self.duration = duration;
        self
    }

    /// Set to true to produce fines with Fluvial::sim().
    pub fn fines(mut self, fines: bool) -> Self {
        self.fines = fines;
        self
    }

    /// Set the index range for transit_times().
    pub fn index(mut self, index: std::ops::Range<i32>) -> Self {
        self.index = index;
        self
    }

    /// Set period length of model runs.
    pub fn obs(mut self, obs: &[f64]) -> Self {
        self.obs = obs.to_owned();
        self
    }

    /// Set number of observations to collect for representative sample.
    pub fn obs_len(mut self, obs: &[f64]) -> Self {
        self.obs_len = obs.len();
        self
    }

    /// Set period length of model runs.
    pub fn period(mut self, period: f64) -> Self {
        self.period = period;
        self
    }

    /// Set the period range for Reservoir::fit_rates().
    pub fn periods(mut self, periods: std::ops::Range<f64>) -> Self {
        self.periods = periods;
        self
    }

    /// Set model seed.
    pub fn range(mut self, seed: u64) -> Self {
        self.range = rand::SeedableRng::seed_from_u64(seed);
        self
    }

    /// Set the range of input/output rates for Reservoir::fit_rates().
    pub fn rates(mut self, rates: std::ops::Range<f64>) -> Self {
        self.rates = rates;
        self
    }

    /// Set number of model runs.
    pub fn runs(mut self, runs: usize) -> Self {
        self.runs = runs;
        self
    }

    /// Return a vector of struct clones with unique but reproducible seeds.
    pub fn seed_clones(mut self) -> Vec<ModelManager> {
        let mut res = Vec::with_capacity(self.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        for _ in 0..self.runs {
            res.push(self.clone().range(seeder.sample(&mut self.range)));
        }
        res
    }

    /// Set the number of runs when estimating source mass (per run of Fluvial::sim()).
    pub fn source_runs(mut self, source_runs: usize) -> Self {
        self.source_runs = source_runs;
        self
    }

    /// Set the range of storage rates for fines.
    pub fn storage_fines(mut self, storage_fines: std::ops::Range<f64>) -> Self {
        self.storage_fines = storage_fines;
        self
    }

    /// Set the range of storage rates for gravels.
    pub fn storage_gravels(mut self, storage_gravels: std::ops::Range<f64>) -> Self {
        self.storage_gravels = storage_gravels;
        self
    }

    /// Set statistical thresholds for goodness-of-fit tests.
    pub fn thresholds(mut self, ad: f64, ch: f64, kp: f64, ks: f64) -> Self {
        self.thresholds = Thresholds::new(ad, ch, kp, ks);
        self
    }

    /// Set the range of turnover periods for sediment storage.
    pub fn turnover(mut self, turnover: std::ops::Range<f64>) -> Self {
        self.turnover = turnover;
        self
    }
}

impl Default for ModelManager {
    fn default() -> ModelManager {
        ModelManager::new()
    }
}

/// Struct for model run statistics.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct FluvialFit {
    capture_rate_fines: f64,
    capture_rate_gravels: f64,
    storage_rate_fines: f64,
    storage_rate_gravels: f64,
    turnover: f64,
    ad1: f64,
    ad2: f64,
    ch: f64,
    kp: f64,
    ks1: f64,
    ks2: f64,
}

impl FluvialFit {
    /// Convert csv record to Gof struct.
    pub fn read(path: &str) -> Result<Vec<FluvialFit>, errors::ResError> {
        let mut fit = Vec::new();
        let var = std::fs::File::open(path)?;
        let mut rdr = csv::Reader::from_reader(var);
        for result in rdr.records() {
            let row = result?;
            let row: FluvialFit = row.deserialize(None)?;
            fit.push(row);
        }
        Ok(fit)
    }

    /// Write statistical results to csv file.
    pub fn record(rec: &mut Vec<FluvialFit>, title: &str) -> Result<(), errors::ResError> {
        let mut wtr = csv::Writer::from_path(title)?;
        for i in rec {
            wtr.serialize(i)?;
        }
        wtr.flush()?;
        Ok(())
    }
}

/// Fluvial traversal struct for gravels and fines.
#[derive(Clone, Debug)]
pub struct Fluvial {
    capture_rate_fines: f64,
    capture_rate_gravels: f64,
    manager: ModelManager,
    mass: Vec<f64>,
    source: Vec<f64>,
    storage_rate_fines: f64,
    storage_rate_gravels: f64,
    turnover: f64,
}

impl Fluvial {
    /// New fluvial struct
    pub fn new() -> Fluvial {
        Fluvial {
            capture_rate_fines: 0.0,
            capture_rate_gravels: 0.0,
            manager: ModelManager::new(),
            mass: Vec::new(),
            source: Vec::new(),
            storage_rate_fines: 0.0,
            storage_rate_gravels: 0.0,
            turnover: 0.0,
        }
    }

    /// Cherry pick the closest fit to an observed distribution.
    pub fn cherry_pick(self) -> Vec<f64> {
        let sims = self
            .clone()
            .manager
            .seed_clones()
            .par_iter()
            .map(|x| self.clone().manager(x).sim())
            .collect::<Vec<Fluvial>>();
        let fit = sims
            .par_iter()
            .map(|s| s.clone().gof(&self.manager.obs))
            .collect::<Vec<Vec<f64>>>();
        let min_ks = fit
            .iter()
            .map(|f| f[4])
            .fold(f64::INFINITY, |a, b| a.min(b));
        let min_id = fit
            .par_iter()
            .enumerate()
            .filter(|(_, val)| (val[4] - min_ks) < 0.0001)
            .map(|(i, _)| i)
            .collect::<Vec<usize>>();
        sims[min_id[0]].clone().mass
    }

    /// Hit rate for statistical thresholds.
    pub fn hit_rate(self) -> Vec<f64> {
        let sims = self
            .clone()
            .manager
            .seed_clones()
            .par_iter()
            .map(|x| self.clone().manager(x).sim())
            .collect::<Vec<Fluvial>>();
        let fit = sims
            .par_iter()
            .map(|s| s.clone().gof(&self.manager.obs))
            .collect::<Vec<Vec<f64>>>();
        let hit_ad = fit
            .par_iter()
            .filter(|f| f[0] <= self.manager.thresholds.ad)
            .count() as f64
            / self.manager.runs as f64;
        let hit_ch = fit
            .par_iter()
            .filter(|f| f[2] <= self.manager.thresholds.ch)
            .count() as f64
            / self.manager.runs as f64;
        let hit_kp = fit
            .par_iter()
            .filter(|f| f[3] <= self.manager.thresholds.kp)
            .count() as f64
            / self.manager.runs as f64;
        let hit_ks = fit
            .par_iter()
            .filter(|f| f[4] <= self.manager.thresholds.ks)
            .count() as f64
            / self.manager.runs as f64;
        vec![hit_ad, hit_ch, hit_kp, hit_ks]
    }

    /// Fits a given flux probability and storage rate to an empiric record.
    pub fn hit_rates(mut self) -> FluvialFit {
        let capture_fines = rand::distributions::Uniform::from(self.manager.capture_fines.clone());
        let capture_rate_fines = capture_fines.sample(&mut self.manager.range);
        let capture_gravels =
            rand::distributions::Uniform::from(self.manager.capture_gravels.clone());
        let capture_rate_gravels = capture_gravels.sample(&mut self.manager.range);
        let storage_fines = rand::distributions::Uniform::from(self.manager.storage_fines.clone());
        let storage_rate_fines = storage_fines.sample(&mut self.manager.range);
        let storage_gravels =
            rand::distributions::Uniform::from(self.manager.storage_gravels.clone());
        let storage_rate_gravels = storage_gravels.sample(&mut self.manager.range);
        let turnover_range = rand::distributions::Uniform::from(self.manager.turnover.clone());
        let turnover = turnover_range.sample(&mut self.manager.range);

        let fit = self
            .capture_rate_fines(capture_rate_fines)
            .capture_rate_gravels(capture_rate_gravels)
            .storage_rate_fines(storage_rate_fines)
            .storage_rate_gravels(storage_rate_gravels)
            .turnover(&turnover)
            .hit_rate();
        FluvialFit {
            capture_rate_fines,
            capture_rate_gravels,
            storage_rate_fines,
            storage_rate_gravels,
            turnover,
            ad1: fit[0],
            ad2: 0.0,
            ch: fit[1],
            kp: fit[2],
            ks1: fit[3],
            ks2: 0.0,
        }
    }

    /// Run fit_rates() for a set duration.
    pub fn hit_rates_timed(mut self, path: &str) -> Result<Vec<FluvialFit>, errors::ResError> {
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);

        let dur = std::time::Duration::new(60 * 60 * self.manager.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(path).exists();
        match exists {
            true => {
                let mut stats: Vec<FluvialFit> = FluvialFit::read(path)?;
                rec.append(&mut stats);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let mut fit = self.clone();
            fit.manager = fit
                .manager
                .range(seeder.sample(&mut self.manager.range))
                .clone();
            let new = fit.hit_rates();
            {
                rec.push(new);
            }
            FluvialFit::record(&mut rec, path)?;
        }
        Ok(rec)
    }

    /// Fit number of runs to gof tests and return mean of each.
    pub fn fit(self, other: &[f64]) -> Vec<f64> {
        let mut rng = self.manager.range.clone();
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds: Vec<u64> = seeder
            .sample_iter(&mut rng)
            .take(self.manager.runs)
            .collect();

        let mut ads = Vec::with_capacity(self.manager.runs as usize);
        let mut adas = Vec::with_capacity(self.manager.runs as usize);
        let mut chs = Vec::with_capacity(self.manager.runs as usize);
        let mut kps = Vec::with_capacity(self.manager.runs as usize);
        let mut kss = Vec::with_capacity(self.manager.runs as usize);
        let mut ksas = Vec::with_capacity(self.manager.runs as usize);

        for seed in seeds {
            let fit = self
                .clone()
                .manager(&self.manager.clone().range(seed))
                .gof(other);
            ads.push(fit[0]);
            adas.push(fit[1]);
            chs.push(fit[2]);
            kps.push(fit[3]);
            kss.push(fit[4]);
            ksas.push(fit[5]);
        }

        let res = vec![
            utils::mean(&ads),
            utils::mean(&adas),
            utils::mean(&chs),
            utils::mean(&kps),
            utils::mean(&kss),
            utils::mean(&ksas),
        ];
        res
    }

    /// Fits a given flux probability and storage rate to an empiric record.
    pub fn fit_rate(mut self, other: &[f64]) -> FluvialFit {
        let capture_fines = rand::distributions::Uniform::from(self.manager.capture_fines.clone());
        let capture_rate_fines = capture_fines.sample(&mut self.manager.range);
        let capture_gravels =
            rand::distributions::Uniform::from(self.manager.capture_gravels.clone());
        let capture_rate_gravels = capture_gravels.sample(&mut self.manager.range);
        let storage_fines = rand::distributions::Uniform::from(self.manager.storage_fines.clone());
        let storage_rate_fines = storage_fines.sample(&mut self.manager.range);
        let storage_gravels =
            rand::distributions::Uniform::from(self.manager.storage_gravels.clone());
        let storage_rate_gravels = storage_gravels.sample(&mut self.manager.range);
        let turnover_range = rand::distributions::Uniform::from(self.manager.turnover.clone());
        let turnover = turnover_range.sample(&mut self.manager.range);

        let fit = self
            .capture_rate_fines(capture_rate_fines)
            .capture_rate_gravels(capture_rate_gravels)
            .storage_rate_fines(storage_rate_fines)
            .storage_rate_gravels(storage_rate_gravels)
            .turnover(&turnover)
            .fit(other);
        FluvialFit {
            capture_rate_fines,
            capture_rate_gravels,
            storage_rate_fines,
            storage_rate_gravels,
            turnover,
            ad1: fit[1],
            ad2: fit[0],
            ch: fit[2],
            kp: fit[3],
            ks1: fit[5],
            ks2: fit[4],
        }
    }

    /// Return a stereotypical mass age distribution at the current model parameters.
    pub fn stereotype_rate(self) -> Vec<f64> {
        // let bins = 1000;
        let mut rng = self.manager.range.clone();
        let mut res: Vec<Fluvial> = Vec::with_capacity(self.manager.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds: Vec<u64> = seeder
            .sample_iter(&mut rng)
            .take(self.manager.runs)
            .collect();

        for seed in seeds {
            // make boot number copies of reservoir
            res.push(self.clone().manager(&self.manager.clone().range(seed)));
        }
        res = res
            .par_iter()
            .cloned()
            .map(|x| x.sim())
            .collect::<Vec<Fluvial>>(); // simulate accumulation record for each copy
        let mut ns = res
            .par_iter()
            .cloned()
            .map(|x| x.mass.len() as f64)
            .collect::<Vec<f64>>(); // number of deposits in reservoir
        let mid_n = utils::median(&ns); // median number of deposits
        ns = ns
            .par_iter()
            .map(|x| rand_distr::num_traits::abs((x / mid_n) - 1.0))
            .collect::<Vec<f64>>(); // distance from median length
                                    // collect reservoir masses into single vector and calculate the cdf
        let mut rec = Vec::new(); // vector of mass
        for r in res.clone() {
            rec.append(&mut r.mass.clone()); // add each run to make one long vector
        }
        // let cdf = utils::cdf_bin(&rec, bins); // subsample vector to length bins

        let fit = res
            .par_iter()
            .cloned()
            .map(|x| x.gof(&rec))
            .collect::<Vec<Vec<f64>>>(); // ks and kp values
        let ks = fit.par_iter().cloned().map(|x| x[5]).collect::<Vec<f64>>(); // clip to just ks values
        let mut least = 1.0; // test for lowest fit (set to high value)
        let mut low = Fluvial::new(); // initialize variable to hold lowest fit
        for (i, val) in ns.iter().enumerate() {
            let loss = ks[i] + val; // loss function
            if loss < least {
                // if lowest value
                low = res[i].clone(); // copy to low
                least = loss; // set least to new low value
            }
        }

        low.mass
    }

    /// Fit a range of flux probabilities and storage rates to an empiric record.
    pub fn fit_rates(mut self, other: &[f64]) -> Vec<FluvialFit> {
        // let mut rng = self.manager.range.clone();
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let mut res = Vec::new();
        for _ in 0..self.manager.batch {
            let mut fit = self.clone();
            fit.manager = fit
                .manager
                .range(seeder.sample(&mut self.manager.range))
                .clone();
            res.push(fit.fit_rate(other));
        }
        res
    }

    /// Run fit_rates() for a set duration.
    pub fn fit_rates_timed(
        mut self,
        other: &[f64],
        path: &str,
    ) -> Result<Vec<FluvialFit>, errors::ResError> {
        // let mut rng = self.manager.range.clone();
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);

        let dur = std::time::Duration::new(60 * 60 * self.manager.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(path).exists();
        match exists {
            true => {
                let mut stats: Vec<FluvialFit> = FluvialFit::read(path)?;
                rec.append(&mut stats);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let mut fit = self.clone();
            fit.manager = fit
                .manager
                .range(seeder.sample(&mut self.manager.range))
                .clone();
            let mut new = fit.fit_rates(other);
            {
                rec.append(&mut new);
            }
            FluvialFit::record(&mut rec, path)?;
        }
        Ok(rec)
    }

    /// Represents the probability of fines routing to the storage pool.
    pub fn capture_rate_fines(mut self, capture_rate_fines: f64) -> Self {
        self.capture_rate_fines = capture_rate_fines;
        self
    }

    /// Represents the probability of gravels routing to the storage pool.
    pub fn capture_rate_gravels(mut self, capture_rate_gravels: f64) -> Self {
        self.capture_rate_gravels = capture_rate_gravels;
        self
    }

    /// Goodness-of-fit test suite.
    pub fn gof(self, other: &[f64]) -> Vec<f64> {
        utils::gof(&self.sim().mass, other)
    }

    /// Assign ModelManager to Fluvial struct
    pub fn manager(mut self, manager: &ModelManager) -> Self {
        self.manager = manager.to_owned();
        self
    }

    /// Runs sim() a number of times specified by ModelManager struct.
    pub fn n_sim(self) -> Vec<f64> {
        /*let mut rec = Vec::new();
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);

        for _ in 0..self.manager.source_runs {
            let mut fit = self.clone();
            fit.manager = fit.manager.range(seeder.sample(&mut self.manager.range)).clone();
            let mut new = fit.sim().mass;
            {
                rec.append(&mut new);
            }
        }
        */
        let mut res = Vec::new();
        let _ = self
            .clone()
            .manager
            .seed_clones()
            .iter()
            .map(|x| res.extend(self.clone().manager(x).sim().mass))
            .collect::<()>();
        res
    }

    /// Prints csv of reservoir mass to file path.
    pub fn record_mass(mut self, path: &str) -> Result<(), errors::ResError> {
        utils::record(&mut self.mass, path)?;
        Ok(())
    }

    /// Simulates removal from the reservoir over time based on the rate.
    pub fn sim(mut self) -> Self {
        let wts = rand::distributions::WeightedIndex::new(&self.source).unwrap();
        let mut index = Vec::new();
        for i in self.manager.index.clone() {
            index.push(i as f64);
        }
        let mut source_flux = Vec::new();
        for _ in 0..self.manager.obs_len {
            source_flux.push(index[wts.sample(&mut self.manager.range)]);
        }
        info!("Selection probability for storage.");
        let mut idx = Vec::new();
        let mut ps = Vec::new();
        let mut thresh = 0.0;
        for i in 0..self.manager.period as i32 {
            idx.push(i as f64);
            match self.manager.fines {
                true => {
                    thresh = self.capture_rate_fines;
                    ps.push(utils::fish(
                        self.storage_rate_fines,
                        i as f64 / self.turnover,
                    ));
                }
                false => {
                    thresh = self.capture_rate_gravels;
                    ps.push(utils::fish(
                        self.storage_rate_gravels,
                        i as f64 / self.turnover,
                    ));
                }
            }
        }
        let wts = rand::distributions::WeightedIndex::new(&ps).unwrap();
        for item in &mut source_flux {
            let roll = self.manager.range.gen_range(0.0..1.0);
            if roll < thresh {
                *item += idx[wts.sample(&mut self.manager.range)];
            }
        }

        self.mass = source_flux;
        self
    }

    /// Sets source flux for reservoirs.
    pub fn source(mut self, source: Vec<f64>) -> Self {
        let source_flux = utils::cdf_rng(&source, &self.manager.index);
        // let source_flux = source.clone().transit_times();

        // if self.manager.fines {
        //     let source_gravel = source.model(&model.fines(false)).transit_times();
        //     source_flux.extend(source_gravel);
        // }
        self.source = source_flux;
        self
    }

    /// Sets source flux from csv.
    /// To save and reuse large samples that are computationally intensive.
    pub fn source_from_csv(mut self, path: &str) -> Result<Self, errors::ResError> {
        let source = utils::read_f64(path)?;
        self.source = source;
        Ok(self)
    }

    /// Sets the rate of the Poisson distribution defining the PMF of storage time for fines.
    pub fn storage_rate_fines(mut self, storage_rate_fines: f64) -> Self {
        self.storage_rate_fines = storage_rate_fines;
        self
    }

    /// Sets the rate of the Poisson distribution defining the PMF of storage time for fines.
    pub fn storage_rate_gravels(mut self, storage_rate_gravels: f64) -> Self {
        self.storage_rate_gravels = storage_rate_gravels;
        self
    }

    /// Return CDF of transit times based on model parameters.
    pub fn transit_times(self) -> Vec<f64> {
        let res = self
            .clone()
            .manager
            .seed_clones()
            .par_iter()
            .map(|x| self.clone().manager(x).sim())
            .collect::<Vec<Fluvial>>();

        let index = self.manager.index;
        let mut out = vec![Complex::new(0.0, 0.0); index.len()];
        let mut real_planner = RealFftPlanner::<f64>::new();
        info!("Constructing FftPlanner for fourier transforms.");
        let r2c = real_planner.plan_fft_forward(index.len() + 1);
        for r in res {
            let cdf = utils::cdf_rng(&r.mass, &index);
            let mut pmf = utils::pmf_from_cdf(&cdf);
            let mut spectrum = r2c.make_output_vec();
            info!("Forward transforming the pmf using fft.");
            r2c.process(&mut pmf, &mut spectrum).unwrap();
            info!("Summing transformed pmfs.");
            out = out
                .iter()
                .zip(spectrum.iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<Complex<f64>>>();
        }
        info!("Constructing FftPlanner for inverse fourier transforms.");
        let c2r = real_planner.plan_fft_inverse(index.len() + 1);
        info!("Inverse fourier transform using fft.");
        let mut out_data = c2r.make_output_vec();
        c2r.process(&mut out, &mut out_data).unwrap();
        info!("Normalize output by dividing by length.");
        out_data = out_data
            .par_iter()
            .map(|a| a / index.len() as f64)
            .collect::<Vec<f64>>();
        info!("Normalize output by dividing by sum of pmfs (having added several pmfs together).");
        let sum_out = out_data.iter().fold(0.0, |acc, x| acc + *x);
        out_data = out_data
            .par_iter()
            .map(|a| a / sum_out)
            .collect::<Vec<f64>>();

        out_data
    }

    /// Sets length of the turnover period.
    pub fn turnover(mut self, turnover: &f64) -> Self {
        self.turnover = turnover.to_owned();
        self
    }
}

impl Default for Fluvial {
    fn default() -> Fluvial {
        Fluvial::new()
    }
}

/// Struct for recording reservoir characteristics.
#[derive(Clone, Debug)]
pub struct Reservoir {
    input: Option<Exp<f64>>,
    model: ModelManager,
    /// Ages (in years) of deposits accumulated in reservoir.
    pub mass: Vec<f64>,
    output: Option<Exp<f64>>,
    flux: Vec<f64>,
    inherit: Option<Vec<f64>>,
    range: rand::rngs::StdRng,
}

impl Reservoir {
    /// Fit reservoir to observed record.
    pub fn fit_reservoir(self) -> Vec<f64> {
        let reservoirs = self
            .clone()
            .model
            .seed_clones()
            .par_iter()
            .map(|x| self.clone().model(x).sim())
            .collect::<Vec<Reservoir>>();
        let fits = reservoirs
            .par_iter()
            .map(|x| utils::gof(&x.mass, &self.model.obs))
            .collect::<Vec<Vec<f64>>>();
        let ad = fits.par_iter().map(|x| x[0]).collect::<Vec<f64>>();
        let ch = fits.par_iter().map(|x| x[2]).collect::<Vec<f64>>();
        let kp = fits.par_iter().map(|x| x[3]).collect::<Vec<f64>>();
        let ks = fits.par_iter().map(|x| x[4]).collect::<Vec<f64>>();
        let ns = reservoirs
            .par_iter()
            .map(|x| x.mass.len() as f64)
            .collect::<Vec<f64>>();
        vec![
            utils::mean(&ad),
            utils::mean(&ch),
            utils::mean(&kp),
            utils::mean(&ks),
            utils::mean(&ns),
        ]
    }

    /// Fit reservoirs to observed record.
    pub fn fit_reservoirs(mut self) -> Result<Vec<ReservoirFit>, errors::ResError> {
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds = seeder
            .sample_iter(&mut self.model.range)
            .take(self.model.batch)
            .collect::<Vec<u64>>();
        let rates = rand::distributions::Uniform::from(self.model.rates.clone())
            .sample_iter(&mut self.model.range)
            .take(self.model.batch)
            .collect::<Vec<f64>>();
        let periods = rand::distributions::Uniform::from(self.model.periods.clone())
            .sample_iter(&mut self.model.range)
            .take(self.model.batch)
            .collect::<Vec<f64>>();

        let mut reservoirs = Vec::with_capacity(self.model.batch);
        for i in 0..self.model.batch {
            let model = self.model.clone().period(periods[i]).range(seeds[i]);
            reservoirs.push(
                self.clone()
                    .input(&rates[i])?
                    .output(&rates[i])?
                    .model(&model),
            );
        }

        let fits = reservoirs
            .iter()
            .map(|x| x.clone().fit_reservoir())
            .collect::<Vec<Vec<f64>>>();
        let ad = fits.iter().map(|x| x[0]).collect::<Vec<f64>>();
        let ch = fits.iter().map(|x| x[1]).collect::<Vec<f64>>();
        let kp = fits.iter().map(|x| x[2]).collect::<Vec<f64>>();
        let ks = fits.iter().map(|x| x[3]).collect::<Vec<f64>>();
        let ns = fits.iter().map(|x| x[4]).collect::<Vec<f64>>();

        let mut gofs = Vec::with_capacity(self.model.batch);
        for i in 0..self.model.batch {
            gofs.push(ReservoirFit {
                rate: rates[i],
                period: periods[i],
                ad: ad[i],
                ch: ch[i],
                kp: kp[i],
                ks: ks[i],
                n: ns[i],
            });
        }

        Ok(gofs)
    }

    /// Timed fit of reservoirs to observed record.
    pub fn fit_reservoirs_timed(
        mut self,
        path: &str,
    ) -> Result<Vec<ReservoirFit>, errors::ResError> {
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);

        let dur = std::time::Duration::new(60 * 60 * self.model.duration, 0);
        let now = std::time::SystemTime::now();
        let mut rec = Vec::new();
        let exists = std::path::Path::new(path).exists();
        match exists {
            true => {
                let mut stats: Vec<ReservoirFit> = ReservoirFit::read(path)?;
                rec.append(&mut stats);
            }
            false => {}
        }
        while std::time::SystemTime::now() < now + dur {
            let model = self
                .model
                .clone()
                .range(seeder.sample(&mut self.model.range));
            let fit = self.clone().model(&model);
            let mut new = fit.fit_reservoirs()?;
            {
                rec.append(&mut new);
            }
            ReservoirFit::record(&mut rec, path)?;
        }
        Ok(rec)
    }

    /// Goodness-of-fit test suite.
    pub fn gof(self, other: &[f64]) -> Vec<f64> {
        utils::gof(&self.mass, other)
    }

    /// Compare the accumulated mass in a reservoir to another record.
    /// Produces two goodness-of-fit statistics in a tuple:
    /// the K-S statistic and the Kuiper statistic, respectively.
    /// Called by [fit_rate](#method.fit_rate), you can use it on individual records too.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let ages = vec![10.0, 42.0, 77.7, 12.0, 99.9, 10000.0, 777.7];
    ///     let mut debris_flows = Reservoir::new()
    ///         .input(&0.78)?
    ///         .output(&0.78)?
    ///         .sim();
    ///     let fit = debris_flows.gof(&ages);
    ///     println!("Fit is {:?}.", fit);
    ///     Ok(())
    /// }
    /// ```
    pub fn gof1(&self, other: &[f64]) -> Fit {
        let cdf = utils::cdf_dual(&self.mass, other);
        let ad = utils::ad_dual(&self.mass, other);
        // anderson-darling test
        let lnx = self.mass.len();
        let lny = other.len();
        let k = lnx + lny;
        let ada = (2.492 - 1.0) * (1.0 - (1.55 / k as f64)) + 1.0;
        // chi-squared pearsons test
        let chs = cdf
            .iter()
            .filter(|x| x.0 > 0.0)
            .map(|x| f64::powi(x.1 - x.0, 2) / x.0)
            .sum::<f64>();
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
        Fit {
            ad,
            ada,
            chs,
            kp,
            ks,
            ksa,
        }
    }

    /// Inherited age refers to age of charcoal upon entering the reservoir.
    /// Multiple samples of charcoal from a single deposit produces a vector of inherited ages,
    /// represented by the mean expected age of each charcoal sample in a f64 vector.
    /// The sample age of charcoal is the sum of its inherited age plus transit time through the reservoir.
    /// When simulating a reservoir model, each event entering the reservoir receives
    /// a random amount of inherited age sampled from the vector `ages`.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let start_ages = vec![7.0, 42.0, 401.0, 1234.5, 7777.7, 5.2, 0.1];
    ///     let res = Reservoir::new().inherit(&start_ages);
    ///     Ok(())
    /// }
    /// ```
    pub fn inherit(mut self, ages: &[f64]) -> Self {
        self.inherit = Some(ages.to_vec());
        self
    }

    /// Assign an input rate to a reservoir.
    /// Converts a reference to a float 64 `rate` into an exponential distribution with lamdba `rate` using the rand crate.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     // new reservoirs have no input rate, call input() to set one
    ///     let mut res = Reservoir::new().input(&0.58)?;
    ///     Ok(())
    /// }
    /// ```
    pub fn input(mut self, rate: &f64) -> Result<Self, errors::ResError> {
        let rate = Exp::new(*rate)?;
        self.input = Some(rate);
        Ok(self)
    }

    /// Set model parameters as fields in manager struct.
    pub fn model(mut self, model: &ModelManager) -> Self {
        self.model = model.to_owned();
        self
    }

    /// Create reservoirs using a builder pattern.  Calling new() creates an empty reservoir.
    /// Use the [input](#method.input) and [output](#method.output) methods to set rates, which start at None.
    /// Set inherited age similarly using the method [inherit](#method.inherit).
    ///
    /// # Examples
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     let mut res = Reservoir::new();
    ///     Ok(())
    /// }
    /// ```
    pub fn new() -> Self {
        Reservoir {
            input: None,
            model: ModelManager::new(),
            mass: Vec::new(),
            output: None,
            flux: Vec::new(),
            inherit: None,
            range: rand::rngs::StdRng::from_entropy(),
        }
    }

    /// Runs sim() a number of times specified by ModelManager struct.
    pub fn n_sim(mut self) -> Vec<f64> {
        let mut rec = Vec::new();
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);

        for _ in 0..self.model.source_runs {
            let mut fit = self.clone();
            fit.model = fit
                .model
                .range(seeder.sample(&mut self.model.range))
                .clone();
            let mut new = fit.sim().mass;
            {
                rec.append(&mut new);
            }
        }
        rec
    }

    /// Assign an output rate to a reservoir.
    /// Converts a reference to a float 64 `rate` into an exponential distribution with lamdba `rate` using the rand crate.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     // new reservoirs have no output rate, call input() to set one
    ///     let mut res = Reservoir::new().output(&0.58)?;
    ///     Ok(())
    /// }
    /// ```
    pub fn output(mut self, rate: &f64) -> Result<Self, errors::ResError> {
        let rate = Exp::new(*rate)?;
        self.output = Some(rate);
        Ok(self)
    }

    /// Assign a seed to the random number generator, for reproducibility.
    ///  - `seed` is a number to convert in an RNG seed using the rand crate.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main() -> Result<(), ResError> {
    ///     // new reservoirs have a random seed, call range() to set a specific seed
    ///     let mut res = Reservoir::new().range(10101);
    ///     Ok(())
    /// }
    /// ```
    pub fn range(mut self, seed: u64) -> Self {
        self.range = rand::SeedableRng::seed_from_u64(seed);
        self
    }

    /// Workhorse function for simulating accumulation records in a reservoir.
    /// Runs simulations on reservoir objects created using the builder pattern.
    /// `period` specifies the amount of time to simulate accumulation in years.
    /// While generally this is a function called in series by other functions, you can use
    /// it to simulate a single accumulation record for a reservoir.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// fn main () -> Result<(), ResError> {
    ///     // create reservoirs
    ///     let mut fines = Reservoir::new().input(&0.75)?.output(&0.75)?;
    ///     let mut gravels = Reservoir::new().input(&0.54)?.output(&0.54)?;
    ///
    ///     // simulate accumulation for 30000 years
    ///     fines = fines.sim();
    ///     gravels = gravels.sim();
    ///     Ok(())
    /// }
    /// ```
    pub fn sim(mut self) -> Self {
        // let _ = pretty_env_logger::try_init();
        let mut om = 0f64;
        let mut im = 0f64;
        let mut flux = Vec::new();
        let mut mass = Vec::new(); // time of arrivals in reservoir

        while om < self.model.period {
            info!("Generating a time for removal.");
            match self.output {
                Some(x) => om += x.sample(&mut self.model.range) as f64,
                // TODO: Implement zero rates for input and output
                None => continue,
            }

            while im < om {
                info!("Generating inputs until time for removal.");
                if let Some(x) = self.input {
                    im += x.sample(&mut self.model.range) as f64;
                    mass.push(im);
                }
            }

            if !mass.is_empty() {
                info!("Inputs are available for removal.");
                let mvec: Vec<f64> = mass.par_iter().cloned().filter(|x| x <= &om).collect();
                info!("Selecting from inputs present before time of removal.");
                if !mvec.is_empty() {
                    let rm = rand::distributions::Uniform::from(0..mvec.len())
                        .sample(&mut self.model.range);
                    flux.push(mass[rm]);
                    mass.remove(rm);
                }
            }
        }

        mass = mass.iter().map(|x| self.model.period - x).collect();
        // flux = flux.par_iter().map(|x| period - x).collect();
        if let Some(x) = self.inherit.clone() {
            let ln = x.len();
            mass = mass
                .iter()
                .map(|z| {
                    z + x[rand::distributions::Uniform::from(0..ln).sample(&mut self.model.range)]
                })
                .collect();
            // flux = flux.iter().map(|z| z + x[rand::distributions::Uniform::from(0..ln).sample(&mut self.range)]).collect();
        }

        self.mass = mass;
        self.flux = flux;
        self
    }

    /// Return a stereotypical mass age distribution at the current model parameters.
    pub fn stereotype(&mut self) -> Vec<f64> {
        let mut res: Vec<Reservoir> = Vec::with_capacity(self.model.runs);
        let seeder: rand::distributions::Uniform<u64> =
            rand::distributions::Uniform::new(0, 10000000);
        let seeds = seeder
            .sample_iter(&mut self.model.range)
            .take(self.model.runs)
            .collect::<Vec<u64>>();

        for seed in seeds {
            // make boot number copies of reservoir
            res.push(self.clone().model(&self.model.clone().range(seed)));
        }
        res = res
            .iter()
            .cloned()
            .map(|x| x.sim())
            .collect::<Vec<Reservoir>>(); // simulate accumulation record for each copy
        let mut ns = res
            .iter()
            .cloned()
            .map(|x| x.mass.len() as f64)
            .collect::<Vec<f64>>(); // number of deposits in reservoir
        let mid_n = utils::median(&ns); // median number of deposits
        ns = ns
            .iter()
            .map(|x| rand_distr::num_traits::abs((x / mid_n) - 1.0))
            .collect(); // distance from median length
                        // collect reservoir masses into single vector and calculate the cdf
        let mut rec = Vec::new(); // vector of mass
        for r in res.clone() {
            rec.append(&mut r.mass.clone()); // add each run to make one long vector
        }
        // let cdf = utils::cdf_bin(&rec, bins); // subsample vector to length bins

        let gof = res
            .iter()
            .cloned()
            .map(|x| x.gof1(&rec))
            .collect::<Vec<Fit>>(); // ks and kp values
        let ks = gof.iter().cloned().map(|x| x.ks).collect::<Vec<f64>>(); // clip to just ks values
        let mut least = 1.0; // test for lowest fit (set to high value)
        let mut low = Reservoir::new(); // initialize variable to hold lowest fit
        for (i, val) in ns.iter().enumerate() {
            let loss = ks[i] + val; // loss function
            if loss < least {
                // if lowest value
                low = res[i].clone(); // copy to low
                least = loss; // set least to new low value
            }
        }

        low.mass
    }

    /// Returns raw transit times based on model parameters.
    pub fn bulk_transits(self) -> Vec<f64> {
        let res = self
            .clone()
            .model
            .seed_clones()
            .iter()
            .map(|x| self.clone().model(x).sim().mass)
            .collect::<Vec<Vec<f64>>>();
        let mut mass = Vec::new();
        for mut r in res {
            mass.append(&mut r);
        }
        mass
    }

    /// Return CDF of transit times based on model parameters.
    pub fn transit_times(self) -> Vec<f64> {
        let res = self
            .clone()
            .model
            .seed_clones()
            .iter()
            .map(|x| self.clone().model(x).sim())
            .collect::<Vec<Reservoir>>();

        let index = self.model.index;
        let mut out = vec![Complex::new(0.0, 0.0); index.len()];
        let mut real_planner = RealFftPlanner::<f64>::new();
        info!("Constructing FftPlanner for fourier transforms.");
        let r2c = real_planner.plan_fft_forward(index.len() + 1);
        for r in res {
            let cdf = utils::cdf_rng(&r.mass, &index);
            let mut pmf = utils::pmf_from_cdf(&cdf);
            let mut spectrum = r2c.make_output_vec();
            info!("Forward transforming the pmf using fft.");
            r2c.process(&mut pmf, &mut spectrum).unwrap();
            info!("Summing transformed pmfs.");
            out = out
                .iter()
                .zip(spectrum.iter())
                .map(|(a, b)| a + b)
                .collect::<Vec<Complex<f64>>>();
        }
        info!("Constructing FftPlanner for inverse fourier transforms.");
        let c2r = real_planner.plan_fft_inverse(index.len() + 1);
        info!("Inverse fourier transform using fft.");
        let mut out_data = c2r.make_output_vec();
        c2r.process(&mut out, &mut out_data).unwrap();
        info!("Normalize output by dividing by length.");
        out_data = out_data
            .par_iter()
            .map(|a| a / index.len() as f64)
            .collect::<Vec<f64>>();
        info!("Normalize output by dividing by sum of pmfs (having added several pmfs together).");
        let sum_out = out_data.iter().fold(0.0, |acc, x| acc + *x);
        out_data = out_data
            .par_iter()
            .map(|a| a / sum_out)
            .collect::<Vec<f64>>();

        out_data
    }
}

impl Default for Reservoir {
    fn default() -> Self {
        Self::new()
    }
}

/// Holder struct to read in charcoal sample ages from csv.
#[derive(Debug, Deserialize)]
pub struct Sample {
    id: String,
    /// Mean estimated charcoal age in years before 2000.
    pub age: f64,
    /// Deposit type - either debris flows (DF), fluvial fines (FF) or fluvial gravels (FG).
    pub facies: String,
}

impl Sample {
    /// Converts a csv file of charcoal ages into a Sample struct.
    pub fn read(path: &str) -> Result<Vec<Sample>, errors::ResError> {
        let mut record = Vec::new();
        let var = std::fs::File::open(path)?;
        let mut rdr = csv::Reader::from_reader(var);
        for result in rdr.records() {
            let row = result?;
            let row: Sample = row.deserialize(None)?;
            record.push(row);
        }
        Ok(record)
    }
}

/// Thresholds for statistical tests.
#[derive(Clone, Debug)]
pub struct Thresholds {
    ad: f64,
    ch: f64,
    kp: f64,
    ks: f64,
}

impl Thresholds {
    /// Builder for threshold struct.
    pub fn new(ad: f64, ch: f64, kp: f64, ks: f64) -> Thresholds {
        Thresholds { ad, ch, kp, ks }
    }
}

/// Struct to hold statistics describing transit times of reservoirs.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Transits {
    rate: f64,
    mean: f64,
}

impl Transits {
    /// Convert csv record to Transits struct.
    pub fn read(path: &str) -> Result<Vec<Transits>, errors::ResError> {
        let mut rec = Vec::new();
        let var = std::fs::File::open(path)?;
        let mut rdr = csv::Reader::from_reader(var);
        for result in rdr.records() {
            let row = result?;
            let row: Transits = row.deserialize(None)?;
            rec.push(row);
        }
        Ok(rec)
    }
}

/// The Record trait indicates the data is organized in spreadsheet format with
/// variables by column and observations by row (using struct fields as column names).
/// Trait methods provides convenience functions for wrangling spreadsheet data.
pub trait Record {
    /// Sorts data and divides into bins, returning averages by bin.
    fn bin_ave(&self, bins: usize) -> (Vec<f64>, Vec<f64>);
}

impl Record for Vec<Gof> {
    /// Orders of vector of Gof structs by input rate,
    /// breaks into chunks based on the number of bins desired `bins`,
    /// returns the mean input rate and ks fit per bin.
    fn bin_ave(&self, bins: usize) -> (Vec<f64>, Vec<f64>) {
        let ord = &mut (*self).clone();
        ord.sort_by(|x, y| x.input.partial_cmp(&y.input).unwrap());
        let chunk_ln = (ord.len() as f64 / bins as f64).round() as usize;
        let rates: Vec<f64> = ord.iter().map(|x| x.input).collect();
        let rates: Vec<f64> = rates.chunks(chunk_ln).map(|x| utils::mean(x)).collect();
        let fits: Vec<f64> = ord.iter().map(|x| x.ks).collect();
        let fits: Vec<f64> = fits.chunks(chunk_ln).map(|x| utils::mean(x)).collect();
        (rates, fits)
    }
}
