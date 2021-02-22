//! Structs and methods for Bolin & Rodhes reservoir models.
use crate::utils;
use rand_distr::{Distribution, Exp};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Holder struct for goodness-of-fit statistics.
#[derive(Debug, Deserialize, Serialize)]
pub struct Gof {
    input: f64,
    output: f64,
    ks: f64,
    kp: f64,
    n: f64,
}

impl Gof {
    /// Convert csv record to Gof struct.
    pub fn read(path: &str) -> Result<Vec<Gof>, ResError> {
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
    pub fn record(rec: &mut Vec<Gof>, title: &str) -> Result<(), ResError> {
        let mut wtr = csv::Writer::from_path(title)?;
        for i in rec {
            wtr.serialize(i)?;
        }
        wtr.flush()?;
        Ok(())
    }
}

/// Holds model characteristics associated with a reservoir.
#[derive(Debug)]
pub struct Model {
    reservoir: Reservoir,
    period: f64,
    runs: usize,
    batch: usize,
    duration: u64,
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
    /// ```rust
    /// use reservoirs::prelude::*;
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
    /// model.fit_range(0.01..1.0, 0.01..1.0, &df, "examples/df_fit_1k.csv");
    /// ```
    pub fn fit_range(
        &self,
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
            let mut new = self.reservoir.fit_rng(
                &self.period,
                self.runs,
                self.batch,
                input.clone(),
                output.clone(),
                obs,
            );
            {
                rec.append(&mut new);
            }
            Gof::record(&mut rec, title).unwrap();
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
    /// ```rust
    /// use reservoirs::prelude::*;
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
    /// ```
    pub fn fit_steady(&self, rate: std::ops::Range<f64>, obs: &[f64], title: &str) {
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
            let mut new =
                self.reservoir
                    .steady(&self.period, self.runs, self.batch, rate.clone(), obs);
            {
                rec.append(&mut new);
            }
            Gof::record(&mut rec, title).unwrap();
        }
    }

    /// Creates a new model from a given `reservoir` with default `period`, `runs` and `batch` values.
    pub fn new(reservoir: Reservoir) -> Self {
        Model {
            reservoir,
            period: 30000.0,
            runs: 100,
            batch: 10,
            duration: 1,
        }
    }

    /// Sets the model period to desired time in years.
    pub fn period(mut self, years: &f64) -> Self {
        self.period = years.to_owned();
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

    /// Produces a raw vector of observations of simulated transit times from multiple runs.
    /// Used to produce statistics on the mean, median and statistical moments.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    ///
    /// let period: f64 = 30000.0;  // length of time to simulate accumulation in the reservoir in years
    /// let boot: usize = 1000; // number of simulations to run for estimation
    /// let df_rt: f64 = 0.78; // rate for steady state debris-flow accumulation
    /// let ff_rt: f64 = 0.63; // rate for steady state fluvial fines accumulation
    /// let fg_rt: f64 = 0.58; // rate for steady state fluvial gravels accumulation
    ///
    /// // to estimate transit times, omit inherit age from charcoal
    /// let debris_flows = Reservoir::new().input(&df_rt)?.output(&df_rt)?;
    /// let fines = Reservoir::new().input(&ff_rt)?.output(&ff_rt)?;
    /// let gravels = Reservoir::new().input(&fg_rt)?.output(&fg_rt)?;
    ///
    /// let df_mod = Model::new(debris_flows).period(&30000.0).runs(100);
    /// let ff_mod = Model::new(fines).period(&30000.0).runs(100);
    /// let fg_mod = Model::new(gravels).period(&30000.0).runs(100);
    ///
    /// // vector of transit times
    /// let df_t = df_mod.transit_times();
    /// let ff_t = ff_mod.transit_times();
    /// let fg_t = fg_mod.transit_times();
    ///
    /// println!("Debris-flow transit time quantiles are {:?}", utils::quantiles(&df_t));
    /// println!("Fluvial fines transit time quantiles are {:?}", utils::quantiles(&ff_t));
    /// println!("Fluvial gravels transit time quantiles are {:?}", utils::quantiles(&fg_t));
    ///
    /// ```
    pub fn transit_times(&self) -> Vec<f64> {
        let mut res = Vec::with_capacity(self.runs);
        for _ in 0..self.runs {
            res.push(self.reservoir.clone());
        }
        res = res
            .par_iter()
            .cloned()
            .map(|x| x.sim(&self.period).unwrap())
            .collect();
        let mut rec = Vec::new();
        for r in res {
            let mut mass = r.mass.clone();
            rec.append(&mut mass);
        }
        rec
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
        }
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

/// Struct for recording reservoir characteristics.
#[derive(Debug, Clone)]
pub struct Reservoir {
    input: Option<Exp<f64>>,
    /// Ages (in years) of deposits accumulated in reservoir.
    pub mass: Vec<f64>,
    output: Option<Exp<f64>>,
    flux: Vec<f64>,
    inherit: Option<Vec<f64>>,
}

/// Custom error type for the reservoirs crate.
#[derive(Debug)]
pub enum ResError {
    /// Error type from csv crate.
    CsvError,
    /// Error type from rand crate.
    ExpError,
    /// Error type from std::io.
    IoError,
}

impl std::error::Error for ResError {}

impl std::fmt::Display for ResError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ResError::CsvError => write!(f, "Could not serialize/deserialize csv file."),
            ResError::ExpError => write!(
                f,
                "Could not create exponential distribution from rate provided."
            ),
            ResError::IoError => write!(f, "Could not read file from path provided."),
        }
    }
}

impl From<csv::Error> for ResError {
    fn from(_: csv::Error) -> Self {
        ResError::CsvError
    }
}

impl From<rand_distr::ExpError> for ResError {
    fn from(_: rand_distr::ExpError) -> Self {
        ResError::ExpError
    }
}

impl From<std::io::Error> for ResError {
    fn from(_: std::io::Error) -> Self {
        ResError::IoError
    }
}

impl Reservoir {
    /// Randomly selects rate pairs from ranges `input` and `output`, and simulates `boot` number of accumulation records
    /// in batches of `bat` using [fit_rate](#method.fit_rate).  Returns the selected input/output pair and the mean
    /// goodness-of-fit statistics for each pair from `boot` simulations.  Called by [fit_range](struct.Model.html#method.fit_range).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter().filter(|x| x.facies == "DF").map(|x| x.age).collect();
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter().map(|x| x.age).collect();
    ///
    /// let mut debris_flows = Reservoir::new().input(&0.78)?.output(&0.78)?.inherit(&ia);
    /// // fit 10 randomly selected rate pairs (from range 0.01 to 1.0) to observed debris flows
    /// // by running 1000 simulations for 30000 years for each pair
    /// let gofs = debris_flows.fit_rng(&30000.0, 1000, 10, 0.01..1.0, 0.01..1.0, &df);
    /// ```
    pub fn fit_rng(
        &self,
        period: &f64,
        boot: usize,
        bat: usize,
        input: std::ops::Range<f64>,
        output: std::ops::Range<f64>,
        obs: &[f64],
    ) -> Vec<Gof> {
        let mut roll = rand::thread_rng();
        let mut inputs = Vec::with_capacity(bat);
        let mut outputs = Vec::with_capacity(bat);
        let mut fits = Vec::with_capacity(bat);
        for i in 0..bat {
            inputs.push(rand::distributions::Uniform::from(input.clone()).sample(&mut roll));
            outputs.push(
                rand::distributions::Uniform::from(
                    output.clone().start.max(inputs[i] * 0.975)
                        ..output.clone().end.min(inputs[i] * 1.0125),
                )
                .sample(&mut roll),
            );
            fits.push(
                self.clone()
                    .input(&inputs[i])
                    .unwrap()
                    .output(&outputs[i])
                    .unwrap(),
            );
        }
        let gof: Vec<(f64, f64, f64)> = fits
            .par_iter()
            .map(|x| x.fit_rate(period, &obs, boot))
            .collect();
        let mut gofs = Vec::with_capacity(bat);
        for i in 0..bat {
            gofs.push(Gof {
                input: inputs[i],
                output: outputs[i],
                ks: gof[i].0,
                kp: gof[i].1,
                n: gof[i].2,
            })
        }
        gofs
    }

    /// Runs `boot` number of simulations of length `period` on a reservoir.
    /// Returns the mean goodness-of-fit statistics compared to accumulation record `other`.
    /// Called by [fit_rng](#method.fit_rng) and [steady](#method.steady).  To use,
    /// set characteristics of the reservoir before running.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter().filter(|x| x.facies == "DF").map(|x| x.age).collect();
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter().map(|x| x.age).collect();
    ///
    /// let mut debris_flows = Reservoir::new().input(&0.78)?.output(&0.78)?.inherit(&ia);
    /// // run 1000 simulations for 30000 years and compare the fit against observed debris flows
    /// let (ks, kp, _) = debris_flows.fit_rate(&30000.0, &df, 1000);
    /// println!("K-S fit is {}.", ks);
    /// println!("Kuiper fit is {}.", kp);
    ///
    /// ```
    pub fn fit_rate(&self, period: &f64, other: &[f64], boot: usize) -> (f64, f64, f64) {
        let mut res: Vec<Reservoir> = Vec::with_capacity(boot);
        for _ in 0..boot {
            res.push(self.clone());
        }
        res = res
            .par_iter()
            .cloned()
            .map(|x| x.sim(period).unwrap())
            .collect();
        let fits: Vec<(f64, f64)> = res.par_iter().cloned().map(|x| x.gof(other)).collect();
        let kss: Vec<f64> = fits.par_iter().map(|x| x.0).collect();
        let kps: Vec<f64> = fits.par_iter().map(|x| x.1).collect();
        let ns: Vec<f64> = res
            .iter()
            .map(|x| x.mass.len() as f64 / other.len() as f64)
            .collect();

        (utils::mean(&kss), utils::median(&kps), utils::mean(&ns))
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
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter().filter(|x| x.facies == "DF").map(|x| x.age).collect();
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter().map(|x| x.age).collect();
    ///
    /// let mut debris_flows = Reservoir::new().input(&0.78)?.output(&0.78)?.inherit(&ia);
    /// debris_flows = debris_flows.sim(&30000.0)?;
    /// let (ks, kp) = debris_flows.gof(&df);
    /// println!("K-S fit is {}.", ks);
    /// println!("Kuiper fit is {}.", kp);
    ///
    /// ```
    pub fn gof(&self, other: &[f64]) -> (f64, f64) {
        let mut x = self.mass.clone();
        let mut y = other.to_vec();
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
        let ks = cdf
            .iter()
            .map(|x| rand_distr::num_traits::abs(x.0 - x.1))
            .fold(0.0, f64::max);
        let kp1 = cdf.iter().map(|x| x.0 - x.1).fold(0.0, f64::max);
        let kp2 = cdf.iter().map(|x| x.1 - x.0).fold(0.0, f64::max);
        let kp = kp1 + kp2;
        (ks, kp)
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
    /// // mean expected inherited age by facies
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter().map(|x| x.age).collect();
    ///
    /// let res = Reservoir::new().inherit(&ia);
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
    /// res = Reservoir::new().input(&0.58)?;
    /// ```
    pub fn input(mut self, rate: &f64) -> Result<Self, ResError> {
        let rate = Exp::new(*rate)?;
        self.input = Some(rate);
        Ok(self)
    }

    /// Create reservoirs using a builder pattern.  Calling new() creates an empty reservoir.
    /// Use the [input](#method.input) and [output](#method.output) methods to set rates, which start at None.
    /// Set inherited age similarly using the method [inherit](#method.inherit).
    ///
    /// # Examples
    /// ```
    /// use reservoirs::prelude::*;
    /// let mut res = Reservoir::new();
    /// ```
    pub fn new() -> Self {
        Reservoir {
            input: None,
            mass: Vec::new(),
            output: None,
            flux: Vec::new(),
            inherit: None,
        }
    }

    /// Assign an output rate to a reservoir.
    /// Converts a reference to a float 64 `rate` into an exponential distribution with lamdba `rate` using the rand crate.
    ///
    /// # Examples
    ///
    /// ```
    /// use reservoirs::prelude::*;
    /// res = Reservoir::new().output(&0.58)?;
    /// ```
    pub fn output(mut self, rate: &f64) -> Result<Self, ResError> {
        let rate = Exp::new(*rate)?;
        self.output = Some(rate);
        Ok(self)
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
    ///
    /// // create reservoirs
    /// let mut fines = Reservoir::new().input(&0.75)?.output(&0.75)?;
    /// let mut gravels = Reservoir::new().input(&0.54)?.output(&0.54)?;
    ///
    /// // simulate accumulation for 30000 years
    /// fines = fines.sim(&30000.0)?;
    /// gravels = gravels.sim(&30000.0)?;
    ///
    /// ```
    pub fn sim(mut self, period: &f64) -> Result<Self, ResError> {
        let mut rng = rand::thread_rng();
        let mut om = 0f64;
        let mut im = 0f64;
        let mut mass = Vec::new(); // time of arrivals in reservoir
                                   // let mut flux = Vec::new();  //

        while om < *period {
            // Generate a time for removal
            match self.output {
                Some(x) => om += x.sample(&mut rand::thread_rng()) as f64,
                // TODO: Implement zero rates for input and output
                None => continue,
            }

            while im < om {
                // Generate inputs until time for removal
                if let Some(x) = self.input {
                    im += x.sample(&mut rand::thread_rng()) as f64;
                    mass.push(im);
                }
            }

            if !mass.is_empty() {
                // If there are inputs to remove
                let mvec: Vec<f64> = mass.par_iter().cloned().filter(|x| x <= &om).collect();
                // Only remove inputs younger than the output time
                if !mvec.is_empty() {
                    let rm = rand::distributions::Uniform::from(0..mvec.len()).sample(&mut rng);
                    // flux.push(mass[rm]);
                    mass.remove(rm);
                }
            }
        }

        mass = mass.par_iter().map(|x| period - x).collect();
        // flux = flux.par_iter().map(|x| period - x).collect();
        if let Some(x) = self.inherit.clone() {
            let ln = x.len();
            mass = mass
                .iter()
                .map(|z| z + x[rand::distributions::Uniform::from(0..ln).sample(&mut rng)])
                .collect();
            // flux = flux.iter().map(|z| z + x[Uniform::from(0..ln).sample(&mut rng)]).collect();
        }

        self.mass = mass;
        // self.flux = flux;
        Ok(self)
    }

    /// Selects from among a number of simulated accumulation records, the record most
    /// characteristic of the average number and age distribution of deposits.
    ///  - `self` is a reservoir with characteristics set.
    ///  - `boot` is the number of times to simulate an accumulation record from `self`
    ///  - `bins` is the number of bins with which to construct a cdf.
    ///  - Returns the record with the lowest loss statistic.
    ///
    /// For large numbers of `boot`, the cdf becomes a large vector.  Selecting a smaller
    /// number of `bins` reduces memory consumption and increases speed, but can reduce accuracy.
    ///
    /// The loss statistic is the K-S stat plus a normalized 'distance from average number
    /// of deposits'.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter()
    ///         .filter(|x| x.facies == "DF")
    ///         .map(|x| x.age)
    ///        .collect();
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter()
    ///         .map(|x| x.age)
    ///         .collect();
    ///
    /// // create steady state reservoir with charcoal inherited ages
    /// let res = Reservoir::new().input(&0.78)?
    ///         .output(&0.78)?
    ///         .inherit(&ia);
    /// // sample a stereotypical record from 1000 runs of 30000 years
    /// let eg = res.stereotype(&30000.0, 1000, 200);
    /// // compare the CDF of the synthetic example to the observed debris-flow deposit record
    /// plot::comp_cdf(&eg, &df, "examples/df_cdf.png");
    ///```
    pub fn stereotype(&self, period: &f64, boot: usize, bins: usize) -> Vec<f64> {
        let mut res: Vec<Reservoir> = Vec::with_capacity(boot);
        for _ in 0..boot {
            // make boot number copies of reservoir
            res.push(self.clone());
        }
        res = res
            .par_iter()
            .cloned()
            .map(|x| x.sim(period).unwrap())
            .collect(); // simulate accumulation record for each copy
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
        let gof: Vec<(f64, f64)> = res.par_iter().cloned().map(|x| x.gof(&cdf)).collect(); // ks and kp values
        let ks: Vec<f64> = gof.par_iter().cloned().map(|x| x.0).collect(); // clip to just ks values
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

    /// Randomly selects a rate from ranges `rate` for a steady state reservoir,
    /// and simulates `boot` number of accumulation records
    /// in batches of `bat` using [fit_rate](#method.fit_rate).  Returns the selected input/output pair and the mean
    /// goodness-of-fit statistics compared to `obs` for each pair from `boot` simulations.
    /// Called by [fit_steady](struct.Model.html#method.fit_steady).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use reservoirs::prelude::*;
    ///
    /// // mean expected deposit age and inherited age by facies
    /// let dep = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/dep.csv")?;
    /// let iat = Sample::read("https://github.com/crumplecup/reservoirs/blob/master/examples/iat.csv")?;
    ///
    /// // subset mean ages of debris flows
    /// let df: Vec<f64> = dep.iter().filter(|x| x.facies == "DF").map(|x| x.age).collect();
    /// // subset inherited ages
    /// let ia: Vec<f64> = iat.iter().map(|x| x.age).collect();
    ///
    /// let mut debris_flows = Reservoir::new().input(&0.78)?.output(&0.78)?.inherit(&ia);
    /// // fit 10 steady state reservoirs with randomly selected rates (from range 0.01 to 1.0) to observed debris flows
    /// // by running 1000 simulations for 30000 years for each pair
    /// let gofs = debris_flows.steady(&30000.0, 1000, 10, 0.01..1.0, &df);
    /// ```
    pub fn steady(
        &self,
        period: &f64,
        boot: usize,
        bat: usize,
        rate: std::ops::Range<f64>,
        obs: &[f64],
    ) -> Vec<Gof> {
        let mut roll = rand::thread_rng();
        let mut rates = Vec::with_capacity(bat);
        let mut fits = Vec::with_capacity(bat);
        for i in 0..bat {
            rates.push(rand::distributions::Uniform::from(rate.clone()).sample(&mut roll));
            fits.push(
                self.clone()
                    .input(&rates[i])
                    .unwrap()
                    .output(&rates[i])
                    .unwrap(),
            );
        }
        let gof: Vec<(f64, f64, f64)> = fits
            .par_iter()
            .map(|x| x.fit_rate(period, &obs, boot))
            .collect();
        let mut gofs = Vec::with_capacity(bat);
        for i in 0..bat {
            gofs.push(Gof {
                input: rates[i],
                output: rates[i],
                ks: gof[i].0,
                kp: gof[i].1,
                n: gof[i].2,
            })
        }
        gofs
    }
}

impl Default for Reservoir {
    fn default() -> Self {
        Self::new()
    }
}

impl Sample {
    /// Converts a csv file of charcoal ages into a Sample struct.
    pub fn read(path: &str) -> Result<Vec<Sample>, ResError> {
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
