/// Custom error type for the reservoirs crate.
#[derive(Debug)]
pub enum ResError {
    /// Error type from csv crate.
    CsvError,
    /// Error type from rand crate.
    ExpError,
    /// Error type from std::io.
    IoError,
    /// Error type from Box<dyn StdError>
    BoxError,
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
            ResError::BoxError => write!(f, "Maybe a plot error."),
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

impl From<Box<dyn serde::de::StdError>> for ResError {
    fn from(_: Box<dyn serde::de::StdError>) -> Self {
        ResError::BoxError
    }
}
