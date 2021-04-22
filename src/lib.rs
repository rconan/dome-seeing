use std::env;

pub mod bounding_box;
pub mod pupil;
pub mod temperature;
pub mod ray_tracing;

pub use bounding_box::{BoundingBox, Projection};
pub use pupil::{Pupil, PupilInner};
pub use temperature::TemperatureField;

pub enum InterpolationMethod {
    NearestNeighbor,
    RadialBasis,
}
impl InterpolationMethod {
    pub fn from_env() -> Result<Self> {
        env::var("METHOD").map_or(Ok(Self::NearestNeighbor), |e| match e.as_str() {
            "NEAREST_NEIGHBOR" => Ok(Self::NearestNeighbor),
            "RADIAL_BASIS" => Ok(Self::RadialBasis),
            _ => Err("interpolation method is either NEAREST_NEIGHBOR or RADIAL_BASIS".into()),
        })
    }
}
impl std::fmt::Display for InterpolationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use InterpolationMethod::*;
        match self {
            NearestNeighbor => f.write_str("nearest_neighbor")?,
            RadialBasis => f.write_str("radial_basis")?,
        };
        Ok(())
    }
}
impl From<InterpolationMethod> for String {
    fn from(method: InterpolationMethod) -> Self {
        use InterpolationMethod::*;
        match method {
            NearestNeighbor => String::from("nearest_neighbor"),
            RadialBasis => String::from("radial_basis"),
        }
    }
}
impl From<InterpolationMethod> for &str {
    fn from(method: InterpolationMethod) -> Self {
        use InterpolationMethod::*;
        match method {
            NearestNeighbor => "nearest_neighbor",
            RadialBasis => "radial_basis",
        }
    }
}

pub fn refraction_index(temp: f64) -> f64 {
    let p_ref =  75000.0; // Reference pressure
    let wlm = 0.5;
    7.76e-7*p_ref*(1.+0.00752/(wlm*wlm))/temp
}

pub type BoxError = std::boxed::Box<
    dyn std::error::Error // must implement Error to satisfy ?
        + std::marker::Send // needed for threads
        + std::marker::Sync, // needed for threads
>;
type Result<T> = std::result::Result<T, BoxError>;
pub struct DomeSeeingError {
    pub source: Option<BoxError>,
}
impl std::fmt::Display for DomeSeeingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("An error occurred while processing dome seeing data")?;
        if let Some(error) = &self.source {
            write!(f, "\nCaused by: {}", error)?;
        }
        Ok(())
    }
}
impl std::fmt::Debug for DomeSeeingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <DomeSeeingError as std::fmt::Display>::fmt(self, f)
    }
}
impl std::error::Error for DomeSeeingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.source {
            Some(error) => Some(error.as_ref()),
            None => None,
        }
    }
}
