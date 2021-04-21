pub mod pupil;
pub mod temperature;
pub use pupil::{Pupil, PupilInner};
pub use temperature::{InterpolationMethod, TemperatureField};

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

