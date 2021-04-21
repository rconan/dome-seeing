use complot as plt;
use dome_seeing::{DomeSeeingError, Pupil, TemperatureField};
use std::env;

fn minmax(data: &[f64]) -> (f64, f64) {
    let max = data
        .iter()
        .filter(|x| !x.is_nan())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let min = data
        .iter()
        .filter(|x| !x.is_nan())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    (min, max)
}

fn main() -> std::result::Result<(), DomeSeeingError> {
    let filename = "data/b2019_30z_0az_os_7ms/OPDData_OPD_Data_9.000000e+02.csv";
    let temp_field =
        TemperatureField::load(filename).map_err(|e| DomeSeeingError { source: Some(e) })?;

    let mut pupil = Pupil::default().sample();

    let temperature = temp_field.z_slice(18f64, &mut pupil);

    println!("Temperature min/max : {:7.3?}K", minmax(&temperature));

    env::var("PLOT")
        .and({
            let filename = format!("slice_{}.png", temp_field.method);
            let plot = plt::png_canvas(&filename);
            plt::imagesc(&temperature, &plot);
            Ok(())
        })
        .map_err(|e| DomeSeeingError {
            source: Some(e.into()),
        })?;

    Ok(())
}
