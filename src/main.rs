use complot as plt;
use dome_seeing::{
    refraction_index, BoundingBox, DomeSeeingError, InterpolationMethod, Projection, Pupil,
    TemperatureField,
};
use gmt_kpp::KPP;
use rayon::prelude::*;
use std::time::Instant;
//use std::env;
use indicatif::{ParallelProgressIterator, ProgressBar};

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

    let bbox = BoundingBox::new([0., 0., 2.], 30., 30., 4.);
    let aabb = bbox.aabb();
    println!("Lower: {:?}", aabb.lower());
    println!("Upper: {:?}", aabb.upper());

    let bbox_data = bbox.locate(&temp_field.scattered_data);
    println!("BBox data count: {}", bbox_data.count());

    bbox.show(&temp_field.scattered_data, Some(Projection::SideY), None);

    let now = Instant::now();
    let z_start = 50.;
    let z_end = 4.;
    let dz = 0.5;
    let nz = ((z_start - z_end) / dz) as usize + 1;
    let bbox_side = 30.;
    let pupil = Pupil::default().sample();
    println!("Interpolation ...");
    let interpolated_data: Vec<_> = (0..nz)
        .into_par_iter()
        .progress_with(ProgressBar::new(nz as u64))
        .map(|k| {
            let z = z_start - k as f64 * dz;
            let mut points = pupil.points.clone();
            points.iter_mut().for_each(|p| {
                p[2] = z;
            });
            let bbox = BoundingBox::new([0., 0., z], bbox_side, bbox_side, 0.05);
            bbox.interp(
                &points,
                &temp_field.scattered_data,
                InterpolationMethod::from_env().unwrap(),
            )
            .iter()
            .map(|x| refraction_index(*x))
            .collect::<Vec<f64>>()
        })
        .collect();
    println!("OPL integration ...");
    let mut opl_pupil: Vec<_> = (0..pupil.nnz())
        .into_par_iter()
        .map(|p| {
            (0..nz - 1)
                .into_iter()
                .map(|s| (interpolated_data[s][p] + interpolated_data[s + 1][p]) * 0.5 * dz)
                .sum::<f64>()
        })
        .collect();
    let mean_opd: f64 = opl_pupil.iter().sum::<f64>() / opl_pupil.len() as f64;
    opl_pupil.iter_mut().for_each(|x| *x -= mean_opd);
    let std_opd: f64 =
        (opl_pupil.iter().map(|x| x * x).sum::<f64>() / opl_pupil.len() as f64).sqrt();
    println!("OPD STD: {:.0}nm", 1e9 * std_opd);
    //println!("OPD min/max : {:?}micron", minmax(&opl_pupil));
    let n = pupil.n;
    let mut opl_map: Vec<f64> = vec![0f64; n * n];
    pupil
        .index
        .iter()
        .zip(opl_pupil.iter())
        .for_each(|(k, vals)| {
            opl_map[*k] = *vals;
        });
    println!("OPD min/max : {:7.3?}micron", {
        let v = minmax(&opl_map);
        (v.0 * 1e6, v.1 * 1e6)
    });
    println!("OPD computed in {:#}ms", now.elapsed().as_millis());

    let p = pupil.mask();
    let mut pssn = KPP::new()
        .wavelength(500e-9)
        .pssn(pupil.length, pupil.n, &p);
    let pssn_val = pssn.estimate(&p, Some(&opl_map));
    println!("V PSSn: {:.4}", pssn_val);

    //    let filename = format!("s.png", temp_field.method);
    let plot = plt::png_canvas("opl_map.png");
    plt::imagesc(&opl_map, &plot);

    /*let n = pupil.n;
    let temperature_slices: Vec<_> = interpolated_data
        .into_iter()
        .map(|values| {
            let mut temperature: Vec<f64> = vec![f64::NAN; n * n];
            pupil.index.iter().zip(values.iter()).for_each(|(k, vals)| {
                temperature[*k] = *vals;
            });
            temperature
        })
        .collect();

    for (k, temperature) in temperature_slices.iter().enumerate() {
        let filename = format!("slice_{}.png", k);
        let plot = plt::png_canvas(&filename);
        plt::imagesc(&temperature, &plot);
    }*/

    /*
    let mut pupil = Pupil::default().sample();
    let temperature = temp_field.z_slice(18f64, &mut pupil);

    println!("Temperature min/max : {:7.3?}K", minmax(&temperature));

    match  env::var("PLOT") {
        Ok(_) => {
            let filename = format!("slice_{}.png", temp_field.method);
            let plot = plt::png_canvas(&filename);
            plt::imagesc(&temperature, &plot);
        },
        Err(_) => ()
    };
     */
    Ok(())
}
