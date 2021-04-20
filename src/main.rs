use complot as plt;
use csv;
use geotrans;
use nalgebra::DVector;
use rbf_interp::{Basis, Scatter};
use rstar::{primitives::PointWithData, RTree, AABB};
use serde::Deserialize;
use std::{env, fs::File, io::BufReader, time::Instant};

#[derive(Deserialize)]
pub struct TemperatureData {
    #[serde(rename = "Temperature (K)")]
    temperature: f64,
    #[serde(rename = "X (m)")]
    x: f64,
    #[serde(rename = "Y (m)")]
    y: f64,
    #[serde(rename = "Z (m)")]
    z: f64,
}

type TemperatureField = PointWithData<f64, [f64; 3]>;

fn polywind(x: f64, y: f64, vx: &[f64], vy: &[f64]) -> i32 {
    let n = vx.len();
    let mut p0 = vx[n - 1];
    let mut p1 = vy[n - 1];
    let mut wind = 0i32;
    for i in 0..n {
        let d0 = vx[i];
        let d1 = vy[i];
        let q = (p0 - x) * (d1 - y) - (p1 - y) * (d0 - x);
        if p1 <= y {
            if d1 > y && q > 0.0 {
                wind += 1;
            }
        } else {
            if d1 <= y && q < 0.0 {
                wind -= 1;
            }
        }
        p0 = d0;
        p1 = d1;
    }
    wind
}
fn truss_shadow(x: f64, y: f64) -> bool {
    let vx = vec![
        -3.011774, -2.446105, -3.011774, -2.799304, -2.33903, -1.566412, -1.640648, -1.65,
        -1.640648, -1.566412, -2.347462, -1.597649, -1.725044, -2.392888, -2.799304,
    ];
    let vy = vec![
        -2.902158, 0., 2.902158, 3.107604, 0.07244, 0.518512, 0.175429, 0., -0.175429, -0.518512,
        -0.067572, -3.865336, -3.810188, -0.427592, -3.107604,
    ];
    if (1..4).fold(0, |a, k| {
        let q = geotrans::Quaternion::unit(-120f64.to_radians() * k as f64, geotrans::Vector::k());
        let (vx, vy): (Vec<f64>, Vec<f64>) = vx
            .iter()
            .cloned()
            .zip(vy.iter().cloned())
            .map(|(x, y)| {
                let v = geotrans::Vector::from([x, y, 0.0]);
                let p = geotrans::Quaternion::from(v);
                let pp = &q * p * q.complex_conjugate();
                let u = pp.vector_as_slice();
                (u[0], u[1])
            })
            .unzip();
        a + polywind(x, y, &vx, &vy)
    }) == 0
    {
        false
    } else {
        true
    }
}

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

const L: f64 = 8.71;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let now = Instant::now();
    println!("Loading data into R-Tree ...");
    let filename = "data/b2019_30z_0az_os_7ms/OPDData_OPD_Data_9.000000e+02.csv";
    let file = File::open(filename)?;
    let rdr = BufReader::with_capacity(100_000, file);
    let mut data = csv::Reader::from_reader(rdr);
    let mut tree = RTree::new();
    let mut z_max = f64::NEG_INFINITY;
    let mut z_min = f64::INFINITY;
    for result in data.deserialize() {
        let record: TemperatureData = result?;
        if record.z > z_max {
            z_max = record.z;
        }
        if record.z < z_min {
            z_min = record.z;
        }
        tree.insert(TemperatureField::new(
            record.temperature,
            [record.x, record.y, record.z],
        ));
    }
    println!(" ... in {}ms", now.elapsed().as_millis());
    println!("Tree size: {}", tree.size());
    println!("z minmax [{:8.3},{:8.3}]m", z_min, z_max);

    let z = 20f64;
    let length = 26f64;
    let d = 0.05; //ength / (n - 1) as f64;
    let n = (length / d) as usize + 1;
    println!(
        "Interpolating through slice (z={:.3}m,d={:.3}m,n={}) ...",
        z, d, n
    );

    // Pupil sampling
    let m1_radius = 8.365 * 0.5;
    let m2_radius = 3.3 * 0.5;
    let mut pupil_sample_point = vec![];
    let mut pupil_index = vec![];
    for i in 0..n {
        let x = i as f64 * d - 0.5 * length;
        for j in 0..n {
            let y = j as f64 * d - 0.5 * length;
            let mut sid = 1;
            loop {
                let u = geotrans::Vector::from([x, y, 0.]);
                let v = if sid < 7 {
                    let (s, c) = (90. - 60. * (sid - 1) as f64).to_radians().sin_cos();
                    let t = geotrans::Vector::from([L * c, L * s, 0.]);
                    u - t
                } else {
                    u
                };
                //println!("{:?}",v);
                let r = v[0].hypot(v[1]);
                if r <= m1_radius && !(sid == 7 && (truss_shadow(x, y) || r < m2_radius)) {
                    //println!("sid: {} -> [{},{}]", sid, i, j);
                    let point = [x, y, z];
                    let k = i * n + j;
                    pupil_sample_point.push(point);
                    pupil_index.push(k);
                    break;
                }
                sid += 1;
                if sid > 7 {
                    break;
                }
            }
        }
    }

    let interp_method = env::var("METHOD").unwrap_or("NEAREST_NEIGHBOR".to_string());
    println!("Interpolation method: {}", interp_method);

    // Interpolation
    let mut temperature: Vec<f64> = vec![f64::NAN; n * n];
    match interp_method.as_str() {
        "NEAREST_NEIGHBOR" => pupil_sample_point.iter().zip(pupil_index.iter()).for_each(|(point,k)| {
            temperature[*k] = tree.nearest_neighbor(point).unwrap().data;
        }),
        "RADIAL_BASIS" => {
            let now = Instant::now();
            let bbox_thickness = 0.1;
            let bbox_width = 30.;
            let b_h = bbox_width * 0.5;
            let b_d = bbox_thickness * 0.5;
            let bbox = AABB::from_corners([-b_h, -b_h, z - b_d], [b_h, b_h, z + b_d]);
            let (sampled_points, sampled_data): (Vec<DVector<f64>>, Vec<DVector<f64>>) = tree
                .locate_in_envelope(&bbox)
                .map(|p| {
                    (
                        DVector::from_column_slice(p.position()),
                        DVector::from_column_slice(&[p.data]),
                    )
                })
                .unzip();
            let n_in_aabb = sampled_points.len();
            println!(
                "# of points within bounding box: {} (in {}ms)",
                n_in_aabb,
                now.elapsed().as_millis()
            );

            println!("RBF creation ...");
            let now = Instant::now();
            let scatter = Scatter::create(sampled_points, sampled_data, Basis::PolyHarmonic(2), 2);
            println!(" ... in {}ms", now.elapsed().as_millis());

            println!("RBF interpolation ...");
            let now = Instant::now();
            let values = scatter.eval2(&pupil_sample_point);
            pupil_index
                .iter()
                .zip(values.iter())
                .for_each(|(k, vals)| {
                    temperature[*k] = *vals;
                });
            println!(" ... in {}ms", now.elapsed().as_millis());
        }
        _ => (),
    }

    println!("Temperature Min/Max : {:7.3?}K", minmax(&temperature));

    env::var("PLOT").and({
        let filename = format!("slice_{}.png", interp_method.to_lowercase());
        let plot = plt::png_canvas(&filename);
        plt::imagesc(&temperature, &plot);
        Ok(())
    })?;

    Ok(())
}
