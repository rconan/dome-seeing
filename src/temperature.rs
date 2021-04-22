use super::{PupilInner, Result,InterpolationMethod};
use csv;
use nalgebra::DVector;
use rbf_interp::{Basis, Scatter};
use rstar::{primitives::PointWithData, RTree, AABB};
use serde::Deserialize;
use std::path::Path;
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

type TemperatureNode = PointWithData<f64, [f64; 3]>;
pub struct TemperatureField {
    pub scattered_data: RTree<TemperatureNode>,
    pub z_minmax: (f64, f64),
    pub temp_minmax: (f64, f64),
    pub method: InterpolationMethod,
}
impl TemperatureField {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let now = Instant::now();
        println!("Loading data into R-Tree ...");
        let file = File::open(path)?;
        let rdr = BufReader::with_capacity(100_000, file);
        let mut data = csv::Reader::from_reader(rdr);
        let mut tree = RTree::new();
        let mut z_max = f64::NEG_INFINITY;
        let mut z_min = f64::INFINITY;
        let mut temp_max = f64::NEG_INFINITY;
        let mut temp_min = f64::INFINITY;
        for result in data.deserialize() {
            let record: TemperatureData = result?;
            let z = record.z - 3.;
            if z > z_max {
                z_max = z;
            }
            if z < z_min {
                z_min = z;
            }
            if record.temperature > temp_max {
                temp_max = record.temperature;
            }
            if record.temperature < temp_min {
                temp_min = record.temperature;
            }
            tree.insert(TemperatureNode::new(
                record.temperature,
                [record.x, record.y, z],
            ));
        }
        println!(" ... in {}ms", now.elapsed().as_millis());
        println!("Tree size: {}", tree.size());
        println!("z min/max [{:8.3},{:8.3}]m", z_min, z_max);
        println!("temperature min/max [{:8.3},{:8.3}]m", temp_min, temp_max);
        Ok(Self {
            scattered_data: tree,
            z_minmax: (z_min, z_max),
            temp_minmax: (temp_min, temp_max),
            method: InterpolationMethod::from_env()?,
        })
    }
    pub fn method(self, method: InterpolationMethod) -> Self {
        Self { method, ..self }
    }
    pub fn z_slice(&self, z: f64, pupil: &mut PupilInner) -> Vec<f64> {
        pupil.height(z);
        let n = pupil.n;
        // Interpolation
        let mut temperature: Vec<f64> = vec![f64::NAN; n * n];
        use InterpolationMethod::*;
        match self.method {
            NearestNeighbor => {
                pupil
                    .points
                    .iter()
                    .zip(pupil.index.iter())
                    .for_each(|(point, k)| {
                        temperature[*k] = self.scattered_data.nearest_neighbor(point).unwrap().data;
                    })
            }
            RadialBasis => {
                let now = Instant::now();
                let bbox_thickness = 0.1;
                let bbox_width = 30.;
                let b_h = bbox_width * 0.5;
                let b_d = bbox_thickness * 0.5;
                let bbox = AABB::from_corners([-b_h, -b_h, z - b_d], [b_h, b_h, z + b_d]);
                let (sampled_points, sampled_data): (Vec<DVector<f64>>, Vec<DVector<f64>>) = self
                    .scattered_data
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
                    "# of points within bounding box of {:.2}m thickness: {} (in {}ms)",
                    bbox_thickness,
                    n_in_aabb,
                    now.elapsed().as_millis()
                );

                println!("RBF creation ...");
                let now = Instant::now();
                let scatter =
                    Scatter::create(sampled_points, sampled_data, Basis::PolyHarmonic(2), 2);
                println!(" ... in {}ms", now.elapsed().as_millis());

                println!("RBF interpolation ...");
                let now = Instant::now();
                let values = scatter.eval2(&pupil.points);
                pupil.index.iter().zip(values.iter()).for_each(|(k, vals)| {
                    temperature[*k] = *vals;
                });
                println!(" ... in {}ms", now.elapsed().as_millis());
            }
        }
        temperature
    }
}
