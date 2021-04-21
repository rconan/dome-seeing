use csv;
use nalgebra::DVector;
use rbf_interp::{Basis, Scatter};
use rstar::{primitives::PointWithData, RTree, AABB};
use serde::Deserialize;
use std::path::Path;
use std::{env, fs::File, io::BufReader, time::Instant};

const L: f64 = 8.71;

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
pub struct Pupil {
    length: f64,
    d: f64,
    n: usize,
    z: f64,
}
pub struct PupilInner {
    pub length: f64,
    pub d: f64,
    pub n: usize,
    pub points: Vec<[f64; 3]>,
    pub index: Vec<usize>,
    pub z: f64,
}
impl Default for Pupil {
    fn default() -> Self {
        Self {
            length: 25.5,
            d: 0.05,
            n: (25.5 / 0.05) as usize + 1,
            z: 0f64,
        }
    }
}
impl Pupil {
    pub fn sample(self) -> PupilInner {
        let (length, d, n, z) = (self.length, self.d, self.n, self.z);
        // Pupil sampling
        let m1_radius = 8.365 * 0.5;
        let m2_radius = 3.3 * 0.5;
        let mut points = vec![];
        let mut index = vec![];
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
                        points.push(point);
                        index.push(k);
                        break;
                    }
                    sid += 1;
                    if sid > 7 {
                        break;
                    }
                }
            }
        }
        println!("# of points within the pupil: {}", points.len());
        PupilInner {
            length,
            d,
            n,
            points,
            index,
            z,
        }
    }
}
impl PupilInner {
    pub fn height(&mut self, z: f64) -> &mut Self {
        self.points.iter_mut().for_each(|p| {
            p[2] = z;
        });
        self
    }
}
type TemperatureNode = PointWithData<f64, [f64; 3]>;
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
pub struct TemperatureField {
    scattered_data: RTree<TemperatureNode>,
    pub z_minmax: (f64, f64),
    pub temp_minmax: (f64, f64),
    method: InterpolationMethod,
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
            if record.z > z_max {
                z_max = record.z;
            }
            if record.z < z_min {
                z_min = record.z;
            }
            if record.temperature > temp_max {
                temp_max = record.temperature;
            }
            if record.temperature < temp_min {
                temp_min = record.temperature;
            }
            tree.insert(TemperatureNode::new(
                record.temperature,
                [record.x, record.y, record.z],
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
