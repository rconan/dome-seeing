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
                        let o = -60. * (sid - 1) as f64;
                        let (s, c) = (90. + o).to_radians().sin_cos();
                        let t = geotrans::Vector::from([L * c, L * s, 0.]);
                        let q = geotrans::Quaternion::unit(o.to_radians(), geotrans::Vector::k());
                        let p: geotrans::Quaternion = u.into();
                        geotrans::Vector::from(
                            (q.complex_conjugate() * (p - t.into()) * &q).vector_as_slice(),
                        )
                    } else {
                        u
                    };
                    //println!("{:?}",v);
                    let s = if sid < 7 {
                        13.601685f64.to_radians().cos()
                    } else {
                        1f64
                    };
                    let r = v[0].hypot(v[1] / s);
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
    pub fn nnz(&self) -> usize {
        self.points.len()
    }
    pub fn mask(&self) -> Vec<f64> {
        let mut m = vec![0f64;self.n*self.n];
        self.index.iter().for_each(|k| {
            m[*k] = 1f64;
        });
        m
    }
}
