//!
//! # Analytic Ray Tracing
//!

use nalgebra as na;
use std::fmt;

pub type Vector = [f64; 3];

pub trait Arithmetic {
    fn dot(&self, other: &[f64]) -> f64;
    fn norm_square(&self) -> f64;
    fn norm(&self) -> f64;
    fn normalize(&mut self) -> Self;
    fn add(&self, other: Self) -> Self;
    fn sub(&self, other: Self) -> Self;
}
impl Arithmetic for Vector {
    fn dot(&self, other: &[f64]) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
    fn norm_square(&self) -> f64 {
        self.dot(self)
    }
    fn norm(&self) -> f64 {
        self.norm_square().sqrt()
    }
    fn normalize(&mut self) -> Self {
        let n = self.norm();
        self[0] /= n;
        self[1] /= n;
        self[2] /= n;
        *self
    }
    fn add(&self, other: Self) -> Self {
        [self[0] + other[0], self[1] + other[1], self[2] + other[2]]
    }
    fn sub(&self, other: Self) -> Self {
        [self[0] - other[0], self[1] - other[1], self[2] - other[2]]
    }
}

/// # Ray definition
///
/// A ray is defined with:
///  - a point of origin: $\vec p = [x,y,z]$,
///  - a direction vector: $\vec u = [k,l,m]$ such as $\| \vec u \|=1$.
///
/// The ray tracing equation is given by: $$\vec{p^\prime} = \vec p + s \vec u,$$ where $s$ is the optical path length.
pub struct Ray {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
/// # Ray builder
///
/// Build a new [`Ray`](crate::analytic::Ray)
pub struct NewRay {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
impl Default for NewRay {
    fn default() -> Self {
        Self {
            p: [0f64; 3],
            u: [0f64, 0f64, -1f64],
        }
    }
}
impl NewRay {
    /// Build the [`Ray`](crate::analytic::Ray)
    pub fn build(self) -> Ray {
        Ray {
            p: self.p,
            u: self.u,
        }
    }
    /// Set the [`Ray`](crate::analytic::Ray) point of origin
    pub fn point_of_origin(self, p: Vector) -> Self {
        Self { p, ..self }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector
    pub fn direction_vector(self, u: Vector) -> Self {
        Self { u, ..self }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector from polar coordinates
    pub fn polar_direction_vector(self, z: f64, a: f64) -> Self {
        let ca = a.cos();
        let sa = a.sin();
        let sz = z.sin();
        let cz = z.cos();
        let u = [sz * ca, sz * sa, -cz].normalize();
        self.direction_vector(u)
    }
}
/// Create a [`NewRay`](crate::analytic::NewRay) at the origin propagate downward (z<0)
pub fn new_ray() -> NewRay {
    NewRay::default()
}
impl Ray {
    /// Compute the distance $s$ from the ray current location to [`Conic`](crate::analytic::Conic)
    /// We find the distance $s$ from:
    /// $$s=\frac{mR-\vec\alpha\cdot\vec\beta + \sqrt{(\vec\alpha\cdot\vec\beta-mR)^2-\|\vec\beta\|^2(\|\vec\alpha\|^2-2zR)}}{\|\vec\beta\|^2}$$
    /// $$\vec\alpha = [x,y,z\sqrt{\kappa+1}]$$
    /// $$\vec\beta = [k,l,m\sqrt{\kappa+1}]$$
    pub fn distance_to(&self, conic: &Conic) -> f64 {
        let q = (conic.constant + 1f64).sqrt();
        let p: Vector = [
            self.p[0] - conic.origin[0],
            self.p[1] - conic.origin[1],
            self.p[2] - conic.origin[2],
        ];
        let alpha: Vector = [p[0], p[1], p[2] * q];
        let beta: Vector = [self.u[0], self.u[1], self.u[2] * q];
        let a = beta.norm_square();
        let b = 2f64 * (alpha.dot(&beta) - self.u[2] * conic.radius);
        let c = alpha.norm_square() - 2f64 * p[2] * conic.radius;
        0.5 * (-b + (b * b - 4f64 * a * c).sqrt()) / a
    }
    /// Trace ray from ray current position to [`Conic`](crate::analytic::Conic)
    pub fn trace_to(&mut self, conic: &Conic) {
        let s = self.distance_to(conic);
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
    }
    pub fn trace(&mut self, s: f64) {
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
    }
    /// Solve ray tracing equation for $z$ given $x$ and $y$
    pub fn solve_for_z(&self, x: f64, y: f64) -> f64 {
        let x = x - self.p[0];
        let y = y - self.p[1];
        let num = x * x + y * y;
        let denom = self.u[0] * self.u[0] + self.u[1] * self.u[1];
        if denom < 1e-30 {
            std::f64::INFINITY
        } else {
            self.p[2] + self.u[2] * (num / denom).sqrt()
        }
    }
}
impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "P: [{:+15.9},{:+15.9},{:+15.9}] ; U: [{:+.9},{:+.9},{:+.9}]",
            self.p[0], self.p[1], self.p[2], self.u[0], self.u[1], self.u[2],
        )
    }
}

/// # Conic definition
///
/// A conic surface is defined by the set of coordinates $(x,y,z)$ that satisfies
/// $$ F(x,y,z) = r^2 - 2zR + z^2(\kappa+1)=0,$$
/// where $r^2=x^2+y^2$, $R$ is the radius of curvature and $\kappa$ is the conic constant.
pub struct Conic {
    /// Conic constant $\kappa$
    pub constant: f64,
    /// Radius of curvature
    pub radius: f64,
    /// Origin vector
    pub origin: Vector,
}
impl Conic {
    /// Creates a new `Conic`
    pub fn new(constant: f64, radius: f64) -> Self {
        Self {
            constant,
            radius,
            origin: [0f64; 3],
        }
    }
    /// Creates a new `Conic` with GMT M1 prescription:
    ///  - $\kappa=-0.9982857$
    ///  - $R=36$
    pub fn gmt_m1() -> Self {
        Self::new(-0.9982857, 36.0)
    }
    /// Creates a new `Conic` with GMT M2 prescription:
    ///  - $\kappa=-0.71692784$
    ///  - $R=-4.1639009$
    pub fn gmt_m2() -> Self {
        Self {
            constant: -0.71692784,
            radius: -4.1639009,
            origin: [0f64, 0f64, 20.26247614],
        }
    }
    /// Solve conic surface $F(x,y,z)$ for z :
    pub fn height_at(&self, v: Vector) -> Vector {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        _v[2] = self.c() * r2 / (1f64 + self.sqrt_(r2));
        _v
    }
    /// Conic $x$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial x} = -\frac{x}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn x_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[0] / self.sqrt_(r2)
    }
    /// Conic $y$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial y} = -\frac{y}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn y_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[1] / self.sqrt_(r2)
    }
    /// Conic $z$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial z} = 1
    /// $$
    ///
    pub fn z_partial_at(&self, _v: Vector) -> f64 {
        1f64
    }
    /// Normal vector, $\vec n = [\partial_x F,\partial_y F,\partial_z F]$, to conic surface
    pub fn normal_at(&self, v: Vector) -> Vector {
        [
            self.x_partial_at(v),
            self.y_partial_at(v),
            self.z_partial_at(v),
        ]
        .normalize()
    }
    /// Reflect ray from conic surface
    /// The ray direction vector after reflection is given by: $\vec{u^\prime} = \vec u - 2 (\vec u \cdot \vec n)\vec n$
    pub fn reflect(&self, ray: &mut Ray) {
        let n = self.normal_at(ray.p);
        let q = 2f64 * ray.u.dot(&n);
        ray.u[0] -= q * n[0];
        ray.u[1] -= q * n[1];
        ray.u[2] -= q * n[2];
        ray.u = ray.u.normalize();
    }
    fn kp1(&self) -> f64 {
        self.constant + 1f64
    }
    fn c(&self) -> f64 {
        1f64 / self.radius
    }
    fn sqrt_(&self, r2: f64) -> f64 {
        (1f64 - self.kp1() * self.c() * self.c() * r2).sqrt()
    }
}
pub struct Gmt {
    pub m1: Conic,
    pub m2: Conic,
}
impl Gmt {
    pub fn new() -> Self {
        Self {
            m1: Conic::gmt_m1(),
            m2: Conic::gmt_m2(),
        }
    }
    pub fn trace(&self, rays: &mut [Ray]) {
        rays.iter_mut().for_each(|mut r| {
            self.m1.reflect(&mut r);
            r.trace_to(&self.m2);
            self.m2.reflect(&mut r);
        })
    }
    pub fn focal_point(&self, marginals: Vec<Vector>, z: f64, a: f64) -> Vec<Ray> {
        // Chief ray at M1 vertex from field angle (z,a)
        let chief_ray = new_ray().polar_direction_vector(z, a).build();
        // Marginal ray on M1 surface from field angle (z,a)
        let marginal_rays: Vec<Ray> = marginals
            .into_iter()
            .map(|m| {
                new_ray()
                    .point_of_origin(self.m1.height_at(m))
                    .polar_direction_vector(z, a)
                    .build()
            })
            .collect();
        // Collecting the rays
        let mut rays = vec![chief_ray];
        marginal_rays.into_iter().for_each(|m| rays.push(m));
        // Ray tracing to and reflecting from M2
        self.trace(&mut rays);
        // De-structuring chief ray
        let chief_ray = rays.remove(0);
        let p0 = chief_ray.p;
        let u0 = chief_ray.u;
        // De-structuring marginal rays
        let p: Vec<Vector> = rays.iter().map(|m| m.p).collect();
        let u: Vec<Vector> = rays.iter().map(|m| m.u).collect();
        // Build A matrix
        let n_u = u.len();
        let z: Vector = [0f64; 3];
        let mut cols: Vec<Vec<f64>> = vec![u0.to_vec();n_u];
        for i_row in 0..n_u {
            let mut el = vec![];
            for i_col in 0..n_u {
                if i_row == i_col {
                    el.push(z.sub(u[i_col]).to_vec());
                } else {
                    el.push(z.clone().to_vec());
                }
            }
            cols.push(el.into_iter().flatten().collect());
        }
        let el: Vec<f64> = cols.into_iter().flatten().collect();
        let a = na::DMatrix::from_column_slice(n_u * 3, n_u + 1, &el);
        // Building b vector
        let b = na::DVector::from_vec(
            p.iter()
                .map(|x| x.sub(p0).to_vec())
                .flatten()
                .collect::<Vec<f64>>(),
        );
        // Solving As=b
        let s = a.svd(true, true).solve(&b, std::f64::EPSILON).unwrap();
        // Ray tracing to focal plane
        rays.insert(0, chief_ray);
        rays.iter_mut().zip(s.into_iter()).for_each(|x| x.0.trace(*x.1));
        rays
    }
}
