use super::InterpolationMethod;
use nalgebra::DVector;
use plotters::prelude::*;
use rbf_interp::{Basis, Scatter};
use rstar::{primitives::PointWithData, RTree, AABB};

pub enum Projection {
    Top,
    SideX,
    SideY,
}

pub struct BoundingBox {
    pub width: f64,
    pub height: f64,
    pub depth: f64,
    center: [f64; 3],
}

impl BoundingBox {
    pub fn new(center: [f64; 3], width: f64, depth: f64, height: f64) -> Self {
        Self {
            width,
            height,
            depth,
            center,
        }
    }
    pub fn aabb(&self) -> AABB<[f64; 3]> {
        let lower = [
            self.center[0] - self.width * 0.5,
            self.center[1] - self.depth * 0.5,
            self.center[2] - self.height * 0.5,
        ];
        let upper = [
            self.center[0] + self.width * 0.5,
            self.center[1] + self.depth * 0.5,
            self.center[2] + self.height * 0.5,
        ];
        AABB::from_corners(lower, upper)
    }
    pub fn locate<'a>(
        &self,
        data: &'a RTree<PointWithData<f64, [f64; 3]>>,
    ) -> impl Iterator<Item = &'a PointWithData<f64, [f64; 3]>> {
        data.locate_in_envelope(&self.aabb())
    }
    pub fn locate_as_vector(
        &self,
        data: &RTree<PointWithData<f64, [f64; 3]>>,
    ) -> (Vec<DVector<f64>>, Vec<DVector<f64>>) {
        self.locate(data)
            .map(|p| {
                (
                    DVector::from_column_slice(p.position()),
                    DVector::from_column_slice(&[p.data]),
                )
            })
            .unzip()
    }
    pub fn interp(
        &self,
        points: &[[f64; 3]],
        data: &RTree<PointWithData<f64, [f64; 3]>>,
        method: InterpolationMethod,
    ) -> Vec<f64> {
        match method {
            InterpolationMethod::NearestNeighbor => points
                .iter()
                .map(|point| data.nearest_neighbor(point).unwrap().data)
                .collect(),
            InterpolationMethod::RadialBasis => {
                let (sampled_points, sampled_data): (Vec<DVector<f64>>, Vec<DVector<f64>>) =
                    self.locate_as_vector(data);

                //println!("RBF creation ...");
                let scatter =
                    Scatter::create(sampled_points, sampled_data, Basis::PolyHarmonic(2), 2);
                //println!("RBF interpolation ...");
                scatter.eval2(points)
            }
        }
    }
    pub fn show(
        &self,
        data: &RTree<PointWithData<f64, [f64; 3]>>,
        projection: Option<Projection>,
        decimation: Option<usize>,
    ) {
        let aabb = self.aabb();
        let l = aabb.lower();
        let u = aabb.upper();

        let root = SVGBackend::new("point_in_bbox.svg", (1024, 1024)).into_drawing_area();

        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .margin(40)
            //.caption("Empty 3D Figure", ("sans-serif", 40))
            .build_cartesian_3d(l[0]..u[0], l[2]..u[2], l[1]..u[1])
            .unwrap();
        match projection {
            Some(Projection::Top) => {
                chart.with_projection(|mut pb| {
                    pb.pitch = std::f64::consts::FRAC_PI_2;
                    pb.yaw = 0.;
                    pb.scale = 1.;
                    pb.into_matrix()
                });
            }
            Some(Projection::SideX) => {
                chart.with_projection(|mut pb| {
                    pb.pitch = 0.;
                    pb.yaw = 0.;
                    pb.scale = 1.;
                    pb.into_matrix()
                });
            }
            Some(Projection::SideY) => {
                chart.with_projection(|mut pb| {
                    pb.pitch = 0.;
                    pb.yaw = std::f64::consts::FRAC_PI_2;
                    pb.scale = 1.;
                    pb.into_matrix()
                });
            }
            None => (),
        };
        chart.configure_axes().draw().unwrap();

        chart
            .draw_series(
                data.locate_in_envelope(&aabb)
                    .step_by(decimation.unwrap_or(1))
                    .map(|p| {
                        let pos = p.position();
                        Circle::new((pos[0], pos[2], pos[1]), 1, BLACK.mix(0.5).filled())
                    }),
            )
            .unwrap();
    }
}
