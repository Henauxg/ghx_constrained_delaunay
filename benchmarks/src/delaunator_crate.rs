use delaunator::Point;

use crate::Distribution;

#[derive(Default)]
pub struct DelaunatorCrate {
    vertices: Vec<Point>,
}

impl crate::DelaunayCrate for DelaunatorCrate {
    type ResultType = delaunator::Triangulation;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>) {
        self.vertices = vertices
            .map(|vertex| Point {
                x: vertex[0],
                y: vertex[1],
            })
            .collect();
    }

    fn run_triangulation(&self, _distribution: Distribution) -> Self::ResultType {
        delaunator::triangulate(&self.vertices)
    }
}
