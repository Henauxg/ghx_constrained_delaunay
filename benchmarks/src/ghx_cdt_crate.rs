use ghx_constrained_delaunay::{triangulation::TriangulationConfiguration, types::Vertex};

use crate::Distribution;

#[derive(Default)]
pub struct GhxCDTCrate {
    vertices: Vec<Vertex>,
}

impl crate::DelaunayCrate for GhxCDTCrate {
    type ResultType = ghx_constrained_delaunay::Triangulation;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>) {
        self.vertices = vertices
            .map(|vertex| Vertex::new(vertex[0], vertex[1]))
            .collect();
    }

    fn run_triangulation(&self, _distribution: Distribution) -> Self::ResultType {
        ghx_constrained_delaunay::triangulation_from_2d_vertices(
            &self.vertices,
            TriangulationConfiguration::default(),
        )
        .unwrap()
    }
}
