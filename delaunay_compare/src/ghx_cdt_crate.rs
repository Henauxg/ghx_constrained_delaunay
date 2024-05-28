use ghx_constrained_delaunay::{triangulation::TriangulationConfiguration, types::Vertice};

use crate::Distribution;

fn bin_density_from_distribution(dis: Distribution) -> f64 {
    match dis {
        Distribution::Uniform => 0.85,
        Distribution::Local => 0.5,
    }
}

#[derive(Default)]
pub struct GhxCDTCrate {
    vertices: Vec<Vertice>,
}

impl crate::DelaunayCrate for GhxCDTCrate {
    type ResultType = ghx_constrained_delaunay::Triangulation;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>) {
        self.vertices = vertices
            .map(|vertex| Vertice::new(vertex[0], vertex[1]))
            .collect();
    }

    fn run_creation(&self, distribution: Distribution) -> Self::ResultType {
        ghx_constrained_delaunay::triangulation_from_2d_vertices(
            &self.vertices,
            TriangulationConfiguration {
                bin_vertex_density_power: bin_density_from_distribution(distribution),
            },
        )
    }
}
