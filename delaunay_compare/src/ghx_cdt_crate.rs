use ghx_constrained_delaunay::{triangulation::TriangulationConfiguration, types::Vertice};

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

    fn run_creation(&self) -> Self::ResultType {
        ghx_constrained_delaunay::triangulation_from_2d_vertices(
            &self.vertices,
            TriangulationConfiguration::default(),
        )
    }
}
