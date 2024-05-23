use ghx_constrained_delaunay::glam::Vec2;

#[derive(Default)]
pub struct GhxCDTCrate {
    vertices: Vec<Vec2>,
}

impl crate::DelaunayCrate for GhxCDTCrate {
    type ResultType = ghx_constrained_delaunay::Triangulation;

    fn init(&mut self, vertices: impl Iterator<Item = [f64; 2]>) {
        self.vertices = vertices
            .map(|vertex| Vec2::new(vertex[0] as f32, vertex[1] as f32))
            .collect();
    }

    fn run_creation(&self) -> Self::ResultType {
        ghx_constrained_delaunay::triangulation_from_2d_vertices(&self.vertices)
    }
}
