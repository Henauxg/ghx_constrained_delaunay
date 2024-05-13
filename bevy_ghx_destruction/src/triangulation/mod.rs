pub mod constrained_triangulation;
pub mod triangulation;

pub type VertexId = usize;
pub type TriangleId = usize;

#[derive(Debug, Clone)]
pub struct TriangleData {
    // Vertices ids
    pub v1: VertexId,
    pub v2: VertexId,
    pub v3: VertexId,
    // Neighbours
    pub edge12: Option<TriangleId>,
    pub edge23: Option<TriangleId>,
    pub edge31: Option<TriangleId>,
    //pub edges: [TriangleId;3],
}
