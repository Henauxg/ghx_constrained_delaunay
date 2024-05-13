pub mod constrained_triangulation;
pub mod triangulation;

pub type VertexId = usize;
pub type Neighbor = Option<VertexId>;
pub type TriangleId = usize;
pub type Edge = (VertexId, VertexId);

pub type TriangleVertexIndex = usize;
pub const VERT_1: TriangleVertexIndex = 0;
pub const VERT_2: TriangleVertexIndex = 1;
pub const VERT_3: TriangleVertexIndex = 2;

pub type TriangleEdgeIndex = usize;
pub const EDGE_12: TriangleEdgeIndex = 0;
pub const EDGE_23: TriangleEdgeIndex = 1;
pub const EDGE_31: TriangleEdgeIndex = 2;

#[derive(Debug, Clone)]
pub struct TriangleData {
    /// Triangle vertices indexes
    pub verts: [VertexId; 3],
    /// Adjacent neighbors on each edge.
    ///
    /// Can be [`None`] for the container triangle added by the triangulation process which has no neighbors,
    /// and contain all the points to triangulate
    pub neighbors: [Neighbor; 3],
}

impl TriangleData {
    #[inline]
    pub fn v1(&self) -> VertexId {
        self.verts[VERT_1]
    }

    #[inline]
    pub fn v2(&self) -> VertexId {
        self.verts[VERT_2]
    }

    #[inline]
    pub fn v3(&self) -> VertexId {
        self.verts[VERT_3]
    }

    #[inline]
    pub fn edge12(&self) -> Edge {
        (self.verts[VERT_1], self.verts[VERT_2])
    }

    #[inline]
    pub fn edge23(&self) -> Edge {
        (self.verts[VERT_2], self.verts[VERT_3])
    }

    #[inline]
    pub fn edge31(&self) -> Edge {
        (self.verts[VERT_3], self.verts[VERT_1])
    }

    #[inline]
    pub fn neighbor12(&self) -> Neighbor {
        self.neighbors[EDGE_12]
    }

    #[inline]
    pub fn neighbor23(&self) -> Neighbor {
        self.neighbors[EDGE_23]
    }

    #[inline]
    pub fn neighbor31(&self) -> Neighbor {
        self.neighbors[EDGE_31]
    }

    #[inline]
    pub fn edges(&self) -> [Edge; 3] {
        [
            (self.verts[VERT_1], self.verts[VERT_2]),
            (self.verts[VERT_2], self.verts[VERT_3]),
            (self.verts[VERT_3], self.verts[VERT_1]),
        ]
    }
}

#[derive(Debug)]
pub struct Quad {
    pub v1: VertexId,
    pub v2: VertexId,
    pub v3: VertexId,
    pub v4: VertexId,
}
