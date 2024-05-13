pub mod constrained_triangulation;
pub mod triangulation;

pub type VertexId = usize;
pub type TriangleId = usize;
pub type Edge = (VertexId, VertexId);

pub type VertexLocalIndex = usize;
pub const VERT_1: VertexLocalIndex = 0;
pub const VERT_2: VertexLocalIndex = 1;
pub const VERT_3: VertexLocalIndex = 2;

pub type NeighborIndex = usize;
pub const NEIGHBOR_12: NeighborIndex = 0;
pub const NEIGHBOR_23: NeighborIndex = 1;
pub const NEIGHBOR_31: NeighborIndex = 2;

pub type TriangleEdge = (VertexLocalIndex, VertexLocalIndex);
pub const TRIANGLE_EDGES: [TriangleEdge; 3] =
    [(VERT_1, VERT_2), (VERT_2, VERT_3), (VERT_3, VERT_1)];

pub type Neighbor = Option<VertexId>;

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
    pub fn set_v1(&mut self, val: VertexId) {
        self.verts[VERT_1] = val;
    }

    #[inline]
    pub fn set_v2(&mut self, val: VertexId) {
        self.verts[VERT_2] = val;
    }

    #[inline]
    pub fn set_v3(&mut self, val: VertexId) {
        self.verts[VERT_3] = val;
    }

    #[inline]
    pub fn set_verts(&mut self, val: [VertexId; 3]) {
        self.verts = val;
    }

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
    pub fn v1_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_1]
    }

    #[inline]
    pub fn v2_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_2]
    }

    #[inline]
    pub fn v3_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_3]
    }

    #[inline]
    pub fn neighbor12(&self) -> Neighbor {
        self.neighbors[NEIGHBOR_12]
    }

    #[inline]
    pub fn neighbor23(&self) -> Neighbor {
        self.neighbors[NEIGHBOR_23]
    }

    #[inline]
    pub fn neighbor31(&self) -> Neighbor {
        self.neighbors[NEIGHBOR_31]
    }

    #[inline]
    pub fn set_neighbor12(&mut self, neighbor: Neighbor) {
        self.neighbors[NEIGHBOR_12] = neighbor;
    }

    #[inline]
    pub fn set_neighbor23(&mut self, neighbor: Neighbor) {
        self.neighbors[NEIGHBOR_23] = neighbor;
    }

    #[inline]
    pub fn set_neighbor31(&mut self, neighbor: Neighbor) {
        self.neighbors[NEIGHBOR_31] = neighbor;
    }

    #[inline]
    pub fn set_neighbors(&mut self, val: [Neighbor; 3]) {
        self.neighbors = val;
    }

    #[inline]
    pub fn edge_from_triangle_edge(&self, edge: &TriangleEdge) -> Edge {
        (self.verts[edge.0], self.verts[edge.1])
    }
}

#[derive(Debug)]
pub struct Quad {
    pub v1: VertexId,
    pub v2: VertexId,
    pub v3: VertexId,
    pub v4: VertexId,
}
