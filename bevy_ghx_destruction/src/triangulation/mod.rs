use bevy::math::Vec2;

use crate::utils::{egdes_intersect, EdgesIntersectionResult};

pub mod constrained_triangulation;
pub mod triangulation;

pub type VertexId = usize;
pub type Neighbor = Option<VertexId>;
pub type TriangleId = usize;

// TODO New types to avoid some bugs
pub type TriangleVertexIndex = usize;
pub const VERT_1: TriangleVertexIndex = 0;
pub const VERT_2: TriangleVertexIndex = 1;
pub const VERT_3: TriangleVertexIndex = 2;

pub type TriangleEdgeIndex = usize;
pub const EDGE_12: TriangleEdgeIndex = 0;
pub const EDGE_23: TriangleEdgeIndex = 1;
pub const EDGE_31: TriangleEdgeIndex = 2;

#[derive(Debug, Copy, Clone, Eq, Hash, PartialEq)]
pub struct Edge {
    pub from: VertexId,
    pub to: VertexId,
}
impl Edge {
    #[inline]
    pub fn new(from: VertexId, to: VertexId) -> Self {
        Self { from, to }
    }

    pub fn undirected_equals(&self, other: &Edge) -> bool {
        self == other || (self.from == other.to && self.to == other.from)
    }

    pub fn to_vertices(&self, vertices: &Vec<Vec2>) -> EdgeVertices {
        (vertices[self.from], vertices[self.to])
    }

    pub fn contains(&self, vert: VertexId) -> bool {
        self.from == vert || self.to == vert
    }
}
impl From<(VertexId, VertexId)> for Edge {
    fn from(vertices: (VertexId, VertexId)) -> Edge {
        Edge::new(vertices.0, vertices.1)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
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
        Edge::new(self.verts[VERT_1], self.verts[VERT_2])
    }
    #[inline]
    pub fn edge23(&self) -> Edge {
        Edge::new(self.verts[VERT_2], self.verts[VERT_3])
    }
    #[inline]
    pub fn edge31(&self) -> Edge {
        Edge::new(self.verts[VERT_3], self.verts[VERT_1])
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
            Edge::new(self.verts[VERT_1], self.verts[VERT_2]),
            Edge::new(self.verts[VERT_2], self.verts[VERT_3]),
            Edge::new(self.verts[VERT_3], self.verts[VERT_1]),
        ]
    }

    #[inline]
    pub fn edge(&self, edge: TriangleEdgeIndex) -> Edge {
        Edge::new(self.verts[edge], self.verts[(edge + 1) % 3])
    }

    #[inline]
    pub fn other_edges(&self, edge: TriangleEdgeIndex) -> [Edge; 2] {
        [self.edge((edge + 1) % 3), self.edge((edge + 2) % 3)]
    }

    #[inline]
    pub fn to_vertices(&self, vertices: &Vec<Vec2>) -> TriangleVertices {
        (
            vertices[self.verts[VERT_1]],
            vertices[self.verts[VERT_2]],
            vertices[self.verts[VERT_3]],
        )
    }

    #[inline]
    pub fn vertex_index(&self, vertex_id: VertexId) -> Option<TriangleVertexIndex> {
        for (v_index, v) in self.verts.iter().enumerate() {
            if *v == vertex_id {
                return Some(v_index);
            }
        }
        None
    }

    /// `vertex` MUST be a vertex of the triangle
    #[inline]
    pub fn get_opposite_edge_index_from_vertex(&self, vertex: VertexId) -> TriangleEdgeIndex {
        TriangleData::get_opposite_edge_index(self.vertex_index(vertex).unwrap())
    }

    /// `edge` MUST be an edge of the triangle
    #[inline]
    pub fn get_opposite_vertex_index(&self, edge: &Edge) -> TriangleVertexIndex {
        if !edge.contains(self.verts[VERT_1]) {
            self.verts[VERT_1]
        } else if !edge.contains(self.verts[VERT_2]) {
            self.verts[VERT_2]
        } else {
            self.verts[VERT_3]
        }
    }
}

impl TriangleData {
    // vert_index: index of the constrained edge start vertex in the current triangle
    // Left -> second_vertex_index = (vert_index + 1) % 3 => edge_index
    // 0 -> 1 => 0
    // 1 -> 2 => 1
    // 2 -> 0 => 2
    // Right -> second_vertex_index = (vert_index + 2) % 3  => edge_index
    // 0 -> 2 => 2
    // 1 -> 0 => 0
    // 2 -> 1 => 1
    #[inline]
    pub fn get_left_edge_index_around(vert_index: TriangleVertexIndex) -> TriangleEdgeIndex {
        vert_index
    }
    #[inline]
    pub fn get_right_edge_index_around(vert_index: TriangleVertexIndex) -> TriangleEdgeIndex {
        (vert_index + 2) % 3
    }

    // Opposite -> edge_index = (vert_index + 1) % 3
    // 0 -> 1
    // 1 -> 2
    // 2 -> 0
    #[inline]
    pub fn get_opposite_edge_index(vert_index: TriangleVertexIndex) -> TriangleEdgeIndex {
        (vert_index + 1) % 3
    }
}

pub type QuadVertexIndex = usize;
pub const QUAD_1: QuadVertexIndex = 0;
pub const QUAD_2: QuadVertexIndex = 1;
pub const QUAD_3: QuadVertexIndex = 2;
pub const QUAD_4: QuadVertexIndex = 3;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Quad {
    pub verts: [VertexId; 4],
}
impl Quad {
    #[inline]
    pub fn new(verts: [VertexId; 4]) -> Self {
        Self { verts }
    }

    #[inline]
    pub fn v1(&self) -> VertexId {
        self.verts[QUAD_1]
    }
    #[inline]
    pub fn v2(&self) -> VertexId {
        self.verts[QUAD_2]
    }
    #[inline]
    pub fn v3(&self) -> VertexId {
        self.verts[QUAD_3]
    }
    #[inline]
    pub fn v4(&self) -> VertexId {
        self.verts[QUAD_4]
    }

    #[inline]
    pub fn to_vertices(&self, vertices: &Vec<Vec2>) -> QuadVertices {
        QuadVertices([
            vertices[self.verts[QUAD_1]],
            vertices[self.verts[QUAD_2]],
            vertices[self.verts[QUAD_3]],
            vertices[self.verts[QUAD_4]],
        ])
    }
}

// TODO may change to structured arrays [Vec2;n]
pub type EdgeVertices = (Vec2, Vec2);
pub type TriangleVertices = (Vec2, Vec2, Vec2);

pub struct QuadVertices(pub [Vec2; 4]);
impl QuadVertices {
    #[inline]
    pub fn q1(&self) -> Vec2 {
        self.0[QUAD_1]
    }
    #[inline]
    pub fn q2(&self) -> Vec2 {
        self.0[QUAD_2]
    }
    #[inline]
    pub fn q3(&self) -> Vec2 {
        self.0[QUAD_3]
    }
    #[inline]
    pub fn q4(&self) -> Vec2 {
        self.0[QUAD_4]
    }

    pub fn diagonals_intersection_test(&self) -> EdgesIntersectionResult {
        egdes_intersect(
            &(self.0[QUAD_1], self.0[QUAD_3]),
            &(self.0[QUAD_2], self.0[QUAD_4]),
        )
    }
}
