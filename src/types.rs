use crate::utils::{egdes_intersect, EdgesIntersectionResult};

#[cfg(not(feature = "f64"))]
pub mod f32;
#[cfg(not(feature = "f64"))]
pub use f32::*;

#[cfg(feature = "f64")]
pub mod f64;
#[cfg(feature = "f64")]
pub use f64::*;

#[cfg(not(feature = "u64_indexes"))]
pub mod u32;
#[cfg(not(feature = "u64_indexes"))]
pub use u32::IndexType;

#[cfg(feature = "u64_indexes")]
pub mod u64;
#[cfg(feature = "u64_indexes")]
pub use u64::IndexType;

pub type VertexId = IndexType;
pub type TriangleId = IndexType;

// TODO Could use new types to avoid possible unwanted misuses

/// Local index of a vertex in a triangle
pub type TriangleVertexIndex = u8;
pub const VERT_1: TriangleVertexIndex = 0;
pub const VERT_2: TriangleVertexIndex = 1;
pub const VERT_3: TriangleVertexIndex = 2;

/// Local index of an edge in a triangle
pub type TriangleEdgeIndex = u8;
pub const EDGE_12: TriangleEdgeIndex = 0;
pub const EDGE_23: TriangleEdgeIndex = 1;
pub const EDGE_31: TriangleEdgeIndex = 2;

/// Local index of a vertex in a quad
pub type QuadVertexIndex = u8;
pub const QUAD_1: QuadVertexIndex = 0;
pub const QUAD_2: QuadVertexIndex = 1;
pub const QUAD_3: QuadVertexIndex = 2;
pub const QUAD_4: QuadVertexIndex = 3;

/// From a TriangleEdgeIndex, gives the next edge index in a clockwise order
pub const NEXT_CLOCKWISE_EDGE_INDEX: [TriangleEdgeIndex; 3] = [EDGE_23, EDGE_31, EDGE_12];
/// From a TriangleEdgeIndex, gives the next edge index in a counter-clockwise order
pub const NEXT_COUNTER_CLOCKWISE_EDGE_INDEX: [TriangleEdgeIndex; 3] = [EDGE_31, EDGE_12, EDGE_23];
/// From a TriangleEdgeIndex, gives the corresponding pair of TriangleVertexIndex
pub const EDGE_TO_VERTS: [[TriangleVertexIndex; 2]; 3] =
    [[VERT_1, VERT_2], [VERT_2, VERT_3], [VERT_3, VERT_1]];
/// From a TriangleEdgeIndex, gives the opposite vertex index
pub const OPPOSITE_VERTEX_INDEX: [TriangleVertexIndex; 3] = [VERT_3, VERT_1, VERT_2];

/// From a TriangleVertexIndex, gives the next vertex index in a clockwise order
pub const NEXT_CW_VERTEX_INDEX: [TriangleVertexIndex; 3] = [VERT_2, VERT_3, VERT_1];
/// From a TriangleVertexIndex, gives the next vertex index in a counter-clockwise order
pub const NEXT_CCW_VERTEX_INDEX: [TriangleVertexIndex; 3] = [VERT_3, VERT_1, VERT_2];
/// From a TriangleVertexIndex, gives the opposite edge index
pub const OPPOSITE_EDGE_INDEX: [TriangleEdgeIndex; 3] = [EDGE_23, EDGE_31, EDGE_12];
/// From a TriangleVertexIndex, gives the next edge index in a CCW order
pub const NEXT_CCW_EDGE_INDEX: [TriangleEdgeIndex; 3] = [EDGE_31, EDGE_12, EDGE_23];

/// From a TriangleVertexIndex, gives the outermost TriangleEdgeIndex in a clockwise order
///
/// Used for looping around a vertex
pub const OUTERMOST_CW_EDGE_INDEX_AROUND_VERTEX: [TriangleEdgeIndex; 3] =
    [EDGE_31, EDGE_12, EDGE_23];
/// From a TriangleVertexIndex, gives the outermost TriangleEdgeIndex in a counter-clockwise order
///
/// Used for looping around a vertex
pub const OUTERMOST_CCW_EDGE_INDEX_AROUND_VERTEX: [TriangleEdgeIndex; 3] =
    [EDGE_12, EDGE_23, EDGE_31];

pub const ADJACENT_QUAD_VERTICES_INDEXES: [[TriangleVertexIndex; 2]; 4] = [
    [QUAD_3, QUAD_4],
    [QUAD_4, QUAD_3],
    [QUAD_2, QUAD_1],
    [QUAD_1, QUAD_2],
];

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

    #[inline]
    pub fn undirected_equals(&self, other: &Edge) -> bool {
        self == other || (self.from == other.to && self.to == other.from)
    }

    #[inline]
    pub fn to_vertices(&self, vertices: &Vec<Vertex>) -> EdgeVertices {
        (vertices[self.from as usize], vertices[self.to as usize])
    }

    #[inline]
    pub fn contains(&self, vert: VertexId) -> bool {
        self.from == vert || self.to == vert
    }

    #[inline]
    pub fn opposite(&self) -> Self {
        Self {
            from: self.to,
            to: self.from,
        }
    }
}
impl From<(VertexId, VertexId)> for Edge {
    fn from(vertices: (VertexId, VertexId)) -> Edge {
        Edge::new(vertices.0, vertices.1)
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
/// Type used to represent an optionnal [TriangleId].
///
/// More packed in memory than an [Option], for efficiency.
pub struct Neighbor {
    pub id: TriangleId,
}

impl Neighbor {
    pub const NONE: Neighbor = Neighbor {
        id: TriangleId::MAX,
    };

    #[inline]
    pub(crate) fn new(id: TriangleId) -> Self {
        Self { id }
    }

    #[inline]
    pub fn exists(&self) -> bool {
        self.id != TriangleId::MAX
    }
}

impl From<TriangleId> for Neighbor {
    #[inline]
    fn from(id: TriangleId) -> Self {
        Self { id }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TriangleData {
    /// Triangle vertices indexes
    pub verts: [VertexId; 3],
    /// Adjacent neighbors on each edge.
    ///
    /// Can be [`None`] for the container triangle added by the triangulation process which has no neighbors,
    /// and contains all the points to triangulate
    pub neighbors: [Neighbor; 3],
}

impl TriangleData {
    pub(crate) fn new_container_triangle(first_index: VertexId) -> Self {
        TriangleData {
            verts: [first_index, first_index + 1, first_index + 2],
            neighbors: [Neighbor::NONE, Neighbor::NONE, Neighbor::NONE],
        }
    }

    #[inline]
    pub fn v(&self, vertex_index: TriangleVertexIndex) -> VertexId {
        self.verts[vertex_index as usize]
    }
    #[inline]
    pub fn v1(&self) -> VertexId {
        self.verts[VERT_1 as usize]
    }
    #[inline]
    pub fn v2(&self) -> VertexId {
        self.verts[VERT_2 as usize]
    }
    #[inline]
    pub fn v3(&self) -> VertexId {
        self.verts[VERT_3 as usize]
    }

    #[inline]
    pub fn v1_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_1 as usize]
    }
    #[inline]
    pub fn v2_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_2 as usize]
    }
    #[inline]
    pub fn v3_mut(&mut self) -> &mut VertexId {
        &mut self.verts[VERT_3 as usize]
    }

    #[inline]
    pub fn edge12(&self) -> Edge {
        Edge::new(self.verts[VERT_1 as usize], self.verts[VERT_2 as usize])
    }
    #[inline]
    pub fn edge23(&self) -> Edge {
        Edge::new(self.verts[VERT_2 as usize], self.verts[VERT_3 as usize])
    }
    #[inline]
    pub fn edge31(&self) -> Edge {
        Edge::new(self.verts[VERT_3 as usize], self.verts[VERT_1 as usize])
    }

    #[inline]
    pub fn neighbor(&self, egde_index: TriangleEdgeIndex) -> Neighbor {
        self.neighbors[egde_index as usize]
    }

    #[inline]
    pub fn neighbor12(&self) -> Neighbor {
        self.neighbors[EDGE_12 as usize]
    }
    #[inline]
    pub fn neighbor23(&self) -> Neighbor {
        self.neighbors[EDGE_23 as usize]
    }
    #[inline]
    pub fn neighbor31(&self) -> Neighbor {
        self.neighbors[EDGE_31 as usize]
    }

    #[inline]
    pub fn neighbor12_mut(&mut self) -> &mut Neighbor {
        &mut self.neighbors[EDGE_12 as usize]
    }
    #[inline]
    pub fn neighbor23_mut(&mut self) -> &mut Neighbor {
        &mut self.neighbors[EDGE_23 as usize]
    }
    #[inline]
    pub fn neighbor31_mut(&mut self) -> &mut Neighbor {
        &mut self.neighbors[EDGE_31 as usize]
    }

    #[inline]
    pub fn edges(&self) -> [Edge; 3] {
        [self.edge12(), self.edge23(), self.edge31()]
    }

    #[inline]
    pub fn edge(&self, edge_index: TriangleEdgeIndex) -> Edge {
        let vert_indexes = EDGE_TO_VERTS[edge_index as usize];
        Edge::new(
            self.verts[vert_indexes[0] as usize],
            self.verts[vert_indexes[1] as usize],
        )
    }

    #[inline]
    pub fn other_edges(&self, edge_index: TriangleEdgeIndex) -> [Edge; 2] {
        [
            self.edge(NEXT_CLOCKWISE_EDGE_INDEX[edge_index as usize]),
            self.edge(NEXT_COUNTER_CLOCKWISE_EDGE_INDEX[edge_index as usize]),
        ]
    }

    #[inline]
    pub fn to_vertices(&self, vertices: &Vec<Vertex>) -> TriangleVertices {
        (
            vertices[self.verts[VERT_1 as usize] as usize],
            vertices[self.verts[VERT_2 as usize] as usize],
            vertices[self.verts[VERT_3 as usize] as usize],
        )
    }
    #[inline]
    pub fn to_vertices_array(&self, vertices: &Vec<Vertex>) -> [Vertex; 3] {
        [
            vertices[self.verts[VERT_1 as usize] as usize],
            vertices[self.verts[VERT_2 as usize] as usize],
            vertices[self.verts[VERT_3 as usize] as usize],
        ]
    }

    /// `vertex` MUST be a vertex of the triangle
    #[inline]
    pub fn vertex_index(&self, vertex: VertexId) -> TriangleVertexIndex {
        if self.v1() == vertex {
            return VERT_1;
        } else if self.v2() == vertex {
            return VERT_2;
        } else {
            return VERT_3;
        }
    }
    #[inline]
    pub fn get_vertex_index(&self, vertex_id: VertexId) -> Option<TriangleVertexIndex> {
        for (v_index, v) in self.verts.iter().enumerate() {
            if *v == vertex_id {
                return Some(v_index as TriangleVertexIndex);
            }
        }
        None
    }

    /// `vertex` MUST be a vertex of the triangle
    #[inline]
    pub fn opposite_edge_index_from_vertex(&self, vertex: VertexId) -> TriangleEdgeIndex {
        opposite_edge_index(self.vertex_index(vertex))
    }

    /// `edge` MUST be an edge of the triangle
    #[inline]
    pub fn get_opposite_vertex_index(&self, edge: &Edge) -> VertexId {
        if !edge.contains(self.v1()) {
            self.v1()
        } else if !edge.contains(self.v2()) {
            self.v2()
        } else {
            self.v3()
        }
    }

    #[inline]
    pub(crate) fn has_no_container_vertex(&self, min_container_vertex_id: VertexId) -> bool {
        self.v1() < min_container_vertex_id
            && self.v2() < min_container_vertex_id
            && self.v3() < min_container_vertex_id
    }
}

#[derive(Clone)]
pub struct Triangles {
    pub buffer: Vec<TriangleData>,
}
impl Triangles {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
        }
    }

    pub fn new() -> Self {
        Self { buffer: Vec::new() }
    }

    #[inline]
    pub fn get(&self, id: TriangleId) -> &TriangleData {
        &self.buffer[id as usize]
    }
    #[inline]
    pub fn get_mut(&mut self, id: TriangleId) -> &mut TriangleData {
        &mut self.buffer[id as usize]
    }

    #[inline]
    pub fn buffer(&self) -> &Vec<TriangleData> {
        &self.buffer
    }
    #[inline]
    pub fn buffer_mut(&mut self) -> &mut Vec<TriangleData> {
        &mut self.buffer
    }

    #[inline]
    pub fn count(&self) -> usize {
        self.buffer.len()
    }
    #[inline]
    pub fn last_id(&self) -> TriangleId {
        (self.buffer.len() - 1) as TriangleId
    }
    #[inline]
    pub fn next_id(&self) -> TriangleId {
        self.buffer.len() as TriangleId
    }

    #[inline]
    pub fn create(&mut self, verts: [VertexId; 3], neighbors: [Neighbor; 3]) {
        self.buffer.push(TriangleData { verts, neighbors })
    }
    #[inline]
    pub fn push(&mut self, triangle: TriangleData) {
        self.buffer.push(triangle)
    }
}

#[inline]
pub fn is_infinite(vert_id: VertexId, min_container_vertex_id: VertexId) -> bool {
    vert_id >= min_container_vertex_id
}

#[inline]
pub fn next_clockwise_edge_index(edge_index: TriangleEdgeIndex) -> TriangleEdgeIndex {
    NEXT_CLOCKWISE_EDGE_INDEX[edge_index as usize]
}

#[inline]
pub fn next_counter_clockwise_edge_index(edge_index: TriangleEdgeIndex) -> TriangleEdgeIndex {
    NEXT_COUNTER_CLOCKWISE_EDGE_INDEX[edge_index as usize]
}

#[inline]
pub fn outermost_clockwise_edge_index_around(vert_index: TriangleVertexIndex) -> TriangleEdgeIndex {
    OUTERMOST_CW_EDGE_INDEX_AROUND_VERTEX[vert_index as usize]
}

#[inline]
pub fn outermost_counter_clockwise_edge_index_around(
    vert_index: TriangleVertexIndex,
) -> TriangleEdgeIndex {
    OUTERMOST_CCW_EDGE_INDEX_AROUND_VERTEX[vert_index as usize]
}

#[inline]
pub fn opposite_edge_index(vert_index: TriangleVertexIndex) -> TriangleEdgeIndex {
    OPPOSITE_EDGE_INDEX[vert_index as usize]
}

#[inline]
pub fn opposite_vertex_index(edge_index: TriangleEdgeIndex) -> TriangleVertexIndex {
    OPPOSITE_VERTEX_INDEX[edge_index as usize]
}

/// From a TriangleVertexIndex, returns the next edge index in a CW order
#[inline]
pub fn vertex_next_cw_edge_index(vertex_index: TriangleVertexIndex) -> TriangleEdgeIndex {
    vertex_index as TriangleEdgeIndex
}

/// From a TriangleVertexIndex, returns the next edge index in a CCW order
#[inline]
pub fn vertex_next_ccw_edge_index(vertex_index: TriangleVertexIndex) -> TriangleEdgeIndex {
    NEXT_CCW_EDGE_INDEX[vertex_index as usize]
}

#[inline]
pub fn next_clockwise_vertex_index(vertex_index: TriangleVertexIndex) -> TriangleVertexIndex {
    NEXT_CW_VERTEX_INDEX[vertex_index as usize]
}
#[inline]
pub fn next_counter_clockwise_vertex_index(
    vertex_index: TriangleVertexIndex,
) -> TriangleVertexIndex {
    NEXT_CCW_VERTEX_INDEX[vertex_index as usize]
}

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
    pub fn v(&self, vertex_index: TriangleVertexIndex) -> VertexId {
        self.verts[vertex_index as usize]
    }
    #[inline]
    pub fn v1(&self) -> VertexId {
        self.verts[QUAD_1 as usize]
    }
    #[inline]
    pub fn v2(&self) -> VertexId {
        self.verts[QUAD_2 as usize]
    }
    #[inline]
    pub fn v3(&self) -> VertexId {
        self.verts[QUAD_3 as usize]
    }
    #[inline]
    pub fn v4(&self) -> VertexId {
        self.verts[QUAD_4 as usize]
    }

    #[inline]
    pub fn to_vertices(&self, vertices: &Vec<Vertex>) -> QuadVertices {
        QuadVertices {
            verts: [
                vertices[self.v1() as usize],
                vertices[self.v2() as usize],
                vertices[self.v3() as usize],
                vertices[self.v4() as usize],
            ],
        }
    }
}

// TODO may change to structured arrays [Vec2;n]
pub type EdgeVertices = (Vertex, Vertex);
pub type TriangleVertices = (Vertex, Vertex, Vertex);

#[derive(Debug)]
pub struct QuadVertices {
    pub verts: [Vertex; 4],
}
impl QuadVertices {
    #[inline]
    pub fn q1(&self) -> Vertex {
        self.verts[QUAD_1 as usize]
    }
    #[inline]
    pub fn q2(&self) -> Vertex {
        self.verts[QUAD_2 as usize]
    }
    #[inline]
    pub fn q3(&self) -> Vertex {
        self.verts[QUAD_3 as usize]
    }
    #[inline]
    pub fn q4(&self) -> Vertex {
        self.verts[QUAD_4 as usize]
    }

    #[inline]
    pub fn diagonals_intersection_test(&self) -> EdgesIntersectionResult {
        egdes_intersect(
            &(self.verts[QUAD_1 as usize], self.verts[QUAD_2 as usize]),
            &(self.verts[QUAD_3 as usize], self.verts[QUAD_4 as usize]),
        )
    }
}
