use crate::types::{EdgeVertices, Float, Vertex, VertexId};

#[cfg(feature = "progress_log")]
use tracing::info;

#[cfg(feature = "more_profile_traces")]
use tracing::{span, Level};

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum EdgesIntersectionResult {
    None,
    Crossing,
    OnEdgeTip,
    SharedEdges,
}

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum Orientation {
    Colinear,
    Clockwise,
    CounterClockwise,
}

#[derive(thiserror::Error, Debug, Eq, PartialEq)]
/// Returned when an input vertex is NaN or infinity
#[error("an input vertex is NaN or infinity")]
pub struct InvalidVertex;

/// Used to validate a vertices collection
pub fn validate_vertices(vertices: &Vec<Vertex>) -> Result<(), InvalidVertex> {
    for v in vertices.iter() {
        if v.is_nan() || !v.is_finite() {
            return Err(InvalidVertex);
        }
    }
    Ok(())
}

// source: https://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
pub fn egdes_intersect(edge_1: &EdgeVertices, edge_2: &EdgeVertices) -> EdgesIntersectionResult {
    if edge_1.0 == edge_2.0 || edge_1.0 == edge_2.1 || edge_1.1 == edge_2.0 || edge_1.1 == edge_2.1
    {
        return EdgesIntersectionResult::SharedEdges;
    }

    let orientation_1 = triplet_orientation(edge_1.0, edge_1.1, edge_2.0);
    let orientation_2 = triplet_orientation(edge_1.0, edge_1.1, edge_2.1);
    let orientation_3 = triplet_orientation(edge_2.0, edge_2.1, edge_1.0);
    let orientation_4 = triplet_orientation(edge_2.0, edge_2.1, edge_1.1);

    // Special Cases
    // edge_1.0, edge_1.1 and edge_2.0 are collinear and edge_2.0 lies on segment edge_1.0edge_1.1
    if orientation_1 == Orientation::Colinear && on_segment(edge_1.0, edge_2.0, edge_1.1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // edge_1.0, edge_1.1 and edge_2.1 are collinear and edge_2.1 lies on segment edge_1.0edge_1.1
    if orientation_2 == Orientation::Colinear && on_segment(edge_1.0, edge_2.1, edge_1.1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // edge_2.0, edge_2.1 and edge_1.0 are collinear and edge_1.0 lies on segment edge_2.0edge_2.1
    if orientation_3 == Orientation::Colinear && on_segment(edge_2.0, edge_1.0, edge_2.1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // edge_2.0, edge_2.1 and edge_1.1 are collinear and edge_1.1 lies on segment edge_2.0edge_2.1
    if orientation_4 == Orientation::Colinear && on_segment(edge_2.0, edge_1.1, edge_2.1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // General case, should be checked after the OnEdgeTip cases
    if orientation_1 != orientation_2 && orientation_3 != orientation_4 {
        return EdgesIntersectionResult::Crossing;
    }

    EdgesIntersectionResult::None // Doesn't fall in any of the above cases
}

/// Returns the orientation of an ordered triplet (p, q, r).
pub fn triplet_orientation(p: Vertex, q: Vertex, r: Vertex) -> Orientation {
    let val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

    if val == 0. {
        Orientation::Colinear
    } else if val > 0. {
        Orientation::Clockwise
    } else {
        Orientation::CounterClockwise
    }
}

#[derive(Debug)]
pub struct EdgeSideTestResult(pub Float);
impl EdgeSideTestResult {
    #[inline]
    pub fn is_on_right_side(&self) -> bool {
        self.0 < 0.
    }
    #[inline]
    pub fn is_strictly_on_right_side(&self) -> bool {
        self.0 < -Float::EPSILON
    }
    #[inline]
    pub fn is_on_left_side(&self) -> bool {
        self.0 > 0.
    }
    #[inline]
    pub fn is_strictly_on_left_side(&self) -> bool {
        self.0 > Float::EPSILON
    }
    #[inline]
    pub fn is_colinear(&self) -> bool {
        self.0 == 0.
    }
    #[inline]
    pub fn is_near_edge(&self) -> bool {
        self.0.abs() < Float::EPSILON
    }
}

/// Returns the position of the point `p` when compared to the oriented edge `e`
///
/// Uses the perp cross product of vectors e0.e1 and e1.p
#[inline]
pub fn test_point_edge_side(e: EdgeVertices, p: Vertex) -> EdgeSideTestResult {
    EdgeSideTestResult((p.y - e.0.y) * (e.1.x - e.0.x) - (p.x - e.0.x) * (e.1.y - e.0.y))
}

/// Returns `true` if and only if the point `p` is on the right side of the oriented edge `e`
///
/// Uses the perp cross product of vectors e0.e1 and e1.p
#[inline]
pub fn is_point_on_right_side_of_edge(e: EdgeVertices, p: Vertex) -> bool {
    (p.y - e.0.y) * (e.1.x - e.0.x) - (p.x - e.0.x) * (e.1.y - e.0.y) <= 0.
}

/// Returns `true` if and only if the point `p` is strictly on the right side of the oriented edge `e`
///
/// Uses the perp cross product of vectors e0.e1 and e1.p
#[inline]
pub fn is_point_strictly_on_right_side_of_edge(e: EdgeVertices, p: Vertex) -> bool {
    (p.y - e.0.y) * (e.1.x - e.0.x) - (p.x - e.0.x) * (e.1.y - e.0.y) <= -Float::EPSILON
}

/// Returns the slope of the line going from `a` to `b`
///
/// - a and b MUST not share the same x coordinate
#[inline]
pub fn line_slope(a: Vertex, b: Vertex) -> Float {
    (b.y - a.y) / (b.x - a.x)
}

/// Given three collinear points p, q, r, the function checks if point `q` lies on line segment 'pr'
#[inline]
pub fn on_segment(p: Vertex, q: Vertex, r: Vertex) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
}

/// Checks if vertex `p` is inside the circumcircle of the triangle formed by the first three vertices in `triangle`
/// - `triangle` vertices must be in a counter-clockwise order
/// - `p` vertex to check
///
/// ```text
/// v3 --------- v2
/// |          / |
/// |        /   |
/// |      /     |
/// |    /       |
/// |  /         |
/// v1 --------- p
/// ```
///
/// where v1, v2 and v3 are the vertices of the given triangle and p the vertex to check
///
/// See: A. K. Cline and R. Renka,
/// A storage efficient method for construction of a Thiessen triangulation.
/// Rocky Mounfain J. Math. 14, 119-139 (1984)
///
/// Note: Seems to return false for a flat triangle
#[inline(always)]
pub(crate) fn is_vertex_in_triangle_circumcircle(
    t1: Vertex,
    t2: Vertex,
    t3: Vertex,
    p: Vertex,
) -> bool {
    #[cfg(feature = "more_profile_traces")]
    let _span = span!(Level::TRACE, "is_vertex_in_triangle_circumcircle").entered();

    let x13 = t1.x - t3.x;
    let x23 = t2.x - t3.x;
    let y13 = t1.y - t3.y;
    let y23 = t2.y - t3.y;
    let x14 = t1.x - p.x;
    let x24 = t2.x - p.x;
    let y14 = t1.y - p.y;
    let y24 = t2.y - p.y;

    let cos_a = x13 * x23 + y13 * y23;
    let cos_b = x24 * x14 + y24 * y14;

    if cos_a >= 0. && cos_b >= 0. {
        false
    } else if cos_a < 0. && cos_b < 0. {
        true
    } else {
        let sin_a = x13 * y23 - x23 * y13;
        let sin_b = x24 * y14 - x14 * y24;
        let sin_ab = sin_a * cos_b + sin_b * cos_a;
        sin_ab < 0.
    }
}

#[derive(Debug)]
pub struct DegenerateTrianglesInfo {
    pub flat_triangles_count: u64,
    pub ccw_triangles_count: u64,
}
impl DegenerateTrianglesInfo {
    fn new(flat_triangles_count: u64, ccw_triangles_count: u64) -> Self {
        Self {
            flat_triangles_count,
            ccw_triangles_count,
        }
    }
}

/// Check degenerate (flat) triangles
pub fn check_degenerate_triangles(
    triangles: impl IntoIterator<Item = [VertexId; 3]>,
    vertices: &Vec<Vertex>,
) -> DegenerateTrianglesInfo {
    let mut flat_triangles_count = 0;
    let mut ccw_triangles_count = 0;

    for t in triangles.into_iter() {
        let (v1, v2, v3) = (
            vertices[t[0] as usize],
            vertices[t[1] as usize],
            vertices[t[2] as usize],
        );
        let orientation = triplet_orientation(v1, v2, v3);
        match orientation {
            Orientation::Colinear => flat_triangles_count += 1,
            Orientation::CounterClockwise => ccw_triangles_count += 1,
            Orientation::Clockwise => (),
        }
    }
    DegenerateTrianglesInfo::new(flat_triangles_count, ccw_triangles_count)
}

#[derive(Debug)]
pub struct CircumcirclesQualityInfo {
    pub non_optimal_triangles_count: usize,
}
impl CircumcirclesQualityInfo {
    fn new(non_optimal_triangles_count: usize) -> Self {
        Self {
            non_optimal_triangles_count,
        }
    }
}
pub fn get_circumcircle_info(verts: &[Vertex; 3]) -> (Vertex, Float) {
    let (v1, v2, v3) = (verts[0], verts[1], verts[2]);
    let b = v2 - v1;
    let c = v3 - v1;
    let d = 2. * (b.x * c.y - b.y * c.x);
    let u = Vertex::new(
        (c.y * (b.x * b.x + b.y * b.y) - b.y * (c.x * c.x + c.y * c.y)) / d,
        (b.x * (c.x * c.x + c.y * c.y) - c.x * (b.x * b.x + b.y * b.y)) / d,
    );
    let center = u + v1;
    let radius = (u.x * u.x + u.y * u.y).sqrt();
    (center, radius)
}

/// Checks that all triangles in the triangulation have a circumcircle which does not contain any other vertices from the triangulation.
///
/// TODO Results may not be true if used on infinite vertices DURING the triangulation. After the triangulation, infintie vertices are always remove,d so this is not an issue.
pub fn check_circumcircles(
    triangles: impl IntoIterator<Item = [VertexId; 3]> + Clone,
    vertices: &Vec<Vertex>,
    _progress_log: bool,
) -> CircumcirclesQualityInfo {
    let triangles1 = triangles.clone().into_iter();
    let _triangle_count = triangles1.size_hint().0;
    let mut non_optimal_triangles_count = 0;
    for (_t_id, t) in triangles1.enumerate() {
        #[cfg(feature = "progress_log")]
        {
            if _progress_log && _t_id % ((_triangle_count / 50) + 1) == 0 {
                let progress = 100. * _t_id as f32 / _triangle_count as f32;
                info!(
                    "check_circumcircles progress, {}%: {}/{}",
                    progress, _t_id, _triangle_count
                );
            }
        }

        let verts = [
            vertices[t[0] as usize],
            vertices[t[1] as usize],
            vertices[t[2] as usize],
        ];

        let circumcircle = get_circumcircle_info(&verts);

        let mut found_contained_vertex = false;
        for other_t_verts in triangles.clone().into_iter() {
            for &v_id in other_t_verts.iter() {
                if v_id == t[0] || v_id == t[1] || v_id == t[2] {
                    continue;
                }
                let v = vertices[v_id as usize];
                // We cannot use is_vertex_in_triangle_circumcircle here since its results depends highly on the input vertices orientation.
                let dist = (v - circumcircle.0).length();
                if dist < circumcircle.1 {
                    non_optimal_triangles_count += 1;
                    found_contained_vertex = true;
                    break;
                }
            }
            if found_contained_vertex {
                break;
            }
        }
    }

    CircumcirclesQualityInfo::new(non_optimal_triangles_count)
}

#[derive(Debug)]
pub struct DelaunayQualityInfo {
    pub degen_triangles: DegenerateTrianglesInfo,
    pub circumcircles: CircumcirclesQualityInfo,
}
impl DelaunayQualityInfo {
    fn new(
        degen_triangles: DegenerateTrianglesInfo,
        circumcircles: CircumcirclesQualityInfo,
    ) -> Self {
        Self {
            degen_triangles,
            circumcircles,
        }
    }
}

/// Slow and simple check of the Delaunay optimality of a triangulation
pub fn check_delaunay_optimal<T>(
    triangles: T,
    vertices: &Vec<Vertex>,
    progress_log: bool,
) -> DelaunayQualityInfo
where
    T: IntoIterator<Item = [VertexId; 3]> + Clone,
{
    let degenerate_triangles_info = check_degenerate_triangles(triangles.clone(), vertices);
    let circumcircles_info = check_circumcircles(triangles, vertices, progress_log);

    DelaunayQualityInfo::new(degenerate_triangles_info, circumcircles_info)
}

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    use crate::{
        types::Vertex,
        utils::{
            egdes_intersect, is_point_on_right_side_of_edge, on_segment, EdgesIntersectionResult,
        },
    };

    #[test]
    fn edge_intersect_none() {
        let edge1 = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let edge2 = (Vertex::new(0., 3.), Vertex::new(3., 3.));

        let intersect = egdes_intersect(&edge1, &edge2);

        assert_eq!(EdgesIntersectionResult::None, intersect);
    }

    #[test]
    fn edge_intersect_crossing() {
        let edge1 = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let edge2 = (Vertex::new(1., -2.), Vertex::new(1., 3.));

        let intersect = egdes_intersect(&edge1, &edge2);

        assert_eq!(EdgesIntersectionResult::Crossing, intersect);
    }

    #[test]
    fn edge_intersect_on_edge_tip() {
        let edge1 = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let edge2 = (Vertex::new(1., 0.), Vertex::new(3., 3.));

        let intersect = egdes_intersect(&edge1, &edge2);

        assert_eq!(EdgesIntersectionResult::OnEdgeTip, intersect);
    }

    #[test]
    fn edge_intersect_shared_edge() {
        let edge1 = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let edge2 = (Vertex::new(0., 0.), Vertex::new(3., 3.));

        let intersect = egdes_intersect(&edge1, &edge2);

        assert_eq!(EdgesIntersectionResult::SharedEdges, intersect);
    }

    #[test]
    fn point_edge_orientatiion_left() {
        let edge = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let p = Vertex::new(0., 3.);

        let orientation = is_point_on_right_side_of_edge(edge, p);

        assert_eq!(false, orientation);
    }

    #[test]
    fn point_edge_orientation_right() {
        let edge = (Vertex::new(0., 0.), Vertex::new(3., 0.));
        let p = Vertex::new(0., -3.);

        let orientation = is_point_on_right_side_of_edge(edge, p);

        assert_eq!(true, orientation);
    }

    #[test]
    fn point_inside_segment() {
        let p = Vertex::new(0., 3.);
        let r = Vertex::new(0., -3.);
        let q = Vertex::new(0., 0.);

        let orientation = on_segment(p, q, r);

        assert_eq!(true, orientation);
    }

    #[test]
    fn point_outside_segment() {
        let p = Vertex::new(0., 3.);
        let r = Vertex::new(0., -3.);
        let q = Vertex::new(8., 0.);

        let orientation = on_segment(p, q, r);

        assert_eq!(false, orientation);
    }
}
