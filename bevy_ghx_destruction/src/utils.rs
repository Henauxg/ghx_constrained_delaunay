use bevy::math::{Vec2, Vec3A};
use rand::Rng;

use crate::triangulation::EdgeVertices;

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

    // TODO General case should be the fastest one to be compted and thus should be at the top.
    // General case
    if orientation_1 != orientation_2 && orientation_3 != orientation_4 {
        return EdgesIntersectionResult::Crossing;
    }

    EdgesIntersectionResult::None // Doesn't fall in any of the above cases
}

/// Returns the orientation of an ordered triplet (p, q, r).
pub fn triplet_orientation(p: Vec2, q: Vec2, r: Vec2) -> Orientation {
    let val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

    if val == 0. {
        Orientation::Colinear
    } else if val > 0. {
        Orientation::Clockwise
    } else {
        Orientation::CounterClockwise
    }
}

#[inline]
pub fn get_random_normalized_vec() -> Vec3A {
    let mut rng = rand::thread_rng();
    Vec3A::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()).normalize()
}

/// Cross product of vectors e0e1 and e1p
#[inline]
pub fn is_point_on_right_side_of_edge(e0: Vec2, e1: Vec2, p: Vec2) -> bool {
    ((p.x - e0.x) * (e1.y - e0.y) - (p.y - e0.y) * (e1.x - e0.x)) >= 0.
}

/// Given three collinear points p, q, r, the function checks if point `q` lies on line segment 'pr'
#[inline]
pub fn on_segment(p: Vec2, q: Vec2, r: Vec2) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
}

/// Cheks if vertex `p` is inside the circumcircle of the triangle formed by the first three vertices in `triangle`
/// - `triangle` are the vertices of the triangle.
///     - length of `triangle` **MUST** be >= 3.
///     - `triangle` vertices must be in a counter-clockwise order
/// - `p` vertex to check
///
/// v3 --------- v2
/// |          / |
/// |        /   |
/// |      /     |
/// |    /       |
/// |  /         |
/// v1 --------- p
///
/// where v1, v2 and v3 are the vertices of the given triangle and p the vertex to check
///
/// See: A. K. Cline and R. Renka,
/// A storage efficient method for construction of a Thiessen triangulation.
/// Rocky Mounfain J. Math. 14, 119-139 (1984)
///
pub fn is_vertex_in_triangle_circumcircle(triangle: &[Vec2], p: Vec2) -> bool {
    let x13 = triangle[0].x - triangle[2].x;
    let x23 = triangle[1].x - triangle[2].x;
    let y13 = triangle[0].y - triangle[2].y;
    let y23 = triangle[1].y - triangle[2].y;
    let x14 = triangle[0].x - p.x;
    let x24 = triangle[1].x - p.x;
    let y14 = triangle[0].y - p.y;
    let y24 = triangle[1].y - p.y;

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
