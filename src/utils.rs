use crate::types::{EdgeVertices, Vertex};

#[cfg(feature = "profile_traces")]
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

    // General case
    if orientation_1 != orientation_2 && orientation_3 != orientation_4 {
        return EdgesIntersectionResult::Crossing;
    }

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

/// Returns `true` if and only if the point `p` is on the right side of the oriented edge `e`
///
/// Uses the ross product of vectors e0.e1 and e1.p
#[inline]
pub fn is_point_on_right_side_of_edge(e: EdgeVertices, p: Vertex) -> bool {
    ((p.x - e.0.x) * (e.1.y - e.0.y) - (p.y - e.0.y) * (e.1.x - e.0.x)) >= 0.
}

/// Given three collinear points p, q, r, the function checks if point `q` lies on line segment 'pr'
#[inline]
pub fn on_segment(p: Vertex, q: Vertex, r: Vertex) -> bool {
    q.x <= p.x.max(r.x) && q.x >= p.x.min(r.x) && q.y <= p.y.max(r.y) && q.y >= p.y.min(r.y)
}

/// Checks if vertex `p` is inside the circumcircle of the triangle formed by the first three vertices in `triangle`
/// - `triangle` contains the vertices of the triangle.
///     - length of `triangle` **MUST** be >= 3.
///     - `triangle` vertices must be in a counter-clockwise order
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
#[inline(always)]
pub(crate) fn is_vertex_in_triangle_circumcircle(triangle: &[Vertex], p: Vertex) -> bool {
    #[cfg(feature = "profile_traces")]
    let _span = span!(Level::TRACE, "is_vertex_in_triangle_circumcircle").entered();

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

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    use crate::{
        types::{Float, Vertex},
        utils::{
            egdes_intersect, is_point_on_right_side_of_edge, is_vertex_in_triangle_circumcircle,
            on_segment, EdgesIntersectionResult,
        },
    };

    #[test]
    fn vertex_in_triangle_circumcircle() {
        let unit_circle = [
            Vertex::new(-1., 0.),
            Vertex::new(1., 0.),
            Vertex::new(0., 1.),
        ];

        let step = 100;
        for i in -step..step {
            for j in -step..step {
                let p = Vertex::new(i as Float / step as Float, j as Float / step as Float);
                let p_length = p.length();
                let p_in_circle = is_vertex_in_triangle_circumcircle(&unit_circle, p);
                if p_length < 1. {
                    assert_eq!(true, p_in_circle, "p_length < 1, p should be in the circle");
                } else if p_length > 1. {
                    assert_eq!(
                        false, p_in_circle,
                        "p_length > 1, p should be out of the circle"
                    );
                }
            }
        }
    }

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
