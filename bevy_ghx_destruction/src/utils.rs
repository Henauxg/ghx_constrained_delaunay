use bevy::math::{Vec2, Vec3A};
use rand::Rng;

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

// source: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
pub fn egdes_intersect(p1: Vec2, q1: Vec2, p2: Vec2, q2: Vec2) -> EdgesIntersectionResult {
    if p1 == p2 || p1 == q2 || q1 == p2 || q1 == q2 {
        return EdgesIntersectionResult::SharedEdges;
    }

    let orientation_1 = triplet_orientation(p1, q1, p2);
    let orientation_2 = triplet_orientation(p1, q1, q2);
    let orientation_3 = triplet_orientation(p2, q2, p1);
    let orientation_4 = triplet_orientation(p2, q2, q1);

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if orientation_1 == Orientation::Colinear && on_segment(p1, p2, q1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if orientation_2 == Orientation::Colinear && on_segment(p1, q2, q1) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if orientation_3 == Orientation::Colinear && on_segment(p2, p1, q2) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if orientation_4 == Orientation::Colinear && on_segment(p2, q1, q2) {
        return EdgesIntersectionResult::OnEdgeTip;
    }

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
