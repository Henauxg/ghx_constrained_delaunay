pub mod constrained_triangulation;
pub mod infinite;
pub mod triangulation;
pub mod types;
pub mod utils;

#[cfg(feature = "debug_context")]
pub mod debug;

pub use glam;
pub use hashbrown;

pub use constrained_triangulation::{
    constrained_triangulation_from_2d_vertices, constrained_triangulation_from_3d_planar_vertices,
};
pub use triangulation::{
    triangulation_from_2d_vertices, triangulation_from_3d_planar_vertices, Triangulation,
};

///////////////////////////////////////////////////////////
///                                                     ///
///                        Tests                        ///
///                                                     ///
///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    use glam::DVec3;

    use crate::{
        constrained_triangulation::ConstrainedTriangulationConfiguration,
        constrained_triangulation_from_2d_vertices,
        triangulation::{TriangleOrientation, TriangulationConfiguration},
        triangulation_from_2d_vertices, triangulation_from_3d_planar_vertices,
        types::{Edge, Float, Vertex, VertexId},
    };

    fn signed_area(vertices: &[Vertex], t: &[VertexId; 3]) -> Float {
        let (a, b, c) = (
            vertices[t[0] as usize],
            vertices[t[1] as usize],
            vertices[t[2] as usize],
        );
        0.5 * (b - a).perp_dot(c - a)
    }

    fn reversed(triangles: &[[VertexId; 3]]) -> Vec<[VertexId; 3]> {
        triangles.iter().map(|t| [t[0], t[2], t[1]]).collect()
    }

    fn dt_config(triangle_orientation: TriangleOrientation) -> TriangulationConfiguration {
        TriangulationConfiguration {
            triangle_orientation,
            // Also covers the parallel filtering path
            #[cfg(feature = "parallel_filtering")]
            filter_parallel_tri_count_threshold: 0,
            #[cfg(feature = "parallel_filtering")]
            filter_parallel_min_batch_len: TriangulationConfiguration::default()
                .filter_parallel_min_batch_len,
            #[cfg(feature = "debug_context")]
            debug_config: Default::default(),
        }
    }

    #[test]
    fn delaunay_counter_clockwise() {
        let vertices: Vec<Vertex> = (0..100)
            .map(|i| {
                Vertex::new(
                    (i % 10) as Float,
                    (i / 10) as Float + 0.01 * (i % 7) as Float,
                )
            })
            .collect();

        let cw =
            triangulation_from_2d_vertices(&vertices, dt_config(TriangleOrientation::Clockwise))
                .expect("Triangulation should succeed");
        let ccw = triangulation_from_2d_vertices(
            &vertices,
            dt_config(TriangleOrientation::CounterClockwise),
        )
        .expect("Triangulation should succeed");

        assert!(cw.triangles.iter().all(|t| signed_area(&vertices, t) < 0.));
        assert!(ccw.triangles.iter().all(|t| signed_area(&vertices, t) > 0.));
        assert_eq!(reversed(&cw.triangles), ccw.triangles);
    }

    #[test]
    fn constrained_delaunay_counter_clockwise() {
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 5.),
            Vertex::new(5., 5.),
            Vertex::new(5., 0.),
            Vertex::new(4., 1.),
            Vertex::new(4., 4.),
            Vertex::new(1., 4.),
            Vertex::new(1., 1.),
        ];
        let constrained_edges = vec![
            Edge::new(0, 1),
            Edge::new(1, 2),
            Edge::new(2, 3),
            Edge::new(3, 0),
            Edge::new(4, 5),
            Edge::new(5, 6),
            Edge::new(6, 7),
            Edge::new(7, 4),
        ];
        let triangulate = |triangle_orientation| {
            constrained_triangulation_from_2d_vertices(
                &vertices,
                &constrained_edges,
                ConstrainedTriangulationConfiguration {
                    triangle_orientation,
                    #[cfg(feature = "debug_context")]
                    debug_config: Default::default(),
                },
            )
            .expect("Triangulation should succeed")
        };

        let cw = triangulate(TriangleOrientation::Clockwise);
        let ccw = triangulate(TriangleOrientation::CounterClockwise);

        assert_eq!(8, ccw.triangles.len());
        assert!(ccw.triangles.iter().all(|t| signed_area(&vertices, t) > 0.));
        assert_eq!(reversed(&cw.triangles), ccw.triangles);
    }

    #[test]
    fn orientation_3d_is_seen_from_the_normal_side() {
        let vertices = [
            DVec3::new(0., 0., 0.),
            DVec3::new(0., 5., 0.),
            DVec3::new(5., 5., 0.),
            DVec3::new(5., 0., 0.),
            DVec3::new(2., 3., 0.),
        ];
        for normal in [DVec3::Z, DVec3::NEG_Z] {
            for orientation in [
                TriangleOrientation::Clockwise,
                TriangleOrientation::CounterClockwise,
            ] {
                let triangulation = triangulation_from_3d_planar_vertices(
                    &vertices,
                    normal,
                    dt_config(orientation),
                )
                .expect("Triangulation should succeed");

                assert_eq!(4, triangulation.triangles.len());
                for t in &triangulation.triangles {
                    let (a, b, c) = (
                        vertices[t[0] as usize],
                        vertices[t[1] as usize],
                        vertices[t[2] as usize],
                    );
                    let side = (b - a).cross(c - a).dot(normal);
                    match orientation {
                        TriangleOrientation::Clockwise => assert!(side < 0., "{normal} {t:?}"),
                        TriangleOrientation::CounterClockwise => {
                            assert!(side > 0., "{normal} {t:?}")
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn delaunay_level_1() {
        // 1-------------2
        // |  \          |
        // |     \       |
        // |        \    |
        // |           \ |
        // 0-------------3
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 5.),
            Vertex::new(5., 5.),
            Vertex::new(5., 0.),
        ];

        let triangulation =
            triangulation_from_2d_vertices(&vertices, TriangulationConfiguration::default())
                .expect("Triangulation should succeed");

        assert_eq!(vec![[2, 3, 1], [0, 1, 3]], triangulation.triangles);
    }

    #[test]
    fn delaunay_level_2() {
        // 1-------------5
        // |           / |
        // |        /    |
        // |     /       |
        // |  /          |
        // 2-------------4
        // |  \          |
        // |     \       |
        // |        \    |
        // |           \ |
        // 0-------------3
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 10.),
            Vertex::new(0., 5.),
            Vertex::new(5., 0.),
            Vertex::new(5., 5.),
            Vertex::new(5., 10.),
        ];

        let triangulation =
            triangulation_from_2d_vertices(&vertices, TriangulationConfiguration::default())
                .expect("Triangulation should succeed");

        assert_eq!(
            vec![[2, 5, 4], [1, 5, 2], [0, 2, 4], [3, 0, 4]],
            triangulation.triangles
        );
    }

    #[test]
    fn delaunay_hull_is_convex() {
        // A small island, then a far point above it. Vertex 4 is the highest point of the island, but ends up
        // strictly inside the convex hull once vertex 5 is added.
        let vertices = vec![
            Vertex::new(-128.3, -24.41259765625),
            Vertex::new(-128.290087890625, -24.39736328125001),
            Vertex::new(-128.2908203125, -24.36464843750001),
            Vertex::new(-128.30361328125, -24.333593750000006),
            Vertex::new(-128.330126953125, -24.3232421875),
            Vertex::new(-177.3541015625, 28.217041015625),
        ];

        let triangulation =
            triangulation_from_2d_vertices(&vertices, TriangulationConfiguration::default())
                .expect("Triangulation should succeed");

        // 6 vertices, 5 of them on the convex hull: 2 * 6 - 2 - 5 triangles
        assert_eq!(5, triangulation.triangles.len());
        for v in 0..6 {
            assert!(triangulation.triangles.iter().any(|t| t.contains(&v)));
        }
    }

    #[test]
    fn delaunay_merges_close_vertices() {
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 1.),
            Vertex::new(1., 1.),
            Vertex::new(1., 0.),
            Vertex::new(0.5, 0.5),
            // Too close to (0.5, 0.5), and not on a horizontal or vertical line through it
            Vertex::new(0.5 + Float::EPSILON / 2., 0.5 + Float::EPSILON / 2.),
        ];

        let triangulation =
            triangulation_from_2d_vertices(&vertices, TriangulationConfiguration::default())
                .expect("Triangulation should succeed");

        assert_eq!(4, triangulation.triangles.len());
        let uses = |v| triangulation.triangles.iter().any(|t| t.contains(&v));
        assert!(
            uses(4) != uses(5),
            "Exactly one of the two close vertices is kept"
        );
    }

    #[test]
    fn constrained_delaunay_level_1() {
        // 1-------------2
        // |          /  |
        // |       /     |
        // |    /        |
        // | /           |
        // 0-------------3
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 5.),
            Vertex::new(5., 5.),
            Vertex::new(5., 0.),
        ];
        let edges = vec![Edge::new(0, 2)];

        let triangulation = constrained_triangulation_from_2d_vertices(
            &vertices,
            &edges,
            ConstrainedTriangulationConfiguration::default(),
        )
        .expect("Triangulation should succeed");

        assert_eq!(vec![[2, 3, 0], [2, 0, 1]], triangulation.triangles);
    }

    #[test]
    fn constrained_delaunay_level_2() {
        // 1-------------2
        // |          /  |
        // |       /     |
        // |    /        |
        // | /           |
        // 0-------------3
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 5.),
            Vertex::new(5., 5.),
            Vertex::new(5., 0.),
        ];
        let edges = vec![Edge::new(0, 2)];

        let triangulation = constrained_triangulation_from_2d_vertices(
            &vertices,
            &edges,
            ConstrainedTriangulationConfiguration::default(),
        )
        .expect("Triangulation should succeed");

        assert_eq!(vec![[2, 3, 0], [2, 0, 1]], triangulation.triangles);
    }

    #[test]
    fn constrained_delaunay_level_3() {
        // ```text
        // 1------>-------------2
        // |    keep (CW)       |
        // |  6----<---------5  |
        // |  | remove (CCW) |  |
        // ^  v              ^  v
        // |  7-->-----------4  |
        // 0--------<-----------3
        // ```
        let vertices = vec![
            Vertex::new(0., 0.),
            Vertex::new(0., 5.),
            Vertex::new(5., 5.),
            Vertex::new(5., 0.),
            Vertex::new(4., 1.),
            Vertex::new(4., 4.),
            Vertex::new(1., 4.),
            Vertex::new(1., 1.),
        ];
        let constrained_edges = vec![
            // CCW
            Edge::new(0, 1),
            Edge::new(1, 2),
            Edge::new(2, 3),
            Edge::new(3, 0),
            // CW
            Edge::new(4, 5),
            Edge::new(5, 6),
            Edge::new(6, 7),
            Edge::new(7, 4),
        ];

        let triangulation = constrained_triangulation_from_2d_vertices(
            &vertices,
            &constrained_edges,
            ConstrainedTriangulationConfiguration::default(),
        )
        .expect("Triangulation should succeed");

        assert_eq!(
            vec![
                [2, 6, 1],
                [1, 6, 7],
                [0, 1, 7],
                [0, 7, 3],
                [4, 3, 7],
                [4, 2, 3],
                [4, 5, 2],
                [2, 5, 6]
            ],
            triangulation.triangles
        );
    }

    #[test]
    fn constrained_delaunay_near_colinear_edges() {
        let vertices = vec![
            Vertex::new(-0.0782132158260142, 0.0),
            Vertex::new(-0.04346227742441594, 0.10349781560157462),
            Vertex::new(0.0388118494777121, 0.2521774205252516),
            Vertex::new(0.10736492278480181, 0.3364695064078743),
            Vertex::new(-0.03490517918394376, 0.12472698060408025),
            Vertex::new(0.021085609120375826, 0.22442779757191253),
            Vertex::new(0.12874905992983615, 0.3600245316864931),
            Vertex::new(0.1901120866083647, 0.4212073366884178),
            Vertex::new(-0.1901120866083661, 0.0),
            Vertex::new(-0.1901120866083661, 2.0),
            Vertex::new(0.19011208660836504, 0.0),
            Vertex::new(0.19011208660836504, 2.0),
        ];

        let edges = vec![
            Edge::new(11, 7),
            Edge::new(8, 9),
            Edge::new(7, 10),
            Edge::new(9, 11),
            Edge::new(4, 1),
            Edge::new(1, 0),
            Edge::new(0, 8),
            Edge::new(6, 3),
            Edge::new(2, 5),
            Edge::new(3, 2),
            Edge::new(5, 4),
            Edge::new(6, 7),
            Edge::new(10, 0),
        ];

        let bounding_box_area: Float = (0.19011208660836504 + 0.1901120866083661) * 2.0;

        // Flips and rotations move the flat triangles next to other infinite vertices.
        // The last 4 are mirrors.
        let transforms: [fn(Vertex) -> Vertex; 8] = [
            |v| v,
            |v| Vertex::new(-v.y, v.x),
            |v| Vertex::new(-v.x, -v.y),
            |v| Vertex::new(v.y, -v.x),
            |v| Vertex::new(-v.x, v.y),
            |v| Vertex::new(v.x, -v.y),
            |v| Vertex::new(v.y, v.x),
            |v| Vertex::new(-v.y, -v.x),
        ];
        for (index, transform) in transforms.iter().enumerate() {
            let vertices: Vec<Vertex> = vertices.iter().map(|v| transform(*v)).collect();
            // A mirror turns CW contours into CCW ones, which would discard the domain.
            let edges: Vec<Edge> = edges
                .iter()
                .map(|e| {
                    if index >= 4 {
                        Edge::new(e.to, e.from)
                    } else {
                        *e
                    }
                })
                .collect();
            assert_valid_cdt(&vertices, &edges, bounding_box_area, index);
        }
    }

    fn assert_valid_cdt(vertices: &[Vertex], edges: &[Edge], expected_area: Float, case: usize) {
        let triangulation = constrained_triangulation_from_2d_vertices(
            vertices,
            edges,
            ConstrainedTriangulationConfiguration::default(),
        )
        .unwrap_or_else(|e| panic!("Case {case}: triangulation failed: {e:?}"));

        for edge in edges {
            assert!(
                triangulation
                    .triangles
                    .iter()
                    .any(|t| t.contains(&edge.from) && t.contains(&edge.to)),
                "Case {case}: constrained edge {edge:?} is missing"
            );
        }
        let mut area = 0.;
        for t in &triangulation.triangles {
            let (a, b, c) = (
                vertices[t[0] as usize],
                vertices[t[1] as usize],
                vertices[t[2] as usize],
            );
            let signed_area = 0.5 * (b - a).perp_dot(c - a);
            assert!(
                signed_area < 0.,
                "Case {case}: triangle {t:?} is not clockwise"
            );
            area -= signed_area;
        }
        assert!(
            (area - expected_area).abs() < 1e-6 * expected_area,
            "Case {case}: area {area} instead of {expected_area}"
        );
    }
}
