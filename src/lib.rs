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

    use crate::{
        constrained_triangulation::ConstrainedTriangulationConfiguration,
        constrained_triangulation_from_2d_vertices,
        triangulation::TriangulationConfiguration,
        triangulation_from_2d_vertices,
        types::{Edge, Vertex},
    };

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

        assert_eq!(vec![[0, 1, 3], [1, 2, 3]], triangulation.triangles);
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
            vec![[0, 2, 3], [2, 4, 3], [1, 5, 2], [5, 4, 2]],
            triangulation.triangles
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
                [1, 7, 0],
                [3, 0, 7],
                [3, 7, 4],
                [2, 3, 4],
                [2, 4, 5],
                [6, 2, 5],
                [1, 2, 6],
                [1, 6, 7]
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

        constrained_triangulation_from_2d_vertices(
            &vertices,
            &edges,
            ConstrainedTriangulationConfiguration::default(),
        )
        .expect("Triangulation should succeed for near-colinear constrained edges");
    }
}
