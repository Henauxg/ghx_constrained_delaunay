pub mod constrained_triangulation;
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
    use glam::Vec3;

    use crate::{
        constrained_triangulation::{
            constrained_triangulation_from_3d_planar_vertices,
            ConstrainedTriangulationConfiguration,
        },
        triangulation::{triangulation_from_3d_planar_vertices, TriangulationConfiguration},
        types::{Edge, Float},
    };

    #[test]
    fn delaunay_level_1() {
        // ---------------
        // |  \          |
        // |     \       |
        // |        \    |
        // |           \ |
        // ---------------
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([0., 0., 0.]);
        vertices.push([0., 5., 0.]);
        vertices.push([5., 5., 0.]);
        vertices.push([5., 0., 0.]);

        let plane_normal = Vec3::Z;
        let triangulation = triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            TriangulationConfiguration::default(),
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in triangulation.vert_indices {
            constrained_vertices.push(vertices[vertexid as usize]);
        }

        assert_eq!(6, constrained_vertices.len());
        assert_eq!(
            Vec::from([
                [0.0, 5.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.0, 5.0, 0.0],
                [0.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [0.0, 5.0, 0.0]
            ]),
            constrained_vertices
        );
    }

    #[test]
    fn delaunay_level_2() {
        // ---------------
        // |           / |
        // |        /    |
        // |     /       |
        // |  /          |
        // ---------------
        // |  \          |
        // |     \       |
        // |        \    |
        // |           \ |
        // ---------------
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([0., 0., 0.]);
        vertices.push([0., 10., 0.]);
        vertices.push([0., 5., 0.]);
        vertices.push([5., 0., 0.]);
        vertices.push([5., 5., 0.]);
        vertices.push([5., 10., 0.]);

        let plane_normal = Vec3::Z;
        let triangulation = triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            TriangulationConfiguration::default(),
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in triangulation.vert_indices {
            constrained_vertices.push(vertices[vertexid as usize]);
        }

        assert_eq!(12, constrained_vertices.len());
        assert_eq!(
            Vec::from([
                [5.0, 0.0, 0.0],
                [5.0, 5.0, 0.0],
                [0.0, 5.0, 0.0],
                [0.0, 5.0, 0.0],
                [5.0, 5.0, 0.0],
                [5.0, 10.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 5.0, 0.0],
                [5.0, 10.0, 0.0],
                [0.0, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [0.0, 5.0, 0.0]
            ]),
            constrained_vertices
        );
    }

    #[test]
    fn constrained_delaunay_level_3() {
        // ---------------<-------------------------------------
        // |                        keep                       |
        // |   ----------->---------------------------------   |
        // |   |                   remove                  |   |
        // |   |   -------<--------     -------<--------   |   |
        // |   |   |              |     |              |   |   |
        // v   ∧   v     keep     ∧     v     keep     ∧   v   ∧
        // |   |   |              |     |              |   |   |
        // |   |   ------->--------     ------->--------   |   |
        // |   |                                           |   |
        // |   -----------<---------------------------------   |
        // |                                                   |
        // --------------->-------------------------------------
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([-4., 3., 0.]);
        vertices.push([3., 3., 0.]);
        vertices.push([3., -2., 0.]);
        vertices.push([-4., -2., 0.]);

        vertices.push([-3., 2., 0.]);
        vertices.push([2., 2., 0.]);
        vertices.push([2., -1., 0.]);
        vertices.push([-3., -1., 0.]);

        vertices.push([-2., 1., 0.]);
        vertices.push([-1., 1., 0.]);
        vertices.push([-1., 0., 0.]);
        vertices.push([-2., 0., 0.]);

        vertices.push([0., 1., 0.]);
        vertices.push([1., 1., 0.]);
        vertices.push([1., 0., 0.]);
        vertices.push([0., 0., 0.]);

        let mut constrained_edges = Vec::new();

        // mesh frontier MUST be counter clockwise
        // then other domains MUST alternate between clockwise and counter clockwise

        constrained_edges.extend([
            // counter clockwise
            Edge::new(0, 3),
            Edge::new(3, 2),
            Edge::new(2, 1),
            Edge::new(1, 0),
            //clockwise
            Edge::new(4, 5),
            Edge::new(5, 6),
            Edge::new(6, 7),
            Edge::new(7, 4),
            //counter clockwise
            Edge::new(8, 11),
            Edge::new(11, 10),
            Edge::new(10, 9),
            Edge::new(9, 8),
            //counter clockwise
            Edge::new(12, 15),
            Edge::new(15, 14),
            Edge::new(14, 13),
            Edge::new(13, 12),
        ]);

        let plane_normal = Vec3::Z;
        let triangulation = constrained_triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            &constrained_edges,
            ConstrainedTriangulationConfiguration::default(),
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in triangulation.vert_indices {
            constrained_vertices.push(vertices[vertexid as usize]);
        }

        assert_eq!(36, constrained_vertices.len());
        assert_eq!(
            Vec::from([
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [3.0, 3.0, 0.0],
                [2.0, -1.0, 0.0],
                [3.0, -2.0, 0.0],
                [-4.0, -2.0, 0.0],
                [3.0, -2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-4.0, -2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-4.0, -2.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-3.0, 2.0, 0.0],
                [2.0, 2.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-3.0, 2.0, 0.0],
                [3.0, 3.0, 0.0],
                [-4.0, 3.0, 0.0],
                [2.0, 2.0, 0.0],
                [3.0, 3.0, 0.0],
                [2.0, 2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-1.0, 1.0, 0.0],
                [-2.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [-2.0, 1.0, 0.0],
                [-2.0, 0.0, 0.0],
                [-1.0, 1.0, 0.0],
            ]),
            constrained_vertices
        );
    }

    #[test]
    fn constrained_delaunay_level_4() {
        let mut vertices = Vec::<[Float; 3]>::new();

        vertices.push([-4., 3., 0.]);
        vertices.push([3., 3., 0.]);
        vertices.push([3., -2., 0.]);
        vertices.push([-4., -2., 0.]);

        vertices.push([-3., 2., 0.]);
        vertices.push([2., 2., 0.]);
        vertices.push([2., -1., 0.]);
        vertices.push([-3., -1., 0.]);

        vertices.push([-2., 1., 0.]);
        vertices.push([-1., 1., 0.]);
        vertices.push([-1., 0., 0.]);
        vertices.push([-2., 0., 0.]);

        vertices.push([0., 1., 0.]);
        vertices.push([1., 1., 0.]);
        vertices.push([1., 0., 0.]);
        vertices.push([0., 0., 0.]);

        let mut constrained_edges = Vec::new();

        // mesh frontier MUST be counter clockwise
        // then other domains MUST alternate between clockwise and counter clockwise

        constrained_edges.extend([
            // counter clockwise
            Edge::new(0, 3),
            Edge::new(3, 2),
            Edge::new(2, 1),
            Edge::new(1, 0),
            //clockwise
            Edge::new(4, 5),
            Edge::new(5, 6),
            Edge::new(6, 7),
            Edge::new(7, 4),
            //counter clockwise
            Edge::new(8, 11),
            Edge::new(11, 10),
            Edge::new(10, 9),
            Edge::new(9, 8),
            //counter clockwise
            Edge::new(12, 15),
            Edge::new(15, 14),
            Edge::new(14, 13),
            Edge::new(13, 12),
        ]);

        let plane_normal = Vec3::Z;
        let triangulation = constrained_triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            &constrained_edges,
            ConstrainedTriangulationConfiguration::default(),
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in triangulation.vert_indices {
            constrained_vertices.push(vertices[vertexid as usize]);
        }

        assert_eq!(36, constrained_vertices.len());
        assert_eq!(
            Vec::from([
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [3.0, 3.0, 0.0],
                [2.0, -1.0, 0.0],
                [3.0, -2.0, 0.0],
                [-4.0, -2.0, 0.0],
                [3.0, -2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-4.0, -2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-4.0, -2.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-3.0, -1.0, 0.0],
                [-3.0, 2.0, 0.0],
                [2.0, 2.0, 0.0],
                [-4.0, 3.0, 0.0],
                [-3.0, 2.0, 0.0],
                [3.0, 3.0, 0.0],
                [-4.0, 3.0, 0.0],
                [2.0, 2.0, 0.0],
                [3.0, 3.0, 0.0],
                [2.0, 2.0, 0.0],
                [2.0, -1.0, 0.0],
                [-1.0, 1.0, 0.0],
                [-2.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [-2.0, 1.0, 0.0],
                [-2.0, 0.0, 0.0],
                [-1.0, 1.0, 0.0],
            ]),
            constrained_vertices
        );
    }
}
