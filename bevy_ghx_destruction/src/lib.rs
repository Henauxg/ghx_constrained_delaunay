pub mod slicing;
pub mod triangulation;
pub mod utils;
#[cfg(test)]
mod tests {
    use bevy::{math::Vec3, utils::hashbrown::HashSet};

    use crate::triangulation::{
        constrained_triangulation::constrained_triangulation_from_3d_planar_vertices,
        triangulation::triangulation_from_3d_planar_vertices, Edge,
    };

    #[test]
    fn delaunay_level_1() {
        // ---------------
        // |  \          |
        // |     \       |
        // |        \    |
        // |           \ |
        // ---------------
        let mut vertices = Vec::<[f32; 3]>::new();

        vertices.push([0., 0., 0.]);
        vertices.push([0., 5., 0.]);
        vertices.push([5., 5., 0.]);
        vertices.push([5., 0., 0.]);

        let plane_normal = Vec3::Z;
        let (vertex_id, _) = triangulation_from_3d_planar_vertices(&vertices, plane_normal.into());

        let mut constrained_vertices = Vec::new();

        for vertexid in vertex_id {
            constrained_vertices.push(vertices[vertexid]);
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
        let mut vertices = Vec::<[f32; 3]>::new();

        vertices.push([0., 0., 0.]);
        vertices.push([0., 10., 0.]);
        vertices.push([0., 5., 0.]);
        vertices.push([5., 0., 0.]);
        vertices.push([5., 5., 0.]);
        vertices.push([5., 10., 0.]);

        let plane_normal = Vec3::Z;
        let (vertex_id, _) = triangulation_from_3d_planar_vertices(&vertices, plane_normal.into());

        let mut constrained_vertices = Vec::new();

        for vertexid in vertex_id {
            constrained_vertices.push(vertices[vertexid]);
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
        let mut vertices = Vec::<[f32; 3]>::new();

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

        let mut constrained_edges: HashSet<Edge> = HashSet::new();

        // mesh frontier MUST be counter clockwise
        // then other domains MUST alternate between clockwise and counter clockwise

        // counter clockwise
        constrained_edges.insert(Edge::new(0, 3));
        constrained_edges.insert(Edge::new(3, 2));
        constrained_edges.insert(Edge::new(2, 1));
        constrained_edges.insert(Edge::new(1, 0));

        //clockwise
        constrained_edges.insert(Edge::new(4, 5));
        constrained_edges.insert(Edge::new(5, 6));
        constrained_edges.insert(Edge::new(6, 7));
        constrained_edges.insert(Edge::new(7, 4));

        //counter clockwise
        constrained_edges.insert(Edge::new(8, 11));
        constrained_edges.insert(Edge::new(11, 10));
        constrained_edges.insert(Edge::new(10, 9));
        constrained_edges.insert(Edge::new(9, 8));

        //counter clockwise
        constrained_edges.insert(Edge::new(12, 15));
        constrained_edges.insert(Edge::new(15, 14));
        constrained_edges.insert(Edge::new(14, 13));
        constrained_edges.insert(Edge::new(13, 12));

        let plane_normal = Vec3::Z;
        let (vertex_id, _) = constrained_triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            &constrained_edges,
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in vertex_id {
            constrained_vertices.push(vertices[vertexid]);
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
        
        let mut vertices = Vec::<[f32; 3]>::new();

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

        let mut constrained_edges: HashSet<Edge> = HashSet::new();

        // mesh frontier MUST be counter clockwise
        // then other domains MUST alternate between clockwise and counter clockwise

        // counter clockwise
        constrained_edges.insert(Edge::new(0, 3));
        constrained_edges.insert(Edge::new(3, 2));
        constrained_edges.insert(Edge::new(2, 1));
        constrained_edges.insert(Edge::new(1, 0));

        //clockwise
        constrained_edges.insert(Edge::new(4, 5));
        constrained_edges.insert(Edge::new(5, 6));
        constrained_edges.insert(Edge::new(6, 7));
        constrained_edges.insert(Edge::new(7, 4));

        //counter clockwise
        constrained_edges.insert(Edge::new(8, 11));
        constrained_edges.insert(Edge::new(11, 10));
        constrained_edges.insert(Edge::new(10, 9));
        constrained_edges.insert(Edge::new(9, 8));

        //counter clockwise
        constrained_edges.insert(Edge::new(12, 15));
        constrained_edges.insert(Edge::new(15, 14));
        constrained_edges.insert(Edge::new(14, 13));
        constrained_edges.insert(Edge::new(13, 12));

        let plane_normal = Vec3::Z;
        let (vertex_id, _) = constrained_triangulation_from_3d_planar_vertices(
            &vertices,
            plane_normal.into(),
            &constrained_edges,
        );

        let mut constrained_vertices = Vec::new();

        for vertexid in vertex_id {
            constrained_vertices.push(vertices[vertexid]);
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
