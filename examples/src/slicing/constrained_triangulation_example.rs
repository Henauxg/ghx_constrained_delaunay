use bevy::prelude::*;

use bevy_mod_billboard::plugin::BillboardPlugin;
use examples::plugin::ExamplesPlugin;
use utils::{
    constrained_triangulation::triangulate_3d_planar_vertices_constrained,
    debug_utils::{
        create_displayed_vertices, draw_triangles_debug_data, update_triangles_debug_index,
        update_triangles_debugs_labels, TriangleDebugDataUpdated, TrianglesDebugData,
    },
};

pub mod utils;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin, BillboardPlugin))
        .add_event::<TriangleDebugDataUpdated>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                update_triangles_debug_index,
                update_triangles_debugs_labels,
                draw_triangles_debug_data,
            )
                .chain(),
        )
        .run();
}

fn setup(mut commands: Commands) {
    let mut vertices = Vec::<[f32; 3]>::new();
    // vertices.push([0., 0., 0.]);
    // vertices.push([0., 10., 0.]);
    // vertices.push([0., 5., 0.]);
    // vertices.push([5., 0., 0.]);
    // vertices.push([5., 5., 0.]);
    // vertices.push([5., 10., 0.]);

    // vertices.push([5., 5., 0.]);
    // vertices.push([0., 6., 0.]);
    // vertices.push([-4., 4., 0.]);
    // vertices.push([-5., 0., 0.]);
    // vertices.push([-3., -2., 0.]);
    // vertices.push([0., -4., 0.]);
    // vertices.push([4., -5., 0.]);
    // vertices.push([6., -3., 0.]);
    // vertices.push([7., 0., 0.]);
    // vertices.push([5., 1., 0.]);
    // vertices.push([0., -6., 0.]);
    // vertices.push([4., -6., 0.]);
    // vertices.push([0., -8., 0.]);
    // vertices.push([4., -8., 0.]);

    vertices.push([0., 0., 0.]); //0
    vertices.push([6., 0., 0.]); //1
    vertices.push([7., 4., 0.]); //2
    vertices.push([8., 2., 0.]); //3
    vertices.push([10., 2., 0.]); //4
    vertices.push([10., 8., 0.]); //5
    vertices.push([9., 10., 0.]); //6
    vertices.push([7., 10., 0.]); //7
    vertices.push([5., 8., 0.]); //8
    vertices.push([4., 10., 0.]); //9
    vertices.push([0., 10., 0.]); //10
    vertices.push([0., 3., 0.]); //11 A shape begin
    vertices.push([1., 5., 0.]); //12
    vertices.push([1., 1., 0.]); //13
    vertices.push([1.5, 1., 0.]); //14
    vertices.push([1.5, 6., 0.]); //15
    vertices.push([0., 4., 0.]); //16 A shape begin

    let mut constrained_edges = Vec::new();
    constrained_edges.push((0,1));
    constrained_edges.push((1,2));
    constrained_edges.push((2,3));
    constrained_edges.push((3,4));
    constrained_edges.push((4,5));
    constrained_edges.push((5,6));
    constrained_edges.push((6,7));
    constrained_edges.push((7,8));
    constrained_edges.push((8,9));
    constrained_edges.push((9,10));
    //constrained_edges.push((10,11));
    // constrained_edges.push((11,12));
    // constrained_edges.push((12,13));
    // constrained_edges.push((13,14));
    // constrained_edges.push((14,15));
    constrained_edges.push((13,14));
    constrained_edges.push((14,15));
    constrained_edges.push((15,16));

    constrained_edges.push((10,16));
    // constrained_edges.push((12,13));
    // constrained_edges.push((11,12));

    let plane_normal = Vec3::Z;
    let (_, debug_data) = triangulate_3d_planar_vertices_constrained(
        &vertices,
        plane_normal.into(),
        &constrained_edges,
    );

    let displayed_vertices = create_displayed_vertices(vertices, plane_normal);

    commands.insert_resource(TrianglesDebugData {
        vertices: displayed_vertices,
        triangles_buffers: debug_data,
        current_buffer_index: 0,
    })
}
