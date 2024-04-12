use std::f32::consts::PI;

use bevy::{
    app::{App, Startup, Update},
    asset::{Assets, Handle},
    ecs::{
        system::{Commands, NonSend, ResMut},
        world::World,
    },
    gizmos::gizmos::Gizmos,
    math::{primitives::Plane3d, EulerRot, Quat, Vec3},
    pbr::{
        AmbientLight, CascadeShadowConfigBuilder, DirectionalLight, DirectionalLightBundle,
        PbrBundle, StandardMaterial,
    },
    render::{
        color::Color,
        mesh::{Indices, Mesh, Meshable, PrimitiveTopology},
        render_asset::RenderAssetUsages,
    },
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use examples::{plugin::ExamplesPlugin, COLORS};
use tritet::Tetgen;

fn main() {
    let mut app = App::new();

    app.insert_resource(AmbientLight {
        color: Color::WHITE,
        brightness: 2000.,
    })
    .add_plugins((DefaultPlugins, ExamplesPlugin));

    app.add_systems(Startup, (setup_sandbox, setup_procedural_mesh))
        .add_systems(Update, draw_delaunay);

    app.run();
}

pub fn setup_sandbox(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // Plane
    commands.spawn(PbrBundle {
        mesh: meshes.add(Plane3d::default().mesh().size(500000.0, 500000.0)),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3)),
        transform: Transform::from_xyz(0., -0.5, 0.),
        ..default()
    });

    // Light
    commands.spawn(DirectionalLightBundle {
        transform: Transform::from_rotation(Quat::from_euler(EulerRot::ZYX, 0.0, 1.0, -PI / 4.)),
        directional_light: DirectionalLight {
            shadows_enabled: true,
            ..default()
        },
        cascade_shadow_config: CascadeShadowConfigBuilder {
            first_cascade_far_bound: 200.0,
            maximum_distance: 400.0,
            ..default()
        }
        .into(),
        ..default()
    });
}

pub struct DestructionMesh {
    delaunay: Tetgen,
}

pub fn setup_procedural_mesh(world: &mut World) {
    // allocate data for 8 points
    let mut delaunay = Tetgen::new(8, None, None, None).unwrap();

    // set points
    delaunay
        .set_point(0, 0, 0.0, 0.0, 0.0)
        .unwrap()
        .set_point(1, 0, 1.0, 0.0, 0.0)
        .unwrap()
        .set_point(2, 0, 1.0, 1.0, 0.0)
        .unwrap()
        .set_point(3, 0, 0.0, 1.0, 0.0)
        .unwrap()
        .set_point(4, 0, 0.0, 0.0, 1.0)
        .unwrap()
        .set_point(5, 0, 1.0, 0.0, 1.0)
        .unwrap()
        .set_point(6, 0, 1.0, 1.0, 1.0)
        .unwrap()
        .set_point(7, 0, 0.0, 1.0, 1.0)
        .unwrap();
    delaunay.generate_delaunay(false).unwrap();
    delaunay.generate_mesh(false, false, None, None).unwrap();

    world.insert_non_send_resource(DestructionMesh { delaunay });
}

pub fn draw_delaunay(destr_mesh: NonSend<DestructionMesh>, mut gizmos: Gizmos) {
    const TETRA_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

    let delaunay = &destr_mesh.delaunay;

    let tetras_count = delaunay.out_ncell();
    // let attributes = Vec::with_capacity(tetras_count);
    for tetra_id in 0..tetras_count {
        // attributes.push(delaunay.out_cell_attribute(tetra_idx));
        let mut tetra_vertices = Vec::with_capacity(4);
        for vertex_local_id in 0..=3 {
            let vert_id = delaunay.out_cell_point(tetra_id, vertex_local_id);
            tetra_vertices.push(Vec3::new(
                delaunay.out_point(vert_id, 0) as f32,
                delaunay.out_point(vert_id, 1) as f32,
                delaunay.out_point(vert_id, 2) as f32,
            ));
        }
        let color = COLORS[tetra_id % COLORS.len()];
        for (vertex_local_id_a, vertex_local_id_b) in TETRA_EDGES {
            gizmos.line(
                tetra_vertices[vertex_local_id_a],
                tetra_vertices[vertex_local_id_b],
                color,
            );
        }
    }
}
