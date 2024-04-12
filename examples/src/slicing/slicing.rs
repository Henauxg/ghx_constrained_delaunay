use bevy_inspector_egui::quick::WorldInspectorPlugin;
use bevy_rapier3d::prelude::*;
use examples::plugin::ExamplesPlugin;
use std::f32::consts::*;

use bevy::{
    ecs::query,
    pbr::CascadeShadowConfigBuilder,
    prelude::*,
    render::{
        mesh::{shape::Quad, MeshVertexAttribute, PrimitiveTopology},
        render_asset::RenderAssetUsages,
    },
};

pub mod utils;

pub const CUBE_FRAC_ASSET_PATH: &str = "cube_frac.glb#Scene0";
pub const CUBE_ASSET_PATH: &str = "cube.glb#Scene0";

fn main() {
    App::new()
        .insert_resource(AmbientLight {
            color: Color::WHITE,
            brightness: 1.0 / 5.0f32,
        })
        .insert_resource(Time::<Fixed>::from_seconds(10.5))
        .add_plugins((DefaultPlugins, ExamplesPlugin, WorldInspectorPlugin::new()))
        .add_plugins(RapierPhysicsPlugin::<NoUserData>::default())
        // .add_plugins(RapierDebugRenderPlugin::default())
        .add_systems(Startup, setup_scene)
        .add_systems(
            Update,
            (respawn_cube, attach_physics_components_to_cells).chain(),
        )
        .run();
}

#[derive(Component, Debug)]
struct FragmentedCubeRoot {
    physics_applied: bool,
}

impl Default for FragmentedCubeRoot {
    fn default() -> Self {
        Self {
            physics_applied: false,
        }
    }
}

#[derive(Component)]
struct ExampleCube;

fn setup_scene(
    mut commands: Commands,
    asset_server: Res<AssetServer>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // // Plane
    let radius = 2000.;
    let height = 200.;
    commands.spawn((
        PbrBundle {
            mesh: meshes.add(Cylinder::new(radius, height)),
            material: materials.add(Color::WHITE),
            transform: Transform::from_xyz(0.0, -height / 2., 0.0),
            ..default()
        },
        Collider::cylinder(height / 2., radius),
        (ActiveCollisionTypes::default()),
        Friction::coefficient(0.7),
        Restitution::coefficient(0.3),
    ));

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

    spawn_frac_cube(&asset_server, &mut commands);
}

fn spawn_cube(asset_server: &Res<AssetServer>, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_ASSET_PATH),
            transform: Transform::from_xyz(0., 2., 0.),
            ..default()
        },
        ExampleCube,
    ));
}

fn spawn_frac_cube(asset_server: &Res<AssetServer>, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_FRAC_ASSET_PATH),
            transform: Transform::from_xyz(0., 2., 0.),
            ..default()
        },
        FragmentedCubeRoot::default(),
        ExampleCube,
    ));
}

fn respawn_cube(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    query_cubes: Query<Entity, With<ExampleCube>>,
    asset_server: Res<AssetServer>,
    mut commands: Commands,
) {
    if keyboard_input.pressed(KeyCode::KeyF) {
        for entity in query_cubes.iter() {
            commands.entity(entity).despawn_recursive();
        }
        spawn_cube(&asset_server, &mut commands);
    }
    if keyboard_input.pressed(KeyCode::KeyG) {
        for entity in query_cubes.iter() {
            commands.entity(entity).despawn_recursive();
        }
        spawn_frac_cube(&asset_server, &mut commands);
    }
}

fn attach_physics_components_to_cells(
    mut commands: Commands,
    mut fractured_scene: Query<(Entity, &mut FragmentedCubeRoot)>,
    children: Query<&Children>,
    meshes_handles: Query<&Handle<Mesh>>,
    meshes: ResMut<Assets<Mesh>>,
) {
    for (fractured_scene_entity, mut fragments_root) in fractured_scene.iter_mut() {
        // We only want to attach the physics components once
        if !fragments_root.physics_applied {
            for entity in children.iter_descendants(fractured_scene_entity) {
                // Attach the physics components only to the meshes entities
                if let Ok(mesh_handle) = meshes_handles.get(entity) {
                    let mesh = meshes.get(mesh_handle).unwrap();
                    info!("Attaching physics components to entity {:?}", entity);
                    let collider =
                        Collider::from_bevy_mesh(mesh, &ComputedColliderShape::ConvexHull).unwrap();
                    commands.entity(entity).insert((
                        RigidBody::Dynamic,
                        collider,
                        ActiveCollisionTypes::default(),
                        Friction::coefficient(0.7),
                        Restitution::coefficient(0.05),
                        ColliderMassProperties::Density(2.0),
                    ));
                    // Children do not seem to be created at the same time as the cube root entity, so we only update the flag once children are present
                    fragments_root.physics_applied = true;
                }
            }
        }
    }
}

// fn slice(
//     source_object: Entity,
//     fragment_template: Entity,
//     slice_normal: Vec3,
//     slice_origin: Vec3,
//     parent: Transform,
// ) {

//     // 1 mesh cube

//     // slice en 2

//     // spawn frags
// }

// fn create_fragment(
//     mut commands: Commands,
//     source_object: Entity,
//     meshes_handles: Query<&Handle<Mesh>>,
//     query_rigid_bodies: Query<&RigidBody>,
//     query_trasnform: Query<&Transform>,
//     fragment_template: Entity,
//     parent: Transform,
//     fragment_mesh_data: Vec<f32>,
//     meshes_assets: ResMut<Assets<Mesh>>,
// ) {
//     // Fragment mesh creation
//     let mut fragment_mesh = Mesh::new(
//         PrimitiveTopology::TriangleList,
//         RenderAssetUsages::MAIN_WORLD,
//     );

//     // Add vertices to mesh
//     fragment_mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, fragment_mesh_data);

//     // List of all fragmented mesh
//     let mut meshes = Vec::new();
//     meshes.push(fragment_mesh);

//     //Get BB size and Rigidbody from the source object
//     if let Ok(source_object_mesh_handle) = meshes_handles.get(source_object) {
//         let source_object_mesh = meshes_assets.get(source_object_mesh_handle).unwrap();

//         //BB
//         let parent_bounding_box_size = source_object_mesh.compute_aabb().unwrap().max()
//             - source_object_mesh.compute_aabb().unwrap().min();
//         //Rigidbody
//         let parent_rigid_body = query_rigid_bodies.get(source_object).unwrap();
//     }

//     for i in 0..meshes.len() {
//         // Create a new entity based on fragment template entity
//         let fragment = fragment_template.clone();

//         commands.entity(parent).push_children(fragment);

//         // Set all the transform parameters
//         let mut fragment_transform = *query_trasnform.get(fragment).unwrap();
//         let source_object_transform = *query_trasnform.get(source_object).unwrap();
//         fragment_transform.translation = Vec3::new(0., 0., 0.);
//         fragment_transform.rotation = Quat::IDENTITY;
//         fragment_transform.scale = source_object_transform.scale;
//     }
// }
