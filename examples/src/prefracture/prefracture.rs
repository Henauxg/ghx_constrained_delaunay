use bevy::{
    app::{App, Startup, Update},
    asset::{AssetServer, Assets, Handle},
    ecs::{
        component::Component,
        entity::Entity,
        event::EventReader,
        query::With,
        schedule::IntoSystemConfigs,
        system::{Commands, Query, Res, ResMut},
    },
    hierarchy::{Children, DespawnRecursiveExt, HierarchyQueryExt, Parent},
    input::{keyboard::KeyCode, ButtonInput},
    render::mesh::Mesh,
    scene::SceneBundle,
    time::{Fixed, Time},
    transform::components::Transform,
    utils::default,
    DefaultPlugins,
};
use bevy_rapier3d::{
    dynamics::{FixedJointBuilder, ImpulseJoint, RigidBody},
    geometry::{
        ActiveCollisionTypes, Collider, ColliderMassProperties, ComputedColliderShape, Friction,
        Restitution,
    },
    pipeline::{CollisionEvent, ContactForceEvent},
    plugin::{NoUserData, RapierPhysicsPlugin},
};
use examples::plugin::ExamplesPlugin;

pub const CUBE_FRAC_ASSET_PATH: &str = "cube_frac.glb#Scene0";
pub const CUBE_ASSET_PATH: &str = "cube.glb#Scene0";

fn main() {
    App::new()
        .insert_resource(Time::<Fixed>::from_seconds(10.5))
        .add_plugins((DefaultPlugins, ExamplesPlugin))
        .add_plugins((
            RapierPhysicsPlugin::<NoUserData>::default(),
            // RapierDebugRenderPlugin::default(),
        ))
        .add_systems(Startup, setup_scene)
        .add_systems(
            Update,
            (
                (respawn_cube, attach_physics_components_to_cells).chain(),
                handle_collisions,
            ),
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

fn setup_scene(mut commands: Commands, asset_server: Res<AssetServer>) {
    spawn_frac_cube(&asset_server, &mut commands);
}

fn spawn_cube(asset_server: &Res<AssetServer>, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_ASSET_PATH),
            transform: Transform::from_xyz(0., 1., 0.),
            ..default()
        },
        ExampleCube,
    ));
}

fn spawn_frac_cube(asset_server: &Res<AssetServer>, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_FRAC_ASSET_PATH),
            transform: Transform::from_xyz(0., 1., 0.),
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
#[derive(Component)]
pub struct Fragment;

fn attach_physics_components_to_cells(
    mut commands: Commands,
    mut fractured_scene: Query<(Entity, &mut FragmentedCubeRoot)>,
    children: Query<&Children>,
    fragments: Query<(&Parent, &Handle<Mesh>)>,
    meshes: ResMut<Assets<Mesh>>,
    transforms: Query<&Transform>,
) {
    for (fractured_scene_entity, mut fragments_root) in fractured_scene.iter_mut() {
        commands
            .entity(fractured_scene_entity)
            .insert(RigidBody::Dynamic);

        // We only want to attach the physics components once
        if !fragments_root.physics_applied {
            for entity in children.iter_descendants(fractured_scene_entity) {
                // Attach the physics components only to the meshes entities
                if let Ok((parent, mesh_handle)) = fragments.get(entity) {
                    let fragment_parent_transform = transforms.get(parent.get()).unwrap();

                    let mesh = meshes.get(mesh_handle).unwrap();
                    let collider =
                        Collider::from_bevy_mesh(mesh, &ComputedColliderShape::ConvexHull).unwrap();
                    let joint = FixedJointBuilder::new()
                        .local_anchor1(fragment_parent_transform.translation);
                    commands.entity(entity).insert((
                        Fragment,
                        RigidBody::Dynamic,
                        collider,
                        ActiveCollisionTypes::default(),
                        Friction::coefficient(0.7),
                        Restitution::coefficient(0.05),
                        ColliderMassProperties::Density(2.0),
                        ImpulseJoint::new(fractured_scene_entity, joint),
                    ));
                    // Children do not seem to be created at the same time as the cube root entity, so we only update the flag once children are present
                    fragments_root.physics_applied = true;
                }
            }
        }
    }
}

fn handle_collisions(
    mut commands: Commands,
    mut collision_events: EventReader<CollisionEvent>,
    mut contact_force_events: EventReader<ContactForceEvent>,
    fragments: Query<&Fragment>,
) {
    for _collision_event in collision_events.read() {
        // println!("Received collision event: {:?}", collision_event);
    }

    for contact_force_event in contact_force_events.read() {
        let fragment = if let Ok(_) = fragments.get(contact_force_event.collider1) {
            Some(contact_force_event.collider1)
        } else if let Ok(_) = fragments.get(contact_force_event.collider2) {
            Some(contact_force_event.collider2)
        } else {
            None
        };
        if let Some(fragment_entity) = fragment {
            println!("Handled contact force event: {:?}", contact_force_event);
            commands.entity(fragment_entity).remove::<ImpulseJoint>();
        }
    }
}
