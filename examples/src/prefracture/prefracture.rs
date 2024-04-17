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
    log::info,
    render::mesh::Mesh,
    scene::SceneBundle,
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
use examples::{ball::BallSensor, plugin::ExamplesPlugin};

// "dragon_frac.glb#Scene0";
// "dragon_low_poly.glb#Scene0";
pub const CUBE_FRAC_ASSET_PATH: &str = "cube_frac.glb#Scene0";
pub const CUBE_ASSET_PATH: &str = "cube.glb#Scene0";

const DEFAULT_TRANSFORM: Transform = Transform::from_xyz(0., 1., 0.);

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, ExamplesPlugin))
        .add_plugins((
            RapierPhysicsPlugin::<NoUserData>::default(),
            // RapierDebugRenderPlugin::default(),
        ))
        .add_systems(Startup, setup_scene)
        .add_systems(
            Update,
            (
                (respawn_cube, handle_collisions),
                (
                    attach_physics_components_to_fragments,
                    attach_physics_components_to_full_meshes,
                ),
            )
                .chain(),
        )
        .run();
}

#[derive(Component)]
struct UninitializedFragmentedMesh;

#[derive(Component)]
pub struct UninitializedFullMesh;

#[derive(Component)]
/// Contains a reference to the root Entity
pub struct FullMesh(Entity);

#[derive(Component)]
pub struct FragmentedMesh;

#[derive(Component)]
struct ExampleCube;

fn setup_scene(mut commands: Commands, asset_server: Res<AssetServer>) {
    spawn_frac_cube(&asset_server, DEFAULT_TRANSFORM, &mut commands);
}

fn spawn_cube(asset_server: &Res<AssetServer>, transform: Transform, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_ASSET_PATH),
            transform,
            ..default()
        },
        ExampleCube,
        UninitializedFullMesh,
    ));
}

fn spawn_frac_cube(asset_server: &Res<AssetServer>, transform: Transform, commands: &mut Commands) {
    commands.spawn((
        SceneBundle {
            scene: asset_server.load(CUBE_FRAC_ASSET_PATH),
            transform,
            ..default()
        },
        ExampleCube,
        UninitializedFragmentedMesh,
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
        spawn_cube(&asset_server, DEFAULT_TRANSFORM, &mut commands);
    }
    if keyboard_input.pressed(KeyCode::KeyG) {
        for entity in query_cubes.iter() {
            commands.entity(entity).despawn_recursive();
        }
        spawn_frac_cube(&asset_server, DEFAULT_TRANSFORM, &mut commands);
    }
}

fn attach_physics_components_to_full_meshes(
    mut commands: Commands,
    mut full_scenes: Query<Entity, With<UninitializedFullMesh>>,
    children: Query<&Children>,
    fragments: Query<&Handle<Mesh>>,
    meshes: ResMut<Assets<Mesh>>,
) {
    for full_scene in full_scenes.iter_mut() {
        let mut children_physics_applied = false;
        for entity in children.iter_descendants(full_scene) {
            // Attach the physics components only to the meshes entities
            if let Ok(mesh_handle) = fragments.get(entity) {
                let mesh = meshes.get(mesh_handle).unwrap();
                let collider =
                    Collider::from_bevy_mesh(mesh, &ComputedColliderShape::ConvexHull).unwrap();
                commands.entity(entity).insert((
                    (
                        RigidBody::Dynamic,
                        ActiveCollisionTypes::default(),
                        Friction::coefficient(0.7),
                        Restitution::coefficient(0.05),
                        ColliderMassProperties::Density(20.0),
                        collider,
                    ),
                    FullMesh(full_scene),
                ));
                // Children do not seem to be created at the same time as the cube root entity, so we only update the flag once children are present
                children_physics_applied = true;
            }
        }
        // We only want to attach the physics components once
        if children_physics_applied {
            commands
                .entity(full_scene)
                .remove::<UninitializedFullMesh>();
        }
    }
}

fn attach_physics_components_to_fragments(
    mut commands: Commands,
    mut fractured_scenes: Query<Entity, With<UninitializedFragmentedMesh>>,
    children: Query<&Children>,
    fragments: Query<(&Parent, &Handle<Mesh>)>,
    meshes: ResMut<Assets<Mesh>>,
    transforms: Query<&Transform>,
) {
    for fractured_scene_entity in fractured_scenes.iter_mut() {
        let mut children_physics_applied = false;
        for entity in children.iter_descendants(fractured_scene_entity) {
            // Attach the physics components only to the meshes entities
            if let Ok((parent, mesh_handle)) = fragments.get(entity) {
                let fragment_parent_transform = transforms.get(parent.get()).unwrap();

                let mesh = meshes.get(mesh_handle).unwrap();
                let collider =
                    Collider::from_bevy_mesh(mesh, &ComputedColliderShape::ConvexHull).unwrap();
                let joint =
                    FixedJointBuilder::new().local_anchor1(fragment_parent_transform.translation);
                commands.entity(entity).insert((
                    (
                        RigidBody::Dynamic,
                        ActiveCollisionTypes::default(),
                        Friction::coefficient(0.7),
                        Restitution::coefficient(0.05),
                        ColliderMassProperties::Density(2.0),
                        collider,
                        ImpulseJoint::new(fractured_scene_entity, joint),
                    ),
                    FragmentedMesh,
                ));
                // Children do not seem to be created at the same time as the cube root entity, so we only update the flag once children are present
                children_physics_applied = true;
            }
        }
        // We only want to attach the physics components once
        if children_physics_applied {
            commands
                .entity(fractured_scene_entity)
                .insert(RigidBody::Fixed)
                .remove::<UninitializedFragmentedMesh>();
        }
    }
}

fn handle_collisions(
    mut commands: Commands,
    mut collision_events: EventReader<CollisionEvent>,
    mut contact_force_events: EventReader<ContactForceEvent>,
    asset_server: Res<AssetServer>,
    ball_sensors: Query<(), With<BallSensor>>,
    fragments: Query<&FragmentedMesh>,
    full_meshes: Query<&FullMesh>,
    transforms: Query<&Transform>,
) {
    for collision_event in collision_events.read() {
        match collision_event {
            CollisionEvent::Started(e1, e2, _flags) => {
                let other_entity = if let Ok(_) = ball_sensors.get(*e1) {
                    Some(*e2)
                } else if let Ok(_) = ball_sensors.get(*e2) {
                    Some(*e1)
                } else {
                    None
                };
                if let Some(other_entity) = other_entity {
                    if let Ok((full_mesh_root_entity)) = full_meshes.get(other_entity) {
                        let full_mesh_root_transform =
                            transforms.get(full_mesh_root_entity.0).unwrap();
                        commands.entity(full_mesh_root_entity.0).despawn_recursive();
                        spawn_frac_cube(&asset_server, *full_mesh_root_transform, &mut commands);
                    }
                }
            }
            CollisionEvent::Stopped(e1, e2, flags) => {
                // TODO Could count the sensors in presence for each full mesh and swap back to a full mesh after some time ?
            }
        }
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
            info!("Handled contact force event: {:?}", contact_force_event);
            commands.entity(fragment_entity).remove::<ImpulseJoint>();
        }
    }
}
