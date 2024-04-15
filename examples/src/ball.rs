use std::time::Duration;

use bevy::{
    asset::{Assets, Handle},
    ecs::{
        component::Component,
        entity::Entity,
        query::With,
        system::{Commands, Query, Res, ResMut, Resource},
    },
    hierarchy::DespawnRecursiveExt,
    math::primitives::Sphere,
    pbr::{PbrBundle, StandardMaterial},
    render::{
        camera::Camera,
        color::Color,
        mesh::{Mesh, Meshable},
    },
    time::{Time, Timer, TimerMode},
    transform::components::{GlobalTransform, Transform},
    utils::default,
};
use bevy_rapier3d::{
    dynamics::{ExternalImpulse, RigidBody},
    geometry::{ActiveCollisionTypes, Collider, ColliderMassProperties, Friction, Restitution},
};

pub const BALL_RADIUS: f32 = 0.5;
pub const BALL_DESPAWN_TIMER_S: u64 = 5;
pub const BALL_THROW_FORCE: f32 = 100.0;

#[derive(Resource)]
pub struct BallAssets {
    pub mesh: Handle<Mesh>,
    pub material: Handle<StandardMaterial>,
}

pub fn setup_ball_assets(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let material = materials.add(StandardMaterial {
        emissive: Color::rgb_linear(9000.0, 3000.0, 23000.0),
        ..default()
    });
    let mesh = meshes.add(Sphere::new(BALL_RADIUS).mesh().ico(5).unwrap());

    commands.insert_resource(BallAssets { mesh, material });
}

#[derive(Component)]
pub struct Ball {
    timer: Timer,
}
impl Ball {
    pub fn new() -> Self {
        Self {
            timer: Timer::new(Duration::from_secs(BALL_DESPAWN_TIMER_S), TimerMode::Once),
        }
    }
}

pub fn throw_ball(
    mut commands: Commands,
    ball_assets: Res<BallAssets>,
    camera: Query<&GlobalTransform, With<Camera>>,
) {
    let cam_transform = camera.get_single().unwrap();
    let mut transform = Transform::from(*cam_transform);
    transform.translation += transform.forward() * 4. * BALL_RADIUS;
    commands.spawn((
        PbrBundle {
            mesh: ball_assets.mesh.clone(),
            material: ball_assets.material.clone(),
            transform,
            ..default()
        },
        Ball::new(),
        RigidBody::Dynamic,
        Collider::ball(BALL_RADIUS),
        ActiveCollisionTypes::default(),
        Friction::coefficient(0.7),
        Restitution::coefficient(0.05),
        ColliderMassProperties::Density(2.0),
        ExternalImpulse {
            impulse: transform.forward() * BALL_THROW_FORCE,
            ..default()
        },
    ));
}

pub fn despawn_balls(
    mut commands: Commands,
    time: Res<Time>,
    mut balls: Query<(Entity, &mut Ball)>,
) {
    for (ball_entity, mut ball) in balls.iter_mut() {
        ball.timer.tick(time.delta());
        if ball.timer.finished() {
            commands.entity(ball_entity).despawn_recursive();
        }
    }
}
