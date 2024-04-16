use std::time::Duration;

use bevy::{
    app::{App, Plugin, PostUpdate, Startup, Update},
    asset::{Assets, Handle},
    ecs::{
        component::Component,
        entity::Entity,
        query::{With, Without},
        schedule::IntoSystemConfigs,
        system::{Commands, Query, Res, ResMut, Resource},
    },
    hierarchy::DespawnRecursiveExt,
    input::{common_conditions::input_just_pressed, keyboard::KeyCode},
    math::primitives::Sphere,
    pbr::{AlphaMode, PbrBundle, StandardMaterial},
    render::{
        camera::Camera,
        color::Color,
        mesh::{Mesh, Meshable},
    },
    time::{Time, Timer, TimerMode},
    transform::{
        components::{GlobalTransform, Transform},
        TransformSystem,
    },
    utils::default,
};
use bevy_rapier3d::{
    dynamics::{ExternalImpulse, RigidBody},
    geometry::{
        ActiveCollisionTypes, ActiveEvents, Collider, ColliderMassProperties, Friction, Restitution,
    },
};

pub const BALL_RADIUS: f32 = 0.25;
pub const PREVIS_BALL_RADIUS: f32 = 0.35 / 10.;
pub const BALL_DESPAWN_TIMER_S: u64 = 5;
pub const BALL_THROW_FORCE: f32 = 7500.0;
pub const BALL_CAMERA_DISTANCE: f32 = 5. * BALL_RADIUS;

pub struct BallPlugin;

impl Plugin for BallPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_balls);
        app.add_systems(
            Update,
            (
                throw_ball.run_if(input_just_pressed(KeyCode::Space)),
                despawn_balls,
            ),
        );
        app.add_systems(
            PostUpdate,
            (update_previsualisation_ball.before(TransformSystem::TransformPropagate),),
        );
    }
}

#[derive(Resource)]
pub struct BallAssets {
    pub mesh: Handle<Mesh>,
    pub material: Handle<StandardMaterial>,
}

#[derive(Component)]
pub struct PrevisualisationBall;

pub fn setup_balls(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let material = materials.add(StandardMaterial {
        emissive: Color::rgb_linear(9000.0, 3000.0, 23000.0),
        ..default()
    });
    let mesh = meshes.add(Sphere::new(BALL_RADIUS).mesh().ico(5).unwrap());

    // Spawn previsualisation ball
    let previs_mesh = meshes.add(Sphere::new(PREVIS_BALL_RADIUS).mesh().ico(5).unwrap());
    let previs_material = materials.add(StandardMaterial {
        base_color: Color::WHITE.with_a(0.2),
        alpha_mode: AlphaMode::Blend,
        ..default()
    });
    commands.spawn((
        PbrBundle {
            mesh: previs_mesh.clone(),
            material: previs_material.clone(),
            ..default()
        },
        PrevisualisationBall,
    ));

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
    transform.translation += transform.forward() * BALL_CAMERA_DISTANCE;
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
        ColliderMassProperties::Density(1000.0),
        ActiveEvents::CONTACT_FORCE_EVENTS,
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

pub fn update_previsualisation_ball(
    camera: Query<&Transform, (With<Camera>, Without<PrevisualisationBall>)>,
    mut previs_ball: Query<&mut Transform, With<PrevisualisationBall>>,
) {
    let cam_transform = camera.get_single().unwrap();
    let mut previs_ball_transform = previs_ball.get_single_mut().unwrap();

    *previs_ball_transform = Transform::from(*cam_transform);
    previs_ball_transform.translation += cam_transform.forward() * BALL_CAMERA_DISTANCE;
}
