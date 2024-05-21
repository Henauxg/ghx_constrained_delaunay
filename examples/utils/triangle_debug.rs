use bevy::{
    app::{Plugin, Update},
    ecs::{
        component::Component,
        entity::Entity,
        event::{Event, EventReader, EventWriter},
        query::With,
        schedule::IntoSystemConfigs,
        system::{Commands, Query, Res, ResMut, Resource},
    },
    gizmos::gizmos::Gizmos,
    hierarchy::DespawnRecursiveExt,
    input::{keyboard::KeyCode, ButtonInput},
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
    utils::default,
};
use bevy_mod_billboard::{plugin::BillboardPlugin, BillboardTextBundle};
use ghx_constrained_delaunay::types::TriangleData;
use glam::Vec3;

use super::{COLORS, DEBUG_LABEL_FONT_SIZE};

pub struct TriangleDebugPlugin;
impl Plugin for TriangleDebugPlugin {
    fn build(&self, app: &mut bevy::prelude::App) {
        app.add_plugins(BillboardPlugin);
        app.add_event::<TriangleDebugDataUpdated>().add_systems(
            Update,
            (
                update_triangles_debug_index,
                update_triangles_debugs_labels,
                draw_triangles_debug_data,
            )
                .chain(),
        );
    }
}

pub fn create_displayed_vertices(vertices: Vec<[f32; 3]>, plane_normal: Vec3) -> Vec<Vec3> {
    let mut displayed_vertices: Vec<Vec3> = vertices.iter().map(|v| Vec3::from_slice(v)).collect();
    // Add the container triangle positions to the debug view buffer
    // Calculate container triangle vertices position in world position
    let e1 = (displayed_vertices[0] - displayed_vertices[1]).normalize();
    let e2 = plane_normal;
    let e3 = e1.cross(e2);
    for v in vec![[-100., -100., 0.], [0., 100., 0.], [100., -100., 0.]] {
        let container_vertex = displayed_vertices[1] + v[0] * e1 + v[1] * e3;
        displayed_vertices.push(container_vertex);
    }

    displayed_vertices
}

#[derive(Resource, Debug, Default)]
pub struct TrianglesDebugData {
    pub vertices: Vec<Vec3>,
    pub triangles_buffers: Vec<Vec<TriangleData>>,
    pub current_buffer_index: usize,
}
impl TrianglesDebugData {
    pub fn new(vertices: Vec<Vec3>, triangles_buffers: Vec<Vec<TriangleData>>) -> Self {
        Self {
            vertices,
            triangles_buffers,
            current_buffer_index: 0,
        }
    }
}

impl TrianglesDebugData {
    pub fn advance_cursor(&mut self) {
        self.current_buffer_index = (self.current_buffer_index + 1) % self.triangles_buffers.len();
    }

    pub fn move_back_cursor(&mut self) {
        if self.current_buffer_index == 0 {
            self.current_buffer_index = self.triangles_buffers.len() - 1;
        } else {
            self.current_buffer_index -= 1;
        }
    }

    pub fn current_triangles(&self) -> &Vec<TriangleData> {
        &self.triangles_buffers[self.current_buffer_index]
    }
}

#[derive(Event)]
pub struct TriangleDebugDataUpdated;

pub fn update_triangles_debug_index(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    mut debug_vert_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugDataUpdated>,
) {
    if keyboard_input.just_pressed(KeyCode::ArrowRight) {
        debug_vert_data.advance_cursor();
        debug_data_updates_events.send(TriangleDebugDataUpdated);
    } else if keyboard_input.just_pressed(KeyCode::ArrowLeft) {
        debug_vert_data.move_back_cursor();
        debug_data_updates_events.send(TriangleDebugDataUpdated);
    }
}

#[derive(Component)]
pub struct TriangleDebugLabel;

pub fn update_triangles_debugs_labels(
    mut commands: Commands,
    query_labels: Query<Entity, With<TriangleDebugLabel>>,
    debug_vert_data: Res<TrianglesDebugData>,
    mut debug_data_updates_events: EventReader<TriangleDebugDataUpdated>,
) {
    if !debug_data_updates_events.is_empty() {
        debug_data_updates_events.clear();
    } else {
        return;
    }

    let triangles = debug_vert_data.current_triangles();
    let vertices = &debug_vert_data.vertices;

    for label in query_labels.iter() {
        commands.entity(label).despawn_recursive();
    }

    for (index, t) in triangles.iter().enumerate() {
        let color = COLORS[index % COLORS.len()];
        let (v1, v2, v3) = (vertices[t.v1()], vertices[t.v2()], vertices[t.v3()]);
        let center = 0.3 * (v1 + v2 + v3);

        commands.spawn((
            BillboardTextBundle {
                transform: Transform::from_translation(center).with_scale(Vec3::splat(0.0085)),
                text: Text::from_sections([TextSection {
                    value: index.to_string(),
                    style: TextStyle {
                        font_size: DEBUG_LABEL_FONT_SIZE,
                        color,
                        ..default()
                    },
                }]),
                ..default()
            },
            TriangleDebugLabel,
        ));

        for (v, label) in vec![(v1, "v1"), (v2, "v2"), (v3, "v3")] {
            let v_display_pos = v + (center - v) * 0.1;
            commands.spawn((
                BillboardTextBundle {
                    transform: Transform::from_translation(v_display_pos)
                        .with_scale(Vec3::splat(0.0085)),
                    text: Text::from_sections([TextSection {
                        value: label.to_string(),
                        style: TextStyle {
                            font_size: DEBUG_LABEL_FONT_SIZE,
                            color,
                            ..default()
                        },
                    }]),
                    ..default()
                },
                TriangleDebugLabel,
            ));
        }
    }
}

pub fn draw_triangles_debug_data(mut gizmos: Gizmos, debug_vert_data: Res<TrianglesDebugData>) {
    // Current observed triangles data
    let triangles = debug_vert_data.current_triangles();
    let vertices = &debug_vert_data.vertices;

    for (index, triangle) in triangles.iter().enumerate() {
        let (v1, v2, v3) = (
            vertices[triangle.v1()],
            vertices[triangle.v2()],
            vertices[triangle.v3()],
        );
        let color = COLORS[index % COLORS.len()];
        gizmos.linestrip(vec![v1, v2, v3, v1], color);
    }
}