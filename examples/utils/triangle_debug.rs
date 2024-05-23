use std::cmp::min;

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
    log::info,
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
    utils::default,
};
use bevy_mod_billboard::{plugin::BillboardPlugin, BillboardTextBundle};
use ghx_constrained_delaunay::{
    triangulation::{DebugContext, CONTAINER_TRIANGLE_VERTICES},
    types::TriangleData,
};
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
                update_triangles_debugs_entities,
                draw_triangles_debug_data,
            )
                .chain(),
        );
    }
}

pub fn extend_displayed_vertices_with_container_vertice(
    displayed_vertices: &mut Vec<Vec3>,
    plane_normal: Vec3,
    debug_context: &DebugContext,
    detransform: bool,
) {
    // Add the container triangle positions to the debug view buffer
    // Calculate container triangle vertices position in world position

    let mut container_vertices: Vec<Vec3> = CONTAINER_TRIANGLE_VERTICES
        .iter()
        .map(|v| Vec3::new(v.x, v.y, 0.))
        .collect();
    for v in container_vertices.iter_mut() {
        v.x = (debug_context.scale_factor * v.x) + debug_context.x_min;
        v.y = (debug_context.scale_factor * v.y) + debug_context.y_min;
    }
    if detransform {
        let e1 = (displayed_vertices[0] - displayed_vertices[1]).normalize();
        let e2 = plane_normal;
        let e3 = e1.cross(e2);
        let plane_origin = displayed_vertices[1];
        for v in container_vertices.iter_mut() {
            *v = plane_origin + v.x * e1 + v.y * e3 // + 0. * e2
        }
    }

    displayed_vertices.extend(container_vertices);
}

#[derive(Resource)]
pub struct TrianglesDebugData {
    pub vertices: Vec<Vec3>,
    pub context: DebugContext,
    pub with_labels: bool,
    pub current_buffer_index: usize,
}
impl TrianglesDebugData {
    pub fn new(vertices: Vec<Vec3>, context: DebugContext, with_labels: bool) -> Self {
        Self {
            vertices,
            context,
            with_labels,
            current_buffer_index: 0,
        }
    }
}

impl TrianglesDebugData {
    pub fn advance_cursor(&mut self) {
        self.current_buffer_index =
            (self.current_buffer_index + 1) % self.context.triangles_buffers.len();
    }

    pub fn move_back_cursor(&mut self) {
        if self.current_buffer_index == 0 {
            self.current_buffer_index = self.context.triangles_buffers.len() - 1;
        } else {
            self.current_buffer_index -= 1;
        }
    }

    pub fn current_triangles(&self) -> &Vec<TriangleData> {
        &self.context.triangles_buffers[self.current_buffer_index]
    }

    pub fn cursor(&self) -> usize {
        self.current_buffer_index
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

pub const BILLBOARD_DEFAULT_SCALE: Vec3 = Vec3::splat(0.001);

pub fn update_triangles_debugs_entities(
    mut commands: Commands,
    query_labels: Query<Entity, With<TriangleDebugLabel>>,
    debug_vert_data: Res<TrianglesDebugData>,
    mut debug_data_updates_events: EventReader<TriangleDebugDataUpdated>,
) {
    if !debug_data_updates_events.is_empty() {
        info!(
            "Triangle Debug Cursor set to step {}, {} triangles",
            debug_vert_data.cursor(),
            debug_vert_data.current_triangles().len()
        );
        debug_data_updates_events.clear();
    } else {
        return;
    }

    let triangles = debug_vert_data.current_triangles();
    let vertices = &debug_vert_data.vertices;

    for label in query_labels.iter() {
        commands.entity(label).despawn_recursive();
    }

    // Do not spawn labels entities if they are disabled
    if !debug_vert_data.with_labels {
        return;
    }

    for (index, t) in triangles.iter().enumerate() {
        let color = COLORS[index % COLORS.len()];
        let (v1, v2, v3) = (vertices[t.v1()], vertices[t.v2()], vertices[t.v3()]);
        let center = (v1 + v2 + v3) / 3.;

        let v1v2 = v2 - v1;
        let v1v3 = v3 - v1;
        let v2v3 = v3 - v2;
        // We could also get triangle_area = 0.5 * v1v2.cross(v1v3).length();
        let min_side_length = v1v2.length().min(v1v3.length()).min(v2v3.length());
        let billboard_scale = BILLBOARD_DEFAULT_SCALE * min_side_length * Vec3::ONE;

        commands.spawn((
            BillboardTextBundle {
                transform: Transform::from_translation(center).with_scale(billboard_scale),
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
            let v_display_pos = v + (center - v) * 0.15;
            commands.spawn((
                BillboardTextBundle {
                    transform: Transform::from_translation(v_display_pos)
                        .with_scale(billboard_scale),
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
