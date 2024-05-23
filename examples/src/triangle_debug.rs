use bevy::{
    app::{Plugin, Startup, Update},
    asset::{Assets, Handle},
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
    input::{common_conditions::input_just_pressed, keyboard::KeyCode, ButtonInput},
    log::info,
    pbr::{MaterialMeshBundle, MaterialPlugin},
    prelude::default,
    render::mesh::Mesh,
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
};
use bevy_mod_billboard::{plugin::BillboardPlugin, BillboardTextBundle};
use ghx_constrained_delaunay::{
    triangulation::{DebugContext, CONTAINER_TRIANGLE_VERTICES},
    types::TriangleData,
};
use glam::Vec3;

use crate::lines::{LineList, LineMaterial};

use super::{COLORS, DEBUG_LABEL_FONT_SIZE};

#[derive(Resource, Default, Copy, Clone, Debug, PartialEq, Eq)]
pub enum DrawMode {
    #[default]
    Gizmos,
    MeshBatches {
        batch_size: usize,
    },
}

#[derive(Default)]
pub struct TriangleDebugPlugin {
    pub draw_mode: DrawMode,
}
impl Plugin for TriangleDebugPlugin {
    fn build(&self, app: &mut bevy::prelude::App) {
        app.add_plugins((BillboardPlugin, MaterialPlugin::<LineMaterial>::default()));
        app.insert_resource(self.draw_mode);
        app.add_event::<TriangleDebugDataUpdated>();
        app.add_systems(Startup, setup);
        app.add_systems(
            Update,
            (
                toggle_draw_mode.run_if(input_just_pressed(KeyCode::F4)),
                update_triangles_debug_index,
                update_triangles_debug_entities,
                draw_triangles_debug_data.run_if(gizmo_draw_mode),
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

#[derive(Resource)]
pub struct TriangleDebugAssets {
    pub color_materials: Vec<Handle<LineMaterial>>,
}

pub fn setup(mut commands: Commands, mut materials: ResMut<Assets<LineMaterial>>) {
    let mut color_materials = Vec::new();
    for color in COLORS.iter() {
        color_materials.push(materials.add(LineMaterial { color: *color }));
    }
    commands.insert_resource(TriangleDebugAssets { color_materials });
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
pub struct TriangleDebugEntity;

pub const BILLBOARD_DEFAULT_SCALE: Vec3 = Vec3::splat(0.001);

pub fn update_triangles_debug_entities(
    mut commands: Commands,
    draw_mode: Res<DrawMode>,
    debug_assets: Res<TriangleDebugAssets>,
    query_triangles: Query<Entity, With<TriangleDebugEntity>>,
    debug_vert_data: Res<TrianglesDebugData>,
    mut debug_data_updates_events: EventReader<TriangleDebugDataUpdated>,
    mut meshes: ResMut<Assets<Mesh>>,
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

    // Despawn all previous entities
    for entity in query_triangles.iter() {
        commands.entity(entity).despawn_recursive();
    }

    if let DrawMode::MeshBatches { batch_size } = *draw_mode {
        let mut lines = Some(vec![]);
        let mut batch_index = 0;
        // TODO get diffs from the debugger and spawn the diffs?
        for (index, t) in triangles.iter().enumerate() {
            let (v1, v2, v3) = (vertices[t.v1()], vertices[t.v2()], vertices[t.v3()]);
            lines
                .as_mut()
                .unwrap()
                .extend(vec![(v1, v2), (v2, v3), (v3, v1)]);

            if index % batch_size == batch_size - 1 {
                commands.spawn((
                    MaterialMeshBundle {
                        mesh: meshes.add(LineList {
                            lines: lines.take().unwrap(),
                        }),
                        material: debug_assets.color_materials[batch_index % COLORS.len()].clone(),
                        ..default()
                    },
                    TriangleDebugEntity,
                ));
                lines = Some(vec![]);
                batch_index += 1;
            }
        }
        // TODO Spawn remaining batch

        info!(
            "Spawned all triangle meshes: {} triangles in {} batches",
            triangles.len(),
            batch_index
        );
    }

    // Do not spawn labels entities if they are disabled
    if !debug_vert_data.with_labels {
        return;
    }

    // TODO Only for the diff would be interesting here
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
            TriangleDebugEntity,
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
                TriangleDebugEntity,
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

pub fn gizmo_draw_mode(draw_mode: Option<Res<DrawMode>>) -> bool {
    match draw_mode {
        Some(mode) => *mode == DrawMode::Gizmos,
        None => false,
    }
}

pub const DEFAULT_TRIANGLES_BATCH_SIZE: usize = 15;

pub fn toggle_draw_mode(mut draw_mode: ResMut<DrawMode>) {
    match *draw_mode {
        DrawMode::Gizmos => {
            *draw_mode = DrawMode::MeshBatches {
                batch_size: DEFAULT_TRIANGLES_BATCH_SIZE,
            }
        }
        DrawMode::MeshBatches { batch_size: _ } => *draw_mode = DrawMode::Gizmos,
    };
}
