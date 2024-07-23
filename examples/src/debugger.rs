use bevy::{
    app::{Plugin, Startup, Update},
    asset::{Assets, Handle},
    color::{
        palettes::css::{ALICE_BLUE, ANTIQUE_WHITE},
        Alpha, Color, LinearRgba,
    },
    ecs::{
        entity::Entity,
        event::{EventReader, EventWriter},
        query::With,
        schedule::IntoSystemConfigs,
        system::{Commands, Query, Res, ResMut},
    },
    hierarchy::DespawnRecursiveExt,
    input::{common_conditions::input_just_pressed, keyboard::KeyCode, ButtonInput},
    log::info,
    pbr::{MaterialMeshBundle, MaterialPlugin, StandardMaterial},
    prelude::{default, AlphaMode, Component, Event, Resource},
    render::mesh::Mesh,
    utils::HashSet,
};
use bevy_mod_billboard::plugin::BillboardPlugin;
use ghx_constrained_delaunay::glam::{Vec2, Vec3};
use ghx_constrained_delaunay::{
    debug::{DebugContext, DebugSnapshot},
    infinite::{INFINITE_VERTS_X_DELTAS, INFINITE_VERTS_Y_DELTAS},
    types::{Edge, Float, TriangleId, TriangleVertexIndex, Vector3},
};
use gizmos::draw_triangles_debug_data_gizmos;
use label::spawn_label;

use crate::lines::{LineList, LineMaterial};

use super::COLORS;

pub mod gizmos;
pub mod label;

#[derive(Default)]
pub struct TriangleDebugPlugin;
impl Plugin for TriangleDebugPlugin {
    fn build(&self, app: &mut bevy::prelude::App) {
        app.add_plugins((BillboardPlugin, MaterialPlugin::<LineMaterial>::default()));
        app.init_resource::<TrianglesDebugViewConfig>();
        app.add_event::<TriangleDebugCursorUpdate>();
        app.add_systems(Startup, setup);
        app.add_systems(
            Update,
            (
                switch_vertex_label_mode.run_if(input_just_pressed(KeyCode::F2)),
                switch_label_mode.run_if(input_just_pressed(KeyCode::F3)),
                switch_triangles_draw_mode.run_if(input_just_pressed(KeyCode::F4)),
                update_triangles_debug_index,
                update_triangles_debug_view,
                draw_triangles_debug_data_gizmos,
            )
                .chain(),
        );
    }
}

pub struct Basis {
    plane_origin: Vec3,
    e1: Vec3,
    e2: Vec3,
}
pub fn detransfrom(point: Vec3, basis: &Basis) -> Vec3 {
    basis.plane_origin + point.x * basis.e1 + point.y * basis.e2 // + 0. * e3
}

#[derive(Default, Copy, Clone, Debug, PartialEq, Eq)]
pub enum TrianglesDrawMode {
    #[default]
    AllAsGizmos,
    ChangedAsGizmos,
    AllAsMeshBatches {
        batch_size: usize,
    },
    AllAsContourAndInteriorMeshes,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum LabelMode {
    #[default]
    All,
    Changed,
    None,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum VertexLabelMode {
    #[default]
    LocalIndex,
    GlobalIndex,
}

#[derive(Resource)]
pub struct TrianglesDebugData {
    // Inputs
    pub vertices: Vec<Vector3>,
    pub constraints: HashSet<Edge>,
    pub context: DebugContext,
    // TODO Rework. The debugger needs a collection of vertices transformed in the new basis, so that it can extrapolate from the transformed vertex, and then detransform the result.
    pub basis: Option<Basis>,
    // State
    pub current_buffer_index: usize,
    pub current_changes_bounds: (Vec3, Vec2),
}
impl TrianglesDebugData {
    pub fn new(vertices: Vec<Vector3>, context: DebugContext) -> Self {
        Self {
            vertices,
            constraints: HashSet::new(),
            context,
            basis: None,
            current_buffer_index: 0,
            current_changes_bounds: (Vec3::ZERO, Vec2::ZERO),
        }
    }

    pub fn new_with_constrained_edges(
        vertices: Vec<Vector3>,
        constrained_edges: &Vec<Edge>,
        context: DebugContext,
    ) -> Self {
        Self {
            vertices,
            constraints: HashSet::from_iter(constrained_edges.iter().cloned()),
            context,
            basis: None,
            current_buffer_index: 0,
            current_changes_bounds: (Vec3::ZERO, Vec2::ZERO),
        }
    }

    pub fn with_basis_transformation(&mut self, plane_normal: Vector3) {
        let e1 = (self.vertices[0] - self.vertices[1]).normalize();
        let e2 = e1.cross(-plane_normal);
        let _e3 = plane_normal;
        let plane_origin = self.vertices[1];
        self.basis = Some(Basis {
            plane_origin: plane_origin.as_vec3(),
            e1: e2.as_vec3(),
            e2: e2.as_vec3(),
        });
    }

    pub fn advance_cursor(&mut self) -> usize {
        if self.context.snapshots.len() > 0 {
            self.current_buffer_index =
                (self.current_buffer_index + 1) % self.context.snapshots.len();
        }
        self.current_buffer_index
    }

    pub fn move_back_cursor(&mut self) -> usize {
        if self.current_buffer_index == 0 {
            self.current_buffer_index = self.context.snapshots.len() - 1;
        } else {
            self.current_buffer_index -= 1;
        }
        self.current_buffer_index
    }

    pub fn set_cursor_to_last_snapshot(&mut self) {
        self.current_buffer_index = self.context.snapshots.len() - 1;
    }

    pub fn current_snapshot(&self) -> Option<&DebugSnapshot> {
        if self.context.snapshots.len() > self.current_buffer_index {
            Some(&self.context.snapshots[self.current_buffer_index])
        } else {
            None
        }
    }

    pub fn cursor(&self) -> usize {
        self.current_buffer_index
    }
}

#[derive(Resource, Default, Clone, Copy)]
pub struct TrianglesDebugViewConfig {
    pub label_mode: LabelMode,
    pub vertex_label_mode: VertexLabelMode,
    pub triangles_draw_mode: TrianglesDrawMode,
    pub log_cursor_updates: bool,
}
impl TrianglesDebugViewConfig {
    pub fn new(
        label_mode: LabelMode,
        vertex_label_mode: VertexLabelMode,
        draw_mode: TrianglesDrawMode,
        log_cursor_updates: bool,
    ) -> Self {
        Self {
            label_mode,
            vertex_label_mode,
            triangles_draw_mode: draw_mode,
            log_cursor_updates,
        }
    }
}

#[derive(Resource)]
pub struct TriangleDebugViewAssets {
    pub color_materials: Vec<Handle<LineMaterial>>,
    pub interior_line_materal: Handle<StandardMaterial>,
    pub contour_materal: Handle<StandardMaterial>,
}

pub fn setup(
    mut commands: Commands,
    mut line_materials: ResMut<Assets<LineMaterial>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mut color_materials = Vec::new();
    for color in COLORS.iter() {
        color_materials.push(line_materials.add(LineMaterial {
            color: (*color).into(),
        }));
    }
    commands.insert_resource(TriangleDebugViewAssets {
        color_materials,
        contour_materal: materials.add(StandardMaterial {
            base_color: Color::Srgba(ANTIQUE_WHITE),
            alpha_mode: AlphaMode::Opaque,
            emissive: LinearRgba::rgb(18.54, 17.46, 15.84),
            ..default()
        }),
        interior_line_materal: materials.add(StandardMaterial {
            base_color: Color::Srgba(ALICE_BLUE).with_alpha(0.1),
            alpha_mode: AlphaMode::Blend,
            emissive: LinearRgba::rgb(0.0, 2.25, 3.),
            ..default()
        }),
    });
}

#[derive(Event)]
pub struct TriangleDebugCursorUpdate;

pub fn update_triangles_debug_index(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    mut debug_vert_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    if keyboard_input.just_pressed(KeyCode::ArrowRight) {
        debug_vert_data.advance_cursor();
        debug_data_updates_events.send(TriangleDebugCursorUpdate);
    } else if keyboard_input.just_pressed(KeyCode::ArrowLeft) {
        debug_vert_data.move_back_cursor();
        debug_data_updates_events.send(TriangleDebugCursorUpdate);
    }
}

#[derive(Component)]
pub struct TriangleDebugEntity;

pub fn update_triangles_debug_view(
    mut commands: Commands,
    debug_assets: Res<TriangleDebugViewAssets>,
    view_config: Res<TrianglesDebugViewConfig>,
    mut debug_data: ResMut<TrianglesDebugData>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut debug_data_updates_events: EventReader<TriangleDebugCursorUpdate>,
    triangles_query: Query<Entity, With<TriangleDebugEntity>>,
) {
    if debug_data_updates_events.is_empty() {
        return;
    } else {
        debug_data_updates_events.clear();
    }

    // Despawn all previous entities
    for entity in triangles_query.iter() {
        commands.entity(entity).despawn_recursive();
    }

    let Some(snapshot) = debug_data.current_snapshot() else {
        return;
    };
    if view_config.log_cursor_updates {
        info!(
            "Cursor set to snapshot n°{}, step {}, phase {:?}, {} changes, {} triangles",
            debug_data.cursor(),
            snapshot.step,
            snapshot.triangulation_phase,
            snapshot.changed_ids.len(),
            snapshot.triangles.count()
        );
    }

    let vertices = &debug_data.vertices;

    // Spawn triangle meshes if enabled
    // TODO With a better diff from the debug context, could only spawn the difference.
    match view_config.triangles_draw_mode {
        TrianglesDrawMode::AllAsMeshBatches { batch_size } => {
            let mut lines = Some(vec![]);
            let mut batch_index = 0;
            for (index, t) in snapshot.triangles.buffer().iter().enumerate() {
                let (v1, v2, v3) = (
                    vertices[t.v1() as usize].as_vec3(),
                    vertices[t.v2() as usize].as_vec3(),
                    vertices[t.v3() as usize].as_vec3(),
                );
                lines
                    .as_mut()
                    .unwrap()
                    .extend(vec![(v1, v2), (v2, v3), (v3, v1)]);

                if index % batch_size == batch_size - 1 {
                    spawn_triangle_mesh_batch(
                        &mut commands,
                        &mut meshes,
                        &debug_assets,
                        batch_index,
                        lines.take().unwrap(),
                    );
                    lines = Some(vec![]);
                    batch_index += 1;
                }
            }
            let lines = lines.unwrap();
            if !lines.is_empty() {
                spawn_triangle_mesh_batch(
                    &mut commands,
                    &mut meshes,
                    &debug_assets,
                    batch_index,
                    lines,
                );
            }
            info!("Spawned {} triangle mesh batches", batch_index);
        }
        TrianglesDrawMode::AllAsContourAndInteriorMeshes => {
            let mut contour = vec![];
            let mut interior_lines = vec![];
            for t in snapshot.triangles.buffer().iter() {
                for edge in t.edges() {
                    let collection = if debug_data.constraints.contains(&edge)
                        || debug_data.constraints.contains(&edge.opposite())
                    {
                        &mut contour
                    } else {
                        &mut interior_lines
                    };
                    collection.push((
                        vertices[edge.from as usize].as_vec3(),
                        vertices[edge.to as usize].as_vec3(),
                    ));
                }
            }
            commands.spawn((
                MaterialMeshBundle {
                    mesh: meshes.add(LineList { lines: contour }),
                    material: debug_assets.contour_materal.clone(),
                    ..default()
                },
                TriangleDebugEntity,
            ));
            commands.spawn((
                MaterialMeshBundle {
                    mesh: meshes.add(LineList {
                        lines: interior_lines,
                    }),
                    material: debug_assets.interior_line_materal.clone(),
                    ..default()
                },
                TriangleDebugEntity,
            ));
        }
        _ => (),
    }

    // Spawn label entites if enabled
    match view_config.label_mode {
        LabelMode::All => {
            for (index, triangle_data) in snapshot.triangles.buffer().iter().enumerate() {
                spawn_label(
                    &mut commands,
                    view_config.vertex_label_mode,
                    vertices,
                    &debug_data.basis,
                    index as TriangleId,
                    triangle_data,
                );
            }
        }
        LabelMode::Changed => {
            for index in snapshot.changed_ids.iter() {
                let triangle_data = &snapshot.triangles.get(*index);
                spawn_label(
                    &mut commands,
                    view_config.vertex_label_mode,
                    vertices,
                    &debug_data.basis,
                    *index,
                    triangle_data,
                );
            }
        }
        LabelMode::None => (),
    };

    // Register snapshot changed data
    let (mut x_min, mut x_max, mut y_min, mut y_max) =
        (Float::MAX, Float::MIN, Float::MAX, Float::MIN);
    for index in snapshot.changed_ids.iter() {
        let t = &snapshot.triangles.get(*index);

        for v_id in t.verts.iter() {
            if (*v_id as usize) < vertices.len() {
                let v = vertices[*v_id as usize];
                if v.x < x_min {
                    x_min = v.x;
                }
                if v.x > x_max {
                    x_max = v.x;
                }
                if v.y < y_min {
                    y_min = v.y;
                }
                if v.y > y_max {
                    y_max = v.y;
                }
            } else {
                // TODO Infinite vert approx
            }
        }
    }
    let size = Vec2::new((x_max - x_min) as f32, (y_max - y_min) as f32);
    let pos = Vec3::new(x_min as f32 + size.x / 2., y_min as f32 + size.y / 2., 0.);
    debug_data.current_changes_bounds = (pos, size);
}

pub fn spawn_triangle_mesh_batch(
    commands: &mut Commands,
    meshes: &mut ResMut<Assets<Mesh>>,
    debug_assets: &Res<TriangleDebugViewAssets>,
    batch_index: usize,
    lines: Vec<(Vec3, Vec3)>,
) {
    commands.spawn((
        MaterialMeshBundle {
            mesh: meshes.add(LineList { lines }),
            material: debug_assets.color_materials[batch_index % COLORS.len()].clone(),
            ..default()
        },
        TriangleDebugEntity,
    ));
}

pub const SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR: Float = 50.;

pub fn finite_vertex_from_semi_infinite_edge(
    finite_vertex: Vector3,
    infinite_vert_local_index: TriangleVertexIndex,
) -> Vec3 {
    // TODO Might need a detransformation
    Vec3::new(
        // TODO Change to a factor dependent on the triangle/triangulation
        (finite_vertex.x
            + SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR
                * INFINITE_VERTS_X_DELTAS[infinite_vert_local_index as usize]) as f32,
        (finite_vertex.y
            + SEMI_INFINITE_EDGE_VISUAL_LEN_FACTOR
                * INFINITE_VERTS_Y_DELTAS[infinite_vert_local_index as usize]) as f32,
        0.,
    )
}

pub const DEFAULT_TRIANGLES_BATCH_SIZE: usize = 15;

pub fn switch_triangles_draw_mode(
    mut view_config: ResMut<TrianglesDebugViewConfig>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    match view_config.triangles_draw_mode {
        TrianglesDrawMode::AllAsGizmos => {
            view_config.triangles_draw_mode = TrianglesDrawMode::ChangedAsGizmos
        }
        TrianglesDrawMode::ChangedAsGizmos => {
            view_config.triangles_draw_mode = TrianglesDrawMode::AllAsMeshBatches {
                batch_size: DEFAULT_TRIANGLES_BATCH_SIZE,
            }
        }
        TrianglesDrawMode::AllAsMeshBatches { batch_size: _ } => {
            view_config.triangles_draw_mode = TrianglesDrawMode::AllAsContourAndInteriorMeshes
        }
        TrianglesDrawMode::AllAsContourAndInteriorMeshes => {
            view_config.triangles_draw_mode = TrianglesDrawMode::AllAsGizmos
        }
    };
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
    info!(
        "Triangles draw mode changed to {:?}",
        view_config.triangles_draw_mode
    );
}

pub fn switch_label_mode(
    mut view_config: ResMut<TrianglesDebugViewConfig>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    match view_config.label_mode {
        LabelMode::All => view_config.label_mode = LabelMode::Changed,
        LabelMode::Changed => view_config.label_mode = LabelMode::None,
        LabelMode::None => view_config.label_mode = LabelMode::All,
    };
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
    info!("Label mode changed to {:?}", view_config.label_mode);
}

pub fn switch_vertex_label_mode(
    mut view_config: ResMut<TrianglesDebugViewConfig>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    match view_config.vertex_label_mode {
        VertexLabelMode::GlobalIndex => view_config.vertex_label_mode = VertexLabelMode::LocalIndex,
        VertexLabelMode::LocalIndex => view_config.vertex_label_mode = VertexLabelMode::GlobalIndex,
    };
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
    info!(
        "Vertex label mode changed to {:?}",
        view_config.vertex_label_mode
    );
}
