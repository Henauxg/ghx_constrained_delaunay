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
    math::{Vec2, Vec3},
    pbr::{MaterialPlugin, StandardMaterial},
    prelude::{default, AlphaMode, Component, Event, Resource},
    utils::HashSet,
};
use bevy_mod_billboard::plugin::BillboardPlugin;
use ghx_constrained_delaunay::{
    debug::{DebugContext, DebugSnapshot},
    types::{Edge, TriangleId},
};
use gizmos::draw_triangles_debug_data_gizmos;
use label::spawn_label;
use meshes::{
    spawn_all_as_mesh_batches, spawn_contour_and_interior_meshes, SpawnContourAndInteriorMeshes,
    SpawnMeshBatches,
};

use crate::lines::LineMaterial;

use super::COLORS;

pub mod gizmos;
pub mod infinite;
pub mod label;
pub mod meshes;

pub const DEFAULT_TRIANGLES_BATCH_SIZE: usize = 15;

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
        app.observe(spawn_all_as_mesh_batches);
        app.observe(spawn_contour_and_interior_meshes);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Basis {
    e1: Vec3,
    e2: Vec3,
}

#[derive(Default, Copy, Clone, Debug, PartialEq, Eq)]
pub enum TrianglesDrawMode {
    #[default]
    /// Draws all the triangles, as gizmos.
    AllAsGizmos,
    // Draws only the triangles that changed, as gizmos.
    ChangedAsGizmos,
    /// Draw all the triangles as meshes. Groups `batch_size` triangles per mesh.
    ///
    /// Increase `batch_size` for more performances.
    AllAsMeshBatches {
        batch_size: usize,
    },
    /// Only displays finite triangles.
    ///
    /// 1 mesh for the constrained edges, 1 mesh for all the other lines.
    ContoursAndInteriorAsMeshes,
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
    pub vertices: Vec<Vec3>,
    pub constraints: HashSet<Edge>,
    pub context: DebugContext,
    pub basis: Basis,
    // State
    pub current_buffer_index: usize,
    pub current_changes_bounds: (Vec3, Vec2),
}

impl Default for TrianglesDebugData {
    fn default() -> Self {
        Self {
            vertices: Default::default(),
            constraints: HashSet::new(),
            context: Default::default(),
            basis: Basis {
                e1: Vec3::X,
                e2: Vec3::Y,
            },
            current_buffer_index: 0,
            current_changes_bounds: (Vec3::ZERO, Vec2::ZERO),
        }
    }
}
impl TrianglesDebugData {
    pub fn new(vertices: Vec<Vec3>, context: DebugContext) -> Self {
        Self {
            vertices,
            context,
            ..Default::default()
        }
    }

    pub fn new_with_constrained_edges(
        vertices: Vec<Vec3>,
        constrained_edges: &[Edge],
        context: DebugContext,
    ) -> Self {
        Self {
            vertices,
            constraints: HashSet::from_iter(constrained_edges.iter().cloned()),
            context,
            ..Default::default()
        }
    }

    pub fn with_basis(mut self, basis: Basis) -> Self {
        self.basis = basis;
        self
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

    pub fn set_cursor_to_first_snapshot(&mut self) {
        self.current_buffer_index = 0;
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
    view_config: Res<TrianglesDebugViewConfig>,
    mut debug_data: ResMut<TrianglesDebugData>,
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
            commands.trigger(SpawnMeshBatches(batch_size))
        }
        TrianglesDrawMode::ContoursAndInteriorAsMeshes => {
            commands.trigger(SpawnContourAndInteriorMeshes);
        }
        _ => (),
    }

    // Spawn label entities if enabled
    match view_config.label_mode {
        LabelMode::All => {
            for (triangle_id, triangle_data) in snapshot.triangles.buffer().iter().enumerate() {
                spawn_label(
                    &mut commands,
                    view_config.vertex_label_mode,
                    vertices,
                    triangle_id as TriangleId,
                    triangle_data,
                    &debug_data.basis,
                );
            }
        }
        LabelMode::Changed => {
            for triangle_id in snapshot.changed_ids.iter() {
                let triangle_data = &snapshot.triangles.get(*triangle_id);
                spawn_label(
                    &mut commands,
                    view_config.vertex_label_mode,
                    vertices,
                    *triangle_id,
                    triangle_data,
                    &debug_data.basis,
                );
            }
        }
        LabelMode::None => (),
    };

    // Register snapshot changed data
    let (mut x_min, mut x_max, mut y_min, mut y_max) = (f32::MAX, f32::MIN, f32::MAX, f32::MIN);
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
            view_config.triangles_draw_mode = TrianglesDrawMode::ContoursAndInteriorAsMeshes
        }
        TrianglesDrawMode::ContoursAndInteriorAsMeshes => {
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
