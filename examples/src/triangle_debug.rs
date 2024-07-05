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
    log::{error, info},
    math::primitives::Direction3d,
    pbr::{AlphaMode, MaterialMeshBundle, MaterialPlugin, StandardMaterial},
    prelude::default,
    render::{color::Color, mesh::Mesh},
    text::{Text, TextSection, TextStyle},
    transform::components::Transform,
    utils::HashSet,
};
use bevy_mod_billboard::{plugin::BillboardPlugin, BillboardTextBundle};
use ghx_constrained_delaunay::{
    debug::{DebugContext, DebugSnapshot, EventInfo},
    triangulation::{
        edge_from_semi_infinite_edge, CONTAINER_TRIANGLE_TOP_VERTEX_INDEX,
        CONTAINER_TRIANGLE_VERTICES, DELTA_VALUE, INFINITE_VERTS_X_DELTAS, INFINITE_VERTS_Y_DELTAS,
    },
    types::{
        next_clockwise_vertex_index, next_counter_clockwise_vertex_index, Edge, Float,
        TriangleData, TriangleId, TriangleVertexIndex, Triangles, Vector3, Vertex,
        NEXT_CCW_VERTEX_INDEX, NEXT_CW_VERTEX_INDEX, VERT_1, VERT_2, VERT_3,
    },
};
use glam::{Quat, Vec2, Vec3};

use crate::lines::{LineList, LineMaterial};

use super::{COLORS, DEBUG_LABEL_FONT_SIZE};

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
                update_triangles_debug_entities,
                draw_triangles_debug_data_gizmos,
            )
                .chain(),
        );
    }
}

pub fn extend_displayed_vertices_with_container_vertice(
    displayed_vertices: &mut Vec<Vector3>,
    plane_normal: Vector3,
    debug_context: &DebugContext,
    detransform: bool,
) {
    // Add the container triangle positions to the debug view buffer
    // Calculate container triangle vertices position in world position
    let mut container_vertices: Vec<Vector3> = CONTAINER_TRIANGLE_VERTICES
        .iter()
        .map(|v| Vector3::new(v.x, v.y, 0.))
        .collect();
    for v in container_vertices.iter_mut() {
        v.x = (debug_context.scale_factor * v.x) + debug_context.x_min;
        v.y = (debug_context.scale_factor * v.y) + debug_context.y_min;
    }
    if detransform {
        let e1 = (displayed_vertices[0] - displayed_vertices[1]).normalize();
        let e2 = e1.cross(-plane_normal);
        let _e3 = plane_normal;
        let plane_origin = displayed_vertices[1];
        let basis = Basis {
            plane_origin: plane_origin.as_vec3(),
            e1: e2.as_vec3(),
            e2: e2.as_vec3(),
        };
        for v in container_vertices.iter_mut() {
            *v = plane_origin + v.x * e1 + v.y * e2 // + 0. * e3
        }
    }

    displayed_vertices.extend(container_vertices);
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
    // TODO Rework
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

    pub fn new_with_constraintss(
        vertices: Vec<Vector3>,
        constraints: &Vec<Edge>,
        context: DebugContext,
    ) -> Self {
        Self {
            vertices,
            constraints: HashSet::from_iter(constraints.iter().cloned()),
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
}

#[derive(Resource, Default, Clone, Copy)]
pub struct TrianglesDebugViewConfig {
    pub label_mode: LabelMode,
    pub vertex_label_mode: VertexLabelMode,
    pub triangles_draw_mode: TrianglesDrawMode,
}
impl TrianglesDebugViewConfig {
    pub fn new(
        label_mode: LabelMode,
        vertex_label_mode: VertexLabelMode,
        draw_mode: TrianglesDrawMode,
    ) -> Self {
        Self {
            label_mode,
            vertex_label_mode,
            triangles_draw_mode: draw_mode,
        }
    }
}

impl TrianglesDebugData {
    pub fn advance_cursor(&mut self) -> usize {
        self.current_buffer_index = (self.current_buffer_index + 1) % self.context.snapshots.len();
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

    pub fn current_triangles(&self) -> &Triangles {
        &self.context.snapshots[self.current_buffer_index].triangles
    }

    pub fn current_snapshot(&self) -> &DebugSnapshot {
        &self.context.snapshots[self.current_buffer_index]
    }

    pub fn cursor(&self) -> usize {
        self.current_buffer_index
    }
}

#[derive(Resource)]
pub struct TriangleDebugAssets {
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
        color_materials.push(line_materials.add(LineMaterial { color: *color }));
    }
    commands.insert_resource(TriangleDebugAssets {
        color_materials,
        contour_materal: materials.add(StandardMaterial {
            base_color: Color::ANTIQUE_WHITE,
            alpha_mode: AlphaMode::Opaque,
            emissive: Color::rgb_linear(9800., 9200.0, 8400.0),
            ..default()
        }),
        interior_line_materal: materials.add(StandardMaterial {
            base_color: Color::ALICE_BLUE.with_a(0.1),
            alpha_mode: AlphaMode::Blend,
            emissive: Color::rgb_linear(0.0, 8000.0, 10000.0),
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

pub const BILLBOARD_DEFAULT_SCALE: Vec3 = Vec3::splat(0.0005);

pub fn update_triangles_debug_entities(
    mut commands: Commands,
    debug_assets: Res<TriangleDebugAssets>,
    view_config: Res<TrianglesDebugViewConfig>,
    mut debug_data: ResMut<TrianglesDebugData>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut debug_data_updates_events: EventReader<TriangleDebugCursorUpdate>,
    triangles_query: Query<Entity, With<TriangleDebugEntity>>,
) {
    if !debug_data_updates_events.is_empty() {
        debug_data_updates_events.clear();
    } else {
        return;
    }
    let snapshot = debug_data.current_snapshot();
    info!(
        "Cursor set to snapshot n°{}, step {}, phase {:?}, {} changes, {} triangles",
        debug_data.cursor(),
        snapshot.step,
        snapshot.triangulation_phase,
        snapshot.changed_ids.len(),
        snapshot.triangles.count()
    );
    // info!("triangles {:?}", snapshot.triangles.buffer);

    let vertices = &debug_data.vertices;

    // Despawn all previous entities
    for entity in triangles_query.iter() {
        commands.entity(entity).despawn_recursive();
    }

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
    debug_assets: &Res<TriangleDebugAssets>,
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

pub fn spawn_label(
    commands: &mut Commands,
    vertex_label_mode: VertexLabelMode,
    vertices: &Vec<Vector3>,
    basis: &Option<Basis>,
    triangle_id: TriangleId,
    triangle: &TriangleData,
) {
    let color = COLORS[(triangle_id as usize) % COLORS.len()];

    let mut infinite_verts = Vec::new();
    if triangle.v1() as usize >= vertices.len() {
        infinite_verts.push(VERT_1);
    }
    if triangle.v2() as usize >= vertices.len() {
        infinite_verts.push(VERT_2);
    }
    if triangle.v3() as usize >= vertices.len() {
        infinite_verts.push(VERT_3);
    }

    let (v1, v2, v3) = if infinite_verts.is_empty() {
        (
            vertices[triangle.v1() as usize].as_vec3(),
            vertices[triangle.v2() as usize].as_vec3(),
            vertices[triangle.v3() as usize].as_vec3(),
        )
    } else if infinite_verts.len() == 1 {
        let infinite_vert_local_index =
            (triangle.v(infinite_verts[0]) as usize - vertices.len()) as TriangleVertexIndex;

        let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_a_id = triangle.v(finite_vert_a_index);
        let finite_vertex_a = vertices[finite_vert_a_id as usize];
        let mut finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex_a, infinite_vert_local_index);

        let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_b_id = triangle.v(finite_vert_b_index);
        let finite_vertex_b = vertices[finite_vert_b_id as usize];
        let mut finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex_b, infinite_vert_local_index);

        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }

        if infinite_verts[0] == VERT_1 {
            (
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v2() as usize].as_vec3(),
                vertices[triangle.v3() as usize].as_vec3(),
            )
        } else if infinite_verts[0] == VERT_2 {
            (
                vertices[triangle.v1() as usize].as_vec3(),
                (finite_tip_a + finite_tip_b) / 2.,
                vertices[triangle.v3() as usize].as_vec3(),
            )
        } else {
            (
                vertices[triangle.v1() as usize].as_vec3(),
                vertices[triangle.v2() as usize].as_vec3(),
                (finite_tip_a + finite_tip_b) / 2.,
            )
        }
    } else if infinite_verts.len() == 2 {
        let finite_v_index = 3 - (infinite_verts[0] + infinite_verts[1]);
        let finite_vert_id = triangle.v(finite_v_index);
        let finite_vertex = vertices[finite_vert_id as usize];

        let infinite_vert_a_local_index = (triangle.v(next_clockwise_vertex_index(finite_v_index))
            as usize
            - vertices.len()) as TriangleVertexIndex;
        let mut finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);

        let infinite_vert_b_local_index =
            (triangle.v(next_counter_clockwise_vertex_index(finite_v_index)) as usize
                - vertices.len()) as TriangleVertexIndex;
        let mut finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_b_local_index);

        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }

        if finite_v_index == VERT_1 {
            (finite_vertex.as_vec3(), finite_tip_a, finite_tip_b)
        } else if finite_v_index == VERT_2 {
            (finite_tip_b, finite_vertex.as_vec3(), finite_tip_a)
        } else {
            (finite_tip_a, finite_tip_b, finite_vertex.as_vec3())
        }
    } else {
        error!("Handle 2 infinite vertices for labels");
        return;
    };

    let center = (v1 + v2 + v3) / 3.;

    let v1v2 = v2 - v1;
    let v1v3 = v3 - v1;
    let v2v3 = v3 - v2;
    // We could also use triangle_area = 0.5 * v1v2.cross(v1v3).length();
    let max_side_length = v1v2.length().max(v1v3.length()).max(v2v3.length()) as f32;
    let billboard_scale = BILLBOARD_DEFAULT_SCALE * max_side_length * Vec3::ONE;

    commands.spawn((
        BillboardTextBundle {
            transform: Transform::from_translation(center).with_scale(billboard_scale),
            text: Text::from_sections([TextSection {
                value: format!("t{}", triangle_id.to_string()),
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

    for (v, v_id, label) in vec![
        (v1, triangle.v1(), "v1"),
        (v2, triangle.v2(), "v2"),
        (v3, triangle.v3(), "v3"),
    ] {
        let v_display_pos = v + (center - v) * 0.15;
        let vertex_label = match vertex_label_mode {
            VertexLabelMode::LocalIndex => label.to_string(),
            VertexLabelMode::GlobalIndex => {
                if (v_id as usize) < vertices.len() {
                    v_id.to_string()
                } else {
                    format!("inf_{}", v_id.to_string())
                }
            }
        };
        commands.spawn((
            BillboardTextBundle {
                transform: Transform::from_translation(v_display_pos).with_scale(billboard_scale),
                text: Text::from_sections([TextSection {
                    value: vertex_label,
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

pub fn draw_triangles_debug_data_gizmos(
    mut gizmos: Gizmos,
    debug_data: Res<TrianglesDebugData>,
    view_config: Res<TrianglesDebugViewConfig>,
) {
    let snapshot = &debug_data.current_snapshot();
    let vertices = &debug_data.vertices;

    // Draw triangles
    match view_config.triangles_draw_mode {
        TrianglesDrawMode::AllAsGizmos => {
            for (triangle_id, triangle) in snapshot.triangles.buffer().iter().enumerate() {
                draw_triangle_gizmo(
                    &mut gizmos,
                    vertices,
                    &debug_data.basis,
                    triangle_id as TriangleId,
                    triangle,
                );
            }
        }
        TrianglesDrawMode::ChangedAsGizmos => {
            for triangle_id in snapshot.changed_ids.iter() {
                let triangle = &snapshot.triangles.get(*triangle_id);
                draw_triangle_gizmo(
                    &mut gizmos,
                    vertices,
                    &debug_data.basis,
                    *triangle_id,
                    triangle,
                );
            }
        }
        _ => (),
    }

    // Draw last changes boundaries
    gizmos.rect(
        debug_data.current_changes_bounds.0,
        Quat::IDENTITY,
        debug_data.current_changes_bounds.1,
        Color::WHITE,
    );

    // Draw specific event info
    if let EventInfo::Split(vertex_id) = snapshot.event {
        // Set circle radius as proportional to the changes bounding box
        let circle_radius = debug_data
            .current_changes_bounds
            .1
            .x
            .min(debug_data.current_changes_bounds.1.y)
            / 10.;
        gizmos.circle(
            vertices[vertex_id as usize].as_vec3(),
            Direction3d::Z,
            circle_radius,
            Color::RED,
        );
    }
}

pub fn draw_triangle_gizmo(
    gizmos: &mut Gizmos,
    vertices: &Vec<Vector3>,
    basis: &Option<Basis>,
    triangle_id: TriangleId,
    triangle: &TriangleData,
) {
    let color = COLORS[triangle_id as usize % COLORS.len()];

    let mut infinite_verts = Vec::new();
    if triangle.v1() as usize >= vertices.len() {
        infinite_verts.push(VERT_1);
    }
    if triangle.v2() as usize >= vertices.len() {
        infinite_verts.push(VERT_2);
    }
    if triangle.v3() as usize >= vertices.len() {
        infinite_verts.push(VERT_3);
    }

    let linestrip = if infinite_verts.is_empty() {
        let (v1, v2, v3) = (
            vertices[triangle.v1() as usize].as_vec3(),
            vertices[triangle.v2() as usize].as_vec3(),
            vertices[triangle.v3() as usize].as_vec3(),
        );
        vec![v1, v2, v3, v1]
    } else if infinite_verts.len() == 1 {
        let infinite_vert_local_index =
            (triangle.v(infinite_verts[0]) as usize - vertices.len()) as TriangleVertexIndex;

        let finite_vert_a_index = NEXT_CW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_a_id = triangle.v(finite_vert_a_index);
        let finite_vertex_a = vertices[finite_vert_a_id as usize];
        let mut finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex_a, infinite_vert_local_index);

        let finite_vert_b_index = NEXT_CCW_VERTEX_INDEX[infinite_verts[0] as usize];
        let finite_vert_b_id = triangle.v(finite_vert_b_index);
        let finite_vertex_b = vertices[finite_vert_b_id as usize];
        let mut finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex_b, infinite_vert_local_index);

        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }
        // info!("draw gizmos, infinite_verts.len() == 1");

        vec![
            finite_tip_a,
            finite_vertex_a.as_vec3(),
            finite_vertex_b.as_vec3(),
            finite_tip_b,
        ]
    } else if infinite_verts.len() == 2 {
        //TODO
        let finite_v_index = 3 - (infinite_verts[0] + infinite_verts[1]);
        let finite_vert_id = triangle.v(finite_v_index);
        let finite_vertex = vertices[finite_vert_id as usize];

        let infinite_vert_a_local_index = (triangle.v(next_clockwise_vertex_index(finite_v_index))
            as usize
            - vertices.len()) as TriangleVertexIndex;
        let mut finite_tip_a =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_a_local_index);

        let infinite_vert_b_local_index =
            (triangle.v(next_counter_clockwise_vertex_index(finite_v_index)) as usize
                - vertices.len()) as TriangleVertexIndex;
        let mut finite_tip_b =
            finite_vertex_from_semi_infinite_edge(finite_vertex, infinite_vert_b_local_index);

        // // TODO Really quick & dirty. Does not work, finite_vertex is not transformed. Would need transformed finite_vertex -> point extrapolation -> detransform
        // if let Some(basis) = basis {
        //     finite_tip_a = detransfrom(finite_tip_a, basis);
        //     finite_tip_b = detransfrom(finite_tip_b, basis);
        // }

        vec![finite_tip_a, finite_vertex.as_vec3(), finite_tip_b]
    } else {
        vec![]
    };
    gizmos.linestrip(linestrip, color);
    // TODO Could draw something for 3 infinite vertices
}

pub fn finite_vertex_from_semi_infinite_edge(
    finite_vertex: Vector3,
    infinite_vert_local_index: TriangleVertexIndex,
) -> Vec3 {
    // TODO Might need a detransformation
    Vec3::new(
        // TODO 100.0 Chane to a factor dependent on the triangle/triangulation
        (finite_vertex.x + 20. * INFINITE_VERTS_X_DELTAS[infinite_vert_local_index as usize])
            as f32,
        (finite_vertex.y + 20. * INFINITE_VERTS_Y_DELTAS[infinite_vert_local_index as usize])
            as f32,
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
