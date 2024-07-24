use bevy::{
    asset::Assets,
    log::info,
    math::Vec3,
    pbr::MaterialMeshBundle,
    prelude::{Commands, Event, Mesh, Res, ResMut, Trigger},
    utils::default,
};
use ghx_constrained_delaunay::{
    infinite::collect_infinite_triangle_vertices, types::opposite_vertex_index_from_edge,
};

use crate::{
    infinite::{
        get_2_extrapolated_vertices_from_1_finite_vertex,
        get_2_extrapolated_vertices_from_2_finite_vertices,
    },
    lines::LineList,
    COLORS,
};

use super::{TriangleDebugEntity, TriangleDebugViewAssets, TrianglesDebugData};

#[derive(Event)]
pub struct SpawnContourAndInteriorMeshes;

#[derive(Event)]
pub struct SpawnMeshBatches(pub usize);

pub fn spawn_all_as_mesh_batches(
    spawn_event: Trigger<SpawnMeshBatches>,
    mut commands: Commands,
    debug_data: Res<TrianglesDebugData>,
    debug_assets: Res<TriangleDebugViewAssets>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let Some(snapshot) = debug_data.current_snapshot() else {
        return;
    };
    let vertices = &debug_data.vertices;
    let batch_size = spawn_event.event().0;

    let mut lines = Some(vec![]);
    let mut batch_index = 0;
    for (index, triangle) in snapshot.triangles.buffer().iter().enumerate() {
        let infinite_verts = collect_infinite_triangle_vertices(&triangle.verts);

        let triangle_lines = if infinite_verts.is_empty() {
            let (v1, v2, v3) = (
                vertices[triangle.v1() as usize].as_vec3(),
                vertices[triangle.v2() as usize].as_vec3(),
                vertices[triangle.v3() as usize].as_vec3(),
            );
            vec![(v1, v2), (v2, v3), (v3, v1)]
        } else if infinite_verts.len() == 1 {
            let (finite_vertex_a, finite_tip_a, finite_vertex_b, finite_tip_b) =
                get_2_extrapolated_vertices_from_2_finite_vertices(
                    vertices,
                    triangle,
                    &debug_data.basis,
                    infinite_verts[0],
                );
            vec![
                (finite_tip_a, finite_vertex_a),
                (finite_vertex_a, finite_vertex_b),
                (finite_vertex_b, finite_tip_b),
            ]
        } else if infinite_verts.len() == 2 {
            let (finite_tip_a, finite_vertex, finite_tip_b) =
                get_2_extrapolated_vertices_from_1_finite_vertex(
                    vertices,
                    triangle,
                    &debug_data.basis,
                    opposite_vertex_index_from_edge(infinite_verts[0], infinite_verts[1]),
                );
            vec![(finite_tip_a, finite_vertex), (finite_vertex, finite_tip_b)]
        } else {
            vec![]
            // TODO Could draw something for 3 infinite vertices
        };

        lines.as_mut().unwrap().extend(triangle_lines);

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

pub fn spawn_contour_and_interior_meshes(
    _spawn_event: Trigger<SpawnContourAndInteriorMeshes>,
    mut commands: Commands,
    debug_data: Res<TrianglesDebugData>,
    debug_assets: Res<TriangleDebugViewAssets>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let Some(snapshot) = debug_data.current_snapshot() else {
        return;
    };
    let vertices = &debug_data.vertices;

    let mut contour = vec![];
    let mut interior_lines = vec![];
    for t in snapshot.triangles.buffer().iter() {
        if t.is_pseudo_infinite() {
            continue;
        }
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
