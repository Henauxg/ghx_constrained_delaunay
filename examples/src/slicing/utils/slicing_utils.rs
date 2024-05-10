use std::{
    collections::{HashMap, HashSet},
    iter::zip,
};

use bevy_rapier3d::prelude::*;

use bevy::{
    math::Vec3A,
    prelude::*,
    reflect::Map,
    render::{
        mesh::{Indices, PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
    },
};
use ordered_float::OrderedFloat;
use rand::Rng;

use super::{mesh_mapping::MeshMapping, triangulation::VertexId};

fn get_mesh_center(
    meshes_handles: Query<&Handle<Mesh>>,
    entity: Entity,
    meshes_assets: ResMut<Assets<Mesh>>,
) -> Option<Vec3A> {
    if let Ok(source_object_mesh_handle) = meshes_handles.get(entity) {
        // TODO Fix unwraps
        let source_object_mesh = meshes_assets.get(source_object_mesh_handle).unwrap();
        Some(source_object_mesh.compute_aabb().unwrap().center)
    } else {
        None
    }
}

pub(crate) fn get_random_normalized_vec() -> Vec3A {
    let mut rng = rand::thread_rng();
    Vec3A::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()).normalize()
}

/// Operate the slice on a mesh
pub(crate) fn slice(
    commands: &mut Commands,
    meshes_assets: &mut ResMut<Assets<Mesh>>,
    mesh: &Mesh,
    materials: ResMut<Assets<StandardMaterial>>,
) {
    // get mesh center
    let mesh_center = mesh.compute_aabb().unwrap().center;

    // get random normalized vector for the cut plane
    let normal_vec = get_random_normalized_vec();

    // execute the slice
    internal_slice(
        commands,
        mesh_center,
        normal_vec,
        mesh,
        materials,
        meshes_assets,
    )
}

fn internal_slice(
    mut commands: &mut Commands,
    origin_point: Vec3A,
    normal_vec: Vec3A,
    main_mesh: &Mesh,
    materials: ResMut<Assets<StandardMaterial>>,
    meshes_assets: &mut ResMut<Assets<Mesh>>,
) {
    // get the top and bottom meshes from the slice operation and the respective pos
    let (top_mesh, bottom_mesh, pos) =
        compute_slice(&origin_point, &normal_vec, main_mesh).unwrap();

    // spawn the fragments from the meshes at the pos
    for mesh in vec![top_mesh, bottom_mesh] {
        let mesh_handle = meshes_assets.add(mesh.clone());
        spawn_fragments(&mesh, &mesh_handle, &materials, &mut commands, pos.into());
    }
}

pub fn compute_slice(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    main_mesh: &Mesh,
) -> Option<(Mesh, Mesh, Vec3A)> {
    // vertices of the main mesh
    let Some(VertexAttributeValues::Float32x3(mesh_vertices)) =
        main_mesh.attribute(Mesh::ATTRIBUTE_POSITION)
    else {
        return None;
    };

    let mut vertices = mesh_vertices.clone();

    info!("mesh vertices {:?}", mesh_vertices);
    info!("normal_vec {:?}", normal_vec);
    info!("origin_point {:?}", origin_point);

    // Split in half the main mesh into mesh mappings
    let (top_slice_splited, bottom_slice_splited) =
        cut_mesh(&origin_point, &normal_vec, &mut vertices, main_mesh);

    // create the top mesh
    let mut fragment_top = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );

    //set vertices
    fragment_top.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        MeshMapping::into_vertices(top_slice_splited.get_vertex_buffer()),
    );

    // //set normals
    // fragment_top.insert_attribute(Mesh::ATTRIBUTE_NORMAL);

    // //set uv
    // fragment_top.insert_attribute(Mesh::ATTRIBUTE_UV_0);

    // create the bottom mesh
    let mut fragment_bottom = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );

    //set vertices
    fragment_bottom.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        MeshMapping::into_vertices(bottom_slice_splited.get_vertex_buffer()),
    );

    // //set normals
    // fragment_top.insert_attribute(Mesh::ATTRIBUTE_NORMAL);

    // //set uv
    // fragment_top.insert_attribute(Mesh::ATTRIBUTE_UV_0);

    Some((fragment_top, fragment_bottom, *origin_point))
}

fn is_above_plane(point: Vec3A, plane_normal: &Vec3A, plane_point: &Vec3A) -> bool {
    let vector_to_plane = (point - *plane_point).normalize();
    let distance = -vector_to_plane.dot(*plane_normal);
    return distance < 0.;
}

fn cut_mesh(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    vertices: &mut Vec<[f32; 3]>,
    mesh: &Mesh,
) -> (MeshMapping, MeshMapping) {
    // flaf buffer to track the position of the vertices along the slice plane
    let mut sides = Vec::new();

    // Get the vertices indexes of the current mesh
    let mesh_indices = mesh.indices().unwrap();

    for vertex in &mut *vertices {
        let vertex_above_plane =
            is_above_plane(Vec3A::from_array(vertex.clone()), normal_vec, origin_point);

        sides.push(vertex_above_plane);
    }

    // create the two index buffers for each future meshes
    let (index_buffer_top_slice, index_buffer_bottom_slice) =
        split_triangles(vertices, origin_point, normal_vec, mesh_indices, &sides);

    todo!()
}

trait Indexable {
    fn at(&self, idx: usize) -> usize;
}

impl Indexable for Indices {
    fn at(&self, idx: usize) -> usize {
        match self {
            Indices::U16(vec) => vec[idx] as usize,
            Indices::U32(vec) => vec[idx] as usize,
        }
    }
}

fn split_triangles(
    vertices: &mut Vec<[f32; 3]>,
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    main_mesh_vertices_index_buffer: &Indices,
    side: &Vec<bool>,
) -> (Vec<usize>, Vec<usize>) {
    // Creating the index buffers for the top/bottom slices
    let mut index_buffer_top_slice = Vec::new();
    let mut index_buffer_bottom_slice = Vec::new();

    let mut vertices_added = HashMap::new();

    for (index, vertex) in vertices.iter().enumerate() {
        vertices_added.insert(
            [
                OrderedFloat(vertex[0]),
                OrderedFloat(vertex[1]),
                OrderedFloat(vertex[2]),
            ],
            main_mesh_vertices_index_buffer.,
        )
    }

    for index in (0..main_mesh_vertices_index_buffer.len()).step_by(3) {
        let vertex_indexe_1 = main_mesh_vertices_index_buffer.at(index);
        let vertex_indexe_2 = main_mesh_vertices_index_buffer.at(index + 1);
        let vertex_indexe_3 = main_mesh_vertices_index_buffer.at(index + 2);

        // if triangle is completely above plane
        if side[vertex_indexe_1] && side[vertex_indexe_2] && side[vertex_indexe_3] {
            index_buffer_top_slice.push(vertex_indexe_1);
            index_buffer_top_slice.push(vertex_indexe_2);
            index_buffer_top_slice.push(vertex_indexe_3);
        }
        // if triangle is bellow plane
        else if !side[vertex_indexe_1] && !side[vertex_indexe_2] && !side[vertex_indexe_3] {
            index_buffer_bottom_slice.push(vertex_indexe_1);
            index_buffer_bottom_slice.push(vertex_indexe_2);
            index_buffer_bottom_slice.push(vertex_indexe_3);
        }
        // the triangle is beeing intercept by the slide plane:
        else {
            //two vertices above and one below:
            if side[vertex_indexe_1] && side[vertex_indexe_2] && !side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    true,
                    &mut vertices_added,
                );
            }

            if side[vertex_indexe_1] && !side[vertex_indexe_2] && side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_3,
                    vertex_indexe_2,
                    true,
                    &mut vertices_added,
                );
            }

            if !side[vertex_indexe_1] && side[vertex_indexe_2] && side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    vertex_indexe_1,
                    true,
                    &mut vertices_added,
                );
            }

            // two vertices below and one above:
            if side[vertex_indexe_1] && !side[vertex_indexe_2] && !side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    vertex_indexe_1,
                    false,
                    &mut vertices_added,
                );
            }

            if !side[vertex_indexe_1] && side[vertex_indexe_2] && !side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_3,
                    vertex_indexe_2,
                    false,
                    &mut vertices_added,
                );
            }

            if !side[vertex_indexe_1] && !side[vertex_indexe_2] && side[vertex_indexe_3] {
                split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    &mut index_buffer_top_slice,
                    &mut index_buffer_bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    false,
                    &mut vertices_added,
                );
            }
        }
    }

    info!("index_buffer_top_slice {:?}", index_buffer_top_slice);
    info!("index_buffer_bottom_slice {:?}", index_buffer_bottom_slice);
    info!("vertices len {:?}", vertices);

    (index_buffer_top_slice, index_buffer_bottom_slice)
}

/// if two edges on top:
///  1-----------------2
///   |              |
///     |          |
///   --13|-------| 23---- plane
///         |   |
///           |
///           3
/// and the top quad is divided into 1 2 13 and 2 23 13, and we get the bottom triangle 13 23 3 (indices might be differents)
///
/// if two edges bellow:
///           3
///           |
///         |   |
///   --13|-------| 23---- plane
///     |           |
///   |               |
///  1------------------2
/// and the bottom quad is divided into 3 1 13 and 2 13 23, and we get the top triangle 23 13 3 (indices might be differents)
fn split_small_triangles(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    vertices: &mut Vec<[f32; 3]>,
    index_buffer_top_slice: &mut Vec<usize>,
    index_buffer_bottom_slice: &mut Vec<usize>,
    vertex_indexe_1: usize,
    vertex_indexe_2: usize,
    vertex_indexe_3: usize,
    two_edges_on_top: bool,
    vertices_added: &mut HashMap<[OrderedFloat<f32>; 3], VertexId>,
) {
    // vertices of the current triangle crossed by the slice plane
    let v1 = Vec3A::from_array(*vertices.get(vertex_indexe_1).unwrap());
    let v2 = Vec3A::from_array(*vertices.get(vertex_indexe_2).unwrap());
    let v3 = Vec3A::from_array(*vertices.get(vertex_indexe_3).unwrap());

    // check if the two edges of the triangle are crossed
    let (intersec_13, v13) = find_intersection_line_plane(v1, v1 - v3, *origin_point, *normal_vec);
    let (intersec_23, v23) = find_intersection_line_plane(v2, v2 - v3, *origin_point, *normal_vec);

    // Check if the cut plane do intersect the triangle
    if intersec_13 && intersec_23 {
        // /!\ Interpolate normals and UV coordinates

        // /!\ Add vertices/normals/uv for the intersection points to each mesh

        let ordered_v13: [OrderedFloat<f32>; 3] = [
            OrderedFloat(v13.x),
            OrderedFloat(v13.y),
            OrderedFloat(v13.z),
        ];
        let ordered_v23 = [
            OrderedFloat(v23.x),
            OrderedFloat(v23.y),
            OrderedFloat(v23.z),
        ];

        // create the indices
        let mut index13 = vertices.len() - 2;
        let mut index23 = vertices.len() - 1;

        if !vertices_added.contains_key(&ordered_v13) {
            // add the new vertices in the hash map
            vertices_added.insert(ordered_v13, index13);

            // add the new vertices in the list
            vertices.push(v13.to_array());
        } else {
            index13 = *vertices_added.get(&ordered_v13).unwrap();
        }

        if !vertices_added.contains_key(&ordered_v23) {
            // add the new vertices in the hash map
            vertices_added.insert(ordered_v23, index23);

            // add the new vertices in the list
            vertices.push(v23.to_array());
        } else {
            index23 = *vertices_added.get(&ordered_v23).unwrap();
        }

        if two_edges_on_top {
            // add two triangles on top
            index_buffer_top_slice.push(index23);
            index_buffer_top_slice.push(index13);
            index_buffer_top_slice.push(vertex_indexe_2);

            index_buffer_top_slice.push(index13);
            index_buffer_top_slice.push(vertex_indexe_1);
            index_buffer_top_slice.push(vertex_indexe_2);

            // and one bellow
            index_buffer_bottom_slice.push(vertex_indexe_3);
            index_buffer_bottom_slice.push(index13);
            index_buffer_bottom_slice.push(index23);
        } else {
            // add two triangles bellow
            index_buffer_bottom_slice.push(vertex_indexe_1);
            index_buffer_bottom_slice.push(vertex_indexe_2);
            index_buffer_bottom_slice.push(index13);

            index_buffer_bottom_slice.push(vertex_indexe_2);
            index_buffer_bottom_slice.push(index23);
            index_buffer_bottom_slice.push(index13);

            // and one on top
            index_buffer_top_slice.push(index13);
            index_buffer_top_slice.push(index23);
            index_buffer_top_slice.push(vertex_indexe_3);
        }
    }
}

fn find_intersection_line_plane(
    line_point: Vec3A,
    line_direction: Vec3A,
    origin_point: Vec3A,
    normal_vec: Vec3A,
) -> (bool, Vec3A) {
    let output = Vec3A::ZERO;

    let d = normal_vec.dot(origin_point);
    if normal_vec.dot(line_direction) == f32::NAN {
        (false, output)
    } else {
        let x = (d - normal_vec.dot(line_point)) / normal_vec.dot(line_direction);
        (true, line_point + line_direction.normalize() * x)
    }
}

fn spawn_fragments(
    fragment_mesh: &Mesh,
    fragment_mesh_handle: &Handle<Mesh>,
    mut materials: &ResMut<Assets<StandardMaterial>>,
    commands: &mut Commands,
    pos: Vec3,
) {
    commands.spawn((
        PbrBundle {
            mesh: fragment_mesh_handle.clone(),
            transform: Transform::from_translation(pos),
            ..default()
        },
        RigidBody::Dynamic,
        Collider::from_bevy_mesh(fragment_mesh, &ComputedColliderShape::ConvexHull).unwrap(),
        ActiveCollisionTypes::default(),
        Friction::coefficient(0.7),
        Restitution::coefficient(0.05),
        ColliderMassProperties::Density(2.0),
    ));
}
