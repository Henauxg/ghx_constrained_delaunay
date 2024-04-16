use std::slice::Chunks;

use bevy_rapier3d::{
    parry::query::point::PointCompositeShapeProjWithFeatureBestFirstVisitor, prelude::*,
};

use bevy::{
    math::Vec3A,
    prelude::*,
    render::{
        mesh::{Indices, PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
    },
};
use rand::Rng;

use super::mesh_mapping::MeshMapping;

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

fn get_random_normalized_vec() -> Vec3A {
    let mut rng = rand::thread_rng();
    Vec3A::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()).normalize()
}

pub fn slice(
    commands: &mut Commands,
    meshes_assets: &mut ResMut<Assets<Mesh>>,
    mesh: &Mesh,
    material: &Handle<StandardMaterial>,
) {
    let mesh_center = mesh.compute_aabb().unwrap().center;
    let normal_vec = get_random_normalized_vec();
    internal_slice(
        commands,
        mesh_center,
        normal_vec,
        mesh,
        material,
        meshes_assets,
    )
}

fn internal_slice(
    mut commands: &mut Commands,
    origin_point: Vec3A,
    normal_vec: Vec3A,
    main_mesh: &Mesh,
    material: &Handle<StandardMaterial>,
    meshes_assets: &mut ResMut<Assets<Mesh>>,
) {
    let (top_mesh, bottom_mesh, pos) =
        compute_slice(&origin_point, &normal_vec, main_mesh).unwrap();

    for mesh in vec![top_mesh, bottom_mesh] {
        let mesh_handle = meshes_assets.add(mesh.clone());
        spawn_fragments(&mesh, &mesh_handle, material, &mut commands, pos.into());
    }
}

fn compute_slice(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    main_mesh: &Mesh,
) -> Option<(Mesh, Mesh, Vec3A)> {
    let Some(VertexAttributeValues::Float32x3(values)) =
        main_mesh.attribute(Mesh::ATTRIBUTE_POSITION)
    else {
        return None;
    };

    // Split in half all triangles that are inside the slice plane:
    let (top_slice_splited, bottom_slice_splited) =
        separate_triangles_in_half(&origin_point, &normal_vec, values, main_mesh);

    let mut fragment_top = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );
    fragment_top.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        top_slice_splited.get_vertex_buffer().into(), // a voir pourquoi ???
    );

    let mut fragment_bottom = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );
    fragment_bottom.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        bottom_slice_splited.get_vertex_buffer().into(), // a voir pourquoi ???
    );

    Some((fragment_top, fragment_bottom, *origin_point))
}

fn is_above_plane(point: Vec3A, plane_normal: &Vec3A, plane_point: &Vec3A) -> bool {
    let vector_to_plane = (point - *plane_point).normalize();
    let distance = -vector_to_plane.dot(*plane_normal);
    return distance < 0.;
}

fn separate_triangles_in_half(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    vertices: &Vec<[f32; 3]>,
    mesh: &Mesh,
) -> (MeshMapping, MeshMapping) {
    // let mut top_mesh_vertices = Vec::new(); // contains top vertices for top mesh
    // let mut bottom_mesh_vertices = Vec::new(); // contains bottom vertices for bottom mesh

    // let mut top_mesh_indexes = Vec::new(); // contains top vertices indexes for top mesh
    // let mut bottom_mesh_indexes = Vec::new(); // contains bottom vertices indexes for bottom mesh

    let mut sides = Vec::new();

    // Get the vertices indexes of the current mesh
    let mesh_indices = mesh.indices().unwrap();
    // let mut index_buffer = Vec::new();

    // for index in mesh_indices.iter() {
    //     index_buffer.push(index);
    // }

    // Create the two slices buffers
    let mut top_slice = MeshMapping::new();
    let mut bottom_slice = MeshMapping::new();

    for vertex_id in 0..vertices.len() {
        let vertex = &vertices[vertex_id];
        let vertex_above_plane =
            is_above_plane(Vec3A::from_array(*vertex), normal_vec, origin_point);

        // Put vertices and indexes in two buffers for the two areas (top and bottom)
        match vertex_above_plane {
            true => top_slice.add_mapped_vertex(vertex, vertex_id),
            false => bottom_slice.add_mapped_vertex(vertex, vertex_id),
        };
        sides.push(vertex_above_plane);
    }

    // // Create the two slices buffers
    // let mut top_slice = MeshMapping::new(top_mesh_vertices, top_mesh_indexes);
    // let mut bottom_slice = MeshMapping::new(bottom_mesh_vertices, bottom_mesh_indexes);

    let (mut top_slice_splited, mut bottom_slice_splited) = split_triangles(
        vertices,
        origin_point,
        normal_vec,
        mesh_indices,
        &sides,
        top_slice,
        bottom_slice,
    );

    (top_slice_splited, bottom_slice_splited)
}

// trait NewTrait {
//     /// Returns an iterator over the indices.
//     fn chunks(&self, chunk_size: usize) -> Chunks<'_, u32>;
// }

// impl NewTrait for Indices {
//     /// Returns an iterator over the indices.
//     fn chunks(&self, chunk_size: usize) -> Chunks<'_, u32> {
//         match self {
//             Indices::U16(vec) => vec.chunks(chunk_size),
//             Indices::U32(vec) => vec.chunks(chunk_size),
//         }
//     }
// }

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
    vertices: &Vec<[f32; 3]>,
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    index_buffer: &Indices,
    side: &Vec<bool>,
    mut top_slice: MeshMapping,
    mut bottom_slice: MeshMapping,
) -> (MeshMapping, MeshMapping) {
    // Creating the index buffers for the top/bottom slices
    let mut index_buffer_top_slice = top_slice.get_index_buffer().clone();
    let mut index_buffer_bottom_slice = bottom_slice.get_index_buffer().clone();

    // for index in (0..index_buffer.len()).step_by(3) {
    //     let vertex_indexe_1 = index_buffer[index];
    //     let vertex_indexe_2 = index_buffer[index + 1];
    //     let vertex_indexe_3 = index_buffer[index + 3];

    for index in index_buffer.iter().step_by(3) {
        let vertex_indexe_1 = index_buffer.at(index);
        let vertex_indexe_2 = index_buffer.at(index + 1);
        let vertex_indexe_3 = index_buffer.at(index + 3);

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
            // two vertices above and one below:
            if side[vertex_indexe_1] && side[vertex_indexe_2] && !side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    true,
                );
            }

            if side[vertex_indexe_1] && !side[vertex_indexe_2] && side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_3,
                    vertex_indexe_2,
                    true,
                );
            }

            if !side[vertex_indexe_1] && side[vertex_indexe_2] && side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    vertex_indexe_1,
                    true,
                );
            }

            // two vertices below and one above:
            if side[vertex_indexe_1] && !side[vertex_indexe_2] && !side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    vertex_indexe_1,
                    false,
                );
            }

            if !side[vertex_indexe_1] && side[vertex_indexe_2] && !side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_3,
                    vertex_indexe_2,
                    false,
                );
            }

            if !side[vertex_indexe_1] && !side[vertex_indexe_2] && side[vertex_indexe_3] {
                (top_slice, bottom_slice) = split_small_triangles(
                    origin_point,
                    normal_vec,
                    vertices,
                    top_slice,
                    bottom_slice,
                    vertex_indexe_1,
                    vertex_indexe_2,
                    vertex_indexe_3,
                    false,
                );
            }
        }
    }

    (top_slice, bottom_slice)
}

fn split_small_triangles(
    origin_point: &Vec3A,
    normal_vec: &Vec3A,
    vertices: &Vec<[f32; 3]>,
    mut top_slice_splited: MeshMapping,
    mut bottom_slice_splited: MeshMapping,
    vertex_indexe_1: usize,
    vertex_indexe_2: usize,
    vertex_indexe_3: usize,
    is_bellow_cut_plane: bool,
) -> (MeshMapping, MeshMapping) {
    let v1 = Vec3A::from_array(*vertices.get(vertex_indexe_1).unwrap());
    let v2 = Vec3A::from_array(*vertices.get(vertex_indexe_2).unwrap());
    let v3 = Vec3A::from_array(*vertices.get(vertex_indexe_3).unwrap());

    let (intersec_13, _v13) = find_intersection_line_plane(v1, v1 - v3, *origin_point, *normal_vec);
    let (intersec_23, _v23) = find_intersection_line_plane(v2, v2 - v3, *origin_point, *normal_vec);

    // Check if the cut plane do intersect the triangle
    if intersec_13 && intersec_23 {
        // /!\ Interpolate normals and UV coordinates => demander car WTF

        // /!\ Add vertices/normals/uv for the intersection points to each mesh => ALED

        // Indices for the intersection vertices => A voir ensemble car bizzre ??
        let index13_a = top_slice_splited.get_vertex_buffer().len() - 2;
        let index23_a = top_slice_splited.get_vertex_buffer().len() - 1;
        let index13_b = bottom_slice_splited.get_vertex_buffer().len() - 2;
        let index23_b = bottom_slice_splited.get_vertex_buffer().len() - 1;

        if is_bellow_cut_plane {
            // add two triangles on top
            top_slice_splited.add_triangle(index23_a, index13_a, vertex_indexe_2);
            top_slice_splited.add_triangle(index13_a, vertex_indexe_1, vertex_indexe_2);

            // and one bellow
            bottom_slice_splited.add_triangle(vertex_indexe_3, index13_b, index23_b);

            // /!\ When looking at the cut-face, the edges should wind counter-clockwise => vient du délire précédent, demander
        } else {
            // add two triangles bellow
            bottom_slice_splited.add_triangle(vertex_indexe_1, vertex_indexe_2, index13_b);
            bottom_slice_splited.add_triangle(vertex_indexe_2, index23_b, index13_b);

            // and one on top
            top_slice_splited.add_triangle(index13_a, index23_a, vertex_indexe_3);

            // /!\ When looking at the cut-face, the edges should wind counter-clockwise => vient du délire précédent, demander
        }
    }

    (top_slice_splited, bottom_slice_splited)
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
    material: &Handle<StandardMaterial>,
    commands: &mut Commands,
    pos: Vec3,
) {
    commands.spawn((
        PbrBundle {
            mesh: fragment_mesh_handle.clone(),
            material: material.clone(),
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

pub fn substract_two_vertices(vertexe_1: [f32; 3], vertexe_2: [f32; 3]) -> [f32; 3] {
    let x = vertexe_1[0] - vertexe_2[0];
    let y = vertexe_1[1] - vertexe_2[1];
    let z = vertexe_1[2] - vertexe_2[2];
    [x, y, z]
}

pub fn normalize_vertexe(vertexe: [f32; 3]) -> [f32; 3] {
    let vec = Vec3A::from_array(vertexe);
    vec.normalize();
    vec.to_array()
}

pub fn cross_from_vec3a(vertexe_1: [f32; 3], vertexe_2: Vec3A) -> [f32; 3] {
    let vec = Vec3A::from_array(vertexe_1);
    let dot_product = vec.cross(vertexe_2);
    dot_product.to_array()
}

pub fn dot_from_array(vertexe_1: [f32; 3], vertexe_2: [f32; 3]) -> f32 {
    let vec1 = Vec3A::from_array(vertexe_1);
    let vec2 = Vec3A::from_array(vertexe_2);
    let dot_product = vec1.dot(vec2);
    dot_product
}
