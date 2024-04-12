use bevy_inspector_egui::quick::WorldInspectorPlugin;
use bevy_rapier3d::prelude::*;
use examples::plugin::ExamplesPlugin;
use rand::Rng;
use std::f32::consts::*;

use bevy::{
    ecs::query,
    math::Vec3A,
    pbr::CascadeShadowConfigBuilder,
    prelude::*,
    render::{
        mesh::{shape::Quad, MeshVertexAttribute, PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
        render_resource::encase::vector,
    },
};

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

fn get_random_vec() -> Vec3 {
    let mut rng = rand::thread_rng();
    Vec3::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>())
}

fn compute_slice(
    mut commands: Commands,
    origin_point: Vec3A,
    normal_vec: Vec3A,
    main_mesh: &Mesh,
    material: &Handle<StandardMaterial>,
    mut meshes_assets: ResMut<Assets<Mesh>>,
) {
    let (top_mesh, bottom_mesh, pos) = slice(&origin_point, &normal_vec, main_mesh).unwrap();

    for mesh in vec![top_mesh, bottom_mesh] {
        let mesh_handle = meshes_assets.add(mesh.clone());
        create_fragments(&mesh, &mesh_handle, material, &mut commands, pos.into());
    }
}

fn slice(
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
    let (top_mesh_vertices, bottom_mesh_vertices) =
        separate_triangles_in_half(&origin_point, &normal_vec, values);

    let mut fragment_top = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );
    fragment_top.insert_attribute(Mesh::ATTRIBUTE_POSITION, top_mesh_vertices);

    let mut fragment_bottom = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD,
    );
    fragment_bottom.insert_attribute(Mesh::ATTRIBUTE_POSITION, bottom_mesh_vertices);

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
) -> (Vec<[f32; 3]>, Vec<[f32; 3]>) {
    let mut top_mesh_vertices = Vec::new(); // contains top vertices for top mesh
    let mut bottom_mesh_vertices = Vec::new(); // contains bottom vertices for bottom mesh

    for vertex_id in 0..vertices.len() {
        let vertex_1 = vertices[vertex_id];
        let v1_side = is_above_plane(Vec3A::from_array(vertex_1), normal_vec, origin_point);

        // TODO Proper triangle splitting
        // TODO Vert indices
        match v1_side {
            true => top_mesh_vertices.push(vertex_1),
            false => bottom_mesh_vertices.push(vertex_1),
        };
    }
    (top_mesh_vertices, bottom_mesh_vertices)
}

fn create_fragments(
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

// let v2_side = is_above_plane(Vec3A::from_array(vertex_2), normal_vec, origin_point);
// let v3_side = is_above_plane(Vec3A::from_array(vertex_3), normal_vec, origin_point);

// // if triangle is completely above plane
// if v1_side && v2_side && v3_side {
//     top_mesh_vertices.push(vertex_1);
//     top_mesh_vertices.push(vertex_2);
//     top_mesh_vertices.push(vertex_3);
// }

// // if triangle is not above plane
// if !v1_side && !v2_side && !v3_side {
//     bottom_mesh_vertices.push(vertex_1);
//     bottom_mesh_vertices.push(vertex_2);
//     bottom_mesh_vertices.push(vertex_3);
// }
// // the triangle is beeing intercept by the slide plane:
// else {
//     // two vertices above and one below:
//     if !v1_side && v2_side && v3_side {
//         bottom_mesh_vertices.push(vertex_1);
//         top_mesh_vertices.push(vertex_2);
//         top_mesh_vertices.push(vertex_3);
//     }
//     if v1_side && !v2_side && v3_side {
//         top_mesh_vertices.push(vertex_1);
//         bottom_mesh_vertices.push(vertex_2);
//         top_mesh_vertices.push(vertex_3);
//     }
//     if v1_side && v2_side && !v3_side {
//         top_mesh_vertices.push(vertex_1);
//         top_mesh_vertices.push(vertex_2);
//         bottom_mesh_vertices.push(vertex_3);
//     }
//     // two vertices below and one above:
//     if !v1_side && !v2_side && v3_side {
//         bottom_mesh_vertices.push(vertex_1);
//         bottom_mesh_vertices.push(vertex_2);
//         top_mesh_vertices.push(vertex_3);
//     }
//     if !v1_side && v2_side && !v3_side {
//         bottom_mesh_vertices.push(vertex_1);
//         top_mesh_vertices.push(vertex_2);
//         bottom_mesh_vertices.push(vertex_3);
//     }
//     if v1_side && !v2_side && !v3_side {
//         top_mesh_vertices.push(vertex_1);
//         bottom_mesh_vertices.push(vertex_2);
//         bottom_mesh_vertices.push(vertex_3);
//     }
// }
