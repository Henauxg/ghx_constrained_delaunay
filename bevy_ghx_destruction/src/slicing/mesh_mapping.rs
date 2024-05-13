#[derive(Debug)]
pub struct MeshMapping {
    vertex_buffer: Vec<[f32; 3]>,
    index_buffer: Vec<u32>,
    index_map: Vec<usize>,
}

impl MeshMapping {
    pub fn new(nbr_vertices: usize) -> MeshMapping {
        MeshMapping {
            vertex_buffer: Vec::new(),
            index_buffer: Vec::new(),
            index_map: Vec::new(),
        }
    }

    pub fn get_vertex_buffer(&self) -> &Vec<[f32; 3]> {
        &self.vertex_buffer
    }

    pub fn into_vertices(vertex_buffer: &Vec<[f32; 3]>) -> Vec<[f32; 3]> {
        let mut vertices: Vec<[f32; 3]> = vec![[0., 0., 0.]; vertex_buffer.len()];
        for (index, vertex) in vertex_buffer.iter().enumerate() {
            let x = vertex[0];
            let y = vertex[1];
            let z = vertex[2];
            vertices[index] = [x, y, z];
        }

        vertices
    }

    pub fn get_index_buffer(&self) -> &Vec<u32> {
        &self.index_buffer
    }

    pub fn get_mapped_index_buffer(&self) -> &Vec<usize> {
        &self.index_map
    }

    pub fn add_triangle(&mut self, i1: u32, i2: u32, i3: u32) {
        let _ = &self.index_buffer.push(i1);
        let _ = &self.index_buffer.push(i2);
        let _ = &self.index_buffer.push(i3);
    }

    pub fn add_mapped_vertex(&mut self, vertex: &[f32; 3], vertex_id: usize) {
        let _ = &self.vertex_buffer.push(*vertex);
        self.index_map.push(vertex_id);
        // TODO Mapping
    }
}
