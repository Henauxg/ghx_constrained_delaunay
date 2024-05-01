pub struct MeshMapping {
    vertex_buffer: Vec<[f32; 3]>,
    index_buffer: Vec<u32>,
    index_map: Vec<u32>,
}

impl MeshMapping {
    // pub fn new(vertex_buffer: Vec<[f32; 3]>, index_buffer: Vec<usize>) -> MeshMapping {
    //     MeshMapping {
    //         vertex_buffer: vertex_buffer,
    //         index_buffer: index_buffer,
    //     }
    // }

    pub fn new() -> MeshMapping {
        MeshMapping {
            vertex_buffer: Vec::new(),
            index_buffer: Vec::new(),
            index_map: Vec::new(),
        }
    }

    pub fn get_vertex_buffer(&self) -> &Vec<[f32; 3]> {
        &self.vertex_buffer
    }

    pub fn get_index_buffer(&self) -> &Vec<u32> {
        &self.index_buffer
    }

    pub fn add_triangle(&mut self, i1: u32, i2: u32, i3: u32) {
        let _ = &self.index_buffer.push(i1);
        let _ = &self.index_buffer.push(i2);
        let _ = &self.index_buffer.push(i3);
    }

    pub fn add_mapped_vertex(&self, vertex: &[f32; 3], vertex_id: u32) {
        // self.vertex_buffer.push(*vertex);
        // self.index_map[vertex_id as usize] = (self.vertex_buffer.len() - 1) as u32;
        // // TODO Mapping
    }
}
