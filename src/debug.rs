use crate::types::{Float, Neighbor, TriangleData, TriangleId, VertexId};

pub struct DebugSnapshot {
    pub step: usize,
    pub triangulation_phase: TriangulationPhase,
    pub changed_ids: Vec<TriangleId>,
    pub triangles: Vec<TriangleData>,
}
impl DebugSnapshot {
    pub(crate) fn new(
        step: usize,
        triangulation_phase: TriangulationPhase,
        triangles: Vec<TriangleData>,
        changed_ids: Vec<TriangleId>,
    ) -> Self {
        Self {
            step,
            triangulation_phase,
            triangles,
            changed_ids,
        }
    }
}

pub struct DebugContext {
    // TODO Could return container triangle vertices as world coordinates
    pub scale_factor: Float,
    pub x_min: Float,
    pub y_min: Float,

    pub snapshots: Vec<DebugSnapshot>,
    pub curent_step: usize,
}

impl DebugContext {
    pub(crate) fn new(scale_factor: Float, x_min: Float, y_min: Float) -> Self {
        Self {
            snapshots: Vec::new(),
            scale_factor,
            x_min,
            y_min,
            curent_step: 0,
        }
    }

    pub(crate) fn push_snapshot(
        &mut self,
        phase: TriangulationPhase,
        triangles: &Vec<TriangleData>,
        triangle_ids: &[TriangleId],
        opt_neighbor_ids: &[Neighbor],
    ) {
        let mut step_changes = Vec::with_capacity(triangle_ids.len() + opt_neighbor_ids.len());
        step_changes.extend(triangle_ids);
        // TODO May consider neighbors differently
        step_changes.extend(
            opt_neighbor_ids
                .iter()
                .filter(|o| o.is_some())
                .map(|o| o.unwrap())
                .collect::<Vec<TriangleId>>(),
        );

        self.snapshots.push(DebugSnapshot::new(
            self.curent_step,
            phase,
            triangles.clone(),
            step_changes,
        ));
    }

    pub(crate) fn set_step(&mut self, step: usize) {
        self.curent_step = step;
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum TriangulationPhase {
    ContainerVerticesInsertion,
    SplitTriangle(VertexId),
    SwapQuadDiagonal,
    RemoveWrapping,
    //
    AfterConstraints,
}
