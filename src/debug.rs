use hashbrown::HashSet;

use crate::types::{Float, Neighbor, TriangleId, Triangles, VertexId};

/// Defines debug recording of the Triangulation
#[derive(Debug, Clone)]
pub enum PhaseRecord {
    /// Records all the phases
    All,
    /// Records the steps during the specified phase
    In(TriangulationPhase),
    /// Records the steps during the specified phases
    InAny(HashSet<TriangulationPhase>),
}

/// Defines debug recording of the Triangulation
#[derive(Debug, Clone)]
pub enum StepsRecord {
    /// Records all the steps
    All,
    /// Records all the steps after the specified one (inclusive)
    From(usize),
    /// Records all the steps until the specified one (inclusive)
    Until(usize),
    /// Records all the steps between the specified ones (inclusive)
    Between(usize, usize),
}

#[derive(Debug, Clone)]
pub struct DebugConfiguration {
    pub phase_record: PhaseRecord,
    pub steps_record: StepsRecord,
    /// [None] means that it won't force an early exit
    pub force_end_at_step: Option<usize>,
}
impl Default for DebugConfiguration {
    fn default() -> Self {
        Self {
            phase_record: PhaseRecord::All,
            steps_record: StepsRecord::All,
            force_end_at_step: Default::default(),
        }
    }
}

pub struct DebugSnapshot {
    pub step: usize,
    pub triangulation_phase: TriangulationPhase,
    pub changed_ids: Vec<TriangleId>,
    pub triangles: Triangles,
}
impl DebugSnapshot {
    pub(crate) fn new(
        step: usize,
        triangulation_phase: TriangulationPhase,
        triangles: Triangles,
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
    pub config: DebugConfiguration,

    // TODO Could return container triangle vertices as world coordinates
    pub scale_factor: Float,
    pub x_min: Float,
    pub y_min: Float,

    pub snapshots: Vec<DebugSnapshot>,
    pub current_step: usize,
}

impl DebugContext {
    pub(crate) fn new(
        config: DebugConfiguration,
        scale_factor: Float,
        x_min: Float,
        y_min: Float,
    ) -> Self {
        Self {
            config,
            scale_factor,
            x_min,
            y_min,
            snapshots: Vec::new(),
            current_step: 0,
        }
    }

    pub(crate) fn push_snapshot(
        &mut self,
        phase: TriangulationPhase,
        triangles: &Triangles,
        triangle_ids: &[TriangleId],
        opt_neighbor_ids: &[Neighbor],
    ) {
        let record = match &self.config.phase_record {
            PhaseRecord::All => true,
            PhaseRecord::InAny(phases) => phases.contains(&phase),
            PhaseRecord::In(rec_phase) => phase == *rec_phase,
        };
        if !record {
            return;
        }
        if match self.config.steps_record {
            StepsRecord::All => true,
            StepsRecord::From(from) => self.current_step >= from,
            StepsRecord::Until(to) => self.current_step <= to,
            StepsRecord::Between(from, to) => self.current_step >= from && self.current_step <= to,
        } {
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
                self.current_step,
                phase,
                triangles.clone(),
                step_changes,
            ));
        }
    }

    /// Returns true if the algorithm should stop
    pub(crate) fn set_step(&mut self, step: usize) -> bool {
        self.current_step = step;
        match self.config.force_end_at_step {
            Some(end_step) => step >= end_step,
            None => false,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum TriangulationPhase {
    ContainerVerticesInsertion,
    SplitTriangle(VertexId),
    SwapQuadDiagonal,
    RemoveWrapping,
    //
    BeforeConstraints,
    AfterConstraints,
}
