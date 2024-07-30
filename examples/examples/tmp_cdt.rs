use bevy::{
    app::{App, Startup},
    ecs::system::Commands,
    log::{error, info},
    math::{DVec2, Vec3},
    prelude::{EventWriter, IntoSystemConfigs, Query, Res, ResMut},
    DefaultPlugins,
};

use examples::{
    camera::PanOrbitCamera, setup_camera, ExamplesPlugin, LabelMode, TriangleDebugCursorUpdate,
    TriangleDebugPlugin, TrianglesDebugData, TrianglesDebugViewConfig, TrianglesDrawMode,
    VertexLabelMode,
};
use ghx_constrained_delaunay::{
    constrained_triangulation::ConstrainedTriangulationConfiguration,
    constrained_triangulation_from_2d_vertices,
    debug::DebugConfiguration,
    types::{Edge, Vertex},
    utils::check_delaunay_optimal,
};

// Modify to change the display scale
const DISPLAY_SCALE: f32 = 500000.;

fn main() {
    App::new()
        .add_plugins((
            DefaultPlugins,
            ExamplesPlugin,
            TriangleDebugPlugin::default(),
        ))
        .add_systems(
            Startup,
            (
                setup,
                display_first_snapshot,
                set_camera_position.after(setup_camera),
            )
                .chain(),
        )
        .run();
}

fn setup(mut commands: Commands) {
    let vertices = vec![
        DVec2::new(31.462459564208984, 4.182472229003906),
        DVec2::new(29.765398025512695, 4.182472229003906),
        DVec2::new(22.388437271118164, 3.416529893875122),
        DVec2::new(20.982324600219727, 4.377071857452393),
        DVec2::new(32.80760192871094, 6.361194133758545),
        DVec2::new(31.462459564208984, 4.182472229003906),
        DVec2::new(33.678627014160156, 7.771997451782227),
        DVec2::new(26.823562622070313, 3.1497762203216553),
        DVec2::new(26.25342559814453, 3.076498508453369),
        DVec2::new(27.03416633605957, 3.1768453121185303),
        DVec2::new(29.045753479003906, 3.6033806800842285),
        DVec2::new(29.634510040283203, 3.7282204627990723),
        DVec2::new(29.803268432617188, 3.7640037536621094),
        DVec2::new(25.360563278198242, 3.8880200386047363),
        DVec2::new(24.0087947845459, 3.812800645828247),
        DVec2::new(25.7857608795166, 3.911680221557617),
        DVec2::new(22.485071182250977, 3.4445977210998535),
        DVec2::new(22.388437271118164, 3.416529893875122),
        DVec2::new(23.02713966369629, 3.280472755432129),
        DVec2::new(29.834842681884766, 3.881032943725586),
        DVec2::new(29.803268432617188, 3.7640037536621094),
        DVec2::new(29.77996253967285, 4.119258880615234),
        DVec2::new(29.769569396972656, 4.164374828338623),
        DVec2::new(29.765398025512695, 4.182472229003906),
        DVec2::new(26.198028564453125, 3.171741008758545),
        DVec2::new(25.84834861755371, 3.7729573249816895),
        DVec2::new(26.25342559814453, 3.076498508453369),
        DVec2::new(25.7857608795166, 3.911680221557617),
        DVec2::new(34.85775375366211, 8.384631156921387),
        DVec2::new(33.678627014160156, 7.771997451782227),
        DVec2::new(23.83258628845215, 3.67209529876709),
        DVec2::new(23.02713966369629, 3.280472755432129),
        DVec2::new(24.0087947845459, 3.812800645828247),
        DVec2::new(30.783798217773438, -10.525899887084961),
        DVec2::new(28.80564308166504, -11.123200416564941),
        DVec2::new(34.610294342041016, -9.428776741027832),
        DVec2::new(32.81171798706055, -9.942837715148926),
        DVec2::new(32.69080352783203, -9.977603912353516),
        DVec2::new(30.783798217773438, -10.525899887084961),
        DVec2::new(42.71263122558594, -1.885225772857666),
        DVec2::new(39.20940399169922, -4.862123966217041),
        DVec2::new(49.4381103515625, 3.879741668701172),
        DVec2::new(44.691070556640625, 10.19693660736084),
        DVec2::new(49.77289962768555, 5.225648403167725),
        DVec2::new(49.89583969116211, 4.2704925537109375),
        DVec2::new(49.4381103515625, 3.879741668701172),
        DVec2::new(50.3620719909668, 4.654744625091553),
        DVec2::new(50.045562744140625, 4.961440563201904),
        DVec2::new(50.3620719909668, 4.654744625091553),
        DVec2::new(49.93797302246094, 5.065690040588379),
        DVec2::new(49.77289962768555, 5.225648403167725),
        DVec2::new(10.620545387268066, -1.8670755624771118),
        DVec2::new(10.711962699890137, -0.5022733211517334),
        DVec2::new(37.098388671875, -6.952149868011475),
        DVec2::new(34.610294342041016, -9.428776741027832),
        DVec2::new(37.496925354003906, -6.555449485778809),
        DVec2::new(39.20940399169922, -4.862123966217041),
        DVec2::new(30.925596237182617, 7.6599249839782715),
        DVec2::new(33.92063522338867, 8.211918830871582),
        DVec2::new(28.791208267211914, 7.266552925109863),
        DVec2::new(33.96668243408203, 8.220404624938965),
        DVec2::new(33.983856201171875, 8.223569869995117),
        DVec2::new(34.36891174316406, 8.294536590576172),
        DVec2::new(34.85775375366211, 8.384631156921387),
        DVec2::new(27.158231735229492, 6.9655914306640625),
        DVec2::new(26.948719024658203, 6.926977157592773),
        DVec2::new(25.96977424621582, 6.746554851531982),
        DVec2::new(43.9307975769043, 10.056818008422852),
        DVec2::new(44.691070556640625, 10.19693660736084),
        DVec2::new(43.71799087524414, 10.017596244812012),
        DVec2::new(37.700401306152344, 8.908539772033691),
        DVec2::new(35.5172119140625, 8.506172180175781),
        DVec2::new(35.47850799560547, 8.499039649963379),
        DVec2::new(35.44377899169922, 8.492637634277344),
        DVec2::new(35.297508239746094, 8.465679168701172),
        DVec2::new(44.841468811035156, 10.224654197692871),
        DVec2::new(44.87495422363281, 10.230826377868652),
        DVec2::new(44.6927604675293, 10.197247505187988),
        DVec2::new(25.908000946044922, 6.717207431793213),
        DVec2::new(25.96977424621582, 6.746554851531982),
        DVec2::new(25.891643524169922, 6.709437370300293),
        DVec2::new(24.027339935302734, 5.823724746704102),
        DVec2::new(23.07302474975586, 5.3703389167785645),
        DVec2::new(22.55023193359375, 5.121965408325195),
        DVec2::new(22.38591766357422, 5.0439019203186035),
        DVec2::new(21.875999450683594, 4.801645755767822),
        DVec2::new(21.103647232055664, 4.43471097946167),
        DVec2::new(20.982324600219727, 4.377071857452393),
        DVec2::new(25.96978187561035, 6.746557235717773),
        DVec2::new(25.969783782958984, 6.74655818939209),
        DVec2::new(25.96978187561035, 6.746557235717773),
        DVec2::new(25.96978187561035, 6.746557235717773),
        DVec2::new(25.969778060913086, 6.746557235717773),
        DVec2::new(25.969772338867188, 6.746555328369141),
        DVec2::new(9.277804374694824, -1.1836272478103638),
        DVec2::new(10.711962699890137, -0.5022733211517334),
        DVec2::new(10.75268840789795, -0.48292315006256104),
        DVec2::new(17.675247192382813, 2.8059115409851074),
        DVec2::new(13.090802192687988, -3.12442684173584),
        DVec2::new(10.682821273803711, -1.898775339126587),
        DVec2::new(13.798712730407715, -3.4847490787506104),
        DVec2::new(15.23499870300293, -4.215813159942627),
        DVec2::new(18.231229782104492, -5.740880012512207),
        DVec2::new(18.980133056640625, -6.122066497802734),
        DVec2::new(20.083602905273438, -6.6837263107299805),
        DVec2::new(23.733049392700195, -8.541275978088379),
        DVec2::new(25.90675926208496, -9.647684097290039),
        DVec2::new(28.80564308166504, -11.123200416564941),
        DVec2::new(10.620545387268066, -1.8670755624771118),
        DVec2::new(29.011558532714844, -11.228011131286621),
        DVec2::new(29.387332916259766, -11.419279098510742),
        DVec2::new(29.419221878051758, -11.435511589050293),
        DVec2::new(30.117847442626953, -11.791106224060059),
        DVec2::new(9.277801513671875, -1.183626413345337),
        DVec2::new(9.277800559997559, -1.18362557888031),
        DVec2::new(9.277804374694824, -1.1836270093917847),
        DVec2::new(9.277803421020508, -1.1836270093917847),
        DVec2::new(9.277804374694824, -1.1836270093917847),
        DVec2::new(9.277804374694824, -1.1836272478103638),
    ];

    let constrained_edges = vec![
        Edge { from: 0, to: 1 },
        Edge { from: 2, to: 3 },
        Edge { from: 4, to: 5 },
        Edge { from: 6, to: 4 },
        Edge { from: 7, to: 8 },
        Edge { from: 9, to: 7 },
        Edge { from: 10, to: 9 },
        Edge { from: 11, to: 10 },
        Edge { from: 12, to: 11 },
        Edge { from: 13, to: 14 },
        Edge { from: 15, to: 13 },
        Edge { from: 16, to: 17 },
        Edge { from: 18, to: 16 },
        Edge { from: 19, to: 20 },
        Edge { from: 21, to: 19 },
        Edge { from: 22, to: 21 },
        Edge { from: 23, to: 22 },
        Edge { from: 24, to: 25 },
        Edge { from: 26, to: 24 },
        Edge { from: 25, to: 27 },
        Edge { from: 28, to: 29 },
        Edge { from: 30, to: 31 },
        Edge { from: 32, to: 30 },
        Edge { from: 33, to: 34 },
        Edge { from: 35, to: 36 },
        Edge { from: 37, to: 38 },
        Edge { from: 36, to: 37 },
        Edge { from: 39, to: 40 },
        Edge { from: 41, to: 39 },
        Edge { from: 42, to: 43 },
        Edge { from: 44, to: 45 },
        Edge { from: 46, to: 44 },
        Edge { from: 47, to: 48 },
        Edge { from: 49, to: 47 },
        Edge { from: 50, to: 49 },
        Edge { from: 51, to: 52 },
        Edge { from: 53, to: 54 },
        Edge { from: 55, to: 53 },
        Edge { from: 56, to: 55 },
        Edge { from: 57, to: 58 },
        Edge { from: 59, to: 57 },
        Edge { from: 58, to: 60 },
        Edge { from: 60, to: 61 },
        Edge { from: 61, to: 62 },
        Edge { from: 62, to: 63 },
        Edge { from: 64, to: 59 },
        Edge { from: 65, to: 64 },
        Edge { from: 66, to: 65 },
        Edge { from: 67, to: 68 },
        Edge { from: 69, to: 67 },
        Edge { from: 70, to: 69 },
        Edge { from: 71, to: 70 },
        Edge { from: 72, to: 71 },
        Edge { from: 73, to: 72 },
        Edge { from: 74, to: 73 },
        Edge { from: 63, to: 74 },
        Edge { from: 75, to: 76 },
        Edge { from: 77, to: 75 },
        Edge { from: 68, to: 77 },
        Edge { from: 78, to: 79 },
        Edge { from: 80, to: 78 },
        Edge { from: 81, to: 80 },
        Edge { from: 82, to: 81 },
        Edge { from: 83, to: 82 },
        Edge { from: 84, to: 83 },
        Edge { from: 85, to: 84 },
        Edge { from: 86, to: 85 },
        Edge { from: 87, to: 86 },
        Edge { from: 88, to: 89 },
        Edge { from: 90, to: 88 },
        Edge { from: 91, to: 90 },
        Edge { from: 92, to: 91 },
        Edge { from: 93, to: 92 },
        Edge { from: 79, to: 93 },
        Edge { from: 94, to: 95 },
        Edge { from: 96, to: 97 },
        Edge { from: 95, to: 96 },
        Edge { from: 97, to: 87 },
        Edge { from: 98, to: 99 },
        Edge { from: 100, to: 98 },
        Edge { from: 101, to: 100 },
        Edge { from: 102, to: 101 },
        Edge { from: 103, to: 102 },
        Edge { from: 104, to: 103 },
        Edge { from: 105, to: 104 },
        Edge { from: 106, to: 105 },
        Edge { from: 107, to: 106 },
        Edge { from: 99, to: 108 },
        Edge { from: 109, to: 107 },
        Edge { from: 110, to: 109 },
        Edge { from: 111, to: 110 },
        Edge { from: 112, to: 111 },
        Edge { from: 113, to: 114 },
        Edge { from: 115, to: 113 },
        Edge { from: 116, to: 115 },
        Edge { from: 117, to: 116 },
        Edge { from: 118, to: 117 },
        Edge { from: 108, to: 118 },
    ];

    let triangulation_result = constrained_triangulation_from_2d_vertices(
        &vertices,
        &constrained_edges,
        ConstrainedTriangulationConfiguration {
            debug_config: DebugConfiguration {
                force_end_at_step: Some(10),
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let debug_context = match triangulation_result {
        Ok(triangulation) => {
            let delaunay_quality = check_delaunay_optimal(
                triangulation.triangles,
                &vertices.iter().map(|v| Vertex::new(v[0], v[1])).collect(),
                false,
            );
            info!("CDT quality info: {:?}", delaunay_quality);
            triangulation.debug_context
        }
        Err(err) => {
            error!("Failed triangulation: {:?}", err.msg);
            err.debug_context
        }
    };

    // Error ? Found vertex place for vertex 116, OnTriangleEdge(4, 0) should be OnTriangleEdge(4, 2)

    let displayed_vertices = vertices
        .iter()
        .map(|v| DISPLAY_SCALE * Vec3::new(v.x as f32, v.y as f32, 0.))
        .collect();
    commands.insert_resource(TrianglesDebugData::new_with_constrained_edges(
        displayed_vertices,
        &constrained_edges,
        debug_context,
    ));
    commands.insert_resource(TrianglesDebugViewConfig::new(
        LabelMode::All,
        VertexLabelMode::GlobalIndex,
        TrianglesDrawMode::AllAsGizmos,
        true,
    ));
}

pub fn set_camera_position(
    debug_data: Res<TrianglesDebugData>,
    mut camera: Query<&mut PanOrbitCamera>,
) {
    let Ok(mut cam) = camera.get_single_mut() else {
        return;
    };
    let Some(vertex) = debug_data.vertices.get(118) else {
        return;
    };
    cam.focus = *vertex;
    cam.needs_transform_refresh = true;
}

pub fn display_first_snapshot(
    mut debug_data: ResMut<TrianglesDebugData>,
    mut debug_data_updates_events: EventWriter<TriangleDebugCursorUpdate>,
) {
    debug_data.set_cursor_to_first_snapshot();
    debug_data_updates_events.send(TriangleDebugCursorUpdate);
}
