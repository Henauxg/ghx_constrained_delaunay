use bevy::render::color::Color;

mod triangle_debug;
pub use triangle_debug::*;

mod example_setup;
pub use example_setup::*;

pub mod lines;

pub mod camera;

pub const DEBUG_LABEL_FONT_SIZE: f32 = 60.0;

pub const COLORS: &'static [Color] = &[
    Color::GREEN,
    Color::BLUE,
    Color::BLACK,
    Color::RED,
    Color::YELLOW,
    Color::MAROON,
    Color::PURPLE,
    Color::SALMON,
    Color::ORANGE,
    Color::CYAN,
    Color::NAVY,
    Color::OLIVE,
    Color::PINK,
    Color::ALICE_BLUE,
    Color::CRIMSON,
    Color::TURQUOISE,
    Color::YELLOW_GREEN,
    Color::TEAL,
    Color::MIDNIGHT_BLUE,
    Color::VIOLET,
];
