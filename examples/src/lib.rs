use bevy::render::color::Color;

pub mod ball;
pub mod debug_utils;
pub mod fps;
pub mod plugin;

pub const DEFAULT_EXAMPLES_FONT_SIZE: f32 = 17.;

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
];
