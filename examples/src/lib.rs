mod debugger;
use bevy::color::{palettes::css::*, Color};
pub use debugger::*;

mod plugin;
pub use plugin::*;

pub mod lines;

pub mod camera;

pub const COLORS: &'static [Color] = &[
    Color::Srgba(GREEN),
    Color::Srgba(BLUE),
    Color::Srgba(BLACK),
    Color::Srgba(RED),
    Color::Srgba(YELLOW),
    Color::Srgba(MAROON),
    Color::Srgba(PURPLE),
    Color::Srgba(SALMON),
    Color::Srgba(ORANGE),
    Color::Srgba(CADET_BLUE),
    Color::Srgba(NAVY),
    Color::Srgba(OLIVE),
    Color::Srgba(PINK),
    Color::Srgba(ALICE_BLUE),
    Color::Srgba(CRIMSON),
    Color::Srgba(TURQUOISE),
    Color::Srgba(YELLOW_GREEN),
    Color::Srgba(TEAL),
    Color::Srgba(MIDNIGHT_BLUE),
    Color::Srgba(VIOLET),
];
