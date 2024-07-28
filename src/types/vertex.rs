use std::ops::Sub;

use glam::{DVec2, DVec3, Vec2, Vec3};

use super::Float;

/// Defines a simple 2d point
pub trait Vertex2d: Clone + Copy + Sized {
    fn x(self) -> Float;
    fn y(self) -> Float;
}

impl Vertex2d for Vec2 {
    #[inline(always)]
    #[must_use]
    fn x(self) -> Float {
        self.x as Float
    }

    #[inline(always)]
    #[must_use]
    fn y(self) -> Float {
        self.y as Float
    }
}
impl Vertex2d for DVec2 {
    #[inline(always)]
    #[must_use]
    fn x(self) -> Float {
        self.x as Float
    }

    #[inline(always)]
    #[must_use]
    fn y(self) -> Float {
        self.y as Float
    }
}

/// Defines a 3d point with a few basic operations
pub trait Vertex3d:
    Clone + Copy + Sized + Sub<Output = Self> + std::ops::Neg<Output = Self>
{
    fn normalize(self) -> Self;
    fn cross(self, rhs: Self) -> Self;
    fn dot(self, rhs: Self) -> Float;
}
impl Vertex3d for Vec3 {
    #[inline(always)]
    #[must_use]
    fn normalize(self) -> Self {
        self.normalize()
    }

    #[inline(always)]
    #[must_use]
    fn cross(self, rhs: Self) -> Self {
        self.cross(rhs)
    }

    #[inline(always)]
    #[must_use]
    fn dot(self, rhs: Self) -> Float {
        self.dot(rhs) as Float
    }
}
impl Vertex3d for DVec3 {
    #[inline(always)]
    #[must_use]
    fn normalize(self) -> Self {
        self.normalize()
    }

    #[inline(always)]
    #[must_use]
    fn cross(self, rhs: Self) -> Self {
        self.cross(rhs)
    }

    #[inline(always)]
    #[must_use]
    fn dot(self, rhs: Self) -> Float {
        self.dot(rhs) as Float
    }
}
