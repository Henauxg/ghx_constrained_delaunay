use bevy::math::Vec3A;
use rand::Rng;

pub fn get_random_normalized_vec() -> Vec3A {
    let mut rng = rand::thread_rng();
    Vec3A::new(rng.gen::<f32>(), rng.gen::<f32>(), rng.gen::<f32>()).normalize()
}
