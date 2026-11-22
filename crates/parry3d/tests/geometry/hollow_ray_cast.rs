// A hollow ray cast (`solid = false`) from inside a cylinder, cone or convex polyhedron used to
// mix meters and the ray's time unit when the direction wasn't unit length, so the exit time
// came out wrong and could even be negative.

use parry3d::math::{Pose, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::{Cone, ConvexPolyhedron, Cylinder};

fn check_exit(name: &str, shape: &dyn RayCast, normal: Vector) {
    // The exit is 1 m away along +x, whatever the direction's length.
    for dir_length in [0.5, 1.0, 2.0, 10.0] {
        let ray = Ray::new(Vector::ZERO, Vector::new(dir_length, 0.0, 0.0));
        let hit = shape
            .cast_ray_and_get_normal(&Pose::IDENTITY, &ray, 100.0, false)
            .unwrap_or_else(|| panic!("{name}: no exit with a direction of length {dir_length}"));
        assert_relative_eq!(hit.time_of_impact, 1.0 / dir_length, epsilon = 1.0e-3);
        assert_relative_eq!(hit.normal, normal, epsilon = 1.0e-3);
    }
}

#[test]
fn hollow_ray_cast_from_inside_with_a_direction_that_is_not_unit() {
    let x = Vector::new(-1.0, 0.0, 0.0);
    check_exit("cylinder", &Cylinder::new(1.0, 1.0), x);
    // The cone's side leans 45 degrees, so its inward normal points down as well.
    let slant = Vector::new(-1.0, -1.0, 0.0).normalize();
    check_exit("cone", &Cone::new(1.0, 2.0), slant);
    let corners: Vec<_> = (0..8)
        .map(|i| {
            Vector::new(
                if i & 1 == 0 { -1.0 } else { 1.0 },
                if i & 2 == 0 { -1.0 } else { 1.0 },
                if i & 4 == 0 { -1.0 } else { 1.0 },
            )
        })
        .collect();
    let cube = ConvexPolyhedron::from_convex_hull(&corners).unwrap();
    check_exit("cube hull", &cube, x);
}
