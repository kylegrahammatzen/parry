// Cones, round cuboids and round cylinders used the iterative GJK ray cast, which missed some
// rays going deep through small or distant shapes and hollow casts from inside, and stopped a
// few millimeters short of or past the surface. These rays have their exact answer.

use parry3d::math::{Real, Vector};
use parry3d::query::{Ray, RayCast};
use parry3d::shape::{Cone, Cuboid, Cylinder, RoundShape};

#[track_caller]
fn hit(shape: &impl RayCast, o: Vector, d: Vector, solid: bool) -> (Real, Vector) {
    let i = shape
        .cast_local_ray_and_get_normal(&Ray::new(o, d), Real::MAX, solid)
        .expect("expected a hit");
    (i.time_of_impact, i.normal)
}

#[test]
fn cone_exact_cases() {
    let cone = Cone::new(0.5, 0.3);
    // Down the axis onto the apex, up the axis onto the base.
    let (t, n) = hit(&cone, Vector::new(0.0, 10.0, 0.0), -Vector::Y, true);
    assert!((t - 9.5).abs() < 1.0e-6 && (n - Vector::Y).length() < 1.0e-6);
    let (t, n) = hit(&cone, Vector::new(0.0, -10.0, 0.0), Vector::Y, true);
    assert!((t - 9.5).abs() < 1.0e-6 && (n + Vector::Y).length() < 1.0e-6);
    // The side, with an unnormalized direction.
    let (t, n) = hit(
        &cone,
        Vector::new(5.0, 0.0, 0.0),
        Vector::new(-1.0e5, 0.0, 0.0),
        true,
    );
    let slant = Vector::new(1.0, 0.3, 0.0).normalize();
    assert!((t - 4.85e-5).abs() < 1.0e-9 && (n - slant).length() < 1.0e-5);
    // Inside, solid: at once, hollow: the exit, in the ray's own time unit.
    assert_eq!(
        hit(&cone, Vector::ZERO, Vector::new(0.0, 3.0, 0.0), true).0,
        0.0
    );
    for (len, expected) in [(1.0, 0.15), (2.0, 0.075), (0.5, 0.3), (1.0e-5, 15000.0)] {
        let (t, n) = hit(&cone, Vector::ZERO, Vector::new(len, 0.0, 0.0), false);
        assert!((t - expected).abs() <= 1.0e-5 * expected, "{len}: {t}");
        assert!((n + slant).length() < 1.0e-5);
    }
    // Rays that pass just above the apex or point away miss.
    assert!(cone
        .cast_local_ray(
            &Ray::new(Vector::new(-5.0, 0.5001, 0.0), Vector::X),
            Real::MAX,
            true
        )
        .is_none());
    assert!(cone
        .cast_local_ray(
            &Ray::new(Vector::new(0.0, 10.0, 0.0), Vector::Y),
            Real::MAX,
            true
        )
        .is_none());
}

#[test]
fn cone_distant_and_extreme() {
    // A game-sized cone hit through its base from 42 m.
    let cone = Cone::new(0.38772818, 0.2115757);
    let o = Vector::new(-0.98230547, -18.726023, -37.748966);
    let d = Vector::new(0.02429677, 0.43873632, 0.8982873);
    let (t, n) = hit(&cone, o, d, true);
    assert!((t - 41.797986).abs() < 1.0e-4 && (n + Vector::Y).length() < 1.0e-6);
    // A needle and a pancake.
    let (t, _) = hit(
        &Cone::new(300.0, 0.01),
        Vector::new(10.0, -100.0, 0.0),
        -Vector::X,
        true,
    );
    assert!((t - (10.0 - 0.01 * 400.0 / 600.0)).abs() < 1.0e-4);
    let (t, _) = hit(
        &Cone::new(0.01, 300.0),
        Vector::new(100.0, 50.0, 0.0),
        -Vector::Y,
        true,
    );
    assert!((t - (50.0 - 0.01 / 3.0)).abs() < 1.0e-4);
}

#[test]
fn round_cuboid() {
    let shape = RoundShape {
        inner_shape: Cuboid::new(Vector::new(1.0882847, 0.6951978, 0.95831984)),
        border_radius: 0.17048365,
    };
    // From 3.6 m onto a rounded edge (the time of impact was 3 mm short).
    let o = Vector::new(-3.5890112, 1.0028619, -0.59255254);
    let d = Vector::new(0.95109946, -0.04488127, 0.30560666);
    let (t, _) = hit(&shape, o, d, true);
    assert!((t - 3.0565175).abs() < 1.0e-4, "{t}");

    // A hollow cast from inside a small round cuboid (it was a miss).
    let shape = RoundShape {
        inner_shape: Cuboid::new(Vector::new(0.02544053, 0.02435734, 0.026508704)),
        border_radius: 0.010296405,
    };
    let o = Vector::new(0.024400692, 0.010736465, 0.003901688);
    let d = Vector::new(-0.01640589, -0.8264535, 0.56276596);
    let (t, _) = hit(&shape, o, d, false);
    assert!((t - 0.05197953).abs() < 1.0e-5, "{t}");
}

#[test]
fn round_cylinder() {
    let shape = RoundShape {
        inner_shape: Cylinder::new(0.07200076, 0.2429958),
        border_radius: 0.039718576,
    };
    // From 25 m onto a rounded rim (the time of impact was 1.4 mm past it).
    let o = Vector::new(11.938063, 20.482256, 6.9343233);
    let d = Vector::new(-0.4881502, -0.8293023, -0.27196878);
    let (t, _) = hit(&shape, o, d, true);
    assert!((t - 24.567307).abs() < 2.0e-4, "{t}");

    // A hollow cast from inside a small round cylinder (it was a miss).
    let shape = RoundShape {
        inner_shape: Cylinder::new(0.042399637, 0.02439157),
        border_radius: 0.012544327,
    };
    let o = Vector::new(-0.020912038, 0.0021028484, -0.010534);
    let d = Vector::new(0.778064, 0.47697848, 0.40878835);
    let (t, _) = hit(&shape, o, d, false);
    assert!((t - 0.06865924).abs() < 1.0e-5, "{t}");
}
