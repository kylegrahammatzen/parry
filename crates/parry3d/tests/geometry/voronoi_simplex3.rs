// Reducing a tetrahedron to its face abd, acd or bcd used to overwrite the dropped vertex, so
// `prev_point` returned the new vertex in place of one of the previous triangle's. GJK builds its
// witness points from that previous triangle when it stalls, and they could end up far apart.

use parry3d::math::Vector;
use parry3d::query::gjk::{CsoPoint, VoronoiSimplex};

#[test]
fn tetrahedron_reduced_to_a_face_keeps_the_previous_triangle() {
    // Three vertices on the plane z = 1 around the origin's projection, and one far vertex.
    let near = [
        Vector::new(-1.0, -1.0, 1.0),
        Vector::new(1.0, -1.0, 1.0),
        Vector::new(0.0, 1.0, 1.0),
    ];
    let far = Vector::new(0.0, 0.0, 5.0);

    // The far vertex goes in slot c, b then a, so the origin projects on face abd, acd then bcd.
    for far_slot in [2, 1, 0] {
        let mut triangle = [near[0], near[1], near[2]];
        triangle[far_slot] = far;
        let new_vertex = near[far_slot];

        let mut simplex = VoronoiSimplex::new();
        simplex.reset(CsoPoint::single_point(triangle[0]));
        assert!(simplex.add_point(CsoPoint::single_point(triangle[1])));
        assert!(simplex.add_point(CsoPoint::single_point(triangle[2])));
        assert!(simplex.add_point(CsoPoint::single_point(new_vertex)));

        let proj = simplex.project_origin_and_reduce();
        assert_relative_eq!(proj, Vector::new(0.0, 0.0, 1.0), epsilon = 1.0e-6);
        assert_eq!(simplex.dimension(), 2);

        assert_eq!(simplex.prev_dimension(), 2);
        for (i, expected) in triangle.iter().enumerate() {
            assert_eq!(
                simplex.prev_point(i).point,
                *expected,
                "far vertex in slot {far_slot}: previous point {i} changed"
            );
        }
    }
}
