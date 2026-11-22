use crate::math::{ComplexField, Real, Vector};
use crate::query::gjk::VoronoiSimplex;
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Cone, FeatureId};

impl RayCast for Cone {
    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        if self.half_height <= 0.0 || self.half_height.is_nan() {
            // A flat cone (a disk): no slant to solve for.
            return crate::query::details::local_ray_intersection_with_support_map_with_params(
                self,
                &mut VoronoiSimplex::new(),
                ray,
                max_time_of_impact,
                solid,
            );
        }
        ray_toi_and_normal_with_cone(self.half_height, self.radius, ray, solid)
            .filter(|inter| inter.time_of_impact <= max_time_of_impact)
    }
}

/// Computes the intersection of a ray with a cone centered at the origin, with its axis along Y,
/// its apex at `y = half_height` and its base of the given radius at `y = -half_height`.
///
/// The solid cone is the part of the infinite cone below the apex that is above the base plane,
/// so the ray is inside it for the times it is both inside the infinite cone and above the base
/// plane: the hit is the later of the two entry times. This is exact, unlike the generic GJK ray
/// cast, which could stall and miss rays going deep through small cones. Works with an
/// unnormalized ray direction.
pub fn ray_toi_and_normal_with_cone(
    half_height: Real,
    radius: Real,
    ray: &Ray,
    solid: bool,
) -> Option<RayIntersection> {
    let d = ray.dir;
    let dd = d.length_squared();

    if dd == 0.0 {
        // Not a ray: it has no direction (as for the cylinder, and for the GJK ray cast
        // cones used before).
        return None;
    }

    // Everything is solved from a point of the ray inside the cone's bounding box (`shift`
    // along it), so the terms stay the size of the cone: from far away, or from far above or
    // below a flat cone, the terms of the quadratic would nearly cancel.
    let extents = Vector::new(radius, half_height, radius);
    let (mut box_in, mut box_out) = (-Real::MAX, Real::MAX);
    for i in 0..3 {
        if d[i] != 0.0 {
            let t0 = (-extents[i] - ray.origin[i]) / d[i];
            let t1 = (extents[i] - ray.origin[i]) / d[i];
            box_in = box_in.max(t0.min(t1));
            box_out = box_out.min(t0.max(t1));
        } else if ray.origin[i].abs() > extents[i] {
            return None;
        }
    }
    if box_in > box_out || box_out < 0.0 {
        return None;
    }
    // Both are finite: at least one direction component is not zero.
    let mut shift = box_in * 0.5 + box_out * 0.5;
    let mut o = ray.origin + d * shift;

    // The radius of the infinite cone shrinks by `k` per unit of height up to the apex.
    let k = radius / (2.0 * half_height);
    let k2 = k * k;

    // When the ray is inside the double cone x² + z² <= k² w², with w the height below the
    // apex: A t² + 2 B t + C <= 0. As for the ball, B and C are taken at the vertex of this
    // parabola (kept inside the box), where C is its extremum: a ray grazing the cone, or
    // passing near the apex, then gets its roots without cancellation.
    let a = d.x * d.x + d.z * d.z - k2 * d.y * d.y;
    if a != 0.0 {
        let b = o.x * d.x + o.z * d.z + k2 * (half_height - o.y) * d.y;
        let vertex = (-b / a).max(box_in - shift).min(box_out - shift);
        o += d * vertex;
        shift += vertex;
    }
    let w0 = half_height - o.y;

    // When the ray is below the apex (w = w0 - d.y t >= 0).
    let (below_in, below_out) = if d.y > 0.0 {
        (-Real::MAX, w0 / d.y)
    } else if d.y < 0.0 {
        (w0 / d.y, Real::MAX)
    } else if w0 >= 0.0 {
        (-Real::MAX, Real::MAX)
    } else {
        return None;
    };

    let b = o.x * d.x + o.z * d.z + k2 * w0 * d.y;
    let c = o.x * o.x + o.z * o.z - k2 * w0 * w0;

    // The two roots, in increasing order, in the numerically stable form.
    let roots = |delta: Real| {
        let root = <Real as ComplexField>::sqrt(delta);
        let q = if b >= 0.0 { -(b + root) } else { root - b };
        if q == 0.0 {
            // b = 0 and delta = 0, so c = 0: a double root at 0.
            (0.0, 0.0)
        } else {
            let (t0, t1) = (q / a, c / q);
            (t0.min(t1), t0.max(t1))
        }
    };

    let (dbl_in, dbl_out) = if a > 0.0 {
        // The ray crosses the cone's side less steeply than its slant: inside between the roots.
        let delta = b * b - a * c;
        if delta < 0.0 {
            return None;
        }
        roots(delta)
    } else if a < 0.0 {
        // Steeper than the slant: the line crosses both nappes, and is inside the lower one
        // (below the apex) before the first root or after the second. It always meets the
        // double cone, so rounding must not make the discriminant negative.
        let (t0, t1) = roots((b * b - a * c).max(0.0));
        if d.y > 0.0 {
            (-Real::MAX, t0)
        } else {
            (t1, Real::MAX)
        }
    } else if b != 0.0 {
        // Parallel to a generatrix of the cone: one crossing.
        let t = -c / (2.0 * b);
        if b > 0.0 {
            (-Real::MAX, t)
        } else {
            (t, Real::MAX)
        }
    } else if c <= 0.0 {
        (-Real::MAX, Real::MAX)
    } else {
        return None;
    };

    let side_in = dbl_in.max(below_in);
    let side_out = dbl_out.min(below_out);

    // When the ray is above the base plane.
    let (base_in, base_out) = if d.y > 0.0 {
        ((-half_height - o.y) / d.y, Real::MAX)
    } else if d.y < 0.0 {
        (-Real::MAX, (-half_height - o.y) / d.y)
    } else if o.y >= -half_height {
        (-Real::MAX, Real::MAX)
    } else {
        return None;
    };

    let enter = side_in.max(base_in);
    let exit = side_out.min(base_out);
    if enter > exit || shift + exit < 0.0 {
        return None;
    }

    // The outward normal where the ray crosses the surface at (relative) time `t`.
    let normal_at = |t: Real, through_base: bool| {
        if through_base {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            let p = o + d * t;
            let rxz = <Real as ComplexField>::sqrt(p.x * p.x + p.z * p.z);
            // At the apex any direction of its normal cone is right, and +Y is one.
            Vector::new(p.x, k * rxz, p.z)
                .try_normalize()
                .unwrap_or(Vector::new(0.0, 1.0, 0.0))
        }
    };

    if shift + enter >= 0.0 {
        // The ray starts outside.
        let normal = normal_at(enter, base_in > side_in);
        Some(RayIntersection::new(
            shift + enter,
            normal,
            FeatureId::Unknown,
        ))
    } else if solid {
        // The ray starts inside a solid cone: it hits at once. The normal is the one the GJK
        // ray cast used to return here, facing the ray.
        let normal = (-d).try_normalize().unwrap_or(Vector::ZERO);
        Some(RayIntersection::new(0.0, normal, FeatureId::Unknown))
    } else {
        // The ray starts inside a hollow cone: it hits the surface on its way out, which faces
        // the ray's origin.
        let normal = -normal_at(exit, base_out < side_out);
        Some(RayIntersection::new(
            shift + exit,
            normal,
            FeatureId::Unknown,
        ))
    }
}
