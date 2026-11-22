use crate::math::{ComplexField, Real, Vector};
use crate::query::{Ray, RayIntersection};
use crate::shape::FeatureId;

/// Computes the intersection of a ray with a round cylinder centered at the origin, with its
/// axis along Y: the points within `border` of the cylinder of the given half height and
/// radius.
///
/// The round cylinder fills its bounding cylinder (half height `half_height + border`, radius
/// `radius + border`) except near the two rims, where it is a torus. So the ray is cast on the
/// bounding cylinder exactly, and when the point it enters (or leaves) by lies on a flat cap or
/// on the side, that is the answer. Only when it lies in a rim corner is the torus solved: by
/// Newton's method on the distance to the inner cylinder minus the border, which is convex
/// along the ray, so from outside the shape every step stays short of the surface and the
/// first step that turns uphill proves a miss. Works with an unnormalized ray direction.
pub fn ray_toi_and_normal_with_round_cylinder(
    half_height: Real,
    radius: Real,
    border: Real,
    ray: &Ray,
    solid: bool,
) -> Option<RayIntersection> {
    let (o, d) = (ray.origin, ray.dir);
    if d == Vector::ZERO {
        return None;
    }
    let hh = half_height + border;
    let rr = radius + border;

    // The bounding cylinder, as in `ray_toi_and_normal_with_cylinder`.
    let a = d.x * d.x + d.z * d.z;
    let shift = if a == 0.0 {
        0.0
    } else {
        -(o.x * d.x + o.z * d.z) / a
    };
    let near = o + d * shift;
    let b = near.x * d.x + near.z * d.z;
    let c = near.x * near.x + near.z * near.z - rr * rr;
    let (side_in, side_out) = if a == 0.0 {
        if c > 0.0 {
            return None;
        }
        (-Real::MAX, Real::MAX)
    } else {
        let delta = b * b - a * c;
        if delta < 0.0 {
            return None;
        }
        let root = <Real as ComplexField>::sqrt(delta);
        let q = if b >= 0.0 { -(b + root) } else { root - b };
        if q == 0.0 {
            (shift, shift)
        } else {
            let (t0, t1) = (q / a, c / q);
            (shift + t0.min(t1), shift + t0.max(t1))
        }
    };
    let (cap_in, cap_out) = if d.y == 0.0 {
        if o.y.abs() > hh {
            return None;
        }
        (-Real::MAX, Real::MAX)
    } else {
        let t0 = (-hh - o.y) / d.y;
        let t1 = (hh - o.y) / d.y;
        (t0.min(t1), t0.max(t1))
    };
    let enter = side_in.max(cap_in);
    let exit = side_out.min(cap_out);
    if enter > exit || exit < 0.0 {
        return None;
    }
    // A point of the ray, from the point nearest the axis so it keeps the shape's precision.
    let at = |t: Real| near + d * (t - shift);

    // The horizontal unit vector from the axis towards `p`.
    let radial = |p: Vector, fallback: Vector| {
        Vector::new(p.x, 0.0, p.z)
            .try_normalize()
            .or_else(|| Vector::new(fallback.x, 0.0, fallback.z).try_normalize())
            .unwrap_or(Vector::X)
    };

    // Newton's method on f(s) = dist(p0 + dir * s, inner cylinder) - border from a point `p0`
    // outside the shape (or on it), for s in [0, s_max]. f is convex, so the steps never pass
    // the first root, and a step uphill means the ray misses. Returns the root and the outward
    // normal there.
    let rim = |p0: Vector, dir: Vector, s_max: Real| -> Option<(Real, Vector)> {
        let dl = dir.length();
        let tol = (rr + hh) * 1.0e-7;
        let mut s = 0.0;
        let (mut f, mut grad) = rim_eval(p0, half_height, radius, border, dir);
        for _ in 0..64 {
            if f <= tol {
                return Some((s, grad));
            }
            let df = grad.dot(dir);
            if df >= 0.0 {
                return None;
            }
            let step = -f / df;
            s += step;
            if s > s_max {
                return None;
            }
            (f, grad) = rim_eval(p0 + dir * s, half_height, radius, border, dir);
            if step * dl <= tol {
                return Some((s, grad));
            }
        }
        None
    };

    // Is a point of the bounding cylinder's surface on the round cylinder's flat part?
    let flat = |p: Vector, through_cap: bool| {
        if through_cap {
            p.x * p.x + p.z * p.z <= radius * radius
        } else {
            p.y.abs() <= half_height
        }
    };
    let flat_normal = |p: Vector, through_cap: bool| {
        if through_cap {
            Vector::new(0.0, if p.y >= 0.0 { 1.0 } else { -1.0 }, 0.0)
        } else {
            radial(p, d)
        }
    };

    // The exit point, and the outward normal there.
    let leave = || -> Option<(Real, Vector)> {
        let through_cap = cap_out < side_out;
        let p = at(exit);
        if flat(p, through_cap) {
            return Some((exit, flat_normal(p, through_cap)));
        }
        let s_max = exit - enter.max(0.0);
        rim(p, -d, s_max).map(|(s, n)| (exit - s, n))
    };

    if enter >= 0.0 {
        let through_cap = cap_in > side_in;
        let p = at(enter);
        if flat(p, through_cap) {
            return Some(RayIntersection::new(
                enter,
                flat_normal(p, through_cap),
                FeatureId::Unknown,
            ));
        }
        return rim(p, d, exit - enter)
            .map(|(s, n)| RayIntersection::new(enter + s, n, FeatureId::Unknown));
    }

    // The ray starts inside the bounding cylinder.
    let p = at(0.0);
    if rim_eval(p, half_height, radius, border, d).0 > 0.0 {
        // In a rim corner, outside the shape.
        return rim(p, d, exit).map(|(s, n)| RayIntersection::new(s, n, FeatureId::Unknown));
    }
    let (t, n) = leave()?;
    if solid {
        Some(RayIntersection::new(0.0, n, FeatureId::Unknown))
    } else {
        Some(RayIntersection::new(t, -n, FeatureId::Unknown))
    }
}

/// The distance to the inner cylinder minus the border at `p`, and its gradient (the outward
/// normal of the round cylinder there), for a point outside the inner cylinder. Inside it,
/// the value is negative and the gradient meaningless.
#[inline]
fn rim_eval(
    p: Vector,
    half_height: Real,
    radius: Real,
    border: Real,
    fallback: Vector,
) -> (Real, Vector) {
    let rxz = <Real as ComplexField>::sqrt(p.x * p.x + p.z * p.z);
    let dr = rxz - radius;
    let dy = p.y.abs() - half_height;
    if dr <= 0.0 && dy <= 0.0 {
        return (dr.max(dy) - border, Vector::Y);
    }
    let (dr, dy) = (dr.max(0.0), dy.max(0.0));
    let dist = <Real as ComplexField>::sqrt(dr * dr + dy * dy);
    let rad = if rxz > 0.0 {
        Vector::new(p.x / rxz, 0.0, p.z / rxz)
    } else {
        Vector::new(fallback.x, 0.0, fallback.z)
            .try_normalize()
            .unwrap_or(Vector::X)
    };
    let grad = (rad * dr + Vector::new(0.0, p.y.signum() * dy, 0.0)) / dist;
    (dist - border, grad)
}
