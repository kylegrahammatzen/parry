use crate::math::{ComplexField, Real, Vector};
use crate::query::{Ray, RayIntersection};
use crate::shape::FeatureId;

/// Computes the intersection of a ray with a round cuboid centered at the origin: the points
/// within `border_radius` of the box of the given half extents.
///
/// The ray is clipped to the box expanded by the border radius first. Each point of that box
/// lies in the Voronoi region of one feature of the inner box: a face, an edge or a vertex,
/// depending on how many of its coordinates are beyond the half extents. Inside a face region
/// the expanded box and the round cuboid are the same, inside an edge region the round cuboid
/// is a cylinder around the edge and inside a vertex region a ball around the vertex. The ray
/// is walked from region to region (it crosses at most six region boundaries) until it meets
/// the shape, which is exact, unlike the generic GJK ray cast. Works with an unnormalized ray
/// direction.
pub fn ray_toi_and_normal_with_round_cuboid(
    half_extents: Vector,
    border_radius: Real,
    ray: &Ray,
    max_time_of_impact: Real,
    solid: bool,
) -> Option<RayIntersection> {
    if ray.dir == Vector::ZERO {
        // Not a ray: it has no direction.
        return None;
    }
    let h = half_extents.to_array();
    let r = border_radius;
    let o = ray.origin.to_array();
    let d = ray.dir.to_array();
    let e = [h[0] + r, h[1] + r, h[2] + r];

    // Clip the line to the expanded box.
    let mut t_in = -Real::MAX;
    let mut t_out = Real::MAX;
    let (mut ax_in, mut ax_out) = (0, 0);
    for i in 0..3 {
        if d[i] == 0.0 {
            if o[i].abs() > e[i] {
                return None;
            }
        } else {
            // `(h - o) + r` rather than `(h + r) - o`: near the face both differences are
            // exact, where rounding `h + r` first would move the face by half an ulp.
            let a = ((-h[i] - o[i]) - r) / d[i];
            let b = ((h[i] - o[i]) + r) / d[i];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            if lo > t_in {
                t_in = lo;
                ax_in = i;
            }
            if hi < t_out {
                t_out = hi;
                ax_out = i;
            }
        }
    }
    if t_in > t_out || t_out < 0.0 {
        return None;
    }

    // Where the line enters the expanded box, on the face it enters through.
    let mut b_in = [o[0] + d[0] * t_in, o[1] + d[1] * t_in, o[2] + d[2] * t_in];
    b_in[ax_in] = if d[ax_in] > 0.0 { -e[ax_in] } else { e[ax_in] };
    let (s, n) = walk(&h, r, b_in, d, t_out - t_in, ax_in)?;
    let t = t_in + s;
    if t >= 0.0 {
        // The ray starts outside.
        return (t <= max_time_of_impact).then(|| RayIntersection::new(t, n, FeatureId::Unknown));
    }

    // The line enters the shape behind the origin: the origin is inside, or past the exit.
    if solid {
        let q = Vector::from(o) - Vector::from(o).clamp(-half_extents, half_extents);
        return (q.length_squared() <= r * r)
            .then(|| RayIntersection::new(0.0, Vector::ZERO, FeatureId::Unknown));
    }
    let md = [-d[0], -d[1], -d[2]];
    let mut b_out = [
        o[0] + d[0] * t_out,
        o[1] + d[1] * t_out,
        o[2] + d[2] * t_out,
    ];
    b_out[ax_out] = if d[ax_out] > 0.0 {
        e[ax_out]
    } else {
        -e[ax_out]
    };
    let (s, n) = walk(&h, r, b_out, md, t_out - t_in, ax_out)?;
    let t = t_out - s;
    // Hollow: the ray leaves through the surface, whose normal faces the origin.
    (t >= 0.0 && t <= max_time_of_impact).then(|| RayIntersection::new(t, -n, FeatureId::Unknown))
}

/// Walks the line `b + d s` from `s = 0` (a point on face `ax` of the expanded box) to
/// `s_end` (where it leaves the expanded box) and returns the first `s` where it touches the
/// round cuboid, with the outward normal there.
#[inline]
fn walk(
    h: &[Real; 3],
    r: Real,
    b: [Real; 3],
    d: [Real; 3],
    s_end: Real,
    ax: usize,
) -> Option<(Real, Vector)> {
    // Most rays meet the flat part of the face they enter the expanded box through.
    let (j, k) = ((ax + 1) % 3, (ax + 2) % 3);
    if b[j].abs() <= h[j] && b[k].abs() <= h[k] {
        let mut n = Vector::ZERO;
        n[ax] = if b[ax] >= 0.0 { 1.0 } else { -1.0 };
        return Some((0.0, n));
    }
    let mut s0 = 0.0;
    // Each coordinate crosses each of its two region boundaries at most once.
    for _ in 0..7 {
        let mut s1 = s_end;
        for j in 0..3 {
            if d[j] != 0.0 {
                let a = (-h[j] - b[j]) / d[j];
                let c = (h[j] - b[j]) / d[j];
                if a > s0 && a < s1 {
                    s1 = a;
                }
                if c > s0 && c < s1 {
                    s1 = c;
                }
            }
        }
        // The region of this stretch of the line, from its middle.
        let sm = 0.5 * (s0 + s1);
        let m = [b[0] + d[0] * sm, b[1] + d[1] * sm, b[2] + d[2] * sm];
        let out = [m[0].abs() > h[0], m[1].abs() > h[1], m[2].abs() > h[2]];
        let n_out = out[0] as u32 + out[1] as u32 + out[2] as u32;
        if n_out <= 1 {
            // A face region, where the inside of the expanded box is inside the shape.
            let i = if s0 == 0.0 || n_out == 0 {
                ax
            } else if out[0] {
                0
            } else if out[1] {
                1
            } else {
                2
            };
            let p = b[i] + d[i] * s0;
            let mut n = Vector::ZERO;
            n[i] = if p >= 0.0 { 1.0 } else { -1.0 };
            return Some((s0, n));
        }

        // An edge (two coordinates out) or a vertex region (three): the shape there is the
        // points within r of the edge line or of the vertex, so solve |(b + d s - c)_out| = r
        // over the coordinates that are out, from the start of the stretch: its coordinates
        // that are out are within r of the feature, so the terms stay the size of r.
        let mut c = [0.0; 3];
        let (mut qa, mut qb, mut qc) = (0.0, 0.0, -r * r);
        for j in 0..3 {
            if out[j] {
                c[j] = if m[j] >= 0.0 { h[j] } else { -h[j] };
                let p = b[j] + d[j] * s0 - c[j];
                qa += d[j] * d[j];
                qb += p * d[j];
                qc += p * p;
            }
        }
        let hit = if qa == 0.0 {
            (qc <= 0.0).then_some(0.0)
        } else {
            let delta = qb * qb - qa * qc;
            if delta < 0.0 {
                None
            } else {
                let root = <Real as ComplexField>::sqrt(delta);
                let q = if qb >= 0.0 { -(qb + root) } else { root - qb };
                let (r0, r1) = if q == 0.0 {
                    (0.0, 0.0)
                } else {
                    (q / qa, qc / q)
                };
                let (enter, exit) = if r0 < r1 { (r0, r1) } else { (r1, r0) };
                (exit >= 0.0 && enter <= s1 - s0).then(|| enter.max(0.0))
            }
        };
        if let Some(s) = hit {
            let s = s0 + s;
            let mut v = Vector::ZERO;
            let mut f = Vector::ZERO;
            for j in 0..3 {
                if out[j] {
                    v[j] = b[j] + d[j] * s - c[j];
                    f[j] = if m[j] >= 0.0 { 1.0 } else { -1.0 };
                }
            }
            let n = v
                .try_normalize()
                .or_else(|| f.try_normalize())
                .unwrap_or(Vector::X);
            return Some((s, n));
        }
        if s1 >= s_end {
            return None;
        }
        s0 = s1;
    }
    None
}
