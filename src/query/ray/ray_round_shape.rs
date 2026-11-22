use crate::math::Real;
use crate::query::gjk::VoronoiSimplex;
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{RoundShape, SupportMap};

impl<S: SupportMap> RayCast for RoundShape<S> {
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_time_of_impact: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        if let Some(hit) = self.inner_shape.local_round_ray_cast(
            self.border_radius,
            ray,
            max_time_of_impact,
            solid,
        ) {
            return hit;
        }
        crate::query::details::local_ray_intersection_with_support_map_with_params(
            self,
            &mut VoronoiSimplex::new(),
            ray,
            max_time_of_impact,
            solid,
        )
    }
}
