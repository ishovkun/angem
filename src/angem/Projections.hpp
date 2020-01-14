#pragma once
#include "Polygon.hpp"
#include "Collisions.hpp"

namespace angem
{

/* Project polygon 1 onto polygon 2 and compute intersecting polygon */
template<typename Scalar>
Polygon<Scalar> project(const Polygon<Scalar> & poly1,
                        const Polygon<Scalar> & poly2)
{
  const std::vector<Point<3,Scalar>> projected_frac_vertices =
      poly2.plane().project_points(poly1.get_points());
  const Polygon<Scalar> proj_poly(projected_frac_vertices);
  std::vector<Point<3,Scalar>> intersection;
  if (!collision(proj_poly, poly2, intersection))
    throw std::runtime_error("projection does not exist");
  return Polygon<Scalar>(intersection);
}

}  // end namespace angem
