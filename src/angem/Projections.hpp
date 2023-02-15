#pragma once
#include "Polygon.hpp"
#include "Collisions.hpp"

namespace angem::projections
{


/* Project polygon 1 onto polygon 2 and points within intersection.
 * These points might form a polygon, if the projection area is not zero.
 * Therefore, the function does not return a polygon but merely a point list. */
template<typename Scalar>
std::vector<Point<3,Scalar>> project(const Polygon<Scalar> & poly1,
                                     const Polygon<Scalar> & poly2,
                                     const double tol = 1e-6)
{
  if ( std::fabs(poly1.normal().dot( poly2.normal() )) < tol)
    throw std::invalid_argument("Polygons are orthogonal");

  const std::vector<Point<3,Scalar>> projected_frac_vertices =
      poly2.plane().project_points(poly1.get_points());

  const Polygon<Scalar> proj_poly(projected_frac_vertices);
  std::vector<Point<3,Scalar>> intersection;

  if (!collision(proj_poly, poly2, intersection, tol))
    throw std::runtime_error("projection does not exist");
  return intersection;
}


}  // end namespace angem::projections
