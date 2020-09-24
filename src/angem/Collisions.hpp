#pragma once

#include "Point.hpp"
#include "LineSegment.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"
#include "Polyhedron.hpp"
#include "Exceptions.hpp"
#include "angem/utils.hpp"
// #include "PolyGroup.hpp"
#include "CollisionGJK.hpp"
#include "Tensor2.hpp"
#include "utils.hpp"
#include <set>


/* This module contains various algorithms for simple shape intersections.
 * As opposed to CollisionGJK module, here we actually compute the section
 * data. These algorithms can be combined and utilized to get intersection data
 * of complex shape e.g. a mesh with a fracture.
 *
 * All functions in this module are boolean and return true in the case of
 * collision. The resulting section data is written into the last input argument.
 * If the last aregument is a vector, the data is appended.
 */

namespace angem
{

// get intersection of two planes
// returns true if the planes are not parallel
// the result of the intersection is saved into Line intersection.
template <typename Scalar>
bool collision(const Plane<Scalar> & pl1,
               const Plane<Scalar> & pl2,
               Line<3,Scalar>      & intersection)
{
  /* The method selects a third plane P3 with
   * an implicit equation n3 · P = 0
   where n3 = n1 x n2 and d3 = 0 (meaning it passes through the origin).
   This always works since: (1) L is perpendicular to P3 and thus
   intersects it, and (2) the vectors n1, n2, and n3 are linearly independent.
   Thus the planes P1, P2 and P3 intersect in a unique point P0 which must be on L.
   Using the formula for the intersection of 3 planes,
   where d3 = 0 for P3, we get:
        (d2 n1 - d1 n2) × n3
   p3 = --------------------
           (n1 × n2) · n3
   Ref:
   http://geomalgorithms.com/a05-_intersect-1.html
  */

  // direction of intersection line n 3 = n1 x n2
  Point<3,Scalar> n3 = pl1.normal().cross(pl2.normal());
  if (n3.norm() < 1e-16)
    return false;
  n3.normalize();

  intersection.direction = n3;

  const auto & n1 = pl1.normal();
  const auto & d1 = pl1.d;

  const auto & n2 = pl2.normal();
  const auto & d2 = pl2.d;

  Point<3,Scalar> numerator = (d2*n1 - d1*n2).cross(n3);
  Scalar denumenator = (n1.cross(n2)).dot(n3);
  intersection.point = numerator / denumenator;
  return true;
}


// collision of a polygon with a plane
// can be 1 points, two points, or zero points
template <typename Scalar>
bool collision(const Polygon<Scalar>        & poly,
               const Plane<Scalar>          & plane,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  // call collision of all edges
  bool result = false;
  const auto & pts = poly.get_points();
  for (size_t i=0; i<pts.size(); ++i)
  {
    bool loc_collision = false;
    if (i < pts.size() - 1)
      loc_collision = collision(pts[i], pts[i+1], plane, intersection, tol);
    else
      loc_collision = collision(pts[i], pts[0], plane, intersection, tol);
    if (loc_collision)
      result = true;
  }

  intersection = remove_duplicates_slow(intersection, tol);
  return result;
}

/*  Intersetion of two lines in 3D */
template <typename Scalar>
bool collision(const Line<3,Scalar>        & l1,
               const Line<3,Scalar>        & l2,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  /*
  **    Line1                         Line2
  **    -----                         -----
  ** x = a11 * t1 + b11            x = a21 * t2 + b21
  ** y = a12 * t1 + b21            y = a22 * t2 + b22
  ** z = a13 * t1 + b31            z = a23 * t2 + b23
  ** Assuming that these lines actually have an intersection, we need to solve the following
  ** system:
  ** a₁₁ t₁ - a₂₁ t₂ = b₂₁ - b₁₁
  ** a₁₂ t₁ - a₂₂ t₂ = b₂₂ - b₁₂
  ** Since there are actually three equations, and some lines can be aligned with axes,
  ** the determinant
  **     | a₁₁   -a₂₁ |
  ** det | a₁₂   -a₂₂ | must be non-zero.
  ** We must handle such cases carefully.
   */
  using pair = std::pair<size_t, size_t>;
  for (const pair ij : std::vector<pair>{pair(0,1),pair(0,2),pair(1,2)})
  {
    const size_t i = ij.first;
    const size_t j = ij.second;
    const Tensor2<2,Scalar> mat = {l1.direction()[i], -l2.direction()[i],
                                   l1.direction()[j], -l2.direction()[j],};
    if (std::fabs(det(mat)) > tol)
    {
      const Point<2,Scalar> rhs = { l2.point()[i] - l1.point()[i] ,
                                    l2.point()[j] - l1.point()[j] ,};
      const Point<2,Scalar> solution = invert(mat) * rhs;
      intersection.push_back(l1.point() + l1.direction() * solution(0));
      return true;
    }
  }
  return false;
}

// collision of two polygons
template <typename Scalar>
bool collision(const Polygon<Scalar>        & poly1,
               const Polygon<Scalar>        & poly2,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  if (poly1.plane().normal().parallel(poly2.plane().normal(), tol))
  {
    // 1. find vertices of each poly inside another
    // 2. find intersection of edges if any points inside
    PointSet<3,Scalar> pset(tol * 1.5);
    // 1.
    const auto & pts1 = poly1.get_points();
    const auto & pts2 = poly2.get_points();
    bool all_inside1 = true, all_inside2 = true;
    for (const auto & p : pts1)
      if (poly2.point_inside(p, tol))
        pset.insert(p);
      else
        all_inside2 = false;

    for (const auto & p : pts2)
      if (poly1.point_inside(p, tol))
        pset.insert(p);
      else
        all_inside1 = false;

    if (all_inside1)
    {
      pset.points.clear();
      pset.points = pts2;
    }
    if (all_inside2)
    {
      pset.points.clear();
      pset.points = pts1;
    }
    if (all_inside1 && all_inside2)
    {
      for (const auto &p : pts1)
        intersection.push_back(p);
      return true;
    }

    // 2.
    if ( !pset.empty() and !all_inside1 and !all_inside2 )
      for (const auto & edge1 : poly1.get_edges())
      {
        std::vector<Point<3,Scalar>> v_points;
        Plane<Scalar> side = poly1.get_side(edge1);
        for (const auto & edge2 : poly2.get_edges())
        {
          collision(pts2[edge2.first], pts2[edge2.second], side, v_points, tol);
          for (const auto & p : v_points)
            if (poly1.point_inside(p))
              pset.insert(p);
        }
      }

    for (const auto & p: pset.points)
      intersection.push_back(p);

    if (!pset.empty())
      return true;
    else
      return false;
  }
  else // two polygons in non-parallel planes
  {
    std::vector<Point<3,Scalar>> v_section;
    angem::collision(poly1, poly2.plane(), v_section, tol);
    bool result = false;
    for (const auto & p : v_section)
      if (poly2.point_inside(p, tol) and poly1.point_inside(p, tol))
      {
        intersection.push_back(p);
        result = true;
      }
    return result;
  }

  return true;
}


// intersection of a segment with plane
// intersection is appended to!
template <typename Scalar>
bool collision(const Point<3,Scalar>        & l0,
               const Point<3,Scalar>        & l1,
               const Plane<Scalar>          & plane,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  // Plane : (p - p0) · n = 0
  // line p = d l + l0
  // segment : l0, l1
  // intersection: d = (p0 - l0) · n / (l · n)
  // call collision of all edges
  const Scalar d1 = plane.signed_distance(l0);
  const Scalar d2 = plane.signed_distance(l1);

  // both points are on the plane
  if (fabs(d1) + fabs(d2) < tol)
  {
    intersection.push_back(l0);
    intersection.push_back(l1);
    return true;
  }

  if (d1*d2 > 0)  // both points on one side of plane
    return false;

  // compute intersection point
  const Point<3,Scalar> l = l1 - l0;
  const Scalar d = (plane.origin() - l0).dot(plane.normal()) /
                    l.dot(plane.normal());
  intersection.push_back(l0 + d * l);
  return true;
}


// marks polygons above fracture as 1
// polygons below fracture as 0
template <typename Scalar>
void split(const Polygon<Scalar> & poly,
           const Plane<Scalar>   & plane,
           PolyGroup<Scalar>     & result,
           const int               marker_below = 0,
           const int               marker_above = 1,
           const double            tolerance = 1e-6)
{
  std::vector<Point<3,Scalar>> section;
  collision(poly, plane, section, tolerance);

  if (section.empty())  // everythin either above or below
  {
    std::vector<std::size_t> indices;
    bool above = false;
    for (const auto & p : poly.get_points())
    {
      const std::size_t ind = result.vertices.insert(p);
      indices.push_back(ind);

      if (plane.above(p))  // technically need to check only one
        above = true;
    }
    result.polygons.push_back(indices);

    // assign markers
    if (above)
      result.markers.push_back(marker_above);
    else
      result.markers.push_back(marker_below);

    return;
  }

  // if split happened figure outpu which parts are above/below
  std::set<size_t> set_above, set_below;
  for (const auto p : poly.get_points())
  {
    const std::size_t ind = result.vertices.insert(p);
    const Scalar signed_dist = plane.signed_distance(p);
    if (signed_dist > -tolerance)
      set_above.insert(ind);
    if (signed_dist < +tolerance)
      set_below.insert(ind);
  }

  for (Point<3,Scalar> & p : section)
  {
    const std::size_t ind = result.vertices.insert(p);
    set_above.insert(ind);
    set_below.insert(ind);
  }

  std::vector<std::size_t> above(set_above.begin(), set_above.end()),
                           below(set_below.begin(), set_below.end());

  if (above.size() > 2)
  {
    result.polygons.push_back(std::move(above));
    result.markers.push_back(marker_above);
  }
  if (below.size() > 2)
  {
    result.polygons.push_back(std::move(below));
    result.markers.push_back(marker_below);
  }
}

/* splits a polyhedron by plane
 *  marks polygons above fracture as marker_above
 * polygons below fracture as marker_below
 *  marks fracture polygon as marker_split */
template <typename Scalar>
bool split(const Polyhedron<Scalar> & polyhedron,
           const Plane<Scalar>      & plane,
           PolyGroup<Scalar>        & result,
           std::vector<size_t>      & parent_face_indices,
           const int                  marker_below = 0,
           const int                  marker_above = 1,
           const int                  marker_split = 2,
           const double               tolerance = 1e-6)
{
  // first split the faces of the polyhedron
  std::vector<Polygon<Scalar>> face_polygons = polyhedron.get_face_polygons();
  for (std::size_t i = 0; i < face_polygons.size(); ++i)
  {
    const Polygon<Scalar> & face = face_polygons[i];
    const size_t old_n_split_poly = result.polygons.size();
    split<Scalar>(face, plane, result, marker_below, marker_above, tolerance);
    const size_t new_n_split_poly = result.polygons.size();
    if (new_n_split_poly - old_n_split_poly == 2)  // split successful
    {
      parent_face_indices.push_back(i);
      parent_face_indices.push_back(i);
    }
    else // if (new_n_split_poly - old_n_split_poly == 1) // nothing split
      parent_face_indices.push_back(i);
  }

  // sort face vertices so they form polygons
  for (auto & face : result.polygons)
    Polygon<Scalar>::reorder_indices(result.vertices.points, face);

  // add a polygon that represents the intersection
  std::vector<size_t> section_poly_verts;
  // vertices of the polyhedron
  const std::vector<Point<3,double>> & poly_vertices = polyhedron.get_points();
  // number of points in intersecting polygon is equal to the total number of points
  // in the split minus the number of vertices in the polyhedron
  section_poly_verts.reserve(result.vertices.size() - poly_vertices.size());
  // element size
  const double h = poly_vertices[1].distance(poly_vertices[0]);
  // figure out which points lie on the plane
  for (size_t i = 0; i < result.vertices.size(); i++)
    if (std::fabs(plane.signed_distance(result.vertices[i])) < tolerance * h)
      section_poly_verts.push_back(i);
  // postprocessing: create a polygon from the points
  if (!section_poly_verts.empty())
  {
    if (section_poly_verts.size() < 3) // split nothing
    {
      section_poly_verts.clear();
      assert( result.vertices.size() - poly_vertices.size() < 2 );
      return false;
      // std::cout << "original points:" << std::endl;
      // for (auto p : polyhedron.get_points())
      //   std::cout << p << std::endl;
      // std::cout << "split points:" << std::endl;
      // for (size_t i : section_poly_verts)
      //   std::cout << result.vertices[i] << std::endl;
      // throw std::runtime_error("something wrong with the splitting");
    }

    // sort vertices so they form a polygon
    Polygon<Scalar>::reorder_indices(result.vertices.points, section_poly_verts);
    result.polygons.push_back(std::move(section_poly_verts));
    parent_face_indices.push_back(result.polygons.size());
    result.markers.push_back(marker_split);
    return true;
  }
  else return false;
}

template <typename Scalar>
bool coincide(const Line<3,Scalar> & line, const Plane<Scalar> & plane, const double tol = 1e-10)
{
  if (std::fabs(line.direction().dot(plane.normal())) < tol)
    if (std::fabs(plane.signed_distance(line.point())) < tol)
      return true;
  return false;
}

// Compute the intersection  of a line and a plane
//  throws std::runtime_error
template <typename Scalar>
bool collision(const Line<3,Scalar> & line,
               const Plane<Scalar>  & plane,
               Point<3,Scalar>      & intersection,
               const double           tol = 1e-8)
{
  // Plane : (p - p0) · n = 0
  // line p = d l + l0
  // intersection: d = (p0 - l0) · n / (l · n)
  // Note: if (l · n) == 0 then line is parallel to the plane
  if (std::fabs(line.direction().dot(plane.normal())) < tol)
  {
    if (coincide(line, plane, tol))
      throw std::runtime_error("line and plane coinside.");
    return false;
  }

  const Scalar d = (plane.origin() - line.point()).dot(plane.normal()) /
      line.direction().dot(plane.normal());
  intersection = line.point() + d*line.direction();
  return true;
}


// find intersection between a polygon and a line
//  section is a vector cause line can reside on polygon
// appends to vector intersection
// note: polygon should have sorted points
template <typename Scalar>
bool collision(const Line<3,Scalar>         & line,
               const Polygon<Scalar>        & poly,
               std::vector<Point<3,Scalar>> & intersection,
               const double tol = 1e-10)
{
  // find intersection between polygon plane and line
  Point<3,Scalar> p;
  if (coincide(line, poly.plane()))
  {
    const auto & vertices = poly.get_points();
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      const size_t j = (i < vertices.size() - 1) ? (i+1) : 0;
      Line<3,Scalar> edge(vertices[i], vertices[j]-vertices[i]);
      std::vector<Point<3,Scalar>> sec;
      if (collision(line, edge, sec, tol))
      {
        // if between edge vertices
        const auto & p = sec.front();
        if ((p - vertices[i]).dot(p-vertices[j]) < 0)
          intersection.push_back(p);
      }
    }
    assert(intersection.size() == 2);
    return true;
  }
  else if (collision(line, poly.plane(), p))  // case colinear
  {
    return false;
  }
  else
  {
    if (poly.point_inside(p), 1e-4)
    {
      intersection.push_back(p);
      return true;
    }
    else return false;
  }
}


// collision of a line segment with a polygon
template <typename Scalar>
bool collision(const Point<3,Scalar>        & p0,
               const Point<3,Scalar>        & p1,
               const Polygon<Scalar>        & poly,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  std::vector<Point<3,Scalar>> new_section;
  if (poly.point_inside(p0, tol) && fabs(poly.plane().signed_distance(p0)) < tol)
    new_section.push_back(p0);
  if (poly.point_inside(p1, tol) && fabs(poly.plane().signed_distance(p1)) < tol)
    new_section.push_back(p1);

  if (new_section.empty())
  {
    const std::size_t ibegin = new_section.size();

    collision(p0, p1,  poly.plane(), new_section, tol);
    for (std::size_t i=ibegin; i<new_section.size(); ++i)
      if (!poly.point_inside(new_section[i], tol))
        new_section.erase(new_section.begin() + i);
  }

  for (const auto & p : new_section)
    intersection.push_back(p);

  if (new_section.empty())
    return false;
  else
    return true;
}


// collision of a line segment with a polyhedron
// the points l0 and l1 define a line segment
// the output is saved into intersecion
template <typename Scalar>
bool collision(const LineSegment<Scalar>    & segment,
               const Polyhedron<Scalar>     & poly,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  std::vector<Point<3,Scalar>> new_section;
  const auto l0 = segment.first();
  const auto l1 = segment.second();

  // points that are inside poly
  if (poly.point_inside(l0))
    new_section.push_back(l0);
  if (poly.point_inside(l1))
    new_section.push_back(l1);

  if (new_section.size() < 2)
  {
    const auto & points = poly.get_points();
    for (const auto & face : poly.get_faces())
    {
      const size_t ibegin = new_section.size();
      Polygon<Scalar> poly_face(points, face);
      collision(l0, l1,  poly_face.plane(), new_section, tol);

      for (size_t i = ibegin; i < new_section.size();)
        if (!poly_face.point_inside(new_section[i], tol))
          new_section.erase(new_section.begin() + i);
        else i++;
    }

    angem::remove_duplicates_slow(new_section, tol);
  }

  for (const auto & p : new_section)
    intersection.push_back(p);

  if (new_section.empty())
    return false;
  else
    return true;
}


// implements ray casting algorithm
// the surface doesn't have to be convex
// the code doesn't guarantee that the surface is enclosed
// NOTE: runtime is O(surface_vertices.size())
template <typename Scalar>
bool point_inside_surface(const Point<3,Scalar>                       & point,  // point
                          const Point<3,Scalar>                       & external,  // external point
                          const std::vector<angem::Point<3,double>>   & surface_vertices,
                          const std::vector<std::vector<std::size_t>> & surface_elements,
                          const double                                  tol = 1e-4)
{
  const angem::Shape<double> ray({point, external});

  std::vector<angem::Point<3,double>> section_points;
  for (const std::vector<std::size_t> & face_vertices : surface_elements)
  {
    const angem::Polygon<double> tria(surface_vertices,  face_vertices);
    angem::collision(point, external, tria, section_points);
  }

  angem::PointSet<3, double> unique_points(tol);
  for (const auto & p : section_points)
    unique_points.insert(p);

  // if odd return true
  return (unique_points.size() % 2 == 1);
}

}  // end namespace
