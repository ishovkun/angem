// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"
#include "Basis.hpp"
#include "Line.hpp"
#include "utils.hpp"

#include <math.h>    // M_PI
#include <algorithm> // clamp

namespace angem
{

/* 3D plane class. A plane is defined by a point on the plane and
 * a normal vector.
 */
template <typename Scalar>
class Plane
{
 public:
  // Default empty constructor. Creates an invalid plane.
  // use it only if set data after.
  Plane();
  // create plane from a point on the plane and a normal vector
  Plane(const Point<3,Scalar> & point,
        const Point<3,Scalar> & normal);
  // strike and dip in degrees
  Plane(const Point<3,Scalar> & point,
        const Scalar          & dip_angle,      // -90 <= dip <= 90
        const Scalar          & strike_angle);
  // create plane from 3 points
  Plane(const Point<3,Scalar> & p1,
        const Point<3,Scalar> & p2,
        const Point<3,Scalar> & p3);

  // pick some three points from the cloud to initialize a plane
  // convenient when have colinear edges in a polygon but want to initialize
  // a plane
  // Be very careful with this method sine it does not check whether
  // all points belong to the same plane
  Plane(const std::vector<Point<3,Scalar>> & cloud);
  // pick some three points from the cloud to initialize a plane
  // convenient when have colinear edges in a polygon but want to initialize
  // a plane
  // Be very careful with this method sine it does not check whether
  // all points belong to the same plane
  Plane(const std::vector<Point<3,Scalar>> & cloud,
        const std::vector<size_t> & specific_indices);

  // setter
  void set_data(const Point<3,Scalar> & p1,
                const Point<3,Scalar> & p2,
                const Point<3,Scalar> & p3);
  // If you need some custom vectors in basis, you can specify them
  inline void set_basis(const Basis<3,Scalar> & basis) noexcept {_basis = basis;}
  // shift support point in direction p
  void move(const Point<3,Scalar> & p) noexcept;
  // get const reference to the plane normal vector
  const Point<3,Scalar> & normal() const noexcept {return _basis(2);}
  // get non-const reference to the basis object
  Basis<3,Scalar> & get_basis() noexcept {return _basis;}
  // get const reference to the basis object
  const Basis<3,Scalar> & get_basis() const noexcept {return _basis;}
  // return point that is in the plane
  const Point<3,Scalar> & origin() const noexcept {return _origin;}

  // compute strike angle (degrees) from normal
  Scalar strike_angle() const;
  // compute dip angle (degrees) from normal
  Scalar dip_angle() const;

  // get point projection on the plane (in global coordinates)
  Point<3,Scalar> project_point(const Point<3,Scalar> & p) const;
  // get points projection on the plane. convenient wrapper around previous method
  std::vector<Point<3,Scalar>> project_points(const std::vector<Point<3,Scalar>> & points) const;
  /**
   * Get the coordinates of the points in the basis
   */
  Point<3,Scalar> local_coordinates(const Point<3,Scalar> & p) const;
  // project vector (no account for plane location)
  Point<3,Scalar> project_vector(const Point<3,Scalar> & p) const;
  // signed distance from point to plane (> 0 if point is above plane)
  Scalar signed_distance(const Point<3,Scalar> & p) const;

  // true if point is above the plane
  bool above(const Point<3,Scalar> & p) const;
  // set point located in the plane
  inline void set_origin(const Point<3,Scalar> & p) {this->_origin = p;}
  // get d coefficient in plane equation ax + by + cz = d
  inline Scalar get_zero_intercept() const noexcept;
  // change basis that the normal is pointing in the opposite direction
  const angem::Point<3,Scalar> & invert_normal() noexcept;

  // ATTRIBUTES
 protected:
  // return two orthogonal vectors within the plane
  void compute_basis(const Point<3,Scalar> & normal);

  Point<3,Scalar> _origin; // point on the plane
  Basis<3,Scalar> _basis;  // tangent and normal vectors
};


template <typename Scalar>
Plane<Scalar>::Plane()
{}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Point<3,Scalar> & normal)
    : _origin(point)
{
  compute_basis(normal);
}

template <typename Scalar>
Plane<Scalar>::Plane(const std::vector<Point<3,Scalar>> & cloud)
{
  if ( cloud.size() < 3 )
    throw std::invalid_argument("Cannot create a plane from two points only");

  const Point<3,Scalar> p1 = cloud[0];
  size_t v2 = 1;
  for (; v2 < cloud.size(); ++v2)
    if (p1.distance(cloud[v2]) > 1e-8)
      break;
  if ( v2 >= cloud.size() - 1 )
  {
    std::cout << "cloud:" << std::endl;
    for (auto & p : cloud)
      std::cout << p << " (dist " << p1.distance(p) << ")"<< std::endl;
    throw std::invalid_argument("Cannot initialize a plane from a point cloud. "
                                "All vertices are within 1e-8 to each other");
  }
  const auto p2 = cloud[v2];

  Point<3,Scalar> p3;
  bool found = false;
  const Line<3, Scalar> line(p1, p2 - p1);
  for (std::size_t i=v2+1; i < cloud.size(); ++i)
  {
    if (line.distance(cloud[i]) < 1e-6)
      continue;
    found = true;
    p3 = cloud[i];
  }
  if (!found)
  {
    for (auto & p : cloud)
      std::cout << p << std::endl;
    throw std::invalid_argument("Cannot initialize a plane from a point cloud");
  }
  set_data(p1, p2, p3);
}

template <typename Scalar>
Plane<Scalar>::Plane(const std::vector<Point<3,Scalar>> & cloud,
                     const std::vector<size_t> &specific_indices)
{
  assert( cloud.size() > 2 );
  assert( specific_indices.size() > 2 );
  const Point<3,Scalar> p1 = cloud[specific_indices[0]], p2 = cloud[specific_indices[1]];
  Point<3,Scalar> p3;
  bool found = false;
  Line<3, Scalar> line(p1, p2 - p1);
  for (std::size_t i=2; i<specific_indices.size(); ++i)
  {
    if (line.distance(cloud[specific_indices[i]]) < 1e-6)
      continue;
    found = true;
    p3 = cloud[specific_indices[i]];
  }
  if (!found)
  {
    for (const size_t i : specific_indices)
      std::cout << cloud[specific_indices[i]] << std::endl;
    throw std::invalid_argument("Cannot initialize a plane from a point cloud");
  }

  set_data(p1, p2, p3);
}

template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Scalar          & dip_angle,
                     const Scalar          & strike_angle)
    : _origin(point)
{
  assert(dip_angle >= - 90 and dip_angle <= 90);

  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  Point<3,Scalar> normal;
  // dip = polar angle between normal and z
  normal[0] = sin(rdip) * cos(rstrike + M_PI/2.);
  normal[1] = sin(rdip) * sin(rstrike + M_PI/2.);
  normal[2] = cos(rdip);

  compute_basis(normal);
}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & p1,
                     const Point<3,Scalar> & p2,
                     const Point<3,Scalar> & p3)
{
  set_data(p1, p2, p3);
}


template <typename Scalar>
void Plane<Scalar>::set_data(const Point<3,Scalar> & p1,
                             const Point<3,Scalar> & p2,
                             const Point<3,Scalar> & p3)
{
#ifndef NDEBUG
  if (p1 == p2 || p2 == p3 || p1 == p3)
    throw std::invalid_argument("Duplicated points while initializing plane");
  if ( (( p2 - p1 ).cross(p3 - p2)).norm() < 1e-6 )
    throw std::invalid_argument("Initializing plane with colinear vectors");
#endif

  _origin = p1;
  // define two tangent vectors
  const Point<3,Scalar> t1 = p1 - p2;
  const Point<3,Scalar> t2 = p1 - p3;
  Point<3,Scalar> normal = cross(t1, t2);
  normal.normalize();

  compute_basis(normal);
}

template <typename Scalar>
Scalar Plane<Scalar>::get_zero_intercept() const noexcept
{
  return normal().dot(_origin);
}


template <typename Scalar>
Scalar Plane<Scalar>::signed_distance(const Point<3,Scalar> & p) const
{
  /* dot product of point by perpendicular vector:
   * if < 0: point is below the plane
   */
  /* signed distance from point p (vertex) to plane
   * with normal n and containing point x0:
   * d = (p - x0) · n
   * if d<0, point below the plane
   */
  return (p - _origin).dot(normal());
}


template <typename Scalar>
bool Plane<Scalar>::above(const Point<3,Scalar> & p) const
{
  if (signed_distance(p) > 0)
    return true;
  else
    return false;
}


// Determine if a set of points lies within a plane
template <typename Scalar>
bool align_on_plane(const std::vector<Point<3, Scalar>> & points,
                    const Scalar tol = 1e-8)
{
  assert(points.size() > 2);

  if (points.size() == 3)
    return true;

  Plane<Scalar> plane(points);
  for (std::size_t i=2; i<points.size(); ++i)
    if ( fabs(plane.signed_distance(points[i])) > tol )
      return false;

  return true;
}


// return two orthogonal vectors within the plane
template <typename Scalar>
void Plane<Scalar>::compute_basis(const Point<3,Scalar> & normal)
{
  /*
   * Algorithm:
   * chose a random vector rv
   * project it on the plane - first tangent basis vector e1
   * complete the basis with e2 = n x e1
   */
  _basis[2] = normal;
  Point<3,Scalar> rv = normal;
  size_t cnt = 0;
  do
  {
    rv[cnt] += 1;
    cnt = (cnt < 2) ? cnt+1 : 0;
  } while (rv.cross(normal).norm() < 1e-8);

  Point<3, Scalar> e1 = project_vector(rv);
  e1.normalize();

  Point<3,Scalar> e2 = normal.cross(e1);
  e2.normalize();

  _basis[0] = e1;
  _basis[1] = e2;
}


// project the direction of the vector on the plane (no translation)
template <typename Scalar>
inline
Point<3, Scalar>
Plane<Scalar>::project_vector(const Point<3,Scalar> & p) const
{
  Scalar p_n = p.dot(_basis(2));
  return p - p_n * _basis(2);
}


// get point projection on the plane (in global coordinates considering plane translation)
template <typename Scalar>
inline
Point<3, Scalar>
Plane<Scalar>::project_point(const Point<3,Scalar> & p) const
{
  // 1: translate p' = p - s  (s - plane support point)
  // 2. project on normal p'n = p' · n
  // 3. Translate back j = s + j'
  Point<3,Scalar> p_prime = p - _origin;
  Point<3,Scalar> j_prime = project_vector(p_prime);
  return _origin + j_prime;
}


template <typename Scalar>
inline
std::vector<Point<3,Scalar>>
Plane<Scalar>::project_points(const std::vector<Point<3,Scalar>> & points) const
{
  std::vector<Point<3,Scalar>> result;
  result.reserve(points.size());
  for (const auto p : points)
    result.push_back(project_point(p));
  return result;
}

// get point projection on the plane (in both global coordinates and local basis)
template <typename Scalar>
Point<3,Scalar>
Plane<Scalar>::local_coordinates(const Point<3,Scalar> & p) const
{
  // translate
  Point<3,Scalar> p_prime = p - _origin;
  // project on basis vectors
  return _basis.transform(p_prime);
}


template <typename Scalar>
void
Plane<Scalar>::move(const Point<3,Scalar> & p) noexcept
{
  _origin += p;
}


template <typename Scalar>
Scalar
Plane<Scalar>::dip_angle() const
{
  Scalar rdip = static_cast<Scalar>(acos(_basis(2)[2]));
  double dip = 180. * rdip / M_PI;
  if (dip > 90.0)
    dip = 180. - dip;
  return dip;
}


template <typename Scalar>
Scalar
Plane<Scalar>::strike_angle() const
{
  Scalar rdip = static_cast<Scalar>(acos(_basis(2)[2]));

  // avoid taking acos(+- 1) -- causes errors due to roundoff
  // const double v1 = std::clamp( _basis(2)[0] / sin(rdip), -1.0, 1.0);
  // const double v2 = std::clamp( _basis(2)[1] / sin(rdip), -1.0, 1.0);
  const double v1 = std::min(std::max(-1.0, _basis(2)[0]), 1.0);
  const double v2 = std::min(std::max(-1.0, _basis(2)[1]), 1.0);

  Scalar rstrike_from_cos = acos(v1) - M_PI / 2.;
  Scalar rstrike_from_sin = asin(v2) - M_PI / 2.;

  Scalar strike;
  if (rstrike_from_sin >= 0 and rstrike_from_cos >= 0)
  {
    // strike = 180. * rstrike_from_cos / M_PI;
    strike = degrees(rstrike_from_cos);
  }
  else if (rstrike_from_sin < 0 and rstrike_from_cos > 0)
  {
    // strike = - 180. * rstrike_from_cos / M_PI;
    strike = - degrees(rstrike_from_cos);
  }
  else if (rstrike_from_sin > 0 and rstrike_from_cos < 0)
  {
    // strike = 180. * rstrike_from_cos / M_PI;
    strike = degrees(rstrike_from_cos);
    strike = fabs(strike);
  }
  else // if (rstrike_from_sin < 0 and rstrike_from_cos < 0)
  {
    // strike = 180. * rstrike_from_cos / M_PI;
    strike = degrees(rstrike_from_cos);
    strike = fabs(strike);
  }

  return strike;
}

template <typename Scalar>
const angem::Point<3,Scalar> & Plane<Scalar>::invert_normal() noexcept
{
  _basis[2] *= -1;
  const auto tmp = _basis[1];
  _basis[1] = _basis[0];
  _basis[0] = tmp;
  return normal();
}

}  // end namespace
