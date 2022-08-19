#pragma once

#include "angem/Shape.hpp"
#include "angem/Line.hpp"

namespace angem {

template<typename Scalar>
class LineSegment : public Shape<Scalar>
{
 public:
  LineSegment(const angem::Point<3,Scalar> & p1, const angem::Point<3,Scalar> & p2);
  inline angem::Point<3,Scalar> const & first() const noexcept;
  inline angem::Point<3,Scalar> const & second() const noexcept;
  inline angem::Point<3,Scalar> & first()  noexcept { return this->points[0]; }
  inline angem::Point<3,Scalar> & second() noexcept { return this->points[1]; }
  angem::Line<3,Scalar> line() const noexcept;
  Scalar length() const noexcept;
  // check whether point is within line segment
  bool contains(angem::Point<3,Scalar> const &p, Scalar tol = 1e-8) const noexcept;
  // distance between line segment and point
  Scalar distance(angem::Point<3,Scalar> const &p) const noexcept;
  inline angem::Point<3,Scalar> direction() const noexcept;
};

template<typename Scalar>
LineSegment<Scalar>::LineSegment(const angem::Point<3,Scalar> & p1, const angem::Point<3,Scalar> & p2)
{
  if (p1.distance(p2) < 1e-10) throw std::invalid_argument("Cannot create a line segment");
  this->points.push_back(p1);
  this->points.push_back(p2);
}

template<typename Scalar>
const angem::Point<3,Scalar> & LineSegment<Scalar>::first() const noexcept
{
  return this->get_points().front();
}

template<typename Scalar>
const angem::Point<3,Scalar> & LineSegment<Scalar>::second() const noexcept
{
  return this->get_points().back();
}

template<typename Scalar>
angem::Line<3,Scalar> LineSegment<Scalar>::line() const noexcept
{
  return angem::Line<3,Scalar>(first(), second() - first());
}

template<typename Scalar>
Scalar LineSegment<Scalar>::length() const noexcept
{
  return first().distance(second());
}

template<typename Scalar>
bool LineSegment<Scalar>::contains(angem::Point<3,Scalar> const &p, Scalar tol) const noexcept
{
  auto const & ab = second() - first();
  auto const & ac = p - first();
  Scalar const l_ab = ab.norm();
  auto const & ab_n = ab / l_ab;
  Scalar const proj_length = ac.dot(ab_n);
  if      ( proj_length < tol * l_ab ) return false;
  else if ( proj_length > l_ab * (1+tol) ) return false;
  else return line().distance(p) < tol * l_ab;
}

template<typename Scalar>
Scalar LineSegment<Scalar>::distance(angem::Point<3,Scalar> const &p) const noexcept
{
  auto const & ab = second() - first();
  auto const & ac = p - first();
  Scalar const l_ab = ab.norm();
  auto const & ab_n = ab / l_ab;
  Scalar const proj_length = ac.dot(ab_n);
  if      ( proj_length < 0 ) return first().distance(p);
  else if ( proj_length > l_ab ) return second().distance(p);
  else return line().distance( p );
}

template<typename Scalar>
angem::Point<3,Scalar> LineSegment<Scalar>::direction() const noexcept
{
  return (this->points[1] - this->points[0]).normalize();
}


}  // end namespace angem
