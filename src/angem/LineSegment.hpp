#pragma once

#include "angem/Shape.hpp"
#include "angem/Line.hpp"

namespace angem {

template<typename Scalar>
class LineSegment : public Shape<Scalar>
{
 public:
  LineSegment(const angem::Point<3,Scalar> & p1, const angem::Point<3,Scalar> & p2);
  const angem::Point<3,Scalar> & first() const noexcept;
  const angem::Point<3,Scalar> & second() const noexcept;
  angem::Line<3,Scalar> line() const noexcept;
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
  assert( this->get_points().size() == 2 );
  return this->get_points().front();
}

template<typename Scalar>
const angem::Point<3,Scalar> & LineSegment<Scalar>::second() const noexcept
{
  assert( this->get_points().size() == 2 );
  return this->get_points().back();
}

template<typename Scalar>
angem::Line<3,Scalar> LineSegment<Scalar>::line() const noexcept
{
  return angem::Line<3,Scalar>(first(), second() - first());
}

}  // end namespace angem
