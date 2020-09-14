#pragma once

#include "angem/Shape.hpp"
#include "angem/Line.hpp"

namespace angem {

template<typename Scalar>
class LineSegment : public Shape<Scalar>
{
 public:
  LineSegment(const angem::Point<3,Scalar> & p1, const ange::Point<3,Scalar> & p2);
  const angem::Point & first() const noexcept;
  const angem::Point & second() const noexcept;
};

template<typename Scalar>
LineSegment<Scalar>::LineSegment(const angem::Point<3,Scalar> & p1, const ange::Point<3,Scalar> & p2)
{
 this->points.push_back(p1);
 this->points.push_back(p2);
}

template<typename Scalar>
const angem::Point & LineSegment<Scalar>::first() const noexcept
{
  assert( get_points().size() == 2 );
  return get_points().front();
}

template<typename Scalar>
const angem::Point & LineSegment<Scalar>::second() const noexcept
{
  assert( get_points().size() == 2 );
  return get_points().back();
}


}  // end namespace angem
