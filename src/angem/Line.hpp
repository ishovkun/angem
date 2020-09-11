// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"


namespace angem
{

/* Simple utility class defining a Line object.
 * A line is assigned by a direction vector and and point on the line.
 * I think this should rather be a struct.
 */
template <int dim, typename Scalar>
class Line
{
 public:
  // Creates a line from a point and direction vectors
  Line(const Point<dim,Scalar> & point,
       const Point<dim,Scalar> & direction);

  Point<dim,Scalar> point() const noexcept {return _point;}
  Point<dim,Scalar> direction() const noexcept {return _direction;}

  // distance between a line and a point
  Scalar distance(const Point<dim,Scalar> & p) const;

 private:
  Point<dim, Scalar> _point;
  Point<dim, Scalar> _direction;
};

template <int dim, typename Scalar>
Line<dim,Scalar>::Line(const Point<dim,Scalar> & point,
                       const Point<dim,Scalar> & direction)
    :
    _point(point),
    _direction(direction)
{
  _direction.normalize();
}

template <int dim, typename Scalar>
Scalar Line<dim,Scalar>::distance(const Point<dim,Scalar> & p) const
{
  const Point<dim,Scalar> & ab = _direction;
  const Point<dim,Scalar> ac = p - _point;
  const Scalar area = ab.cross(ac).norm();
  const Scalar cd = area / ab.norm();
  return cd;
}

}  // end namespace
