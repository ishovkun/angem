#pragma once
#include "Point.hpp"
#include <array>

namespace angem
{

/* This class implements a orthonormal basis
 */
template <int dim, typename Scalar>
class Basis
{
 public:
  // default constructor. Creates e.g. in 3D three vectors
  // {1,0,0}, {0,1,0}, {0,0,1}
  Basis();
  // create basis from vector of points
  // checks that vecs.size() == dim
  Basis(const std::vector<Point<dim,Scalar>> & vecs);
  // create basis from three vectors
  Basis(Point<dim,Scalar> const & v1,
        Point<dim,Scalar> const & v2,
        Point<dim,Scalar> const & v3);

  // setters
  // assign basis components from the vector of points
  void set_data(const std::vector<Point<dim,Scalar>> & vecs);

  // non-constant getter
  Point<dim,Scalar> & operator[] (int i);
  // constant getter
  const Point<dim,Scalar> & operator() (int i) const;

  // get coordinates of a point in this basis
  Point<dim,Scalar> transform(const Point<dim,Scalar> & p) const;
  // probably do not need this any more since i changed the default constructor
  // bool is_empty() const;

  // invert basis (the direction) while preserving right-handness
  // returns reference to this
  Basis<dim,Scalar> & invert();

  // printout
  template <int d, typename S>
  friend std::ostream &operator<<(std::ostream     & os,
                                  const Basis<d,S> & p);

 private:
  std::array<Point<dim,Scalar>, dim> _vectors;
};


template <int dim, typename Scalar>
Basis<dim,Scalar>::Basis(Point<dim,Scalar> const & v1,
                         Point<dim,Scalar> const & v2,
                         Point<dim,Scalar> const & v3)
{
  _vectors[0] = v1;
  _vectors[1] = v2;
  if (dim > 2) _vectors[2] = v3;
  for (size_t i = 0; i < dim; ++i)
    _vectors[i].normalize();
}

template <int dim, typename Scalar>
Basis<dim,Scalar>::Basis()
{
  for (int i=0; i<dim; ++i)
    _vectors[i][i] = static_cast<Scalar>(1);
}


template <int dim, typename Scalar>
Basis<dim,Scalar>::Basis(const std::vector<Point<dim,Scalar>> & vecs)
{
  assert(vecs.size() == dim);
  std::copy(vecs.begin(), vecs.begin() + dim, _vectors.begin());
}


template <int dim, typename Scalar>
void
Basis<dim,Scalar>::set_data(const std::vector<Point<dim,Scalar>> & vecs)
{
  assert(vecs.size() == dim);
  // vectors = vecs;
  std::copy(vecs.begin(), vecs.begin() + dim, _vectors.begin());
}


template <int dim, typename Scalar>
Point<dim,Scalar> &
Basis<dim,Scalar>::operator[](int i)
{
  assert(i < dim);
  return _vectors[i];
}


template <int dim, typename Scalar>
const Point<dim,Scalar> &
Basis<dim,Scalar>::operator()(int i) const
{
  assert(i < dim);
  return _vectors[i];
}


template <int dim, typename Scalar>
Point<dim,Scalar>
Basis<dim,Scalar>::transform(const Point<dim,Scalar> & p) const
{
  // assert(!is_empty());
  Point<dim,Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p.dot(_vectors[i]);
  return result;
}


// template <int dim, typename Scalar>
// bool
// Basis<dim,Scalar>::is_empty() const
// {
//   Scalar sum = static_cast<Scalar>(0);
//   for (int i=0; i<dim; ++i)
//     sum += _vectors[i].norm();
//   if (fabs(sum - static_cast<Scalar>(3)) < 1e-8)
//     return false;
//   else
//     return true;
// }

template <int dim, typename Scalar>
Basis<dim,Scalar> &
Basis<dim,Scalar>::invert()
{
  for (auto & v : _vectors)
    v *= -1;
  std::swap( _vectors[1] , _vectors[0] );  // preserve right orientation
  return *this;
}


template <int d, typename S>
std::ostream &operator<<(std::ostream & os, const Basis<d,S> & basis)
{
  os << "(" << basis(0) << ",\n " << basis(1) << ",\n " << basis(2) << ")";
  return os;
}

}  // end namespace
