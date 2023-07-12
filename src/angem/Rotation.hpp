#pragma once
#include "Point.hpp"
#include "Tensor2.hpp"

namespace angem {

template <typename T>
class Rotation
{
 public:
  /*
  ** Constructor: defines a rotation.
  ** Input:
  ** \param[in] origin : center of rotation
  ** \param[in] axis   : axis of rotation
  ** \param[in] anble  : rotation angle [radians]
   */
  Rotation(Point<3,T> origin, Point<3,T> axis, double angle);
  // default copy assignment
  Rotation<T> & operator=( const Rotation<T> & other );
  // in-place apply transformation to a range of angem points
  template<typename IteratorType>
  void apply(IteratorType begin, IteratorType const end) const;

private:
  Tensor2<3,T> build_matrix_();

  Point<3,T> _origin;
  Point<3,T> _axis;
  double _angle;
  Tensor2<3,T> _matrix;
};

template<typename T>
Rotation<T>::Rotation(Point<3,T> origin, Point<3,T> axis, double angle)
    : _origin(origin), _axis(axis.normalized()), _angle(angle),
      _matrix(build_matrix_())
{}

template<typename T>
Rotation<T> & Rotation<T>::operator=( const Rotation<T> & other )
{
  _origin = other._origin;
  _axis = other._axis;
  _angle = other._angle;
  _matrix = build_matrix_();
  return *this;
}

template<typename T>
template<typename IteratorType>
void Rotation<T>::apply(IteratorType begin, IteratorType const end) const
{
  for (auto it = begin; it != end; ++it) {
    auto v = *it - _origin;
    v = _matrix * v;
    v += _origin;
    *it = v;
  }
}

template<typename T>
Tensor2<3,T> Rotation<T>::build_matrix_()
{
  T const c = std::cos(_angle);
  T const s = std::sin(_angle);
  T const t = static_cast<T>(1) - c;
  auto const & a = _axis;

  Tensor2<3,T> m;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      m(i, j) = t * a[i]*a[j];

  for (size_t i = 0; i < 3; ++i)
    m(i, i) += c;

  m(0, 1) += - a[2]*s;
  m(0, 2) += + a[1]*s;
  m(1, 0) += + a[2]*s;
  m(2, 0) += - a[1]*s;
  m(2, 1) += - a[0]*s;
  m(1, 2) += + a[0]*s;

  return m;
}

}  // end namespace angem
