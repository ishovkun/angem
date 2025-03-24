#pragma once
#include "Point.hpp"
#include "Tensor2.hpp"
#include "utils.hpp"

namespace angem {
  
template <typename T>
class Transform
{
 public:
  /*
  ** Constructor: defines a rotation.
  ** Input:
  ** \param[in] origin : center of rotation
  ** \param[in] axis   : axis of rotation
  ** \param[in] anble  : rotation angle [radians]
   */
  static Transform<T> createRotation(Point<3,T> axis, double angle) {
      return Transform (angem::Point<3,double>(),
                        build_rotation_matrix(axis, angle));
  }
  
  static Transform<T> createScaling(Point<3,T> scale) {
    Tensor2<3,T> mat{scale[0], 0., 0.,
                      0., scale[1], 0.,
                      0., 0., scale[2]};
    return Transform<T> (angem::Point<3,double>(), mat);
  } 

  static  Transform<T> createTranslation(Point<3,T>const& offset) {
    return Transform<T>(offset, Tensor2<3,T>::make_unit_tensor());
  }

  static Transform<T> createIdentity() {
    return Transform<T>(Point<3,T>(), Tensor2<3,T>::make_unit_tensor());
  }

  /*
   ** Constructor: defines a transform
   */
  Transform(Point<3,T>const & offset, angem::Tensor2<3,T> const &matrix)
    : _offset(offset)
    , _matrix(matrix)

  {}
  // default copy assignment
  Transform<T> & operator=( const Transform<T> & other );
  // in-place apply transformation to a range of angem points
  template<typename IteratorType>
  void apply(IteratorType begin, IteratorType const end) const;
  // apply to a single point
  Point<3,T> apply(Point<3,T> const & p) const;
  
private:

  Point<3,T> _offset;
  Tensor2<3,T> _matrix;
};

template<typename T>
Transform<T> & Transform<T>::operator=( const Transform<T> & other )
{
  _offset = other._offset;
  _matrix = other._matrix;
  return *this;
}

template<typename T>
template<typename IteratorType>
void Transform<T>::apply(IteratorType begin, IteratorType const end) const
{
  for (auto it = begin; it != end; ++it) {
    auto v = *it + _offset;
    v = _matrix * v;
    *it = v;
  }
}

template <typename T>
Point<3,T> Transform<T>::apply(angem::Point<3,T> const & p) const
{
  auto v = p + _offset;
  v = _matrix * v;
  // v += _origin;
  return v;
}


}  // end namespace angem
