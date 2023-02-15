#pragma once
#include "Point.hpp"
#include <limits>

namespace angem {

/*
** Kind of a bounding box.
** Stores minimum and maximum of the point cloud
*/
template<int dim, typename Scalar>
class Envelope
{
 public:
  Envelope(std::vector<Point<dim,Scalar>> const & cloud)
  {
    constexpr Scalar lowest = std::numeric_limits<Scalar>::lowest();
    constexpr Scalar largest = std::numeric_limits<Scalar>::max();
    for (size_t i = 0; i < dim; ++i) {
      _min[i] = largest;
      _max[i] = lowest;
    }

    for ( auto const & p : cloud )
      for (size_t i = 0; i < dim; ++i) {
        _min[i] = std::min( _min[i], p[i] );
        _max[i] = std::max( _max[i], p[i] );
      }
  }

  inline Point<dim,Scalar> min() const { return _min; }
  inline Point<dim,Scalar> max() const { return _max; }
  inline Point<dim,Scalar> size() const { return _max - _min; }

 private:
  Point<dim,Scalar> _min, _max;
};



}  // end namespace angem
