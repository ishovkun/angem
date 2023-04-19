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
  Envelope() {
    for (size_t i = 0; i < dim; ++i) {
      _min[i] = std::numeric_limits<Scalar>::max();
      _max[i] = std::numeric_limits<Scalar>::lowest();
    }
  }

  Envelope(std::vector<Point<dim,Scalar>> const & cloud)
  {
    for (size_t i = 0; i < dim; ++i) {
      _min[i] = std::numeric_limits<Scalar>::max();
      _max[i] = std::numeric_limits<Scalar>::lowest();
    }

    for ( auto const & p : cloud )
      update(p);
  }

  void udpate(Point<3, Scalar> const & p) {
    for (size_t i = 0; i < dim; ++i) {
      _min[i] = std::min( _min[i], p[i] );
      _max[i] = std::max( _max[i], p[i] );
    }
  }

  Envelope ( std::vector<size_t> const & indices,
             std::vector<angem::Point<dim,Scalar>> const & all_verts )
  {
    for (size_t i = 0; i < dim; ++i) {
      _min[i] = std::numeric_limits<Scalar>::max();
      _max[i] = std::numeric_limits<Scalar>::lowest();
    }
    for (size_t k = 0; k < indices.size(); ++k) {
      auto const & p = all_verts[indices[k]];
      for (size_t i = 0; i < dim; ++i) {
        _min[i] = std::min( _min[i], p[i] );
        _max[i] = std::max( _max[i], p[i] );
      }
    }
  }

  inline Point<dim,Scalar> min() const { return _min; }
  inline Point<dim,Scalar> max() const { return _max; }
  inline Point<dim,Scalar> size() const { return _max - _min; }

 private:
  Point<dim,Scalar> _min, _max;
};



}  // end namespace angem
