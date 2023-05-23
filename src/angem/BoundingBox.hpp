#include "Point.hpp"
#include "Hexahedron.hpp"
#include "Tensor2.hpp"

#ifdef WITH_EIGEN
#include <Eigen/Dense>
#include "covariance.hpp"
#endif

namespace angem {

template<typename Scalar>
class BoundingBox {
 public:
  BoundingBox(std::vector<size_t> const & indices, std::vector<Point<3,Scalar>> const & all_verts)
      : _box(compute_(indices, all_verts))
  {}

  operator Hexahedron<Scalar>() const { return _box; }

  angem::Point<3, Scalar> dimensions() const {
    angem::Point<3, Scalar> ans;
    auto const & points = _box.get_points();
    ans[0] = points[0].distance(points[1]);
    ans[1] = points[0].distance(points[3]);
    ans[2] = points[0].distance(points[4]);
    return ans;
  }

  std::array<angem::Point<3,Scalar>, 3> principal_components() const {
    std::array<angem::Point<3,Scalar>, 3> ans;
    auto const & points = _box.get_points();
    ans[0] = points[1] - points[0];
    ans[1] = points[3] - points[0];
    ans[2] = points[4] - points[0];
    return ans;
  }

  bool is_valid() const {
    auto dims = dimensions();
    for (size_t i = 0; i < 3; ++i)
      if ( std::isnan(dims[i]) )
        return false;
    return true;
  }

 private:

#ifdef WITH_EIGEN
  angem::Hexahedron<Scalar> compute_(std::vector<size_t> const & indices,
                                     std::vector<Point<3,Scalar>> const & all_verts)
  {
    auto const c = compute_center_mass( all_verts, indices );
    auto cov = covariance<Scalar>(indices, all_verts);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(cov, Eigen::ComputeEigenvectors);
    Eigen::Matrix3f ev = eigen_solver.eigenvectors();
    // make sure basis is right-handed
    ev.col(2) = ev.col(0).cross( ev.col(1) );

    Tensor2<3,double> transform({
        ev(0, 0), ev(1, 0), ev(2, 0),
        ev(0, 1), ev(1, 1), ev(2, 1),
        ev(0, 2), ev(1, 2), ev(2, 2),
      });

    // compute envelope
    double const upper = std::numeric_limits<double>::max();
    double const lower = std::numeric_limits<double>::lowest();
    angem::Point<3,double> bbox_min = {upper, upper, upper};
    angem::Point<3,double> bbox_max = {lower, lower, lower};
    for (size_t k = 0; k < indices.size(); ++k) {
      auto const & p = all_verts[indices[k]];
      auto const pt = transform * ( p - c );
      for (size_t i = 0; i < 3; ++i) {
        bbox_min[i] = std::min( bbox_min[i], pt[i] );
        bbox_max[i] = std::max( bbox_max[i], pt[i] );
      }
    }

    auto hex = build_hexahedron_(bbox_min, bbox_max);

    auto const transform_inv = invert(transform);
    for ( auto & v : hex.get_points() ) {
      v = c + transform_inv * v;
    }
    return hex;
  }

  angem::Hexahedron<Scalar> build_hexahedron_(angem::Point<3,double> const & min,
                                              angem::Point<3,double> const & max)
  {
    auto const delta = max - min;
    std::vector<angem::Point<3,Scalar>> verts(8);  // hex has 8 vertices
    std::fill( verts.begin(), verts.begin() + 4, min );
    verts[1][0] += delta[0];
    verts[2][0] += delta[0];
    verts[2][1] += delta[1];
    verts[3][1] += delta[1];

    std::copy( verts.begin(), verts.begin() + 4, verts.begin() + 4 );
    std::for_each( verts.begin() + 4, verts.end(), [&delta](auto & v) { v[2] += delta[2];} );

    std::vector<size_t> hex_indices(verts.size());
    std::iota( hex_indices.begin(), hex_indices.end(), 0 );
    return angem::Hexahedron<Scalar>(verts, hex_indices);
  }

#else
  angem::Hexahedron<Scalar> compute_(std::vector<size_t> const & indices,
                                     std::vector<Point<3,Scalar>> const & all_verts)
  {
    static_assert(false, "Cannot compute bounding box without Eigen");
  }
#endif


  angem::Hexahedron<double> _box;
};

}  // end namespace angem
