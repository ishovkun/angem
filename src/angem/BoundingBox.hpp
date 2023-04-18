#include "Point.hpp"
#include "Hexahedron.hpp"
#include "Tensor2.hpp"

#ifdef WITH_EIGEN
#include <Eigen/Dense>
#endif

namespace angem {

template<typename Scalar>
class BoundingBox {
 public:
  BoundingBox(std::vector<Point<3,Scalar>> const & cloud)
      : _box(compute_(cloud))
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
  angem::Hexahedron<Scalar> compute_(std::vector<Point<3,Scalar>> const & cloud)
  {
    auto const c = compute_center_mass(cloud);
    Eigen::Matrix3f covariance = Eigen::Matrix3f::Zero();
    for (const auto& p : cloud) {
        Eigen::Vector3f diff(p[0] - c[0], p[1] - c[1], p[2] - c[2]);
        covariance += diff * diff.transpose();
    }
    covariance /= cloud.size();
    // for (size_t i = 0; i < 3; ++i)
    //   for (size_t j = 0; j < 3; ++j)
    //   {
    //     for (size_t k = 0; k < cloud.size(); ++k)
    //       covariance(i, j) += ( cloud[k][i] - c[i]) * ( cloud[k][j] - c[j] );
    //     // covariance(i, j) /= cloud.size() - 1;
    //     covariance(i, j) /= cloud.size();
    //   }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
    Eigen::Matrix3f ev = eigen_solver.eigenvectors();
    // make sure basis is right-handed
    ev.col(2) = ev.col(0).cross( ev.col(1) );

    Tensor2<3,double> transform({
        ev(0, 0), ev(1, 0), ev(2, 0),
        ev(0, 1), ev(1, 1), ev(2, 1),
        ev(0, 2), ev(1, 2), ev(2, 2),
      });

    // compute limits
    double const upper = std::numeric_limits<double>::max();
    double const lower = std::numeric_limits<double>::lowest();
    angem::Point<3,double> bbox_min = {upper, upper, upper};
    angem::Point<3,double> bbox_max = {lower, lower, lower};

    for (auto const & p : cloud)
    {
      auto const pt = transform * ( p - c );
      for (size_t i = 0; i < 3; ++i) {
        bbox_min[i] = std::min( bbox_min[i], pt[i] );
        bbox_max[i] = std::max( bbox_max[i], pt[i] );
      }
    }
    auto const delta = bbox_max - bbox_min;
    std::vector<angem::Point<3,Scalar>> verts(8);  // hex has 8 vertices
    std::fill( verts.begin(), verts.begin() + 4, bbox_min );
    verts[1][0] += delta[0];
    verts[2][0] += delta[0];
    verts[2][1] += delta[1];
    verts[3][1] += delta[1];

    std::copy( verts.begin(), verts.begin() + 4, verts.begin() + 4 );
    std::for_each( verts.begin() + 4, verts.end(), [&delta](auto & v) { v[2] += delta[2];} );

    // inverse transform for real coordinates
    auto const transform_inv = invert(transform);
    for (auto & v : verts)
      v = c + transform_inv * v;

    // range array
    std::vector<size_t> indices(verts.size());
    std::iota( indices.begin(), indices.end(), 0 );
    return angem::Hexahedron<Scalar>(verts, indices);
  }
#elif
  angem::Hexahedron<Scalar> compute_(std::vector<Point<3,Scalar>> const & cloud)
  {
    static_assert(false, "Cannot compute bounding box without Eigen");
  }
#endif


  angem::Hexahedron<double> _box;
};

}  // end namespace angem
