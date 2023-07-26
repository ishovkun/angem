#pragma once
#include "Point.hpp"
#ifdef WITH_EIGEN
#include <Eigen/Dense>

namespace angem {

template<typename Scalar>
Eigen::Matrix3f covariance(std::vector<Point<3,Scalar>> const & cloud)
{
  auto const c = compute_center_mass(cloud);
  Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();
  for (size_t k = 0; k < cloud.size(); ++k) {
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        cov(i, j) += ( cloud[k][i] - c[i]) * ( cloud[k][j] - c[j] );
  }
  cov /= cloud.size();
  return cov;
}

template<typename Scalar>
Eigen::Matrix3f covariance(std::vector<size_t> const & indices,
                           std::vector<Point<3,Scalar>> const & all_points)
{
  Point<3,Scalar> c;
  compute_center_mass( all_points, indices, c);

  Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();
  for (size_t k = 0; k < indices.size(); ++k) {
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        cov(i, j) += ( all_points[indices[k]][i] - c[i]) * ( all_points[indices[k]][j] - c[j] );
  }
  cov /= indices.size();
  return cov;
}

}  // end namespace angem

#endif
