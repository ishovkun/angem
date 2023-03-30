#pragma once
#include "Polygon.hpp"

namespace angem {


template <typename T>
class OptimizeDistanceToEdges {
 public:
  OptimizeDistanceToEdges(T tol) : _tol(tol) {}
  Point<3,T> get(Polygon<T> const & poly) const;

 private:
  Point<3,T> computeGradient_(Point<3,T> const x, std::vector<Point<3,T>> const & vts) const;
  T computeValue_(Point<3,T> const x, std::vector<Point<3,T>> const & vts) const;

  T _tol;
};

template <typename T>
Point<3,T> OptimizeDistanceToEdges<T>::get(Polygon<T> const & poly) const
{
  auto vts = poly.get_points();
  auto const offset = vts.front();
  std::for_each(vts.begin(), vts.end(), [&offset](auto & p){p -= offset;});

  auto const cm = compute_center_mass(vts);
  auto x = cm;
  auto error = 2*_tol;
  T alpha = 0.9;
  // std::cout << "_tol = " << _tol << " err = " << error << std::endl;
  int iter = 0;
  while ( error > _tol ) {
    auto grad = computeGradient_(x, vts);
    x -= alpha * grad;
    error = grad.norm();
    // std::cout << "iter: " <<  iter << " error = " << error << " value = " << computeValue_(x, vts) << std::endl;
    if (++iter == 1000) {
      // std::cout << "error in optimization" << std::endl;
      break;
      throw std::runtime_error("error in finding pole of inaccessibility");
    }

  }
  // std::cout << "moved by " << x.distance(cm)  << std::endl;
  return x + offset;
}

template <typename T>
Point<3,T> OptimizeDistanceToEdges<T>::computeGradient_(Point<3,T> const p, std::vector<Point<3,T>> const & vts) const
{
  Point<3,T> sum_grad; sum_grad.set_zero();
  T sum_len = static_cast<T>(0.0);
  size_t const n = vts.size();
  for (size_t i = 0; i < n; ++i) {
    auto const & p1 = vts[i];
    auto const & p2 = vts[(i+1)%n];
    auto const p1p2 = p2 - p1;
    auto const len = p1p2.norm();
    auto const av = p1p2 / len;
    auto const xv = p - p1;
    auto const dist_to_line = (xv - xv.dot(av) * av).norm();

    auto const a = av[0], b = av[1], c = av[2];
    auto const x = xv[0], y = xv[1], z = xv[2];
    Point<3,T> grad; grad.set_zero();
    if ( dist_to_line > p.distance(p1) ) {
      grad = 2*(p - p1);
    }
    else if ( dist_to_line > p.distance(p2) ) {
      grad = 2*(p - p2);
    }
    else {
      grad[0] = + 2*(1-a*a)*(x*(1-a*a) - a*b*y - a*c*z)
                - 2*a*b*(y*(1-b*b) - a*b*x - b*c*z)
                - 2*a*c*(z*(1-c*c) - a*c*x - b*c*y);
      grad[1] = - 2*a*b*(x*(1-a*a) - a*b*y - a*c*z)
                + 2*(1-b*b)*(y*(1-b*b) - a*b*x - b*c*z)
                - 2*b*c*(z*(1-c*c) - a*c*x - b*c*y);
      grad[2] = - 2*a*c*(x*(1-a*a) - a*b*y - a*c*z)
                - 2*b*c*(y*(1-b*b) - a*b*x - b*c*z)
                + 2*(1-c*c)*(z*(1-c*c) - a*c*x - b*c*y);
    }
    sum_grad += len*len*grad;
    sum_len += len*len;
  }
  return sum_grad / sum_len;
}

template <typename T>
T OptimizeDistanceToEdges<T>::computeValue_(Point<3,T> const p, std::vector<Point<3,T>> const & vts) const
{
  T sum = static_cast<T>(0.);
  T sum_len = static_cast<T>(0.);
  size_t const n = vts.size();
  for (size_t i = 0; i < n; ++i) {
    auto const & p1 = vts[i];
    auto const & p2 = vts[(i+1)%n];
    auto const p1p2 = p2 - p1;
    auto const l = p1p2.norm();
    auto const av = p1p2 / l;
    auto const xv = p - p1;
    auto dist = (xv - xv.dot(av) * av).norm();
    dist = std::min( dist, p.distance(p1) );
    dist = std::min( dist, p.distance(p2) );
    sum += l*l*dist*dist;
    sum_len += l*l;
  }
  return sum / sum_len;
}



}  // end namespace angem
