#pragma once
#include "Polygon.hpp"
#include "Projections.hpp"
#include "Envelope.hpp"
#include <queue>

namespace angem {

/*
 * This class implements a fast algorithm for finding polygon pole of inaccessibility,
 * i.e. the most distant internal point from the polygon outline (not to be confused with centroid)
 *
 * The code is adapted from https://github.com/mapbox/polylabel to angem primitive types
 */
template <typename T>
class Polylabel {
 public:
  Polylabel(T precision) : _precision(precision)
  {}

  Point<3,T> get(Polygon<T> const & poly) const;

 private:

  Point<2,T> get_polylabel_(std::vector<Point<2,T>> const & pts) const;

  static T point_to_polygon_dist(Point<2,T> const & center, std::vector<Point<2,T>> const & poly);

  struct Quad {
    Quad(const Point<2,T>& center, T half_size, std::vector<Point<2,T>> const & pts)
        : center(center),
          half_size(half_size),
          distance(point_to_polygon_dist(center, pts)),
          potential(distance + half_size * std::sqrt(2))
        {}

    Point<2,T> center; // cell center
    T half_size; // half the cell size
    T distance; // distance from cell center to polygon
    T potential; // max distance to polygon within a cell (defined as distance + cell radius)
  };

  static Quad get_centroid_cell(std::vector<Point<2,T>> const & pts);

  T _precision;
};

template <typename T>
Point<3,T> Polylabel<T>::get(Polygon<T> const & poly) const
{
  auto verts = poly.get_points();

  auto const offset = verts.front();
  std::for_each(verts.begin(), verts.end(), [&offset](auto & p){
     p -= offset;
  });

  angem::Plane<T> plane( verts[0], angem::polygon_average_normal(verts) );
  auto verts2d = plane.get_planar_coordinates(verts);
  auto p2d = get_polylabel_(verts2d);
  auto basis = plane.get_basis();
  auto p3d = offset + plane.origin() + p2d[0]*basis[0] + p2d[1]*basis[1];
  return p3d;
}

template <typename T>
Point<2,T> Polylabel<T>::get_polylabel_(std::vector<Point<2,T>> const & pts) const
{
  // std::cout << "input points:" << std::endl;
  // for ( auto & p : pts ) {
  //   std::cout << p << std::endl;
  // }
  Envelope<2,T> env(pts);
  T const cell_size = std::min( env.size()[0], env.size()[1] );
  if ( std::isnan(1./cell_size) ) return env.min();

  T h = 0.5 * cell_size;

  // a priority queue of cells in order of their "potential" (max distance to polygon)
  using C = Polylabel<T>::Quad;
  auto cmp = [] (C const &a, C const &b) {return a.potential < b.potential;};
  std::priority_queue<C, std::vector<C>, decltype(cmp)> cellq(cmp);

  // cover polygon with initial cells
  for ( T x = env.min()[0]; x < env.max()[0]; x += cell_size )
    for ( T y = env.min()[1]; y < env.max()[1]; y += cell_size )
    {
      C newc({x+h, y+h}, h, pts);
      cellq.push( std::move(newc) );
    }

  // take centroid as the first best guess
  auto best_cell = get_centroid_cell(pts);

  // second guess: bounding box centroid
  C bbox_cell(0.5*(env.min() + env.max()), 0, pts);

  best_cell = (bbox_cell.distance > best_cell.distance) ? bbox_cell : best_cell;

  auto n_prob = cellq.size();
  while ( !cellq.empty() ) {
    // pick the most promising cell from the queue
    auto cell = std::move(cellq.top());
    cellq.pop();

    // update the best cell if we found a better one
    if ( cell.distance > best_cell.distance ) {
      best_cell = cell;
      // std::cout << "found best " << ::round(1e4 * cell.distance) / 1e4 << " after " << n_prob << " probes" << std::endl;
    }

    // do not search down further if there is no chance of a better soultion
    if ( cell.potential - best_cell.distance > _precision ) {

      // split the cell into four cells
      h = 0.5*cell.half_size;
      cellq.push(C({cell.center[0] - h, cell.center[1] - h}, h, pts));
      cellq.push(C({cell.center[0] + h, cell.center[1] - h}, h, pts));
      cellq.push(C({cell.center[0] - h, cell.center[1] + h}, h, pts));
      cellq.push(C({cell.center[0] + h, cell.center[1] + h}, h, pts));
      n_prob += 4;
    }
  }

  // std::cout << "num probes: " << n_prob << std::endl;
  // std::cout << "best distance: " << best_cell.distance << std::endl;
  // std::cout << "best c = " << best_cell.center << std::endl;
  // std::cout << std::endl;

  return best_cell.center;
}

// get squared distance from a point to a segment
template <class T>
T getSegDistSq(Point<2,T> const &p, Point<2,T> const& a, Point<2,T> const& b) {
  auto x = a.x();
  auto y = a.y();
  auto dx = b.x() - x;
  auto dy = b.y() - y;

  if (dx != 0 || dy != 0) {

    auto t = ((p[0] - x) * dx + (p[1] - y) * dy) / (dx * dx + dy * dy);

    if (t > 1) {
      x = b.x();
      y = b.y();

    } else if (t > 0) {
      x += dx * t;
      y += dy * t;
    }
  }

  dx = p.x() - x;
  dy = p.y() - y;

  return dx * dx + dy * dy;
}

template <typename T>
T Polylabel<T>::point_to_polygon_dist(Point<2,T> const & p, std::vector<Point<2,T>> const & pts)
{
  bool inside = false;
  auto midDistSq = std::numeric_limits<T>::max();

  size_t const n = pts.size();
  for (size_t i = 0; i < n; ++i) {
    auto const & p1 = pts[i];
    auto const & p2 = pts[(i+1)%n];

    if ((p1[1] > p[1]) != (p2[1] > p[1]) && (p[0] < (p2[0] - p1[0]) * (p[1] - p1[1]) / (p2[1] - p1[1]) + p1[0]))
      inside = !inside;

    midDistSq = std::min( midDistSq, getSegDistSq(p, p1, p2) );
  }

  return (inside ? 1 : -1) * std::sqrt(midDistSq);
}

template <typename T>
typename Polylabel<T>::Quad Polylabel<T>::get_centroid_cell(std::vector<Point<2,T>> const & pts)
{
  T area = 0;
  Point<2,T> c(0,0);
  for (size_t i = 0; i < pts.size(); ++i) {
    auto const & a = pts[i];
    auto const & b = pts[(i+1)%pts.size()];
    auto f = a.x() * b.y() - b.x() * a.y();
    c[0] += (a[0] + b[0]) * f;
    c[1] += (a[1] + b[1]) * f;
    area += 3*f;
  }
  c = std::isnan(1./area) ? pts[0] : c / area;
  return Polylabel<T>::Quad( c, 0, pts );
}


}  // end namespace angem
