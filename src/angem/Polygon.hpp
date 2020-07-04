#pragma once

#include "Shape.hpp"
#include "Plane.hpp"
#include "PointSet.hpp"

#include "utils.hpp"
#include <iostream>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

namespace angem
{

using Edge = std::pair<std::size_t, std::size_t>;

/* This class implements 2D polygons in 3D space.
 */
template<typename Scalar>
class Polygon: public Shape<Scalar>
{
 public:
  // Default empty constructor. Use only if assigning data later
	Polygon();
  // Create a polygon from a vector of points.
  // Vertices are ordered in a clock-wise manner upon creation.
  Polygon(const std::vector<Point<3,Scalar>> & points_list);
  // Helper constructor. Construct a polygon (face) from some mesh vertices.
  // Vertices are ordered in a clock-wise manner upon creation if reorder_vertices == true.
  // Otherwise the vertices are in the assigned order.
  Polygon(const std::vector<Point<3,Scalar>> & all_mesh_vertices,
          const std::vector<std::size_t>     & indices,
          const bool                           reorder_vertices = true);
  // Helper constructor. Construct a polygon (face) from some mesh vertices.
  // Vertices are ordered in a clock-wise manner upon creation.
  Polygon(const PointSet<3,Scalar>           & all_mesh_vertices,
          const std::vector<std::size_t>     & indices);

  // helper function: get a plane that contains an edge wich the normal in the
  // plane of the polygon
  Plane<Scalar> get_side(const Edge & edge) const;

  // returns true if point is inside a 3D shape formed by
  // lines going through the poly vertices in the direction of
  // poly plane normal
  bool point_inside(const Point<3,Scalar> & p,
                    const Scalar            tol = 1e-10) const;

  // compute the area of the polygon
  Scalar area() const;
  // compute the center of mass of the polygons
  virtual Point<3,Scalar> center() const override;
  // helper function: get a vector of edges represented by pairs of vertex indices
  std::vector<Edge> get_edges() const;

  // setter
  void set_data(const std::vector<Point<3,Scalar>> & point_list,
                const bool                           reorder_vertices = true);
  // shift all points in direction p
  virtual void move(const Point<3,Scalar> & p) override;
  // order vertices in a clockwise fashion
  static  void reorder(std::vector<Point<3,Scalar>> & points);
  // order indices vector so that the corresponding points are in a clockwise fashiok
  static  void reorder_indices(const std::vector<Point<3, Scalar>> &verts,
                               std::vector<std::size_t>            &indices);
  angem::Point<3,double> normal() const { return plane().normal(); }
  inline const Plane<Scalar> & plane() const { return m_plane; }
  inline Plane<Scalar> & plane() { return m_plane; }

 protected:
  // Attributes
  Plane<Scalar> m_plane;
};


template<typename Scalar>
Polygon<Scalar>::Polygon()
    :
    Shape<Scalar>::Shape()
{}


template<typename Scalar>
Polygon<Scalar>::Polygon(const std::vector<Point<3,Scalar>> & point_list)
{
  assert(point_list.size() > 2);
  set_data(point_list);
}


template<typename Scalar>
Polygon<Scalar>::Polygon(const std::vector<Point<3,Scalar>> & all_mesh_vertices,
                         const std::vector<std::size_t>     & indices,
                         const bool                           reorder_vertices)
{
  assert(indices.size() > 2);
  std::vector<Point<3,Scalar>> point_list;
  for (const std::size_t & i : indices)
    point_list.push_back(all_mesh_vertices[i]);
  set_data(point_list, reorder_vertices);
}


template<typename Scalar>
Polygon<Scalar>::Polygon(const PointSet<3,Scalar>           & all_mesh_vertices,
                         const std::vector<std::size_t>     & indices)
{
  assert(indices.size() > 2);
  std::vector<Point<3,Scalar>> point_list;
  for (const std::size_t i : indices)
    point_list.push_back(all_mesh_vertices[i]);
  set_data(point_list);
}


template<typename Scalar>
void
Polygon<Scalar>::set_data(const std::vector<Point<3,Scalar>> & point_list,
                          const bool                           reorder_vertices)
{
  assert(point_list.size() >= 3);
  this->points = point_list;
  if (reorder_vertices)
    reorder(this->points);
  const Point<3, Scalar> cm = compute_center_mass(point_list);
  m_plane = Plane<Scalar>(point_list);
  m_plane.set_origin( cm );
}


template<typename Scalar>
void
Polygon<Scalar>::move(const Point<3,Scalar> & p)
{
  Shape<Scalar>::move(p);
  plane().move(p);
}


template<typename Scalar>
void
Polygon<Scalar>::reorder(std::vector<Point<3, Scalar> > & points)
{
  const std::size_t n_points = points.size();
  assert(n_points > 2);
  if (n_points == 3) return;  // no need to sort a triangle

  const Point<3,Scalar> cm = compute_center_mass(points);
  Plane<Scalar> plane = Plane<Scalar>(points);

  // this is a version of Graham Scan but without enforcing convexity
  // Just sort vertices by the polar angle with respect to the first point
  // select the first point as the support in the polygon first tangent vector direction
  // Scalar max_dist = -std::numeric_limits<Scalar>::max();
  Scalar max_dist = std::numeric_limits<Scalar>::lowest();
  std::size_t first = points.size();
  const auto & x_direction = plane.get_basis()[0];
  for (std::size_t i=0; i<points.size(); ++i)
  {
    const Scalar dist = points[i].dot(x_direction);
    if (dist > max_dist)
    {
      max_dist = dist;
      first = i;
    }
  }
  // select a perpendicular y-direction
  // comput the cosine of the p1-p angle and the y direction
  // since angle ϵ [0, 2·π], and cosine is monotonically decreasing,
  // we can sort point by ther angle
  const auto & y_direction = plane.get_basis()[1];
  std::vector<double>  cosines(points.size(), 0.0);
  for (std::size_t i=0; i<points.size(); ++i)
  {
    if (i == first)
      cosines[i] = - 1.1;  // cos can't be less than one
                           // but we put the first point in the end of the list
    else
    {
      const auto diff = points[i] - points[first];
      cosines[i] = (y_direction).dot( diff ) / diff.norm();
    }
  }

  // initialize original index locations
  std::vector<size_t> idx(n_points);
  iota(idx.begin(), idx.end(), 0);
  // first sort by cosine
  std::sort(idx.begin(), idx.end(), [&cosines](size_t i1, size_t i2)
                                    {return cosines[i1] > cosines[i2];});

  /* Support for hanging nodes
   * If there are several vertices with the same value of cos,
   * we additionally sort them by distance
   * Note that we need to sort only vertices with cos = min(cosines), and
   * cos = max(cosines), since  we do not fully support non-convex polygons */
  // find lowest and highest values skipping -1.1
  Scalar lowest = std::numeric_limits<Scalar>::max(),
         highest = std::numeric_limits<Scalar>::lowest();
  for ( auto & val : cosines )
  {
    if (val != -1.1 && val < lowest)
      lowest = val;
    if (val > highest)
      highest = val;
  }

  // sort vertices with cos = lowest
  size_t ibegin = idx.size(), iend = idx.size();
  for (std::size_t i=0; i<idx.size(); ++i)
    if (std::fabs(cosines[idx[i]] - lowest) < 1e-6)
    {
      if (ibegin == idx.size())
        ibegin = i;
      iend = i;
    }

  if (ibegin  != iend)
    std::sort(idx.begin() + ibegin, idx.begin() + iend + 1,
              [&points, first](size_t i1, size_t i2)
              {
                return (points[first].distance(points[i1]) > points[first].distance(points[i2]));
              });

  // sort vertices with cos = highest
  ibegin = idx.size(), iend = idx.size();
  for (std::size_t i=0; i<idx.size(); ++i)
    if (std::fabs(cosines[idx[i]] - highest) < 1e-6)
    {
      if (ibegin == idx.size())
        ibegin = i;
      iend = i;
    }

  if (ibegin != iend)
    std::sort(idx.begin()+ibegin, idx.begin() + iend + 1,
              [&points, first](size_t i1, size_t i2)
              {
                return (points[first].distance(points[i1]) < points[first].distance(points[i2]));
              });


  const auto tmp = points;
  for (std::size_t i=0; i<n_points; ++i)
    points[i] = tmp[idx[i]];
}


template<typename Scalar>
void
Polygon<Scalar>::reorder_indices(const std::vector<Point<3, Scalar>> &verts,
                                 std::vector<std::size_t>            &indices)
{
  if (indices.size() < 3)
    throw std::invalid_argument("Polygon cannot have less than 3 vertices");
  std::vector<Point<3, Scalar>> points(indices.size());
  for (std::size_t i=0; i<indices.size(); ++i)
    points[i] = verts[indices[i]];

  reorder(points);

  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t ind = find(points[i], verts, 1e-6);
    indices[i] = ind;
  }
}

template<typename Scalar>
Scalar Polygon<Scalar>::area() const
{
  const Point<3, Scalar> cm = compute_center_mass(this->points);
  Scalar total_area = 0;
  for (std::size_t i=0; i<this->points.size(); ++i)
  {
    if (i < this->points.size() - 1)
    {
      const Point<3,Scalar> v1 = this->points[i];
      const Point<3,Scalar> v2 = this->points[i + 1];
      total_area += triangle_area(v1, v2, cm);
    }
    else
    {
      const Point<3,Scalar> v1 = this->points[i];
      const Point<3,Scalar> v2 = this->points[0];
      total_area += triangle_area(v1, v2, cm);
    }
  }

  return total_area;
}


template<typename Scalar>
std::vector<Edge> Polygon<Scalar>::get_edges() const
{
  std::vector<Edge> edges;
  for (std::size_t i=0; i<this->points.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < this->points.size() - 1)
    {
      i1 = i;
      i2 = i + 1;
    }
    else
    {
      i1 = i;
      i2 = 0;
    }

    edges.push_back({i1, i2});
  }

  return edges;
}


template<typename Scalar>
bool Polygon<Scalar>::point_inside(const Point<3, Scalar> & p ,
                                   const Scalar             tol) const
{
  if ( std::fabs(this->plane().signed_distance(p)) > tol )
    return false;

  const Point<3,Scalar> cm = this->center();
  for (const auto & edge : get_edges())
  {
    const Plane<Scalar> side = get_side(edge);
    if (std::fabs( side.signed_distance(p) ) > tol)
      if (side.above(p) != side.above(cm))
        return false;
  }

  return true;
}


template<typename Scalar>
Plane<Scalar> Polygon<Scalar>::get_side(const Edge & edge) const
{
  if (edge.first >= this->points.size() or edge.second >= this->points.size())
    throw std::out_of_range("Edge does not exist");

  const Point<3,Scalar> point3 = this->points[edge.first] +
                                 m_plane.normal() * (this->points[edge.first] -
                                                     this->points[edge.second]).norm();
  Plane<Scalar> side(this->points[edge.first], this->points[edge.second], point3);
  return side;
}


template<typename Scalar>
Point<3,Scalar> Polygon<Scalar>::center() const
{
  Point<3,Scalar> u, v, n, c;
  Scalar poly_area = 0;

  /* Break the poly into triangles with
   * vertices in points[0], points[j], points[j+1] */
  for (std::size_t j=1; j<this->points.size()-1; j++)
  {
    /* compute normal and offset w from first 3 vertices */
    // u and v are tangent vectors
    u = this->points[j] - this->points[0];
    v = this->points[j+1] - this->points[0];
    // normal vector components (not normalized)
    n = cross_product(u, v);
    const Scalar areatmp = 0.5 * n.norm();
    c += areatmp/3 * (this->points[0] + this->points[j] + this->points[j+1]);
    poly_area += areatmp;
  }

  c /= poly_area;
  return c;
}

}
