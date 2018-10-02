#pragma once

#include <Shape.hpp>
#include <Plane.hpp>
#include <utils.hpp>
#include <iostream>

namespace angem
{

template<typename Scalar>
class Polygon: public Shape<Scalar>
{
 public:
	Polygon();
  Polygon(const std::vector<Point<3,Scalar>> & points_list);
  // construct a polygon (face) from some mesh vertices
  Polygon(const std::vector<Point<3,Scalar>> & all_mesh_vertices,
          const std::vector<std::size_t>     & indices);

  // shift all points in direction p
  virtual void set_data(const std::vector<Point<3,Scalar>> & point_list) override;
  virtual void move(const Point<3,Scalar> & p) override;
  static  void reorder(std::vector<Point<3,Scalar>> & points);

  // Attributes
  Plane<Scalar> plane;
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
                         const std::vector<std::size_t>     & indices)
{
  assert(indices.size() > 2);
  std::vector<Point<3,Scalar>> point_list;
  for (const std::size_t & i : indices)
    point_list.push_back(all_mesh_vertices[i]);
  set_data(point_list);
}


template<typename Scalar>
void
Polygon<Scalar>::set_data(const std::vector<Point<3,Scalar>> & point_list)
{
  // TODO: i don't check whether all points aren't on one line
  assert(point_list.size() >= 3);

  // std::cout << "making poly" << std::endl;
  // for (const auto & p : point_list)
  //   std::cout << p << std::endl;

  this->points = point_list;
  reorder(this->points);
  Point<3, Scalar> cm = compute_center_mass(point_list);
  /* i'd like the plane support point (first argument) to be
     the center of the poly
     to create the plane I need to pass three points that are not
     aligned on a line this could happen if I do it like that:
     plane.set_data(cm, point_list[1], point_list[2]);
     Therefore, I need to selects points appropriately.
     only two polygon vertices can potentially be on the same line
     as the center of mass. So i need to do only one check
  */
  if ( ((point_list[0] - cm).cross(point_list[1] - cm)).norm() > 1e-16 )
    plane.set_data(cm, point_list[0], point_list[1]);
  else
    plane.set_data(cm, point_list[0], point_list[2]);

}


template<typename Scalar>
void
Polygon<Scalar>::move(const Point<3,Scalar> & p)
{
  Shape<Scalar>::move(p);
  plane.move(p);
}


template<typename Scalar>
void
Polygon<Scalar>::reorder(std::vector<Point<3, Scalar> > &points)
{
  const std::size_t n_points = points.size();
  assert(n_points > 2);
  if (n_points == 3)
    return;

  Plane<Scalar> plane(points[0], points[1], points[2]);
  Point<3,Scalar> normal = plane.normal();

  std::vector<Point<3, Scalar> > v_points;
  std::vector<Point<3,Scalar>> copy = points;
  v_points.push_back(copy.front());
  copy.erase(copy.begin());

  while (!copy.empty())
  {
    if (copy.size() == 1)
    {
      v_points.push_back(copy[0]);
      break;
    }
    // find such vertex that all other vertices are on one side of the edge
    for (std::size_t i=0; i<copy.size(); ++i)
    {
      // make plane object that we use to check on which side of the plane
      // any point is
      Scalar len = (copy[i] - v_points.back()).norm();
      Point<3, Scalar> p_perp = v_points.back() + normal * len;
      Plane<Scalar> pln(v_points.back(), p_perp, copy[i]);

      bool all_above = true;
      bool orientation;
      bool orientation_set = false;  // set after first assignment
      for (std::size_t j=0; j<points.size(); ++j)
      {
        if (points[j] == copy[i] || points[j] == v_points.back())
          continue;
        const bool above = pln.above(points[j]);
        if (!orientation_set)
        {
          orientation = above;
          orientation_set = true;
        }
        if (above != orientation)
        {
          all_above = false;
          break;
        }
      }
      if (all_above)
      {
        v_points.push_back(copy[i]);
        copy.erase(copy.begin() + i);
        break;
      }

    }
  }

  points = v_points;

}

}