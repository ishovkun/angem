#pragma once

#include "Point.hpp"
#include <algorithm>
#include <limits>  // std::numeric_limits

namespace angem
{

template<int dim, typename Scalar>
Point<dim,Scalar>
compute_center_mass(const std::vector<Point<dim,Scalar>> & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += p;
  center /= static_cast<Scalar>(points.size());
  return center;
}

template<int dim, typename Scalar, typename Iterable>
Point<dim,Scalar> compute_center_mass(const Iterable & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += p;
  center /= static_cast<Scalar>(points.size());
  return center;
}


template<int dim, typename Scalar>
Point<dim,Scalar>
compute_center_mass(const std::vector<Point<dim,Scalar> *> & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += *p;
  center /= static_cast<Scalar>(points.size());
  return center;
}


// this function is O(n²)
template<int dim, typename Scalar>
void
remove_duplicates_slow(const std::vector<Point<dim,Scalar>> & points,
                       std::vector<Point<dim,Scalar>>       & result,
                       const double                           tolerance = 0)
{
  // returns result vector that contains only unique entries of points vector
  // two points are considered duplicate if the distance between them is
  // less than tolerance

  if (result.empty())
    result.clear();

  for (const auto & p : points)
    if (find(p, result, tolerance) == result.size())
      result.push_back(p);
}


template<int dim, typename Scalar>
std::vector<Point<dim,Scalar>>
remove_duplicates_slow(const std::vector<Point<dim,Scalar>> & points,
                       const double                           tolerance = 0)
{
  // returns result vector that contains only unique entries of points vector
  // two points are considered duplicate if the distance between them is
  // less than tolerance
  std::vector<Point<dim,Scalar>> result;
  remove_duplicates_slow(points, result, tolerance);
  return result;
}


// remove duplicates (with tolerance) from the vector of points
// this function is O(n)
// template<int dim, typename Scalar>
// void remove_duplicates(std::vector<Point<dim,Scalar>> & points,
                       // const double tolerance = 1e-6)
// {
  // PointSet<dim,Scalar> pset;
  // for (const auto & p : points)
    // pset.insert(p);
  // points = std::move(pset.points);
// }


template<int dim, typename Scalar>
Scalar triangle_area(const Point<dim,Scalar> & p1,
                     const Point<dim,Scalar> & p2,
                     const Point<dim,Scalar> & p3)
{
  return static_cast<Scalar>(0.5) * (cross( p2 - p1, p3 - p1 )).norm();
}


// convert radians to degrees
template<typename Scalar>
inline
Scalar degrees(const Scalar angle)
{
  return angle / static_cast<Scalar>(M_PI) * static_cast<Scalar>(180.);
}


// convert degrees to degrees
template<typename Scalar>
inline
Scalar radians(const Scalar angle)
{
  return angle * static_cast<Scalar>(M_PI) / static_cast<Scalar>(180.);
}


// template<int dim, typename Scalar>
template<int dim, typename Scalar, typename Iterable>
std::size_t
find_closest_index(const Point<dim,Scalar> & point,
                   const std::vector<Point<dim,Scalar>> & points)
{
  Scalar min_dist = std::numeric_limits<Scalar>::max();
  std::size_t closest = 0;
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    const Scalar current = point.distance(points[i]);
    if (current < min_dist)
    {
      min_dist = current;
      closest = i;
    }
  }
  return closest;
}


template<int dim, typename Scalar, typename Iterable>
Point<dim,Scalar>
find_closest(const Point<dim,Scalar>              & point,
             const std::vector<Point<dim,Scalar>> & points)
{
  Scalar min_dist = std::numeric_limits<Scalar>::max();
  Point<dim,Scalar> closest;
  for (auto it = points.begin(); it != points.end(); ++it)
  {
    const Scalar current_dist = point.distance(*it);
    if (current_dist < min_dist)
    {
      min_dist = current_dist;
      closest = *it;
    }
  }
  return closest;
}


template<int dim, typename Scalar>
std::size_t
find_closest_vertex(const Point<dim,Scalar>              & point,
                    const std::vector<Point<dim,Scalar>> & all_vertices,
                    const std::vector<std::size_t>       & subset)
{
  Scalar min_dist = std::numeric_limits<Scalar>::max();
  std::size_t closest_index = 0;
  for (const std::size_t i : subset)
  {
    const Scalar current_dist = point.distance(all_vertices[i]);
    if (current_dist < min_dist)
    {
      min_dist = current_dist;
      closest_index = i;
    }
  }
  return closest_index;
}

/* Reorder array so that elements are placed in index location.
 * Value from arr[i] will be placed into arr[index[i]].
 * Input:
 * \param[in,out] arr : array to be reordered
 * \param[in,out] idx : mapping
 *
 * Caution: after reordering, the array of indices becomes range(0, n)
 * NOTE: reorder_to(reorder_from(arr, idx), idx) = arr
 * Complexity: A O(n) time and O(1) extra space
 *
 * Example:
 * Input: arr=[0 1 2 3 4 5], idx=[3 2 1 4 5 0]
 * Result: arr=[0 2 1 0 3 4], idx=[0 1 2 3 4 5]
*/
template<typename Scalar, typename IdxType>
void reorder_to(std::vector<Scalar> &arr, std::vector<IdxType> & index)
{
  assert( arr.size() == index.size() );
  for (IdxType i = 0; i < arr.size(); i++)
  {
    IdxType src = i;
    IdxType dst = index[i];
    Scalar pending = arr[src];

    while (src != dst) {
      Scalar const tmp = arr[dst];
      arr[dst] = pending;
      src = dst;
      dst = index[dst];
      pending = tmp;
      index[src] = src;
    }
  }
}

/* Reorder array so that elements are taken from index location.
 * Value from arr[index[i]] will be placed into arr[i].
 * Input:
 * \param[in,out] arr : array to be reordered
 * \param[in,out] idx : mapping
 *
 * Caution: after reordering, the array of indices becomes range(0, n)
 * NOTE: reorder_to(reorder_from(arr, idx), idx) = arr
 * Complexity: A O(n) time and O(1) extra space
 *
 * Example:
 * Input: arr=[0 1 2 3 4 5], idx=[3 2 1 4 5 0]
 * Result: arr=[3 2 1 4 5 0], idx=[0 1 2 3 4 5]
 *
*/
template<typename Scalar, typename IdxType>
void reorder_from(std::vector<Scalar> &arr, std::vector<IdxType> & index)
{
  assert( arr.size() == index.size() );
  for (int i = 0; i < arr.size(); i++)
  {
    Scalar tmp = arr[i];
    IdxType src = index[i];
    IdxType dst = i;
    while (src != i)
    {
      arr[dst] = arr[src];
      index[dst] = dst;
      dst = src;
      src = index[src];
    }
    arr[dst] = tmp;
    index[dst] = dst;
  }
}

// // Same as reorder_from
// template<typename Scalar, typename IdxType>
// void reorder(std::vector<Scalar> &arr, std::vector<IdxType> & index)
// {
//   reorder_from(arr, index);
// }

}
