#pragma once

#include <Point.hpp>
#include <Polygon.hpp>
#include <PointSet.hpp>
// #include <PolyGroup.hpp>

#include <unordered_set>


namespace angem
{

using Edge = std::pair<std::size_t, std::size_t>;

template <typename Scalar>
class SurfaceMesh // : PolyGroup<Scalar>
{
 public:
  SurfaceMesh(const double tol = 1e-6,
              const std::size_t max_edges = 1e8);
  // add new polygon (element) to the mesh. no check for duplicates
  void insert(const Polygon<Scalar> & poly);
  // true if set vertices is not empty
  bool empty() const;

  // GETTERS
  // get vector of neighbor indices
  std::vector<std::size_t> get_neighbors( const std::size_t ielement ) const;
  // vector of indices of cells neighboring an edge
  std::vector<std::size_t> get_neighbors( const Edge & edge ) const;
  // get vector of index pairs representing edges
  std::vector<Edge> get_edges( const std::size_t ielement ) const;
  //
  void delete_element(const std::size_t ielement);


  // merge polygon with its largest neighbor
  // complexity linear to number of faces (edges)
  // does not remove vertices when no connected neighbors
  std::size_t merge_element(const std::size_t element);

  PointSet<3,Scalar>                    vertices;
  std::vector<std::vector<std::size_t>> polygons;  // indices
  // std::vector<std::size_t> neighbors;

  // hash of two vert indices -> vector polygons
  // essentially edge -> neighbors
  std::unordered_map<std::size_t, std::vector<std::size_t>> map_edge_neighbors;

 private:
  // merge two elements haaving a comon edge
  // the edge is vertex indices
  // private cause no check if they are neighbors
  // jelement is the element to be merged into ielement
  // does not remove vertices when no connected neighbors
  void merge_elements(const std::size_t ielement,
                      const std::size_t jelement,
                      const Edge      & edge);

  std::size_t hash_value(const std::size_t ind1,
                         const std::size_t ind2) const
  {
    if (ind1 < ind2)
      return max_edges*ind1 + ind2;
    else
      return max_edges*ind2 + ind1;
  }


  Scalar tol;
  std::size_t max_edges;  // used to hash a int-int pair
};


template <typename Scalar>
SurfaceMesh<Scalar>::SurfaceMesh(const double tol,
                                 const std::size_t max_edges)
    :
    tol(tol),
    max_edges(max_edges)
{
  assert(tol > 0);
}


template <typename Scalar>
void SurfaceMesh<Scalar>::insert(const Polygon<Scalar> & poly)
{
  // this part is a copy of PolyGroup add method
  // append vertices, compute vertex indices, and append these indices
  std::vector<std::size_t> indices;
  const std::vector<Point<3,Scalar>> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }
  polygons.push_back(indices);

  // now we need to compute edges
  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < points.size() - 1)
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[i+1]);
    }
    else
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[0]);
    }

    if (i1 == i2)
      throw std::runtime_error("what kind of polygon is that?");

    const std::size_t hash = hash_value(i1, i2);
    auto iter = map_edge_neighbors.find(hash);
    if (iter != map_edge_neighbors.end())
      (iter->second).push_back(polygons.size() - 1);
    else
      map_edge_neighbors.insert({ {hash, {polygons.size() - 1}} });
  }

}


template <typename Scalar>
std::vector<std::size_t>
SurfaceMesh<Scalar>::get_neighbors( const std::size_t ielement ) const
{
  std::vector<std::size_t> v_neighbors;
  std::vector<Edge> edges = get_edges(ielement);
  for (const Edge & edge : edges)
  {
    for (const auto & ineighbor : get_neighbors(edge))
    {
      if (ineighbor == ielement)
        continue;
      else
        v_neighbors.push_back(ineighbor);
    }
  }
  return std::move(v_neighbors);
}


template <typename Scalar>
std::size_t SurfaceMesh<Scalar>::merge_element( const std::size_t ielement )
{

  // const std::vector<std::size_t> neighbors = get_neighbors(element);
  const std::vector<Edge> edges = get_edges(ielement);

  std::size_t max_iedge = 0;
  std::size_t max_neighbor = 0;
  std::size_t counter = 0;
  Scalar max_area = 0;
  for (const auto & edge : edges)
  {
    const std::vector<std::size_t> edge_neighbors = get_neighbors(edge);
    for (std::size_t element : edge_neighbors)
      if (element != ielement)
      {
        const Polygon<Scalar> poly(vertices.points, polygons[element]);
        const Scalar area = poly.area();
        if (area > max_area)
        {
          max_area = area;
          max_neighbor = element;
          max_iedge = counter;
        }
      }
    counter++;
  }

  merge_elements(max_neighbor, ielement, edges[max_iedge]);
  return max_neighbor;
}


template <typename Scalar>
void SurfaceMesh<Scalar>::delete_element(const std::size_t element)
{
  // delete polygon and reduce all polygons
  // in map_edge_neighbors with higher indices by 1
  polygons.erase(polygons.begin() + element);
  for (auto iter = map_edge_neighbors.begin();
       iter != map_edge_neighbors.end(); ++iter)
  {
    std::size_t counter = 0;
    for (std::size_t i=0; i<iter->second.size(); ++i)
    {
      std::size_t ielement = iter->second[i];
      if (ielement == element)
        iter->second.erase(iter->second.begin() + i);
      else if (ielement < element)
        iter->second[i]--;
    }
  }

}


template <typename Scalar>
void SurfaceMesh<Scalar>::merge_elements(const std::size_t ielement,
                                         const std::size_t jelement,
                                         const Edge      & edge)
{
  std::cout << "doing merge" << std::endl;
  // create set of vertex indices of element vertices
  std::unordered_set<std::size_t> set_vertices;
  for (auto & ivertex : polygons[ielement])
    set_vertices.insert(ivertex);
  for (auto & jvertex : polygons[jelement])
    set_vertices.insert(jvertex);
  // convert it into a vector
  std::vector<std::size_t> v_vertices;
  for (const auto & vertex : set_vertices)
    v_vertices.push_back(vertex);
  std::cout << "reorder" << std::endl;
  std::cout << "v_vertices.size() = " << v_vertices.size() << std::endl;

  Polygon<Scalar>::reorder_indices(vertices.points, v_vertices);
  std::cout << "copy" << std::endl;
  polygons[ielement] = v_vertices;
  std::cout << "done merge" << std::endl;
  delete_element(jelement);
  std::cout << "deleted" << std::endl;
}


template <typename Scalar>
std::vector<Edge>
SurfaceMesh<Scalar>::get_edges( const std::size_t ielement ) const
{
  std::vector<Edge> edges;
  const Polygon<Scalar> poly(vertices.points, polygons[ielement]);
  const std::vector<Point<3,Scalar>> &points = poly.get_points();
  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < points.size() - 1)
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[i+1]);
    }
    else
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[0]);
    }

    edges.push_back({i1, i2});

  }
  return std::move(edges);
}


template <typename Scalar>
std::vector<std::size_t>
SurfaceMesh<Scalar>::get_neighbors( const Edge & edge ) const
{
  std::vector<std::size_t> v_neighbors;
  const std::size_t hash = hash_value(edge.first, edge.second);
  const auto iter = map_edge_neighbors.find(hash);
  if (iter == map_edge_neighbors.end())
    throw std::out_of_range("edge does not exist");

  return iter->second;
}

}
