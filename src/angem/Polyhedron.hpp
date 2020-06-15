#pragma once

#include "Shape.hpp"  // provides Shape<Scalar>
#include "Plane.hpp"  // provides Plane<Scalar>
#include "PointSet.hpp"  // provides PointSet<dim,Scalar>
#include "Polygon.hpp"   // provides Polygon<Scalar>
#include "PolyGroup.hpp" // provides Polygroup
#include "VTK_ID.hpp"    // provides VTK_ID
#include <typeinfo>
#include <exception>
#include <map>

namespace angem
{

using Edge = std::pair<std::size_t, std::size_t>;

template<typename Scalar>
class Polyhedron: public Shape<Scalar>
{
 public:
  // Constructors
  Polyhedron(const int vtk_id = -1);
  Polyhedron(const std::vector<Point<3,Scalar>>          & vertices,
             const std::vector<std::vector<std::size_t>> & faces,
             const int                                     vtk_id = -1);
  // setter
  void set_data(const std::vector<Point<3,Scalar>>          & vertices,
                const std::vector<std::vector<std::size_t>> & faces);
  // getters
  // get vtk id of the polyhedron
  int id() const {return vtk_id;}
  // compute the volume of the polyhedron
  virtual Scalar volume() const;
  // compute the center of mass of the polyhedron
  virtual Point<3,Scalar> center() const override;
  // check whether the point is inside the polyhedron
  bool point_inside(const Point<3,Scalar> & p) const;
  // check whether the point resides on one of the boundaries
  // using the provided tolerance
  bool point_on_boundary(const Point<3,Scalar> & p,
                         const double tolerance = 1e-4) const;

  // get vector of vectors that contain the vertex indices of each polyhedron face
  const std::vector<std::vector<std::size_t>> & get_faces() const;
  std::vector<std::vector<std::size_t>> & get_faces();
  std::vector<Polygon<Scalar>> get_face_polygons() const;

  // return vector of ordered pairs of vertex indices
  // (ordering by index comparison)
  std::vector<Edge> get_edges() const;

 protected:
  std::vector<std::vector<std::size_t>> m_faces;
  const int vtk_id;
};


template<typename Scalar>
Polyhedron<Scalar>::Polyhedron(const int vtk_id)
    :
    vtk_id(vtk_id)
{}


template<typename Scalar>
Polyhedron<Scalar>::Polyhedron(const std::vector<Point<3,Scalar>>          & vertices,
                               const std::vector<std::vector<std::size_t>> & faces,
                               const int vtk_id)
    : vtk_id(vtk_id)
{
  set_data(vertices, faces);
}


template<typename Scalar>
void
Polyhedron<Scalar>::set_data(const std::vector<Point<3,Scalar>>          & vertices,
                             const std::vector<std::vector<std::size_t>> & faces)
{
  assert(vertices.size() > 3);
  m_faces.resize(faces.size());

  std::map<size_t,size_t> global_to_local;  // vertex indices, retains order
  for (const auto & face : faces)
    for (const size_t v : face)
      global_to_local.insert( {v , 0} );
  size_t n_local = 0;                        // number of local vertices
  for (auto & it : global_to_local)
    it.second = n_local++;

  for (size_t iface=0; iface<faces.size(); ++iface)
  {
    const auto & face = faces[iface];
    m_faces[iface].reserve(face.size());
    for(const size_t vert_global : face)
    {
      const size_t local = global_to_local[vert_global];
      m_faces[iface].push_back(local);
    }
  }

  this->points.reserve(global_to_local.size());
  // invert the map
  std::vector<size_t> inv(global_to_local.size());
  for (const auto & it : global_to_local)
    inv[it.second] = it.first;

  for (const size_t vglob: inv )
    this->points.push_back( vertices[ vglob ] );
}


template<typename Scalar>
std::vector<std::vector<std::size_t>> &
Polyhedron<Scalar>::get_faces()
{
  return m_faces;
}


template<typename Scalar>
const std::vector<std::vector<std::size_t>> &
Polyhedron<Scalar>::get_faces() const
{
  return m_faces;
}


template<typename Scalar>
Scalar Polyhedron<Scalar>::volume() const
{
  const Point<3,Scalar> c = Shape<Scalar>::center();
  Scalar vol = 0;
  for (const auto & face_indices : get_faces())
  {
    const auto face_poly = Polygon<Scalar>(this->points, face_indices);
    const Scalar face_area = face_poly.area();
    const Scalar h = std::fabs(face_poly.plane().signed_distance(c));
    vol += 1./3. * h * face_area;
  }
  return vol;
}

template<typename Scalar>
bool Polyhedron<Scalar>::point_inside(const Point<3,Scalar> & p) const
{
  const Point<3,Scalar> c = this->center();
  for (const auto & face : m_faces)
  {
    Plane<Scalar> plane(this->points, face);
    if (plane.above(p) != plane.above(c))
      return false;
  }
  return true;
}

template<typename Scalar>
bool Polyhedron<Scalar>::point_on_boundary(const Point<3,Scalar> & p,
                                           const double tolerance) const
{
  for (const auto & face : m_faces)
  {
    Plane<Scalar> plane(this->points[face[0]], this->points[face[1]], this->points[face[2]]);
    if (plane.distance(p) < tolerance)
      return true;
  }
  return false;
}


template<typename Scalar>
std::ostream &operator<<(std::ostream             & os,
                         const Polyhedron<Scalar> & poly)
{
  const auto & points = poly.get_points();
  const auto & faces = poly.get_faces();
  os << points.size() << " vertices:" << std::endl;
  os << points;
  os << faces.size() << " faces:" << std::endl;
  for (const auto & face : faces)
  {
    for (const auto ivertex: face)
      os << ivertex << "\t";
    os << std::endl;
  }
  return os;
}

template<typename Scalar>
std::vector<Polygon<Scalar>> Polyhedron<Scalar>::get_face_polygons() const
{
  std::vector<Polygon<Scalar>> polys;
  for (const auto & face_indices : get_faces())
  {
    polys.emplace_back(this->points, face_indices,
                       /*reorder_vertices = */ false);
  }

  return polys;
}


template<typename Scalar>
Point<3,Scalar> Polyhedron<Scalar>::center() const
{
  /* Polyhedron center of mass is not the same as the
   * mass center of vertices!!!
   * Algorithm: split polyhedron into pyramids (we know their volume),
   * we also know their center of mass so we're good. */

  const auto faces = get_face_polygons();
  Point<3,Scalar> face_center_mass;
  double total_face_area = 0;
  for (const auto & face : faces)
  {
    Scalar face_area = face.area();
    face_center_mass += face.center() * face_area;
    total_face_area += face_area;
  }

  face_center_mass /= ( total_face_area );

  Point<3,Scalar> c;  // polyhedron center mass
  Scalar vol = 0;
  for (const auto & face : faces)
  {
    const auto face_center = face.center();
    const Scalar h = face.plane().normal().dot( face_center - face_center_mass );
    const Scalar volumetmp = fabs(h * face.area()) / 3.;

    if (std::isnan(volumetmp))
      throw std::runtime_error("invalid polyhedron volume");

    c += (face_center + 0.25*(face_center_mass - face_center))*volumetmp;
    vol += volumetmp;
  }

  c /= vol;
  return c;
}


template<typename Scalar>
std::vector<Edge> Polyhedron<Scalar>::get_edges() const
{
  std::vector<Edge> edges;
  for (const std::vector<size_t> & face : m_faces)
    for (std::size_t i = 0; i < face.size(); ++i)
    {
      std::size_t i1, i2;
      if (i < face.size() - 1)
      {
        i1 = face[i];
        i2 = face[i+1];
      }
      else
      {
        i1 = face[i];
        i2 = face[0];
      }
      auto edge = std::minmax(i1, i2);
      if (std::find_if( edges.begin(), edges.end(),
                     [&edge](const Edge & it)->bool
                     {
                       return it.first == edge.first &&
                           it.second == edge.second;
                     }) == edges.end())
        edges.push_back( std::move(edge) );
    }
  return edges;
}


}  // end namespace
