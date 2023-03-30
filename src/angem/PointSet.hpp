#pragma once
#include "angem/Point.hpp"
#include <vector>
#include <unordered_map>
#include <string>
#include <array>
#include <limits>

namespace angem {

constexpr size_t num_digits(size_t x)
{
  size_t ans = 0;
  while (x > 9) {
    ans++;
    x = x / 10;
  }
  return ans + 1;
}

/// @brief Implements a set of points that are indistinguishable within a specified tolerance
/// @tparam storage_type floating-point presicion of stored points
/// @tparam box_type floating-point presicion used in the cell calculation (implementation detail)
template<typename storage_type = double, typename box_type = long double>
class PointSet {
 public:
  /* Constructor
   * The user can specify the bounds of the point set (cordners of the bounding box).
   * Properly chosen bounds minimize the number of recomputes of underlying cartesian grid.
  */
  PointSet(storage_type tol = 1e-10,
           Point<3,box_type> min = Point<3,box_type>(std::numeric_limits<storage_type>::max(),
                                                     std::numeric_limits<storage_type>::max(),
                                                     std::numeric_limits<storage_type>::max()),
           Point<3,box_type> max = Point<3,box_type>(std::numeric_limits<storage_type>::lowest(),
                                                     std::numeric_limits<storage_type>::lowest(),
                                                     std::numeric_limits<storage_type>::lowest()))
              : _tol(tol), _min(min), _max(max)
  {
    assert( _tol > 0 );
  }

  // get element by index in storage vector
  angem::Point<3,storage_type> const & operator[](size_t i) const { return _storage[i]; }
  
  // get reference to the whole storage vector
  std::vector<angem::Point<3,storage_type>> const & points() const { return _storage; }

  // clear the set (both the mapping and storage).
  void clear() {
    _storage.clear();
    _mapping.clear();
  }

  // if found, returns the index of point in storage
  // else returns the size of the set
  size_t find(Point<3,storage_type> const & p) const
  {
    for (size_t i = 0; i < 3; ++i)
      if ( _min[i] > p[i] )
        return false;

    auto idx = cartesian_indices_(p);
    auto key = compute_key_(idx);
    size_t id = find_(key, p);
    if ( id < size() ) return id;

    // lookup neighbor cells
    for (size_t i = 0; i < 3; ++i)
    {
      idx[i] = idx[i] + 1;
      key = compute_key_(idx);
      id = find_(key, p);
      if ( id < size() ) return id;

      idx[i] = idx[i] - 2;
      key = compute_key_(idx);
      id = find_(key, p);
      if ( id < size() ) return id;

      idx[i]++;
    }

    return size();
  }

  // try to insert a point in the set.
  // if a point within tol from p already exists in the set, return its storage index.
  // otherwise, return the size of the set.
  size_t insert(Point<3,storage_type> const & p)
  {
    // need to adjust the limit somehow
    if ( update_limit_(p) )
      update_mapping_();

    auto idx = cartesian_indices_(p);
    auto key = compute_key_(idx);
    size_t id = find_(key, p);
    if ( id == size() )  {
      _storage.push_back(p);
      _mapping[key] = id;
    }

    return id;
  }

  // return s the number of points in the set
  inline size_t size() const noexcept { return _storage.size(); }

  // true if size() == 0, false otherwise
  inline bool empty() const noexcept { return _storage.empty(); }

 private:

  void update_mapping_()
  {
    _mapping.clear();
    for (size_t i = 0; i < _storage.size(); ++i) {
      auto idx = cartesian_indices_( _storage[i] );
      auto key = compute_key_(idx);
      assert( find_(key, _storage[i]) == size() );
      _mapping[key] = i;
    }
  }

  bool update_limit_(Point<3,storage_type> const &p)
  {
    bool updated = false;
    for (size_t i = 0; i < 3; ++i) {
      _max[i] = std::max(_max[i], static_cast<box_type>(p[i]));

      if ( _min[i] > p[i]) {
        _min[i] = p[i] - 0.5*(_max[i] - p[i]);
        updated = true;
      }
#ifndef NDEBUG
      static constexpr box_type eps = std::numeric_limits<box_type>::epsilon();
      storage_type const diff = _max[i] - _min[i];
      if (diff * eps > _tol) {
        std::cout << "max - min = "<< diff << std::endl;
        std::cout << "eps = " << eps << "\n";
        std::cout << "tol = " << _tol << "\n";
        std::cout << "eps*(max-min) = " << eps*diff << std::endl;
        throw std::runtime_error("Cannot resolve this set with such tolerance. " "Use higher precision");
      }
#endif
    }
    return updated;
  }

  std::string compute_key_(std::array<size_t,3> const & idx) const noexcept
  {
    std::string s;
    for (size_t i = 0; i < 3; ++i) {
      s += std::to_string(idx[i]);
      if ( i < 2 ) s += ",";
    }
    return s;
  }

  std::array<size_t,3> cartesian_indices_(Point<3,storage_type> const & p) const
  {
    std::array<size_t,3> idx;
    for (size_t i = 0; i < 3; ++i)
      idx[i] = static_cast<size_t>( ((box_type)p[i] - _min[i]) / _tol );
    return idx;
  }

  size_t find_(std::string const & key, Point<3,storage_type> const & p) const noexcept
  {
    auto it = _mapping.find(key);
    if (it != _mapping.end())
    {
      auto const & candidate = _storage[it->second];
      if ( candidate.distance(p) < _tol )
        return it->second;
      else return size();
    }
    else return size();
  }

  Point<3,box_type> _min, _max;
  box_type _tol;
  std::vector<Point<3,double>> _storage;
  std::unordered_map<std::string,size_t> _mapping;
};

}  // end namespace aaa
