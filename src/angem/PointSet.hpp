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

class PointSet {
 public:
  PointSet(double tol = 1e-10,
            Point<3,long double> min = Point<3,long double>(std::numeric_limits<double>::max(),
                                                            std::numeric_limits<double>::max(),
                                                            std::numeric_limits<double>::max()))
            : _tol(tol), _min(min)
  {
    assert( _tol > 0 );
#ifndef NDEBUG
    _max = Point<3,long double>(std::numeric_limits<double>::lowest(),
                                std::numeric_limits<double>::lowest(),
                                std::numeric_limits<double>::lowest());
#endif
  }

  angem::Point<3,double> const & operator[](size_t i) const { return _storage[i]; }

  std::vector<angem::Point<3,double>> const & points() const { return _storage; }

  void clear() {
    _storage.clear();
    _mapping.clear();
  }

  // if found, returns the index of point in storage
  // else returns the size of the set
  size_t find(Point<3,double> const & p) const
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

  size_t insert(Point<3,double> const & p)
  {
    // need to adjust the limit somehow
    for (size_t i = 0; i < 3; ++i)
      if ( _min[i] > p[i] )
        {
          update_limit_(p);
          update_mapping_();
        }

    auto idx = cartesian_indices_(p);
    auto key = compute_key_(idx);
    size_t id = find_(key, p);
    if ( id == size() )  {
      _storage.push_back(p);
      _mapping[key] = id;
    }

    return id;
  }

  inline size_t size() const noexcept { return _storage.size(); }

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

  void update_limit_(Point<3,double> const &p)
  {
    for (size_t i = 0; i < 3; ++i) {
      if ( _min[i] > p[i]) {
        _min[i] = ( p[i] >= 0 ) ? 0.5*p[i] : 2*p[i];
      }
#ifndef NDEBUG
      _max[i] = std::max((double)_max[i], p[i]);
      if ((_max[i] - _min[i]) * std::numeric_limits<long double>::epsilon() > _tol) {
        throw std::runtime_error("Cannot resolve this set with such tolerance. " "Use higher precision");
      }
#endif
    }
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

  std::array<size_t,3> cartesian_indices_(Point<3,double> const & p) const
  {
    std::array<size_t,3> idx;
    for (size_t i = 0; i < 3; ++i)
      idx[i] = static_cast<size_t>( ((long double)p[i] - _min[i]) /(long double)_tol );
    return idx;
  }

  size_t find_(std::string const & key, Point<3,double> const & p) const noexcept
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

  Point<3,long double> _min;
  double _tol;
  std::vector<Point<3,double>> _storage;
  std::unordered_map<std::string,size_t> _mapping;
#ifndef NDEBUG
  Point<3,long double> _max;
#endif
};

}  // end namespace aaa
