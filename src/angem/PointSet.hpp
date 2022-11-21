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
  PointSet(double tol = 1e-6)
      : _tol(tol)
      , _upper(_tol*_dmax / 2)
      , _lower(- _upper)
  {
    assert( tol > 0 );
    assert( (_upper - _lower) / _tol > _dmax );
    // std::cout << "_lower = " << _lower << std::endl;
  }

  PointSet & operator=(PointSet const & other)
  {
    _storage = other._storage;
    _mapping = other._mapping;
    _tol = other._tol;
    _upper = other._upper;
    _lower = other._lower;
    return std::ref(*this);
  }

  PointSet & operator=(PointSet & other)
  {
    _storage.swap(other._storage);
    _mapping.swap(other._mapping);
    _tol = other._tol;
    _upper = other._upper;
    _lower = other._lower;
    return std::ref(*this);
  }

  inline Point<3,double> & operator[](size_t i) { return _storage[i]; }
  inline Point<3,double> const & operator[](size_t i) const { return _storage[i]; }
  inline std::vector<Point<3,double>> & points() noexcept { return _storage; }
  inline std::vector<Point<3,double>> const & points() const noexcept { return _storage; }


  size_t insert(Point<3,double> const & p)
  {
    size_t idx = find(p);
    if ( idx == this->size() )
    {
      _storage.push_back(p);
      auto key = compute_key(compute_idx_(p));
      _mapping[key] = idx;
    }
    return idx;
  }

  size_t find(Point<3,double> const & p) const noexcept
  {
    auto idx = compute_idx_(p);
    auto key = compute_key(idx);
    size_t id = find_(key, p);
    if ( id < size() ) return id;

    for (size_t i = 0; i < 3; ++i)
    {
      idx[i] = idx[i] + 1;
      key = compute_key(idx);
      id = find_(key, p);
      if ( id < size() ) return id;

      idx[i] = idx[i] - 2;
      key = compute_key(idx);
      id = find_(key, p);
      if ( id < size() ) return id;

      idx[i]++;
    }

    return size();
  }

  std::string compute_key(std::array<size_t,3> const & idx) const noexcept
  {
    std::string ans(_hash_len, '0');
    // std::cout << "idx = ";
    // for (size_t i = 0; i < 3; ++i)
    //   std::cout << idx[i] << " ";
    // std::cout << std::endl;

    for (size_t i = 0; i < 3; ++i)
    {
      size_t offset = _dlength*(i+1) - 1;
      size_t x = idx[i];
      while(x > 0) {
        char c = '0' + (x % 10);
        x /= 10;
        ans[offset] = c;
        offset--;
      }
    }
    return ans;
  }

  size_t find_by_hash(std::string const & key);

  size_t size() const noexcept {return _storage.size();}
  bool empty() const noexcept { return _storage.empty(); }

 private:
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

  std::array<size_t,3> compute_idx_(Point<3,double> const & p) const
  {
    std::array<size_t,3> idx;
    for (size_t i = 0; i < 3; ++i)
      idx[i] = static_cast<size_t>( ((long double)(p[i]) - _lower) /_tol );
    return idx;
  }



  std::vector< angem::Point<3,double> > _storage;
  std::unordered_map<std::string, size_t> _mapping;
  long double _tol;

  static constexpr size_t _dmax {std::numeric_limits<size_t>::max()};
  static constexpr size_t _dlength {num_digits(_dmax)};
  static constexpr size_t _hash_len{3*_dlength};
  double _upper, _lower;
};

}  // end namespace aaa
