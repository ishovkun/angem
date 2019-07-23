#pragma once

#include "Point.hpp"

#include <array>
#include <vector>

namespace angem
{

/* Implementation of 2nd order tensor */
template <int dim, typename T>
class Tensor2
{
 public:
  // default constructor
  Tensor2();
  // create point from std::vector
  Tensor2(const std::vector<T> & v);

  // GETTERS
  inline T & operator()(const int i, const int j);
  inline const T & operator()(const int i, const int j) const;

  // fuctions
  template<int d, typename Scalar>
  friend Tensor2<d,Scalar> product(const Tensor2<d,Scalar> & t1,
                                   const Tensor2<d,Scalar> & t2);
  template<int d, typename Scalar>
  friend Point<d,Scalar> product(const Tensor2<d,Scalar> & t,
                                 const Point<d,Scalar>   & p);

 private:
  // storage 2d array (2nd order tensor)
  std::array<std::array<T, dim>, dim> storage;
};


template <int dim, typename T>
Tensor2<dim,T>::Tensor2()
{
  for (auto & row : storage)
    for (auto & item : row)
      item = 0.0;
}


template <int dim, typename T>
Tensor2<dim,T>::Tensor2(const std::vector<T> & v)
{
  assert(v.size() == dim*dim);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      storage[i][j] = v[dim*i + j];
}


template <int dim, typename T>
T & Tensor2<dim,T>::operator()(const int i, const int j)
{
  assert(i > 0); assert(j > 0);
  assert(i < dim); assert(j < dim);
  return storage[i][j];
}


template <int dim, typename T>
const T & Tensor2<dim,T>::operator()(const int i, const int j) const
{
  assert(i > 0); assert(j > 0);
  assert(i < dim); assert(j < dim);
  return storage[i][j];
}


// fuctions
template<int d, typename Scalar>
Tensor2<d,Scalar> product(const Tensor2<d,Scalar> & t1,
                          const Tensor2<d,Scalar> & t2)
{
  Tensor2<d,Scalar> result;
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      for (int k = 0; k < d; k++)
        result(i, k) += t1(i, j) * t2(j, k);
  return result;
}


template<int d, typename Scalar>
Point<d,Scalar> product(const Tensor2<d,Scalar> & t,
                        const Point<d,Scalar>   & p)
{
  Point<d,Scalar> result;

  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      result[i] += t(i, j) * p(j);
}

}  // end namespace
