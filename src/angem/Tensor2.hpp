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
  // returns dim * dim
  static size_t size() { return dim*dim; }
  // create a unit tensor
  static Tensor2<dim,T> make_unit_tensor();
  // GETTERS
  inline T & operator()(const int i, const int j);
  inline const T & operator()(const int i, const int j) const;
  inline T & operator()(const std::size_t i, const std::size_t j);
  inline const T & operator()(const std::size_t i, const std::size_t j) const;

  // dot product with a vector
  Point<dim,T> operator*(const Point<dim,T>   & p) const;
  // inline component multiply by a scalar
  void operator*=(const T & x);
  // inline component divide by a scalar
  void operator/=(const T & x);

  //  FRIEND FUNCTIONS
  // left dot product
  template<int d, typename Scalar>
  friend Point<d,Scalar> operator*(const Point<d,Scalar>   & p,
                                   const Tensor2<d,Scalar> & t);
  // component-wise product with a scalar
  template<int d, typename Scalar>
  friend Tensor2<d,Scalar> operator*(const Scalar            & x,
                                     const Tensor2<d,Scalar> & t);
  // component-wise multiply by a scalar
  template<int d, typename Scalar>
  friend Tensor2<d,Scalar> operator*(const Tensor2<d,Scalar> & t,
                                     const Scalar            & x);
  // component-wise divide by a scalar
  template<int d, typename Scalar>
  friend Tensor2<d,Scalar> operator/(const Tensor2<d,Scalar> & t,
                                     const Scalar            & x);

  // dot product with another tensor
  template<int d, typename Scalar>
  friend Tensor2<d,Scalar> product(const Tensor2<d,Scalar> & t1,
                                   const Tensor2<d,Scalar> & t2);
  // dot product with a vector
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
Tensor2<dim,T>
Tensor2<dim,T>::make_unit_tensor()
{
  Tensor2<dim,T> ut;
  for (int i = 0; i < dim; i++)
    ut(i, i) = 1;
  return ut;
}



template <int dim, typename T>
T & Tensor2<dim,T>::operator()(const int i, const int j)
{
  assert(i >= 0); assert(j >= 0);
  assert(i < dim); assert(j < dim);
  return storage[i][j];
}


template <int dim, typename T>
T & Tensor2<dim,T>::operator()(const std::size_t i, const std::size_t j)
{
  assert(i < dim); assert(j < dim);
  return storage[i][j];
}


template <int dim, typename T>
const T & Tensor2<dim,T>::operator()(const int i, const int j) const
{
  assert(i >= 0); assert(j >= 0);
  assert(i < dim); assert(j < dim);
  return storage[i][j];
}


template <int dim, typename T>
const T & Tensor2<dim,T>::operator()(const std::size_t i, const std::size_t j) const
{
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
  return result;
}


template <int dim, typename T>
Point<dim,T> Tensor2<dim,T>::operator*(const Point<dim,T>   & p) const
{
  Point<dim,T> result;

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      result[i] += storage[i][j] * p(j);
  return result;
}


template <int dim, typename T>
void Tensor2<dim,T>::operator*=(const T & x)
{
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      storage[i][j] *= x;
}


template <int dim, typename T>
void Tensor2<dim,T>::operator/=(const T & x)
{
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      storage[i][j] /= x;
}

//    left dot product
template <int dim, typename T>
Point<dim,T> operator*(const Point<dim,T>   & p,
                       const Tensor2<dim,T> & t)
{
  Point<dim,T> result;

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      result[j] += p(i) * t.storage[i][j];
  return result;

}


template<int d, typename Scalar>
Tensor2<d,Scalar> operator*(const Tensor2<d,Scalar> & t,
                            const Scalar            & x)
{
  auto result = t;
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      result *= x;
  return result;
}

template<int d, typename Scalar>
Tensor2<d,Scalar> operator*(const Scalar            & x,
                            const Tensor2<d,Scalar> & t)
{
  return operator*(t, x);
}


template<int d, typename Scalar>
Tensor2<d,Scalar> operator/(const Tensor2<d,Scalar> & t,
                            const Scalar            & x)
{
  auto result = t;
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
      result /= x;
  return result;
}

}  // end namespace
