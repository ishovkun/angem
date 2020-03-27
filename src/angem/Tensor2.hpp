#pragma once

#include "Point.hpp"
#include "Exceptions.hpp"

#include <array>
#include <vector>
#include <initializer_list>

namespace angem
{

/* Implementation of 2nd order tensor */
template <int dim, typename T>
class Tensor2
{
 public:
  // default constructor. Initializes a tensor filled with zeros.
  Tensor2();
  // construct a tensor from an initializer list
  Tensor2(std::initializer_list<T> l);
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

  // matrix (dot) product with another tensor2
  Tensor2<dim,T> operator*(const Tensor2<dim,T> & other) const;

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

  // printout
  template <int d, typename S>
  friend std::ostream &operator<<(std::ostream     & os,
                                  const Point<d,S> & p);

  // invert tensor for dim = 1
  template<typename S>
  friend Tensor2<1,S> invert(const Tensor2<1,S> & tens);

  // invert tensor for dim = 2
  template<typename S>
  friend Tensor2<2,S> invert(const Tensor2<2,S> & tens);

  // determinant of jacobian 3x3
  // template<typename S>
  // friend Tensor2<3,S> det(const Tensor2<3,S> & tens);

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
Tensor2<dim,T>::Tensor2(std::initializer_list<T> l)
{
  assert( l.size() == dim*dim );
  auto it = l.begin();
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
    {
      storage[i][j] = *it; ++it;
    }

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


template <int dim, typename T>
Tensor2<dim,T> Tensor2<dim,T>::operator*(const Tensor2<dim,T> & other) const
{
  Tensor2<dim,T> result;
  for (size_t i=0; i<dim; ++i)
    for (size_t j=0; j<dim; ++j)
      for (size_t k=0; k<dim; ++k)
        result(i, j) += this->operator()(i, k)  * other(k, j);
  return result;
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

template <typename T>
Tensor2<1,T> invert(const Tensor2<1,T> & tens)
{
  Tensor2<1,T> result;
  result(0, 0) = static_cast<T>(1.0) / tens(0,0);
  return result;
}

template <typename T>
Tensor2<2,T> invert(const Tensor2<2,T> & tens)
{
  Tensor2<2,T> result;
  throw NotImplemented("inverting 2x2 tensor not implemented");
  return result;
}

template <typename T>
inline T det(const Tensor2<3,T> & a)
{
  return a(0,0) * (a(1,1) * a(2,2) - a(1,2) * a(2,1)) -
         a(0,1) * (a(1,0) * a(2,2) - a(1,2) * a(2,0)) +
         a(0,2) * (a(1,0) * a(2,1) - a(1,1) * a(2,0));
}

template <typename T>
Tensor2<3,T> invert(const Tensor2<3,T> & a)
{
  Tensor2<3,T> result;
  const T deta = det<T>(a);
  assert( !std::isnan(1.0/deta));
  result(0,0) =   ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / deta;
  result(0,1) = - ( a(0,1) * a(2,2) - a(0,2) * a(2,1) ) / deta;
  result(0,2) =   ( a(0,1) * a(1,2) - a(0,2) * a(1,1) ) / deta;
  result(1,0) = - ( a(1,0) * a(2,2) - a(1,2) * a(2,0) ) / deta;
  result(1,1) =   ( a(0,0) * a(2,2) - a(0,2) * a(2,0) ) / deta;
  result(1,2) = - ( a(0,0) * a(1,2) - a(0,2) * a(1,0) ) / deta;
  result(2,0) =   ( a(1,0) * a(2,1) - a(1,1) * a(2,0) ) / deta;
  result(2,1) = - ( a(0,0) * a(2,1) - a(0,1) * a(2,0) ) / deta;
  result(2,2) =   ( a(0,0) * a(1,1) - a(0,1) * a(1,0) ) / deta;
  return result;
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

// Printing
template<int dim, typename Scalar>
std::ostream &operator<<(std::ostream            & os,
                         const Tensor2<dim,Scalar> & tensor)
{
  os << "(";
  for (int i=0; i<dim; ++i)
  {
    for (int j=0; j<dim; ++j)
    {
      os << tensor(i, j);
      if (j != dim - 1)
        os << " ";
    }

    if (i != dim - 1)
      os << "; ";
  }
  os << ")";
  return os;
}


}  // end namespace
