#pragma once
#include <list>

#include "Transform.hpp"

namespace angem {

template <typename T>
class TransformSeries
{
 public:
  TransformSeries() = default;
  TransformSeries(std::initializer_list<Transform<T>> transforms);
  angem::Point<3,T> apply(const angem::Point<3,T> & p) const;

 private:
  std::list<Transform<T>> _transforms;
};

template <typename T>
TransformSeries<T>::TransformSeries(std::initializer_list<Transform<T>> transforms)
    : _transforms(transforms)
{}

template <typename T>
angem::Point<3,T> TransformSeries<T>::apply(const angem::Point<3,T> & p) const
{
  angem::Point<3,T> p_new = p;
  for (const auto & transform : _transforms)
    p_new = transform.apply(p_new);
  return p_new;
}

}  // end namespace angem
