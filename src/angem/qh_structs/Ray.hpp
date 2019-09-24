#pragma once

#include "../Point.hpp"

namespace angem {
namespace quickhull {

template <typename T>
struct Ray {
const Point<3,T> m_S;  // start
const Point<3,T> m_V;  // direction
const T m_VInvLengthSquared;

Ray(const Point<3,T>& S,const Point<3,T>& V)
  : m_S(S), m_V(V),
    m_VInvLengthSquared(1.0 / ( m_V.norm() * m_V.norm() ))
{}


inline T distance(const Point<3,T> & p) const
{
  const Point<3,T> s = p - m_S;
  const T t = s.dot(m_V);
  const T s_norm = s.norm();
  return s_norm*s_norm - t * t * m_VInvLengthSquared;
}

};

}


}
