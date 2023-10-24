/*
 * This code was shamelessly stolen from
 * https://github.com/ruddradev/gjk_cplusplus
 * and adapted for the local infrastructure
 * Author: Igor Shovkun <ishovkun@stanford.edu>
 */

#pragma once

#include "angem/Shape.hpp"
#include "angem/Point.hpp"

#include <vector>

namespace angem
{

/* This class implements GJK algorithm for fast collision detection */
template<typename Scalar>
class CollisionGJK
{
public:
  // default empty constructor
	CollisionGJK();
  // returns true if andgem shapes collide
	bool check(const Shape<Scalar> & shape1,
             const Shape<Scalar> & shape2);

 private:
  // get farthest point in the direction
	Point<3,Scalar> support(Point<3,Scalar> & direction);
	Point<3,Scalar> update_direction();

  // attributes
	const Shape<Scalar> * a;
	const Shape<Scalar> * b;
	std::vector<Point<3,Scalar>> simplex;
	Point<3,Scalar> origin;
};


template<typename Scalar>
CollisionGJK<Scalar>::CollisionGJK()
{
	simplex.reserve(4);
}

template<typename Scalar>
Point<3,Scalar>
CollisionGJK<Scalar>::support(Point<3,Scalar> & dir)
{
	Point<3,Scalar> support_point;
	const Point<3,Scalar> negDir = origin-dir;
	const Point<3,Scalar> aDot = a->support(dir);
	const Point<3,Scalar> bDot = b->support(negDir);
	support_point = aDot - bDot;
	return support_point;
}

template<typename Scalar>
Point<3,Scalar>
CollisionGJK<Scalar>::update_direction()
{
  Point<3,Scalar> dir,AB,AC,AD,AO,tPlane;
  std::vector<Point<3,Scalar>> temp;
  temp.reserve(3);
  const size_t n = simplex.size();

  switch(n)
  {
    case 2:
      AB = simplex.at(0) - simplex.at(1);
      AO = origin - simplex.at(1);
      // dir = (AB.cross(AO)).cross(AB);
      // dir = (cross(AB,AO)).cross(AB);
      dir = cross(cross(AB,AO), AB);
      break;

    case 3:
      AB = simplex.at(1) - simplex.at(2);
      AC = simplex.at(0) - simplex.at(2);
      AO = origin - simplex.at(2);
      // tPlane = AC.cross(AB);
      tPlane = cross(AC, AB);
      if ( dot(cross(AC, tPlane), AO) > 0 )
      {
        // dir = AC.cross(AO.cross(AC));
        dir = cross(AC, cross(AO, AC));
        simplex.erase(simplex.begin() + 1);
      }
      else if( dot( cross(tPlane, AB), AO ) <= 0 )
      {
        if (tPlane.dot(AO) > 0)
          dir = tPlane;
        else if(tPlane.dot(AO) == static_cast<Scalar>(0))
          break;
        else
          dir = origin - tPlane;
      }
      else
      {
        dir = cross(AB, cross(AO, AB) );
        simplex.erase(simplex.begin());
      }
      break;

    case 4:
      const Point<3,Scalar> A = simplex.at(3);
      const Point<3,Scalar> B = simplex.at(2);
      const Point<3,Scalar> C = simplex.at(1);
      const Point<3,Scalar> D = simplex.at(0);
      const Point<3,Scalar> AD = D - A;
      const Point<3,Scalar> AC = C - A;
      const Point<3,Scalar> AB = B - A;
      const Point<3,Scalar> normADC = cross(AD, AC);
      const Point<3,Scalar> normACB = cross(AC, AB);
      const Point<3,Scalar> normABD = cross(AB, AD);
      const Point<3,Scalar> AO = origin - A;

      if (normADC.dot(AO) > 0)
      {
        dir = normADC;
        simplex.erase(simplex.begin() + 2);
      }
      else if (normACB.dot(AO) > 0)
      {
        dir = normACB;
        simplex.erase(simplex.begin());
      }
      else if (normABD.dot(AO) > 0)
      {
        dir = normABD;
        simplex.erase(simplex.begin() + 1);
      }
      else break;
  }
  return dir;
}

template<typename Scalar>
bool
CollisionGJK<Scalar>::check(const Shape<Scalar> & shape1,
                            const Shape<Scalar> & shape2)
{
  a = & shape1;
  b = & shape2;
  simplex.clear();

  Point<3,Scalar> d(1, -1, -1);
  Point<3,Scalar> s = support(d);
  simplex.push_back(s);
  d = origin - d;

  int iter = 0;
  while(iter < 50)
  {
    iter++;
    s = support(d);
    if(d.dot(s) < static_cast<Scalar>(0))
    {
      return false;
      break;
    }

    simplex.push_back(s);
    d = update_direction();
    if(d.x() == 0 && d.y() == 0 && d.z() == 0)
      return true;
  }

  return false;
}

}  // end namespace
