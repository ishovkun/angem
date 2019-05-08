# angem/CMakeLists.txt

# boost is optional for 256-bit integers (for mesh library)
# but it works muuuch faster with it
ADD_LIBRARY(angem
  angem/Point.hpp
  angem/PointSet.hpp
  angem/PointSet.cpp
  angem/Basis.hpp
  angem/Plane.hpp
  angem/Shape.hpp
  angem/Polygon.hpp
  angem/Polyhedron.hpp
  angem/Tetrahedron.hpp
  angem/Hexahedron.hpp
  angem/Wedge.hpp
  angem/Rectangle.hpp
  angem/Collisions.hpp
  angem/CollisionGJK.hpp
  angem/Line.hpp
  angem/PolyGroup.hpp
)

SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(angem PRIVATE
	${CMAKE_SOURCE_DIR}/src
	${CMAKE_SOURCE_DIR}/src/angem
)


if(Boost_FOUND)
  TARGET_LINK_LIBRARIES(angem ${Boost_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(angem uint256)
endif()
