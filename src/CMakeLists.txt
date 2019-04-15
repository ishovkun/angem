# angem/CMakeLists.txt

# boost is optional for 256-bit integers (for mesh library)
# but it works muuuch faster with it
find_package(Boost 1.65.1)
if(Boost_FOUND)
  message("Boost found. Using boost libraries for hashing.")
  message("To turn off edit the CMakeLists file.")
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
  add_definitions(-DUSE_BOOST)
else()
  # 256 bit integer lib (for mesh)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/uint256)
  ADD_SUBDIRECTORY(src/uint256)
ENDIF()

ADD_LIBRARY(angem
  Point.hpp
  PointSet.hpp
  PointSet.cpp
  Basis.hpp
  Plane.hpp
  Shape.hpp
  Polygon.hpp
  Polyhedron.hpp
  Tetrahedron.hpp
  Hexahedron.hpp
  Wedge.hpp
  Rectangle.hpp
  Collisions.hpp
  CollisionGJK.hpp
  Line.hpp
  PolyGroup.hpp
)

SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)

TARGET_INCLUDE_DIRECTORIES(angem PRIVATE
	${CMAKE_SOURCE_DIR}/src
)


if(Boost_FOUND)
  TARGET_LINK_LIBRARIES(angem ${Boost_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(angem uint256)
endif()
