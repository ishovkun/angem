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
