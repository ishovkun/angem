# angem_root/CMakeLists.txt

# boost is optional for 256-bit integers (for PointSet)
# but it works muuuch faster with it
if(NOT Boost_FOUND)
  find_package(Boost 1.65.1)
  if(Boost_FOUND)
    message(STATUS "Boost found. Using boost libraries for hashing.")
    message("To turn off edit the CMakeLists file.")
  ENDIF()
ENDIF()

if(Boost_FOUND)
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
  add_definitions(-DUSE_BOOST)
ENDIF()

# angem
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src)
set(angem_include_dirs ${angem_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR}/src  PARENT_SCOPE)

