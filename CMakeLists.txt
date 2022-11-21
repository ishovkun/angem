# angem_root/CMakeLists.txt
project(angem)

# angem
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src)
set(angem_include_dirs ${angem_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR}/src  PARENT_SCOPE)

