cmake_minimum_required (VERSION 3.7)

# set(CMAKE_CXX_STANDARD 17)

# for ycmd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

# boost is optional for 256-bit integers (for mesh library)
# but it works muuuch faster with it
find_package(Boost 1.65.1)
if(Boost_FOUND)
  message("Boost found. Using boost libraries for hashing.")
  message("To turn off edit the CMakeLists file.")
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
  add_definitions(-DUSE_BOOST)
# else()
#   # 256 bit integer lib (for mesh)
#   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/external/uint256)
#   ADD_SUBDIRECTORY(external/uint256)
ENDIF()

# angem
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
ADD_SUBDIRECTORY(src)



# TARGET_INCLUDE_DIRECTORIES(msh2gprs PRIVATE
# 	${CMAKE_SOURCE_DIR}/src
# )

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
    # target_compile_options(msh2gprs PUBLIC "-Wall -Wextra -pedantic")
endif()

# memory leak debugging
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

# for profiling memory leaks
# TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers mesh uint256 -lstdc++fs -lasan)
# last flag for std::filesystem
# TARGET_LINK_LIBRARIES(msh2gprs gprs_data parsers mesh -lstdc++fs)
