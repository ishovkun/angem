# angem_root/src/CMakeLists.txt

# ADD_LIBRARY(angem SHARED
#   angem/QuickHull.cpp
# )

ADD_LIBRARY(angem INTERFACE)


SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)

# if(Boost_FOUND)
#   TARGET_LINK_LIBRARIES(angem ${Boost_LIBRARIES})
#   TARGET_INCLUDE_DIRECTORIES(angem PUBLIC ${Boost_INCLUDE_DIRS} )
# endif()

# if(USE_EIGEN)
#   # find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#   # TARGET_LINK_LIBRARIES(angem Eigen3::Eigen)
#   TARGET_INCLUDE_DIRECTORIES(angem PRIVATE ${EIGEN3_INCLUDE_DIR})
# endif()


install(TARGETS angem LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
