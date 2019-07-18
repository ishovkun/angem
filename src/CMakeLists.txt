# angem_root/src/CMakeLists.txt

# boost is optional for 256-bit integers (for mesh library)
# but it works muuuch faster with it
ADD_LIBRARY(angem SHARED
  angem/PointSet.cpp
)

SET_TARGET_PROPERTIES (
    angem
    PROPERTIES LINKER_LANGUAGE CXX
)

if(Boost_FOUND)
  TARGET_LINK_LIBRARIES(angem ${Boost_LIBRARIES})
endif()
