# not required if this is compiled as part of entropy
project(mc)
cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0042 OLD)
set(CMAKE_VERBOSE_MAKEFILE       OFF)
set(CMAKE_BUILD_TYPE             "Debug")
set(CMAKE_BASE_NAME              "g++")
set(CMAKE_CXX_COMPILER           "g++")
set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
set(CMAKE_CXX_COMPILER_FULLPATH  "g++")
set(EXECUTABLE_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH          ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_CXX_FLAGS_DEBUG        "-g")
set(CMAKE_CXX_FLAGS_RELEASE      "-g")
set(CMAKE_INSTALL_PREFIX         "..")

add_executable(rbo_mc rbo_mc.cpp)
target_link_libraries(rbo_mc Entropy++ boost_filesystem boost_system m glog gflags pthread)

add_executable(dr dr.cpp)
target_link_libraries(dr Entropy++ boost_filesystem boost_system m glog gflags pthread)
 
# install(TARGETS hopping_mc RUNTIME DESTINATION bin)
