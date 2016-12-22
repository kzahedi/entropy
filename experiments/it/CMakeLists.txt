# not required if this is compiled as part of entropy
project(iterativescaling)
cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0042 OLD)
set(CMAKE_VERBOSE_MAKEFILE       OFF)
set(CMAKE_BUILD_TYPE             "Debug")
set(CMAKE_BASE_NAME              "g++")
set(CMAKE_CXX_COMPILER           "g++")
set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
set(CMAKE_CXX_COMPILER_FULLPATH  "g++")
set(EXECUTABLE_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH          ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_CXX_FLAGS_DEBUG        "-g -arch x86_64")
set(CMAKE_CXX_FLAGS_RELEASE      "-g -arch x86_64")
set(CMAKE_INSTALL_PREFIX         "..")

add_executable(random_network random_network.cpp)
target_link_libraries(random_network Entropy++ boost_filesystem boost_system m glog gflags)

# install(TARGETS hopping_mc RUNTIME DESTINATION bin)

