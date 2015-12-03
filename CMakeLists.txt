project(entropy)

cmake_minimum_required(VERSION 3.3)

ENABLE_TESTING()

cmake_policy(SET CMP0042 OLD)

set(EXECUTABLE_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH          ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_VERBOSE_MAKEFILE       OFF)
set(CMAKE_CXX_COMPILER           "g++")
set(CMAKE_CXX_COMPILER_FULLPATH  "g++")

set(CMAKE_CXX_FLAGS_DEBUG        "-g -arch x86_64")
set(CMAKE_CXX_FLAGS_RELEASE      "-g -arch x86_64")
set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")

include_directories(${CMAKE_SOURCE_DIR}/src)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED filesystem system)

add_subdirectory(src/entropy++)
add_subdirectory(tests)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/src/entropy++/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/entropy++
  CONFIGURATIONS Release Debug
  PATTERN "*.cpp" EXCLUDE 
  PATTERN "*.txt" EXCLUDE 
  )

option(COMPILE_HOPPING "Compile MC example with hopping data" OFF)

if(COMPILE_HOPPING)
  add_subdirectory(examples/hopping)
endif(COMPILE_HOPPING)
       
