project(entropy)

cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_COMPILER           "g++")
set(CMAKE_CXX_COMPILER_FULLPATH  "g++")

ENABLE_TESTING()

cmake_policy(SET CMP0042 OLD)

set(EXECUTABLE_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH          ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_VERBOSE_MAKEFILE       OFF)

option(ENTROPY_COMPILE_HOPPING "Compile MC example with hopping data" OFF)
option(ENTROPY_COMPILE_RBO     "Compile MC example for RBOHand data"  OFF)
option(ENTROPY_COMPILE_ITSC    "Compile Iterative Scaling examples"   OFF)
option(ENTROPY_TESTS           "Compile test cases"                   OFF)
option(ENTROPY_OPENMP          "Use OpenMP"                           OFF)
option(ENTROPY_COMPILE_SML     "Compile SML"                          OFF)


if(ENTROPY_OPENMP)
  add_definitions(-DUSE_OPENMP)
endif(ENTROPY_OPENMP)

set(CMAKE_BUILD_TYPE "Debug")

if(APPLE)
  set(CMAKE_CXX_COMPILER           "/usr/local/opt/llvm/bin/clang++")
  set(CMAKE_CXX_COMPILER_FULLPATH  "/usr/local/opt/llvm/bin/clang++")
  link_directories(/usr/local/opt/llvm/lib)
  include_directories(/usr/local/opt/llvm/include)
  if(ENTROPY_OPENMP)
    set(CMAKE_CXX_FLAGS_DEBUG        "-g  -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE      "-g  -fopenmp")
    set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -fopenmp")
  else(ENTROPY_OPENMP)
    set(CMAKE_CXX_FLAGS_DEBUG        "-g")
    set(CMAKE_CXX_FLAGS_RELEASE      "-g")
    set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
  endif(ENTROPY_OPENMP)
else(APPLE)
  if(ENTROPY_OPENMP)
    find_package(OpenMP)
  endif(ENTROPY_OPENMP)
  set(CMAKE_CXX_COMPILER           "g++")
  set(CMAKE_CXX_COMPILER_FULLPATH  "g++")
  set(CMAKE_CXX_FLAGS_DEBUG        "-g   ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE      "-g   ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3  ${OpenMP_CXX_FLAGS}")
endif(APPLE)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

if(UNIX AND NOT APPLE)
  include_directories(~/local/include)
  link_directories(~/local/lib)
endif(UNIX AND NOT APPLE)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED thread filesystem system unit_test_framework)

add_subdirectory(src/entropy++)

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

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  CONFIGURATIONS Release Debug
  PATTERN "*.cpp" EXCLUDE
  PATTERN "*.txt" EXCLUDE
  )

if(ENTROPY_COMPILE_HOPPING)
  add_subdirectory(experiments/hopping)
endif(ENTROPY_COMPILE_HOPPING)

if(ENTROPY_COMPILE_RBO)
  add_subdirectory(experiments/RBOHand)
endif(ENTROPY_COMPILE_RBO)

if(ENTROPY_COMPILE_ITSC)
  add_subdirectory(experiments/itsc)
endif(ENTROPY_COMPILE_ITSC)

if(ENTROPY_COMPILE_SML)
  add_subdirectory(experiments/sml)
endif(ENTROPY_COMPILE_SML)

if(ENTROPY_TESTS)
  add_subdirectory(tests)
endif(ENTROPY_TESTS)


find_package(Doxygen)
if(DOXYGEN_FOUND)
  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc)
    FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc)
  ENDIF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc)

  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc/html)
    FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html)
  ENDIF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc/html)

  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc/latex)
    FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/latex)
  ENDIF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/doc/latex)
  include(${CMAKE_SOURCE_DIR}/cmake/UseDoxygen.cmake)
endif(DOXYGEN_FOUND)

