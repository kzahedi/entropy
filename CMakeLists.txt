project(entropy)

cmake_minimum_required(VERSION 3.0)

ENABLE_TESTING()

cmake_policy(SET CMP0042 OLD)

set(EXECUTABLE_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH          ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_VERBOSE_MAKEFILE       OFF)
#set(CMAKE_CXX_COMPILER           "g++")
#set(CMAKE_CXX_COMPILER_FULLPATH  "g++")

if(APPLE)
set(CMAKE_CXX_FLAGS_DEBUG        "-g")
set(CMAKE_CXX_FLAGS_RELEASE      "-g")
set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
else(APPLE)
set(CMAKE_CXX_FLAGS_DEBUG        "-g  -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE      "-g  -fPIC")
set(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3")
endif(APPLE)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

if(UNIX AND NOT APPLE)
  include_directories(~/local/include)
  link_directories(~/local/lib)
endif(UNIX AND NOT APPLE)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED filesystem system)

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

option(ENTROPY_COMPILE_HOPPING "Compile MC example with hopping data" OFF)
option(ENTROPY_COMPILE_RBO     "Compile MC example for RBOHand data"  OFF)
option(ENTROPY_COMPILE_IT      "Compile IT example"                   OFF)
option(ENTROPY_TESTS           "Compile test cases"                   OFF)
option(ENTROPY_PROGRESS_BAR    "Use progress bar"                     OFF)
option(ENTROPY_MEMORY_EFFICIENT "Compile with memory efficiency. Code will be very slow if this is set to ON." ON)

if(ENTROPY_MEMORY_EFFICIENT)
  add_definitions(-DMEMORY_EFFICIENT)
endif(ENTROPY_MEMORY_EFFICIENT)

if(ENTROPY_PROGRESS_BAR)
  add_definitions(-DUSE_PROGRESS_BAR)
endif(ENTROPY_PROGRESS_BAR)

if(ENTROPY_COMPILE_HOPPING)
  add_subdirectory(experiments/hopping)
endif(ENTROPY_COMPILE_HOPPING)

if(ENTROPY_COMPILE_RBO)
  add_subdirectory(experiments/RBOHand)
endif(ENTROPY_COMPILE_RBO)

if(ENTROPY_COMPILE_IT)
  add_subdirectory(experiments/it)
endif(ENTROPY_COMPILE_IT)

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

