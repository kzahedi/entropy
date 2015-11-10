project(entropy)

cmake_minimum_required(VERSION 3.3)

ENABLE_TESTING()

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
set(CMAKE_CXX_FLAGS_RELEASE      "-arch x86_64")

include_directories(${CMAKE_SOURCE_DIR}/src)

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

