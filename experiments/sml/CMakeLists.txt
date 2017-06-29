add_executable(binary binary.cpp)
target_link_libraries(binary Entropy++ glog gflags ${Boost_LIBRARIES})

add_executable(non-binary non-binary.cpp)
target_link_libraries(non-binary Entropy++ glog gflags ${Boost_LIBRARIES})

INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/
  CONFIGURATIONS Release Debug
  PATTERN "*.cpp" EXCLUDE
  PATTERN "*.txt" EXCLUDE
  )

