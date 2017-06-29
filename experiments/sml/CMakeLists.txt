add_executable(binary binary.cpp)
target_link_libraries(binary Entropy++ glog gflags ${Boost_LIBRARIES})

add_executable(non-binary non-binary.cpp)
target_link_libraries(non-binary Entropy++ glog gflags ${Boost_LIBRARIES})
