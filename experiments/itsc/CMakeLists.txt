add_executable(random_network random_network.cpp)
target_link_libraries(random_network Entropy++ ${Boost_LIBRARIES} gflags glog)
