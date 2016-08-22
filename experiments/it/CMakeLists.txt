add_executable(it it.cpp Feature.cpp FeatureMatrix.cpp)
target_link_libraries(it Entropy++ boost_filesystem boost_system m)
