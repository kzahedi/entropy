add_library(eit Feature.cpp FeatureMatrix.cpp GIS.cpp)
target_link_libraries(eit Entropy++ boost_filesystem boost_system m)

add_executable(it it.cpp)
target_link_libraries(it Entropy++ boost_filesystem boost_system m eit)