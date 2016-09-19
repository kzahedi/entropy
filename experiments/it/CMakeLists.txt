add_library(eit Feature.cpp FeatureMatrix.cpp FeatureMatrixsp.cpp GIS.cpp GISsp.cpp InstanceMatrix.cpp SCGIS.cpp Comp.cpp GISgp.cpp SCGISgp.cpp)
target_link_libraries(eit Entropy++ boost_filesystem boost_system m)

add_executable(it it.cpp)
target_link_libraries(it Entropy++ boost_filesystem boost_system m eit)