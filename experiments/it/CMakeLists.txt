add_library(eit Feature.cpp ITMatrix.cpp FeatureMatrix.cpp IT.cpp GIS.cpp InstanceMatrix.cpp SCGIS.cpp Test.cpp)  
target_link_libraries(eit Entropy++ boost_filesystem boost_system m)

add_executable(it it_main.cpp)
target_link_libraries(it Entropy++ boost_filesystem boost_system m eit)
