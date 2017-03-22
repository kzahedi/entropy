include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(container-tests        container_test.cpp)
target_link_libraries(container-tests Entropy++ ${Boost_LIBRARIES})
ADD_TEST(container "${PROJECT_BINARY_DIR}/bin/container-tests")

add_executable(matrix-tests           matrix_test.cpp)
target_link_libraries(matrix-tests    Entropy++ ${Boost_LIBRARIES})
ADD_TEST(matrix "${PROJECT_BINARY_DIR}/bin/matrix-tests")

add_executable(pi-tests               pi_test.cpp)
target_link_libraries(pi-tests        Entropy++ ${Boost_LIBRARIES})
ADD_TEST(PI "${PROJECT_BINARY_DIR}/bin/pi-tests")

add_executable(mi-tests               mi_test.cpp)
target_link_libraries(mi-tests        Entropy++ ${Boost_LIBRARIES})
ADD_TEST(MI "${PROJECT_BINARY_DIR}/bin/mi-tests")

add_executable(delta-tests            delta_test.cpp)
target_link_libraries(delta-tests     Entropy++ ${Boost_LIBRARIES})
ADD_TEST(Delta "${PROJECT_BINARY_DIR}/bin/delta-tests")

add_definitions(-DTEST_CSV="${PROJECT_SOURCE_DIR}/tests/test.csv")
add_executable(csv-tests              csv_test.cpp)
target_link_libraries(csv-tests       Entropy++ ${Boost_LIBRARIES})
ADD_TEST(csv   "${PROJECT_BINARY_DIR}/bin/csv-tests")

add_executable(entropy-tests          entropy_test.cpp)
target_link_libraries(entropy-tests   Entropy++ ${Boost_LIBRARIES})
ADD_TEST(entropy "${PROJECT_BINARY_DIR}/bin/entropy-tests")

add_executable(cmi-tests              cmi_test.cpp)
target_link_libraries(cmi-tests       Entropy++ ${Boost_LIBRARIES})
ADD_TEST(cmi "${PROJECT_BINARY_DIR}/bin/cmi-tests")


add_executable(original-tests            original_test.cpp)
target_link_libraries(original-tests     Entropy++ cppunit dl)
ADD_TEST(Original       "${PROJECT_BINARY_DIR}/bin/original-tests")

add_executable(gis-tests        gis_test.cpp)
target_link_libraries(gis-tests Entropy++ ${Boost_LIBRARIES})
ADD_TEST(GIS "${PROJECT_BINARY_DIR}/bin/gis-tests")

add_executable(scgis-tests        scgis_test.cpp)
target_link_libraries(scgis-tests Entropy++ ${Boost_LIBRARIES})
ADD_TEST(SCGIS "${PROJECT_BINARY_DIR}/bin/scgis-tests")
