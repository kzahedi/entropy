include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(container-tests        main.cpp  container_test.cpp)
target_link_libraries(container-tests Entropy++ cppunit)

add_executable(matrix-tests           main.cpp  matrix_test.cpp)
target_link_libraries(matrix-tests    Entropy++ cppunit)

add_executable(pi-tests               main.cpp  pi_test.cpp)
target_link_libraries(pi-tests        Entropy++ cppunit)

add_executable(it-tests               main.cpp  it_test.cpp)
target_link_libraries(it-tests        Entropy++ cppunit eit)

add_executable(mi-tests               main.cpp  mi_test.cpp)
target_link_libraries(mi-tests        Entropy++ cppunit)

add_executable(cmi-tests              main.cpp  cmi_test.cpp)
target_link_libraries(cmi-tests       Entropy++ cppunit)

add_definitions(-DTEST_CSV="${PROJECT_SOURCE_DIR}/tests/test.csv")
add_executable(csv-tests              main.cpp  csv_test.cpp)
target_link_libraries(csv-tests       Entropy++ cppunit)

ADD_TEST(container "${PROJECT_BINARY_DIR}/bin/container-tests")
ADD_TEST(matrix    "${PROJECT_BINARY_DIR}/bin/matrix-tests")
ADD_TEST(mi        "${PROJECT_BINARY_DIR}/bin/mi-tests")
ADD_TEST(it        "${PROJECT_BINARY_DIR}/bin/it-tests")
ADD_TEST(cmi       "${PROJECT_BINARY_DIR}/bin/cmi-tests")
ADD_TEST(csv       "${PROJECT_BINARY_DIR}/bin/csv-tests")
