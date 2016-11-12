include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(container-tests        main.cpp  container_test.cpp)
target_link_libraries(container-tests Entropy++ cppunit dl)

add_executable(matrix-tests           main.cpp  matrix_test.cpp)
target_link_libraries(matrix-tests    Entropy++ cppunit dl)

add_executable(pi-tests               main.cpp  pi_test.cpp)
target_link_libraries(pi-tests        Entropy++ cppunit dl)

add_executable(it-tests               main.cpp  it_test.cpp)
target_link_libraries(it-tests        Entropy++ cppunit dl eit)

add_executable(iterative-scaling-tests        main.cpp  iterative_scaling_test.cpp)
target_link_libraries(iterative-scaling-tests Entropy++ cppunit dl eit)

add_definitions(-DPARENT="${CMAKE_SOURCE_DIR}/experiments/hopping/data/")
add_executable(it-cmi-tests           main.cpp  it_cmi_test.cpp)
target_link_libraries(it-cmi-tests    Entropy++ cppunit dl eit)

add_executable(mi-tests               main.cpp  mi_test.cpp)
target_link_libraries(mi-tests        Entropy++ cppunit dl)

add_executable(entropy-tests          main.cpp  entropy_test.cpp)
target_link_libraries(entropy-tests   Entropy++ cppunit dl)

add_executable(cmi-tests              main.cpp  cmi_test.cpp)
target_link_libraries(cmi-tests       Entropy++ cppunit dl)

add_executable(feature-tests          main.cpp  feature_test.cpp)
target_link_libraries(feature-tests   Entropy++ cppunit dl eit)

add_definitions(-DTEST_CSV="${PROJECT_SOURCE_DIR}/tests/test.csv")
add_executable(csv-tests              main.cpp  csv_test.cpp)
target_link_libraries(csv-tests       Entropy++ cppunit dl)

ADD_TEST(container "${PROJECT_BINARY_DIR}/bin/container-tests")
ADD_TEST(matrix    "${PROJECT_BINARY_DIR}/bin/matrix-tests")
ADD_TEST(mi        "${PROJECT_BINARY_DIR}/bin/mi-tests")
ADD_TEST(entropy   "${PROJECT_BINARY_DIR}/bin/entropy-tests")
ADD_TEST(feature   "${PROJECT_BINARY_DIR}/bin/feature-tests")
ADD_TEST(it        "${PROJECT_BINARY_DIR}/bin/it-tests")
ADD_TEST(iterativescaling "${PROJECT_BINARY_DIR}/bin/iterative-scaling-test")
ADD_TEST(itcmi     "${PROJECT_BINARY_DIR}/bin/it-cmi-tests")
ADD_TEST(cmi       "${PROJECT_BINARY_DIR}/bin/cmi-tests")
ADD_TEST(csv       "${PROJECT_BINARY_DIR}/bin/csv-tests")
