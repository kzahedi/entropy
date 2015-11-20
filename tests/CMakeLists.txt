include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_executable(container-tests        main.cpp  container_test.cpp)
target_link_libraries(container-tests Entropy++ cppunit)

add_executable(matrix-tests        main.cpp  matrix_test.cpp)
target_link_libraries(matrix-tests Entropy++ cppunit)

ADD_TEST(container "${PROJECT_BINARY_DIR}/bin/container-tests")
ADD_TEST(matrix    "${PROJECT_BINARY_DIR}/bin/matrix-tests")
