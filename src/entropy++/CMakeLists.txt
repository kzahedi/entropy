include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_library(Entropy++ SHARED
  Container.cpp
  Entropy.cpp
  Matrix.cpp
  SparseMatrix.cpp
  MI.cpp
  MIsd.cpp
  MIs.cpp
  PI.cpp
  PIs.cpp
  CMI.cpp
  CMIs.cpp
  MC_W.cpp
  MC_Ws.cpp
  MC_MI.cpp
  MC_MIsd.cpp
  MC_MIs.cpp
  CsvToContainer.cpp
  )
target_link_libraries(Entropy++ ${Boost_LIBRARIES})

install(TARGETS Entropy++ LIBRARY DESTINATION lib)
