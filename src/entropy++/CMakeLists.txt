include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_library(Entropy++ SHARED
  Container.cpp
  Entropy.cpp
  Matrix.cpp
  SparseMatrix.cpp
  MI.cpp
  MIsd.cpp
  MIssd.cpp
  MIs.cpp
  H.cpp
  Hs.cpp
  PI.cpp
  PIs.cpp
  CMI.cpp
  CMIsd.cpp
  CMIs.cpp
  CMIssd.cpp
  MC_W.cpp
  MC_Wsd.cpp
  MC_Ws.cpp
  MC_Wssd.cpp
  MC_MI.cpp
  MC_MIsd.cpp
  MC_MIssd.cpp
  MC_MIs.cpp
  Csv.cpp
  )
target_link_libraries(Entropy++ ${Boost_LIBRARIES})

install(TARGETS Entropy++ LIBRARY DESTINATION lib)
