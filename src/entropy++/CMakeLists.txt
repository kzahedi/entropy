include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

set(SRC

  CMI.cpp
  H.cpp
  MC_CW.cpp
  MC_MI.cpp
  MC_W.cpp
  MI.cpp
  PI.cpp
  ConditionalEntropy.cpp
  SparseMatrix.cpp
  Matrix.cpp
  Csv.cpp
  powi.cpp
  Random.cpp

  state/CMI.cpp
  state/MC_CW.cpp
  state/MC_MI.cpp
  state/MC_W.cpp
  state/MI.cpp

  sparse/CMI.cpp
  sparse/ConditionalEntropy.cpp
  sparse/H.cpp
  sparse/MC_CW.cpp
  sparse/MC_MI.cpp
  sparse/MC_W.cpp
  sparse/MI.cpp
  sparse/PI.cpp

  sparse/state/ConditionalEntropy.cpp
  sparse/state/CMI.cpp
  sparse/state/MC_CW.cpp
  sparse/state/MC_MI.cpp
  sparse/state/MC_W.cpp
  sparse/state/MI.cpp

  iterativescaling/Delta.cpp
  iterativescaling/DeltaMatcher.cpp
  iterativescaling/Feature.cpp
  iterativescaling/Model.cpp
  iterativescaling/GIS.cpp
  iterativescaling/SCGIS.cpp
  iterativescaling/KL.cpp
  iterativescaling/Original.cpp

)

if(APPLE)
  add_library(Entropy++ SHARED ${SRC})
  install(TARGETS Entropy++ LIBRARY DESTINATION lib)
else(APPLE)
  add_library(Entropy++ STATIC ${SRC})
endif(APPLE)

target_link_libraries(Entropy++ ${Boost_LIBRARIES})
