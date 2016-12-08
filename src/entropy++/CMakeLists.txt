include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

set(SRC sparse/CMI.cpp

  CMI.cpp
  H.cpp
  MC_CW.cpp
  MC_MI.cpp
  MC_W.cpp
  MI.cpp
  PI.cpp
  ConditionalEntropy.cpp

  state/CMI.cpp
  state/MC_CW.cpp
  state/MC_MI.cpp
  state/MC_W.cpp
  state/MI.cpp

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

  iterativescaling/Model.cpp
  iterativescaling/GIS.cpp

  iterativescaling/IterativeScalingBase.cpp
  iterativescaling/Feature.cpp
  iterativescaling/FeatureMatrix.cpp
  # iterativescaling/FeatureMatrixsp.cpp
  iterativescaling/InstanceMatrix.cpp
  iterativescaling/ITMatrix.cpp
  iterativescaling/gis/IterativeScaling.cpp
  iterativescaling/gis/gp/IterativeScaling.cpp
  iterativescaling/scgis/IterativeScaling.cpp
  iterativescaling/scgis/gp/IterativeScaling.cpp

  SparseMatrix.cpp
  # Container.cpp
  Matrix.cpp

  Csv.cpp
  powi.cpp
)

if(ENTROPY_PROGRESS_BAR)
  set(SRC ${SRC} ext/progressbar.c)
endif(ENTROPY_PROGRESS_BAR)

if(APPLE)
  add_library(Entropy++ SHARED ${SRC})
  install(TARGETS Entropy++ LIBRARY DESTINATION lib)
else(APPLE)
  add_library(Entropy++ STATIC ${SRC})
endif(APPLE)

target_link_libraries(Entropy++ ${Boost_LIBRARIES})
if(ENTROPY_PROGRESS_BAR)
  target_link_libraries(Entropy++ ncurses)
endif(ENTROPY_PROGRESS_BAR)

