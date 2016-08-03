include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

set(SRC sparse/CMI.cpp
  sparse/H.cpp
  sparse/MC_MI.cpp
  sparse/MI.cpp
  sparse/MC_W.cpp
  sparse/PI.cpp

  state/CMI.cpp
  state/MC_MI.cpp
  state/MC_W.cpp
  state/MI.cpp

  sparse/state/CMI.cpp
  sparse/state/MC_W.cpp
  sparse/state/MI.cpp
  sparse/state/MC_MI.cpp

  SparseMatrix.cpp
  # Container.cpp
  Entropy.cpp
  Matrix.cpp
  MI.cpp
  H.cpp
  PI.cpp
  CMI.cpp
  MC_W.cpp
  MC_MI.cpp

  Csv.cpp
)

if(APPLE)
	add_library(Entropy++ SHARED ${SRC})
	install(TARGETS Entropy++ LIBRARY DESTINATION lib)
else(APPLE)
	add_library(Entropy++ STATIC ${SRC})
endif(APPLE)

target_link_libraries(Entropy++ ${Boost_LIBRARIES})

