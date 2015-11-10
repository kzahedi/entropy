include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/lib)

add_library(Entropy++ SHARED
  Container.cpp
  MI.cpp
  PI.cpp
  Entropy.cpp
  )
# target_link_libraries(Entropy ${Boost_LIBRARIES} xerces-c glog YarsClientCom dl)

install(TARGETS Entropy++ LIBRARY DESTINATION lib)
