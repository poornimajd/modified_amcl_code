add_library(sensors_lib SHARED
  laser/laser.cpp
  laser/beam_model.cpp
  laser/likelihood_field_model.cpp
  laser/likelihood_field_model_prob.cpp
)
# map_update_cspace
target_link_libraries(sensors_lib pf_lib map_lib)
find_package(tf2_geometry_msgs REQUIRED)
ament_target_dependencies(sensors_lib
  tf2_geometry_msgs
)
install(TARGETS
  sensors_lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
