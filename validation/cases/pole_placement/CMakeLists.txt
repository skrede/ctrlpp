add_executable(pole_placement pole_placement.cpp)
target_link_libraries(pole_placement PRIVATE ctrlpp::ctrlpp)
