add_executable(lqr_closed_loop_settling lqr_closed_loop_settling.cpp)
target_link_libraries(lqr_closed_loop_settling PRIVATE ctrlpp::ctrlpp)
# Reads the board legs' own kernel so the case asks the question they ask; a restated plant could drift from it.
target_include_directories(lqr_closed_loop_settling PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/embedded/shared")
