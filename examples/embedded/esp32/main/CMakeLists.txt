set(CTRLPP_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../..")

idf_component_register(
    SRCS "app_main.cpp"
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart
)

# ctrlpp's public include tree plus the shared board-leg kernel. SYSTEM keeps the
# strict warning surface from tripping on third-party Eigen and the ctrlpp
# template depth.
target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE
    "${CTRLPP_ROOT}/lib/ctrlpp/include"
    "${CTRLPP_ROOT}/examples/embedded/shared")

# Eigen lives outside ESP-IDF's component graph and is not published as a
# managed component; the build expects a local checkout exposed via
# CTRLPP_EIGEN_DIR.
if (DEFINED ENV{CTRLPP_EIGEN_DIR})
    target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE "$ENV{CTRLPP_EIGEN_DIR}")
endif ()

# Belt-and-suspenders: the -fno-exceptions posture already drives
# CTRLPP_HAS_EXCEPTIONS to 0 via __cpp_exceptions, and defining this makes the
# no-throw intent explicit at the build boundary.
target_compile_definitions(${COMPONENT_LIB} PRIVATE CTRLPP_NO_EXCEPTIONS)

# C++20 baseline matches the ctrlpp project floor.
target_compile_features(${COMPONENT_LIB} PRIVATE cxx_std_20)
