set(CTRLPP_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../..")

idf_component_register(
    SRCS "main.cpp"
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart esp_timer
)

# ctrlpp's public include tree. SYSTEM keeps the strict warning surface from
# tripping on third-party Eigen and the ctrlpp template depth.
target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE
    "${CTRLPP_ROOT}/lib/ctrlpp/include")

# Eigen lives outside ESP-IDF's default component graph; it is pulled through
# the IDF Component Manager (see ./idf_component.yml). If a local Eigen checkout
# is preferred, expose it via CTRLPP_EIGEN_DIR.
if (DEFINED ENV{CTRLPP_EIGEN_DIR})
    target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE "$ENV{CTRLPP_EIGEN_DIR}")
endif ()

# Belt-and-suspenders: the -fno-exceptions posture already drives
# CTRLPP_HAS_EXCEPTIONS to 0 via __cpp_exceptions, and defining this makes the
# no-throw intent explicit at the build boundary.
target_compile_definitions(${COMPONENT_LIB} PRIVATE CTRLPP_NO_EXCEPTIONS)

# C++20 baseline matches the ctrlpp project floor.
target_compile_features(${COMPONENT_LIB} PRIVATE cxx_std_20)
