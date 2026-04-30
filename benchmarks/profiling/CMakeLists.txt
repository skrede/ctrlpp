option(CTRLPP_PROFILE_NMPC_PENDULUM "Build pendulum-NMPC profiling harnesses" OFF)

if(CTRLPP_PROFILE_NMPC_PENDULUM)
    add_subdirectory(nmpc_pendulum)
endif()
