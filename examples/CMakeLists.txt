add_subdirectory(pid)
add_subdirectory(lqr)

if (CTRLPP_BUILD_OSQP)
    add_subdirectory(mpc)
endif ()

if (CTRLPP_BUILD_NLOPT)
    add_subdirectory(nmpc)
endif ()
