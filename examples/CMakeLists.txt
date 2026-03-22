add_subdirectory(pid)
add_subdirectory(lqr)
add_subdirectory(ekf)
add_subdirectory(ukf)
add_subdirectory(pf)

if (CTRLPP_BUILD_OSQP)
    add_subdirectory(mpc)
endif ()

if (CTRLPP_BUILD_NLOPT)
    add_subdirectory(nmpc)
endif ()
