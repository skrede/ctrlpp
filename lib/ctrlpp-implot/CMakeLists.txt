if (CTRLPP_CMAKE_FETCH_DEPS)
    FetchContent_Declare(
        imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG        v1.92.6
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(imgui)

    FetchContent_Declare(
        implot
        GIT_REPOSITORY https://github.com/epezent/implot.git
        GIT_TAG        v0.17
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(implot)

    FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG        3.4
        GIT_SHALLOW    TRUE
    )
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_WAYLAND  OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(glfw)

    find_package(OpenGL REQUIRED)

    add_library(ctrlpp_implot_deps STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
        ${implot_SOURCE_DIR}/implot.cpp
        ${implot_SOURCE_DIR}/implot_items.cpp
    )

    target_include_directories(ctrlpp_implot_deps SYSTEM PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
        ${implot_SOURCE_DIR}
    )

    target_link_libraries(ctrlpp_implot_deps PUBLIC glfw OpenGL::GL)

    target_compile_options(ctrlpp_implot_deps PRIVATE
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-w>
        $<$<CXX_COMPILER_ID:MSVC>:/w>
    )

    set(_CTRLPP_IMPLOT_LINK_DEPS ctrlpp_implot_deps)
else ()
    find_package(imgui CONFIG REQUIRED)
    find_package(implot CONFIG REQUIRED)
    find_package(glfw3 CONFIG REQUIRED)
    find_package(OpenGL REQUIRED)

    foreach (_target imgui::imgui implot::implot glfw)
        get_target_property(_aliased ${_target} ALIASED_TARGET)
        if (_aliased)
            set(_real ${_aliased})
        else ()
            set(_real ${_target})
        endif ()
        get_target_property(_inc ${_real} INTERFACE_INCLUDE_DIRECTORIES)
        if (_inc)
            set_target_properties(${_real} PROPERTIES
                INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_inc}"
            )
        endif ()
    endforeach ()

    set(_CTRLPP_IMPLOT_LINK_DEPS imgui::imgui implot::implot glfw OpenGL::GL)
endif ()

add_library(ctrlpp_implot STATIC
    src/app.cpp
    src/sim_harness.cpp
    src/signal_recorder.cpp
    src/scrolling_plot.cpp
    src/static_plot.cpp
)
add_library(ctrlpp::implot ALIAS ctrlpp_implot)

target_include_directories(ctrlpp_implot
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ctrlpp_implot PUBLIC ctrlpp::ctrlpp ${_CTRLPP_IMPLOT_LINK_DEPS})

target_compile_options(ctrlpp_implot PRIVATE ${CTRLPP_WARNING_FLAGS})

install(TARGETS ctrlpp_implot EXPORT ctrlppImplotTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ctrlpp/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ctrlpp
)

install(EXPORT ctrlppImplotTargets
    FILE ctrlppImplotTargets.cmake
    NAMESPACE ctrlpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ctrlpp
)
