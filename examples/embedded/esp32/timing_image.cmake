# The timing image: a second executable linked from the same component libraries
# as the project's own image, with the main component's copy built without the
# allocation sensor in its place. The framework builds one executable per
# project, so this repeats what its project() does for that one: every
# component, whole-archive where the component asks for it, the framework's link
# options, and an elf2image pass that writes the ELF's digest into the image's
# descriptor, so the image has its own fingerprint.

function(ctrlpp_esp32_link_components elf main_lib)
    idf_build_get_property(components BUILD_COMPONENT_ALIASES)
    list(REMOVE_ITEM components idf::main)
    foreach (component ${components})
        idf_component_get_property(whole ${component} WHOLE_ARCHIVE)
        if (whole)
            target_link_libraries(${elf} PRIVATE "-Wl,--whole-archive" ${component} "-Wl,--no-whole-archive")
        else ()
            target_link_libraries(${elf} PRIVATE ${component})
        endif ()
    endforeach ()
    target_link_libraries(${elf} PRIVATE "-Wl,--whole-archive" ${main_lib} "-Wl,--no-whole-archive")
endfunction()

# The project's executable carries every option and linker-script dependency the
# framework version adds, so they are copied rather than restated; only the map
# file is this image's own. Linking after that executable guarantees the
# generated section-placement script both depend on exists.
function(ctrlpp_esp32_copy_link_setup elf map)
    set(project_elf ${CMAKE_PROJECT_NAME}.elf)
    get_target_property(link_options ${project_elf} LINK_OPTIONS)
    get_target_property(link_depends ${project_elf} LINK_DEPENDS)
    list(FILTER link_options EXCLUDE REGEX "^-Wl,--Map=")
    target_link_options(${elf} PRIVATE ${link_options} "-Wl,--Map=${map}")
    set_property(TARGET ${elf} APPEND PROPERTY LINK_DEPENDS "${link_depends}")
    add_dependencies(${elf} ${project_elf})
endfunction()

# 6.0 keeps the tool command and its arguments as component properties; earlier
# releases leave them as variables in the project's scope.
function(ctrlpp_esp32_image_binary elf bin)
    idf_component_get_property(esptool esptool_py ESPTOOLPY_CMD)
    idf_component_get_property(elf2image_args esptool_py ESPTOOL_PY_ELF2IMAGE_ARGS)
    if (NOT esptool)
        set(esptool ${ESPTOOLPY})
        set(elf2image_args ${esptool_elf2image_args})
    endif ()
    add_custom_command(OUTPUT "${bin}"
        COMMAND ${esptool} elf2image ${elf2image_args} -o "${bin}" "$<TARGET_FILE:${elf}>"
        DEPENDS ${elf}
        VERBATIM)
endfunction()

function(ctrlpp_esp32_timing_image)
    set(elf ctrlpp_esp32_timing.elf)
    set(bin "${CMAKE_BINARY_DIR}/ctrlpp_esp32_timing.bin")
    set(elf_src "${CMAKE_BINARY_DIR}/ctrlpp_esp32_timing_src.c")
    if (NOT EXISTS "${elf_src}")
        file(WRITE "${elf_src}" "")
    endif ()
    add_executable(${elf} "${elf_src}")
    ctrlpp_esp32_link_components(${elf} ctrlpp_esp32_timing_main)
    ctrlpp_esp32_copy_link_setup(${elf} "${CMAKE_BINARY_DIR}/ctrlpp_esp32_timing.map")
    ctrlpp_esp32_image_binary(${elf} "${bin}")
    add_custom_target(ctrlpp_esp32_timing_bin ALL DEPENDS "${bin}")
endfunction()
