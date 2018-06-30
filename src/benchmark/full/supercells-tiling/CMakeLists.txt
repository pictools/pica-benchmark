include_directories(
    ../../../../pica/include
    ../../../)

add_executable(full-supercells-tiling
    Main.cpp)

if (APPLE)
    target_link_libraries(full-supercells-tiling
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(full-supercells-tiling
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(full-supercells-tiling
        benchmark-utility)
endif()

add_custom_target(post-full-supercells-tiling ALL DEPENDS full-supercells-tiling
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:full-supercells-tiling>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:full-supercells-tiling>)
set_target_properties(post-full-supercells-tiling PROPERTIES FOLDER BuildInfrastructure)
