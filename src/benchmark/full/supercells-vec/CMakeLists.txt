include_directories(
    ../../../../pica/include
    ../../../)

add_executable(full-supercells-vec
    Main.cpp)

if (APPLE)
    target_link_libraries(full-supercells-vec
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(full-supercells-vec
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(full-supercells-vec
        benchmark-utility)
endif()

add_custom_target(post-full-supercells-vec ALL DEPENDS full-supercells-vec
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:full-supercells-vec>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:full-supercells-vec>)
set_target_properties(post-full-supercells-vec PROPERTIES FOLDER BuildInfrastructure)
