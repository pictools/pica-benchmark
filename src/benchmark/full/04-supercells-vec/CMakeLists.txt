include_directories(
    ../../../../pica/include
    ../../../)

add_executable(04-full-supercells-vec
    Main.cpp)

if (APPLE)
    target_link_libraries(04-full-supercells-vec
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(04-full-supercells-vec
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(04-full-supercells-vec
        benchmark-utility)
endif()

add_custom_target(post-full-supercells-vec ALL DEPENDS 04-full-supercells-vec
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:04-full-supercells-vec>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:04-full-supercells-vec>)
set_target_properties(post-full-supercells-vec PROPERTIES FOLDER BuildInfrastructure)
