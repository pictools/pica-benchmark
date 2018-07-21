include_directories(
    ../../../../pica/include
    ../../../)

add_executable(full-supercells-AoS
    Main.cpp)

if (APPLE)
    target_link_libraries(full-supercells-AoS
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(full-supercells-AoS
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(full-supercells-AoS
        benchmark-utility)
endif()

add_custom_target(post-full-supercells-AoS ALL DEPENDS full-supercells-AoS
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:full-supercells-AoS>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:full-supercells-AoS>)
set_target_properties(post-full-supercells-AoS PROPERTIES FOLDER BuildInfrastructure)
