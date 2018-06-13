include_directories(
    ../../../../pica/include
    ../../../)

add_executable(full-supercells
    Main.cpp)

if (APPLE)
    target_link_libraries(full-supercells
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(full-supercells
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(full-supercells
        benchmark-utility)
endif()

add_custom_target(post-full-supercells ALL DEPENDS full-supercells
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:full-supercells>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:full-supercells>)
set_target_properties(post-full-supercells PROPERTIES FOLDER BuildInfrastructure)
