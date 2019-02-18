include_directories(
    ../../../../pica/include
    ../../../)

add_executable(02-full-supercells
    Main.cpp)

if (APPLE)
    target_link_libraries(02-full-supercells
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(02-full-supercells
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(02-full-supercells
        benchmark-utility)
endif()

add_custom_target(post-full-supercells-tiling ALL DEPENDS 02-full-supercells
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:02-full-supercells>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:02-full-supercells>)
set_target_properties(post-full-supercells-tiling PROPERTIES FOLDER BuildInfrastructure)
