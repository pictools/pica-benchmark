include_directories(
    ../../../../pica/include
    ../../../)

add_executable(full-sorted
    Main.cpp)

if (APPLE)
    target_link_libraries(full-sorted
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(full-sorted
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(full-sorted
        benchmark-utility)
endif()

add_custom_target(post-full-sorted ALL DEPENDS full-sorted
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:full-sorted>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:full-sorted>)
set_target_properties(post-full-sorted PROPERTIES FOLDER BuildInfrastructure)
