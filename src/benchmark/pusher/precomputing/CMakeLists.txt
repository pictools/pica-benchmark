include_directories(
    ../../../../pica/include
    ../../../)

add_executable(pusher-precomputing
    Main.cpp)

if (APPLE)
    target_link_libraries(pusher-precomputing
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(pusher-precomputing
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(pusher-precomputing
        benchmark-utility)
endif()

add_custom_target(post-pusher-precomputing ALL DEPENDS pusher-precomputing
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:pusher-precomputing>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:pusher-precomputing>)
set_target_properties(post-pusher-precomputing PROPERTIES FOLDER BuildInfrastructure)
