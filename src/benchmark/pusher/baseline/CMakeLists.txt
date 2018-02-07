include_directories(
    ../../../../pica/include
    ../../../)

add_executable(pusher-baseline
    Main.cpp)

if (APPLE)
    target_link_libraries(pusher-baseline
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(pusher-baseline
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(pusher-baseline
        benchmark-utility)
endif()

add_custom_target(post-pusher-baseline ALL DEPENDS pusher-baseline
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:pusher-baseline>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:pusher-baseline>)
set_target_properties(post-pusher-baseline PROPERTIES FOLDER BuildInfrastructure)
