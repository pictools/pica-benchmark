include_directories(
    ../../../../pica/include
    ../../../)

add_executable(01-pusher-baseline
    Main.cpp)

if (APPLE)
    target_link_libraries(01-pusher-baseline
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(01-pusher-baseline
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(01-pusher-baseline
        benchmark-utility)
endif()

add_custom_target(post-pusher-baseline ALL DEPENDS 01-pusher-baseline
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:01-pusher-baseline>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:01-pusher-baseline>)
set_target_properties(post-pusher-baseline PROPERTIES FOLDER BuildInfrastructure)
