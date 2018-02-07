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
