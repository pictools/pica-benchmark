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
