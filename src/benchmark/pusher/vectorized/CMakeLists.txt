include_directories(
    ../../../../pica/include
    ../../../)

add_executable(pusher-vectorized
    Main.cpp)

if (APPLE)
    target_link_libraries(pusher-vectorized
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(pusher-vectorized
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(pusher-vectorized
        benchmark-utility)
endif()
