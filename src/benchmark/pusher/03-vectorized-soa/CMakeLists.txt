include_directories(
    ../../../../pica/include
    ../../../)

add_executable(03-pusher-vectorized-soa
    Main.cpp)

if (APPLE)
    target_link_libraries(03-pusher-vectorized-soa
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(03-pusher-vectorized-soa
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(03-pusher-vectorized-soa
        benchmark-utility)
endif()

add_custom_target(post-pusher-vectorized-soa ALL DEPENDS 03-pusher-vectorized-soa
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:03-pusher-vectorized-soa>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:03-pusher-vectorized-soa>)
set_target_properties(post-pusher-vectorized-soa PROPERTIES FOLDER BuildInfrastructure)
