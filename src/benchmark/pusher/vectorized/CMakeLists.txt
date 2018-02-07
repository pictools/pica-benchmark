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

add_custom_target(post-pusher-vectorized ALL DEPENDS pusher-vectorized
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:pusher-vectorized>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:pusher-vectorized>)
set_target_properties(post-pusher-vectorized PROPERTIES FOLDER BuildInfrastructure)
