include_directories(
    ../../../../pica/include
    ../../../)

add_executable(02-pusher-vectorized
    Main.cpp)

if (APPLE)
    target_link_libraries(02-pusher-vectorized
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(02-pusher-vectorized
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(02-pusher-vectorized
        benchmark-utility)
endif()

add_custom_target(post-pusher-vectorized ALL DEPENDS 02-pusher-vectorized
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:02-pusher-vectorized>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:02-pusher-vectorized>)
set_target_properties(post-pusher-vectorized PROPERTIES FOLDER BuildInfrastructure)
