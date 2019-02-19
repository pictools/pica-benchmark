include_directories(
    ../../../../pica/include
    ../../../)

add_executable(04-pusher-vectorized-soa-rcp
    Main.cpp)

if (WIN32)
    target_compile_options(04-pusher-vectorized-soa-rcp PRIVATE /QxCORE-AVX512 /Qno-prec-div)
else()
    target_compile_options(04-pusher-vectorized-soa-rcp PRIVATE -xCORE-AVX512 -no-prec-div)
endif()

if (APPLE)
    target_link_libraries(04-pusher-vectorized-soa-rcp
        -Wl,-all_load
        benchmark-utility
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(04-pusher-vectorized-soa-rcp
        -Wl,-whole-archive
        benchmark-utility
        -Wl,-no-whole-archive)
else()
    target_link_libraries(04-pusher-vectorized-soa-rcp
        benchmark-utility)
endif()

add_custom_target(post-pusher-vectorized-soa-rcp ALL DEPENDS 04-pusher-vectorized-soa-rcp
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:04-pusher-vectorized-soa-rcp>
                               ${ROOT_DIR}/bin/$<TARGET_FILE_NAME:04-pusher-vectorized-soa-rcp>)
set_target_properties(post-pusher-vectorized-soa-rcp PROPERTIES FOLDER BuildInfrastructure)
