cmake_minimum_required(VERSION 3.0.0)

project(pica-benchmark)

set(CMAKE_CXX_FLAGS "-DNDEBUG -O2")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

option(USE_TESTS "Enable unit tests" ON)
option(USE_OPENMP "Use OpenMP" ON)
option(USE_MIC "Use MIC" OFF)
option(USE_AVX "Use AVX" OFF)

if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_EXE_LINKER_FLAGS "/NODEFAULTLIB:libcmt.lib")
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(CMAKE_EXE_LINKER_FLAGS "/machine:x64 ${CMAKE_EXE_LINKER_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "/machine:x64 ${CMAKE_EXE_LINKER_FLAGS}")
        set(CMAKE_MODULE_LINKER_FLAGS "/machine:x64")
        set(CMAKE_STATIC_LINKER_FLAGS "/machine:x64")
    endif()
endif()

set(ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(ARCH "x64")

if (DCMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND)
    message(WARNING "Cannot find ${DCMAKE_CXX_COMPILER} compiler")
endif()

message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER_ID}")

if (USE_MIC)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmic")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
        set(ARCH "mic")
        message(STATUS "Building code for Intel Xeon Phi coprocessor")
    else()
        message(FATAL_ERROR "Use Intel Composer to build code for Xeon Phi coprocessors")
    endif()
endif()

if (USE_AVX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(AVX_FLAG "-xavx")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(AVX_FLAG "-mavx")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        message(WARNING "AVX is not supported by $CMAKE_CXX_COMPILER_ID compiler")
        set(AVX_FLAG "")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AVX_FLAG}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AVX_FLAG}")
    add_definitions(-DUSE_AVX)
endif()

set(LIB_DIR ${ROOT_DIR}/lib)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${LIB_DIR})

add_subdirectory(pica)
add_subdirectory(src)