cmake_minimum_required(VERSION 3.0.0)

project(pica-benchmark)

set(CMAKE_CXX_FLAGS "-g -DNDEBUG -O3 -xCOMMON-AVX512 -qopenmp -qopt-zmm-usage=high")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

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

set(LIB_DIR ${ROOT_DIR}/lib)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${LIB_DIR})

add_subdirectory(pica)
add_subdirectory(src)