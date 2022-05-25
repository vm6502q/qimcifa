cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

include_directories ("include" "include/common")

set (OPENCL_AMDSDK /opt/AMDAPPSDK-3.0 CACHE PATH "Installation path for the installed AMD OpenCL SDK, if used")

# Options used when building the project
find_package (OpenCL)
if (NOT OpenCL_FOUND)
    # Attempt with AMD's OpenCL SDK
    find_library (LIB_OPENCL OpenCL PATHS ${OPENCL_AMDSDK}/lib/x86_64/)
    if (NOT LIB_OPENCL)
        message( FATAL_ERROR "The OpenCL library must be installed, CMake will exit." )
    else ()
        # Found, set the required include path.
        set (OpenCL_INCLUDE_DIRS ${OPENCL_AMDSDK}/include CACHE PATH "AMD OpenCL SDK Header include path")
        set (OpenCL_COMPILATION_OPTIONS
            -Wno-ignored-attributes
            -Wno-deprecated-declarations
            CACHE STRING "AMD OpenCL SDK Compilation Option Requirements")
        message ("OpenCL support found in the AMD SDK")
    endif()
endif ()

message ("Found OpenCL.")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")

add_executable (qimcifa
    src/qimcifa.cpp
    )
target_link_libraries (qimcifa pthread)

