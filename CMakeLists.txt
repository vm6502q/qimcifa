cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

include_directories ("include" "include/common")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/common)

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

# Build the OpenCL command files
find_program (XXD_BIN xxd)
file (GLOB_RECURSE COMPILABLE_RESOURCES "src/common/*.cl")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/common/)
foreach (INPUT_FILE ${COMPILABLE_RESOURCES})
    get_filename_component (INPUT_NAME ${INPUT_FILE} NAME)
    get_filename_component (INPUT_BASENAME ${INPUT_FILE} NAME_WE)
    get_filename_component (INPUT_DIR ${INPUT_FILE} DIRECTORY)

    set (OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/include/common/${INPUT_BASENAME}cl.hpp)

    message (" Creating XXD Rule for ${INPUT_FILE} -> ${OUTPUT_FILE}")
    add_custom_command (
        WORKING_DIRECTORY ${INPUT_DIR}
        OUTPUT ${OUTPUT_FILE}
        COMMAND ${XXD_BIN} -i ${INPUT_NAME} > ${OUTPUT_FILE}
        COMMENT "Building OpenCL Commands in ${INPUT_FILE}"
        )
    list (APPEND COMPILED_RESOURCES ${OUTPUT_FILE})
endforeach ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")

add_executable (qimcifa
    src/qimcifa.cpp
    )
target_link_libraries (qimcifa pthread)

add_executable (qimcifa_cl
    ${COMPILED_RESOURCES}
    src/qimcifa_cl.cpp
    src/common/oclengine.cpp
    )
target_link_libraries (qimcifa_cl OpenCL)

target_compile_definitions (qimcifa_cl PUBLIC CL_HPP_TARGET_OPENCL_VERSION=200)
target_compile_definitions (qimcifa_cl PUBLIC CL_HPP_MINIMUM_OPENCL_VERSION=110)

