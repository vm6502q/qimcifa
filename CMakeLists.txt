cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

include_directories ("include" "include/common")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/common)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3")

add_executable (qimcifa
    src/qimcifa.cpp
    )
target_link_libraries (qimcifa pthread)
# target_link_libraries (qimcifa pthread gmp)

# Options used when building the project
find_package (OpenCL)
if (OpenCL_FOUND)
    # Build the OpenCL command files
    find_program (XXD_BIN xxd)
    if (XXD_BIN)
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

        add_executable (qimcifa_cl
            ${COMPILED_RESOURCES}
            src/qimcifa_cl.cpp
            src/common/oclengine.cpp
            )
        target_link_libraries (qimcifa_cl OpenCL)
        target_compile_definitions (qimcifa_cl PUBLIC CL_HPP_TARGET_OPENCL_VERSION=200)
        target_compile_definitions (qimcifa_cl PUBLIC CL_HPP_MINIMUM_OPENCL_VERSION=110)
    endif()
endif ()
