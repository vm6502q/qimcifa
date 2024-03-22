cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

option (IS_RSA_SEMIPRIME "Optimize for RSA semiprime numbers, only. (Might restrict bounds too far on certain RSA semiprimes.)" ON)
option (IS_PARALLEL "Enable parallelism within nodes." ON)
option (USE_GMP "Use GMP library instead of Boost or pure language for arbitrary precision integers." OFF)
option (USE_BOOST "Use Boost library instead of pure language for arbitrary precision integers." ON)
set(BIG_INT_BITS "64" CACHE STRING "Change the maximum bit width of arbitrary precision 'big integers.'")

message ("Optimize for RSA semiprime numbers: ${IS_RSA_SEMIPRIME}")
message ("Enable parallelism within nodes: ${IS_PARALLEL}")
message ("Enable 'quantum-inspired' randomness: ${IS_RANDOM}")
message ("Enable congruence of squares check: ${IS_SQUARES_CONGRUENCE_CHECK}")
message ("Use GMP library instead of Boost or pure language for arbitrary precision integers: ${USE_GMP}")
message ("Use Boost library instead of pure language for arbitrary precision integers: ${USE_BOOST}")
message ("Maximum bit width of arbitrary precision 'big integers': ${BIG_INT_BITS}")

configure_file(include/common/config.h.in include/common/config.h @ONLY)
include_directories ("include" "include/common")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/common)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3")

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
add_library(prime_gen MODULE src/prime_gen.cpp)
target_link_libraries(prime_gen pybind11::module)
configure_file(src/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py COPYONLY)
configure_file(src/prime_gen.cpp ${CMAKE_CURRENT_BINARY_DIR}/prime_gen.cpp COPYONLY)
configure_file(${CMAKE_CURRENT_BINARY_DIR}/include/common/config.h ${CMAKE_CURRENT_BINARY_DIR}/config.h COPYONLY)
configure_file(include/big_integer.hpp ${CMAKE_CURRENT_BINARY_DIR}/big_integer.hpp COPYONLY)
configure_file(include/dispatchqueue.hpp ${CMAKE_CURRENT_BINARY_DIR}/dispatchqueue.hpp COPYONLY)
configure_file(src/common/big_integer.cpp ${CMAKE_CURRENT_BINARY_DIR}/big_integer.cpp COPYONLY)
configure_file(src/common/dispatchqueue.cpp ${CMAKE_CURRENT_BINARY_DIR}/dispatchqueue.cpp COPYONLY)

if (USE_GMP OR USE_BOOST)
add_executable (qimcifa
    src/qimcifa.cpp
    )
add_executable (prime_generator
    src/prime_generator.cpp
    src/common/dispatchqueue.cpp
    )
else (USE_GMP OR USE_BOOST)
add_executable (qimcifa
    src/qimcifa.cpp
    src/common/big_integer.cpp
    )
add_executable (prime_generator
    src/prime_generator.cpp
    src/common/dispatchqueue.cpp
    src/common/big_integer.cpp
    )
endif (USE_GMP OR USE_BOOST)
if (USE_GMP)
    target_link_libraries (qimcifa pthread gmp)
    target_link_libraries (prime_generator pthread gmp)
else (USE_GMP)
    target_link_libraries (qimcifa pthread)
    target_link_libraries (prime_generator pthread)
endif (USE_GMP)
