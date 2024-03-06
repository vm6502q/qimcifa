cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

option (IS_RSA_SEMIPRIME "Optimize for RSA semiprime numbers, only. (Might restrict bounds too far on certain RSA semiprimes.)" ON)
option (IS_DISTRIBUTED "Enable node distribution configuraiton dialog." ON)
option (IS_SQUARES_CONGRUENCE_CHECK "Optionally additionally check random number generator outputs for factoring via congruence of squares." OFF)
option (USE_GMP "Use GMP library instead of Boost or pure language for arbitrary precision integers." OFF)
option (USE_BOOST "Use Boost library instead of pure language for arbitrary precision integers." OFF)
set(BIG_INT_BITS "128" CACHE STRING "Change the maximum bit width of arbitrary precision 'big integers.'")

configure_file(include/common/config.h.in include/common/config.h @ONLY)
include_directories ("include" "include/common")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/common)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3")

if (USE_GMP OR USE_BOOST)
add_executable (qimcifa
    src/qimcifa.cpp
    )
add_executable (qimcifa_tuner
    src/qimcifa_tuner.cpp
    )
else (USE_GMP OR USE_BOOST)
add_executable (qimcifa
    src/qimcifa.cpp
    src/common/big_integer.cpp
    )
add_executable (qimcifa_tuner
    src/qimcifa_tuner.cpp
    src/common/big_integer.cpp
    )
endif (USE_GMP OR USE_BOOST)
if (USE_GMP)
    target_link_libraries (qimcifa pthread gmp)
    target_link_libraries (qimcifa_tuner pthread gmp)
else (USE_GMP)
    target_link_libraries (qimcifa pthread)
    target_link_libraries (qimcifa_tuner pthread)
endif (USE_GMP)
