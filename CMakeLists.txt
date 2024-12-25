cmake_minimum_required (VERSION 3.9)
project (qimcifa VERSION 0.1.0 DESCRIPTION "Quantum-inspired Monte Carlo integer factoring algorithm")

# Installation commands
include (GNUInstallDirs)

option (IS_RSA_SEMIPRIME "Optimize for RSA semiprime numbers, only. (Might restrict bounds too far on certain RSA semiprimes.)" ON)
option (IS_DISTRIBUTED "Enable node distribution configuration dialog." ON)
option (IS_SQUARES_CONGRUENCE_CHECK "Optionally additionally check random number generator outputs for factoring via congruence of squares." OFF)
option (USE_GMP "Use GMP library instead of Boost or pure language for arbitrary precision integers." OFF)
option (USE_BOOST "Use Boost library instead of pure language for arbitrary precision integers." ON)
set(BIG_INT_BITS "64" CACHE STRING "Change the maximum bit width of arbitrary precision 'big integers.'")

message ("Optimize for RSA semiprime numbers: ${IS_RSA_SEMIPRIME}")
message ("Enable node distribution configuration dialog: ${IS_DISTRIBUTED}")
message ("Enable congruence of squares check: ${IS_SQUARES_CONGRUENCE_CHECK}")
message ("Use GMP library instead of Boost or pure language for arbitrary precision integers: ${USE_GMP}")
message ("Use Boost library instead of pure language for arbitrary precision integers: ${USE_BOOST}")
message ("Maximum bit width of arbitrary precision 'big integers': ${BIG_INT_BITS}")

configure_file(include/common/config.h.in include/common/config.h @ONLY)
include_directories ("include" "include/common")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/common)

if (MSVC)
    set(CMAKE_CXX_FLAGS "/Wall -O3")
else (MSVC)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3")
endif (MSVC)

if (USE_GMP OR USE_BOOST)
add_executable (qimcifa
    src/qimcifa.cpp
    )
add_executable (qimcifa_tuner
    src/qimcifa_tuner.cpp
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
add_executable (qimcifa_tuner
    src/qimcifa_tuner.cpp
    src/common/big_integer.cpp
    )
add_executable (prime_generator
    src/prime_generator.cpp
    src/common/dispatchqueue.cpp
    src/common/big_integer.cpp
    )
endif (USE_GMP OR USE_BOOST)
if (MSVC)
    if (USE_GMP)
        target_link_libraries (qimcifa gmp)
        target_link_libraries (prime_generator gmp)
    else (USE_GMP)
        target_link_libraries (qimcifa)
        target_link_libraries (prime_generator)
    endif (USE_GMP)
else (MSVC)
    if (USE_GMP)
        target_link_libraries (qimcifa pthread gmp)
        target_link_libraries (prime_generator pthread gmp)
    else (USE_GMP)
        target_link_libraries (qimcifa pthread)
        target_link_libraries (prime_generator pthread)
    endif (USE_GMP)
endif (MSVC)
target_compile_features(prime_generator PRIVATE cxx_std_17)
