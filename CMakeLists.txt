cmake_minimum_required(VERSION 3.13)

# CFLAGS/CXXFLAGS
set(CMAKE_C_FLAGS "$ENV{CFLAGS} -std=c89 -pipe")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Wundef")
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -std=c++17 -pipe")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wundef")

# options
option(TIIFHE_BUILD_DEBUG "Build in debug mode" OFF)
option(TIIFHE_BUILD_TESTS "Build tests" ON)
option(TIIFHE_FIX_SEED    "Use fix seeds" OFF)
option(TIIFHE_LOG_ERROR   "Log error" OFF)

if(TIIFHE_BUILD_DEBUG)
	set(CMAKE_C_FLAGS   "-Og -ggdb3 -gdwarf ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-Og -ggdb3 -gdwarf ${CMAKE_CXX_FLAGS}")
else()
	set(CMAKE_C_FLAGS   "-O2 -march=native ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-O2 -march=native ${CMAKE_CXX_FLAGS}")
endif()

# config
set(TIIFHE_MSG   "mpz"  CACHE STRING "fhelib message implementation")
set(TIIFHE_RING  "hexl" CACHE STRING "fhelib ring implementation")

set(TIIFHE_N     ""     CACHE STRING "fhelib ring dimension")
set(TIIFHE_T     ""     CACHE STRING "fhelib plaintext modulus")
set(TIIFHE_Q     ""     CACHE STRING "fhelib ciphertext modulus")
set(TIIFHE_P     ""     CACHE STRING "fhelib key switching modulus")
set(TIIFHE_OMEGA ""     CACHE STRING "fhelib key switching decomposition number")
include(params.cmake)

# project
project(TIIFHE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)
if(TIIFHE_BUILD_TESTS)
	add_subdirectory(tests)
endif()
