cmake_minimum_required(VERSION 3.13)

# fhelib dependency
set(TIIFHE_N "(1 << 13)")
set(TIIFHE_T "\"65537\"")
set(TIIFHE_Q "1125899906826241,1125899906629633,1125899905744897")
set(TIIFHE_P "1125899905351681")
set(TIIFHE_OMEGA "3")
add_subdirectory(fhelib)

# example
project(example)

set(CMAKE_C_FLAGS "$ENV{CFLAGS} -std=c89 -pipe -O2 -march=native")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Wundef")
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -std=c++17 -pipe -O2 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wundef")

add_executable(example example.c)
target_compile_options(example PRIVATE -fsanitize=address)
target_link_options(example PRIVATE -fsanitize=address)
target_link_libraries(example tiifhe)
