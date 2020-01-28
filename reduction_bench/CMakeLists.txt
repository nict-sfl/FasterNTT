cmake_minimum_required(VERSION 3.00)
project(reduction_bench C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Ofast -funroll-loops -fomit-frame-pointer -march=native -mtune=native")

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(bench64 bench64.c)
set_target_properties(bench64 PROPERTIES COMPILE_FLAGS "-save-temps")
target_link_libraries(bench64 gmp)

add_executable(test64 test64.c)
set_target_properties(test64 PROPERTIES COMPILE_FLAGS "-save-temps")
target_link_libraries(test64 gmp)