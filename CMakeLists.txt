cmake_minimum_required(VERSION 3.00)
project(FasterNTT CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Ofast -funroll-loops -fomit-frame-pointer -march=native -mtune=native")

include_directories(
        ./include
        external/NFLlib/include
        external/YELL/include
        /usr/local/include)

link_directories(/usr/local/lib)

set(NFLLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/NFLlib/build)
set(YELL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/YELL/lib)

add_subdirectory(reduction_bench)

add_executable(test test.cpp)
target_link_libraries(test gmp ntl pthread crypto)

add_executable(bench bench.cpp)
target_link_libraries(bench benchmark pthread crypto)
set_target_properties(bench PROPERTIES COMPILE_FLAGS "-save-temps")

add_executable(bench_cmp_ntt_impl bench_cmp_ntt_impl.cpp)
target_link_libraries(bench_cmp_ntt_impl benchmark pthread crypto)

add_executable(bench_cmp_lib bench_cmp_lib.cpp)
target_link_libraries(bench_cmp_lib
        ${NFLLIB_DIR}/libnfllib_static.a
        ${YELL_DIR}/libyell.a
        benchmark
        pthread
        crypto)
