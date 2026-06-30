cmake_minimum_required(VERSION 3.14)
project(implicit_predicates)

add_library(implicit_predicates INTERFACE)

target_include_directories(implicit_predicates SYSTEM INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/implicit_predicates>
)

target_compile_features(implicit_predicates INTERFACE cxx_std_20)

target_compile_options(implicit_predicates INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-frounding-math>
    $<$<CXX_COMPILER_ID:MSVC>:/fp:strict>
)
