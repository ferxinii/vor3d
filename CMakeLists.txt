cmake_minimum_required(VERSION 3.20)
project(vor3d LANGUAGES C)

if (NOT TARGET geometry)
    add_subdirectory(external/geometry)
endif()

if (NOT TARGET gnuplotC)
    add_subdirectory(external/gnuplotC)
endif()

add_library(vor3d STATIC
    src/algebra.c
    src/vor3d.c
    src/bpoly.c
    src/delaunay.c
    src/simplical_complex.c
    src/vdiagram.c
)

target_include_directories(vor3d
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(vor3d 
    PRIVATE
        gnuplotC
    PUBLIC
        geometry  # Users need to know about s_point
)

target_compile_options(vor3d
    PRIVATE
        -Wall -Werror -Wextra -O3
)
