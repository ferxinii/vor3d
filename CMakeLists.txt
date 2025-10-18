cmake_minimum_required(VERSION 3.20)
project(vor3d LANGUAGE C)

add_subdirectory(external/convhull_3d)
add_subdirectory(external/geometry)
add_subdirectory(external/gnuplotC)

add_library(vor3d STATIC
    algebra.c
    vor3d.c
    bpoly.c
    delaunay.c
    simplical_complex.c
    vdiagram.c
    convh.c
)

target_include_directories(vor3d
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(vor3d 
    PRIVATE
        convhull_3d
        geometry
        gnuplotC
)

target_compile_options(vor3d
    PRIVATE
        -Wall -Werror -Wextra -O3
)
