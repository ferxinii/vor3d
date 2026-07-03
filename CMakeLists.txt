cmake_minimum_required(VERSION 3.20)
project(vor3d LANGUAGES C CXX)

if (NOT TARGET geometry)
    add_subdirectory(external/geometry)
endif()

if (NOT TARGET gnuplotC)
    add_subdirectory(external/gnuplotC)
endif()

if (NOT TARGET implicit_predicates)
    add_subdirectory(external/implicit_predicates)  # used by cdt_predicates.cpp (CDT)
endif()

add_library(vor3d STATIC
    src/vor3d.c
    src/bpoly.c
    src/dt.c
    src/cdt.c
    src/cdt_predicates.cpp
    src/scplx.c
    src/acplx.c
    src/vdiagram.c
    src/mirroring.c
    src/ncvx.c
    # src/volsph.c
)

target_include_directories(vor3d
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(vor3d
    PRIVATE
        gnuplotC
        voronoi_predicates
        implicit_predicates  # used by cdt_predicates.cpp (CDT)
    PUBLIC
        geometry  # Users need to know about s_point
)

target_compile_options(vor3d
    PRIVATE
        -Wall -Werror -Wextra
)
