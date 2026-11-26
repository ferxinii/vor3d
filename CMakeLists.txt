cmake_minimum_required(VERSION 3.20)
project(vor3d LANGUAGES C CXX)

if (NOT TARGET geometry)
    # Optional parts of geometry that vor3d needs:
    #  TESTS_POWER    test_orthosphere/test_orthosphere_w (dt.c, acplx.c, ashape.c)
    #  LP             lp2_*/lp3_* LP feasibility, lp.h (vor3d.c, ncvx.c)
    #  TESTS_IMPLICIT itest_* exact implicit predicates (dt.c, cdt.c)
    set(GEOMETRY_COMPILE_TESTS_POWER ON CACHE BOOL "" FORCE)
    set(GEOMETRY_COMPILE_LP ON CACHE BOOL "" FORCE)
    set(GEOMETRY_COMPILE_TESTS_IMPLICIT ON CACHE BOOL "" FORCE)
    add_subdirectory(external/geometry)
endif()

if (NOT TARGET gnuplotC)
    add_subdirectory(external/gnuplotC)
endif()

add_library(vor3d STATIC
    src/vor3d.c
    src/bpoly.c
    src/dt.c
    src/cdt.c
    src/scplx.c
    src/acplx.c
    src/ashape.c
    src/vdiagram.c
    src/ncvx.c
    src/medax.c
    src/volsph.c
)

target_include_directories(vor3d
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(vor3d
    PRIVATE
        gnuplotC
        predicates_basic     # orient3d_dd (dt.c)
        predicates_lp        # lp3_* + power_lp3_* (surf_vkey.h)
        predicates_power     # power_incircle3d (dt.c)
    PUBLIC
        geometry  # Users need to know about s_point, gtests_implicit predicates
)

target_compile_options(vor3d
    PRIVATE
        -Wall -Werror -Wextra
)
