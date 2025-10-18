cmake_minimum_required(VERSION 3.20)
project(convhull_3d LANGUAGES C)

add_library(convhull_3d STATIC
    convhull_3d.c
)

target_include_directories(convhull_3d
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(convhull_3d
    PUBLIC
        -Wall -Wextra -Werror -O3
)

target_link_libraries(convhull_3d
    PRIVATE geometry
)
