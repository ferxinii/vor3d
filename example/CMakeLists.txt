cmake_minimum_required(VERSION 3.20)
project(vor3d_example LANGUAGES C)

add_subdirectory(../vor3d)

add_executable(example example.c)

target_link_libraries(example PRIVATE vor3d)
