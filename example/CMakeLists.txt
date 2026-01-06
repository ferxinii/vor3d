cmake_minimum_required(VERSION 3.20)
project(vor3d_example LANGUAGES C)

if(APPLE)
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(HOMEBREW_LIBOMP)
        set(OpenMP_ROOT "${HOMEBREW_LIBOMP}")
        set(ENV{LDFLAGS} "$ENV{LDFLAGS} -L${OpenMP_ROOT}/lib")
        set(CMAKE_PREFIX_PATH "${OpenMP_ROOT}" ${CMAKE_PREFIX_PATH})
        message(STATUS "Homebrew libomp detected at ${OpenMP_ROOT}")
    else()
        message(FATAL_ERROR "libomp not found via Homebrew; please install it with 'brew install libomp'")
    endif()
endif()

find_package(OpenMP REQUIRED)

# Add the main library (vor3d) as a subdirectory
add_subdirectory(../ vor3d_build)

add_executable(example example.c)

target_link_libraries(example PUBLIC vor3d
        OpenMP::OpenMP_C
)

target_compile_options(example
    PRIVATE
        -Wall -Werror -Wextra -O3
)

# ------ SETUP INPUT / OUTPUT -------
# Create symlink for the lobes directory
add_custom_command(TARGET example POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/lobes
        $<TARGET_FILE_DIR:example>/lobes
    COMMENT "Linking lobes directory into build"
)
# Create empty output directories
foreach(dir plots virtual_population)
    add_custom_command(TARGET example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:example>/${dir}
        COMMENT "Creating empty directory: ${dir}"
    )
endforeach()
