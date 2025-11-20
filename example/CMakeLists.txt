cmake_minimum_required(VERSION 3.20)
project(vor3d_example LANGUAGES C)

# Add the main library (vor3d) as a subdirectory
add_subdirectory(../ vor3d_build)

add_executable(example example.c)

target_link_libraries(example PUBLIC vor3d)

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
