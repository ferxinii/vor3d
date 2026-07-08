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

add_executable(example_basic example_basic.c)
add_executable(example_lungs example_lungs.c)
add_executable(example_nonconvex example_nonconvex.c)
add_executable(example_bug_coplanar example_bug_coplanar.c)
add_executable(example_lobe example_lobe.c)
add_executable(example_medax example_medax.c)
add_executable(example_ashape example_ashape.c)
add_executable(example_merge_orphans example_merge_orphans.c)
# add_executable(example_volsph example_volsph.c)

target_link_libraries(example_basic 
        PUBLIC vor3d OpenMP::OpenMP_C
)
target_link_libraries(example_lungs
        PUBLIC vor3d OpenMP::OpenMP_C
)
target_link_libraries(example_nonconvex
        PUBLIC vor3d OpenMP::OpenMP_C
)
target_link_libraries(example_bug_coplanar
        PUBLIC vor3d OpenMP::OpenMP_C
)
target_link_libraries(example_lobe
        PUBLIC vor3d OpenMP::OpenMP_C
)
target_link_libraries(example_medax
        PUBLIC vor3d m OpenMP::OpenMP_C
)
target_link_libraries(example_ashape
        PUBLIC vor3d m OpenMP::OpenMP_C
)
target_link_libraries(example_merge_orphans
        PUBLIC vor3d m OpenMP::OpenMP_C
)
# target_link_libraries(example_volsph
#         PUBLIC vor3d 
# )

target_compile_options(example_basic 
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_lungs
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_nonconvex
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_lobe
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_medax
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_ashape
    PRIVATE
        -Wall -Werror -Wextra -O3
)
target_compile_options(example_merge_orphans
    PRIVATE
        -Wall -Werror -Wextra -O3
)
# target_compile_options(example_volsph
#     PRIVATE
#         -Wall -Werror -Wextra -O3
# )

# ------ SETUP INPUT / OUTPUT -------
# Create symlink for the lobes directory
# Both example_lungs and example_lobe live in the same build dir and want the
# same `lobes` symlink; `cmake -E create_symlink` errors if it exists and the
# two POST_BUILD steps race under -j.  `ln -sfn` is idempotent and atomic, so
# neither the existing-link nor the race can fail the build.
add_custom_command(TARGET example_lungs POST_BUILD
    COMMAND sh -c "ln -sfn '${CMAKE_CURRENT_SOURCE_DIR}/lobes' '$<TARGET_FILE_DIR:example_lungs>/lobes'"
    COMMENT "Linking lobes directory into build"
)
add_custom_command(TARGET example_lobe POST_BUILD
    COMMAND sh -c "ln -sfn '${CMAKE_CURRENT_SOURCE_DIR}/lobes' '$<TARGET_FILE_DIR:example_lobe>/lobes'"
    COMMENT "Linking lobes directory into build"
)
add_custom_command(TARGET example_medax POST_BUILD
    COMMAND sh -c "ln -sfn '${CMAKE_CURRENT_SOURCE_DIR}/lobes' '$<TARGET_FILE_DIR:example_medax>/lobes'"
    COMMENT "Linking lobes directory into build"
)
# Create empty output directories
foreach(dir plots virtual_population)
    add_custom_command(TARGET example_lungs POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:example_lungs>/${dir}
        COMMENT "Creating empty directory: ${dir}"
    )
endforeach()
