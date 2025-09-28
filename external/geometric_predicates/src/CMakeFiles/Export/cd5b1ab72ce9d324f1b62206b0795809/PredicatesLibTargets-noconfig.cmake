#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Predicates::predicates" for configuration ""
set_property(TARGET Predicates::predicates APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(Predicates::predicates PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libpredicates.dylib"
  IMPORTED_SONAME_NOCONFIG "@rpath/libpredicates.dylib"
  )

list(APPEND _cmake_import_check_targets Predicates::predicates )
list(APPEND _cmake_import_check_files_for_Predicates::predicates "${_IMPORT_PREFIX}/lib/libpredicates.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
