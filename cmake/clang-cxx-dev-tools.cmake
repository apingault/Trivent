# From http://www.labri.fr/perso/fleury/posts/programming/using-clang-tidy-and-clang-format.html
# Additional targets to perform clang-format/clang-tidy
# Get all project files
file(GLOB_RECURSE
     ALL_CXX_SOURCE_FILES
     *.[chi]pp *.[chi]xx *.cc *.hh *.ii *.[CHI]
     )

# Adding clang-format target if executable is found
# Will look for a .clang-format file at the root of the project
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
  add_custom_target(
    clang-format
    COMMAND clang-format
    -i
    -style=file
    ${ALL_CXX_SOURCE_FILES}
    )
endif()

# Adding clang-tidy target if executable is found
# Will look for a .clang-tidy file at the root of the project
find_program(CLANG_TIDY "clang-tidy")
if(CLANG_TIDY)
  add_custom_target(
    clang-tidy
    COMMAND clang-tidy
    ${ALL_CXX_SOURCE_FILES}
    -config=''
    --
    -std=c++14
    ${INCLUDE_DIRECTORIES}
    )
endif()
