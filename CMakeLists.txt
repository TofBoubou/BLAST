cmake_minimum_requiered(VERSION 3.20)
project(BLAST_MODERN
        VERSION 1.0.0
        LANGUAGE CXX
        DESCRIPTION "Modern Boundary Layer Analysis and Solution Tool")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

find_package(Mutation++ REQUIRED)
find_package(Eigen3 REQUIRED)

add_subdirectory(src)
add_subdirectory(tests)