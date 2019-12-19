cmake_minimum_required(VERSION 2.8.12)
project(tetgen)

add_library(tetgen STATIC tetgen.cxx predicates.cxx)
target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)
