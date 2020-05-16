cmake_minimum_required(VERSION 2.8)
project(serreflection)

include_directories (include)
include_directories (lib/boostlog0/include)

add_subdirectory(lib/boostlog0)
add_subdirectory(test)
