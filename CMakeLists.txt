cmake_minimum_required(VERSION 3.25)

project(Tiled_Matrix VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

include(CTest)

add_subdirectory(src)
add_subdirectory(test)

