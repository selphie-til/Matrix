cmake_minimum_required(VERSION 3.25)

add_executable(Matrix_test Matrix_test.cpp)

target_link_libraries(Matrix_test PUBLIC matrix)

target_include_directories(Matrix_test PUBLIC "../src")
