cmake_minimum_required(VERSION 3.25)

add_executable(Matrix_test Matrix_test.cpp)

target_link_libraries(Matrix_test PUBLIC matrix)

target_include_directories(Matrix_test PUBLIC "../src")

include(FetchContent)

FetchContent_Declare(
	"googletest"
	GIT_REPOSITORY "https://github.com/google/googletest.git"
	GIT_TAG "v1.13.0"
	)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(${PROJECT_NAME}_googletest matrix_googletest.cpp)
target_link_libraries(${PROJECT_NAME}_googletest PUBLIC gtest_main)
target_link_libraries(${PROJECT_NAME}_googletest PUBLIC matrix)
target_include_directories(${PROJECT_NAME}_googletest PUBLIC ../src)

add_test(NAME test COMMAND ${PROJECT_NAME}_googletest)
enable_testing()
