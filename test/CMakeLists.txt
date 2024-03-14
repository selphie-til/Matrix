cmake_minimum_required(VERSION 3.25)

set(CommonTargetProps PUBLIC matrix "../src")

add_executable(Matrix_test Matrix_test.cpp)
target_link_libraries(Matrix_test PUBLIC matrix)
target_include_directories(Matrix_test ${CommonTargetProps})

include(FetchContent)
FetchContent_Declare(
	"googletest"
	GIT_REPOSITORY "https://github.com/google/googletest.git"
	GIT_TAG "v1.13.0"
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(gtest_exec ${PROJECT_NAME}_gtest)
add_executable(${gtest_exec} matrix_googletest.cpp)
target_link_libraries(${gtest_exec} PUBLIC gtest_main matrix)
target_include_directories(${gtest_exec} ${CommonTargetProps})

add_test(NAME test COMMAND ${gtest_exec})
enable_testing()
