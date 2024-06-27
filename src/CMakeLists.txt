cmake_minimum_required(VERSION 3.25)

add_library(matrix STATIC matrix.cpp)

target_compile_options(matrix PRIVATE -Wall -Wextra)

target_include_directories(matrix
        PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

include(GNUInstallDirs)

install(TARGETS matrix
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h"
)

