cmake_minimum_required(VERSION 3.15)

project(relxill VERSION 1.4.20210104)


set(PROJECT_ROOT_DIR "-DRELXILL_SOURCE_DIR='\"${CMAKE_SOURCE_DIR}\"'")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic ${PROJECT_ROOT_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-std=c++17")

include_directories($ENV{HEADAS}/include/)
link_directories($ENV{HEADAS}/lib/)

set(SRC_DIR src)

configure_file(${SRC_DIR}/config.h.in config.h) # written into binary tree

list(APPEND EXTRA_LIBS Relxill)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/${SRC_DIR}")

add_subdirectory(${SRC_DIR})
add_subdirectory(test)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(FILES "${PROJECT_BINARY_DIR}/config.h" DESTINATION model)

## testing setup
find_package(Catch2 REQUIRED)
set(test_execname tests)
add_executable(${test_execname} test/unit/${test_execname}.cpp ${SOURCE_FILES_CPP} test/unit/tests-execmodel.cpp test/unit/test-cppspectrum.cpp)
target_link_libraries(${test_execname} m cfitsio)
target_link_libraries(${test_execname} Catch2::Catch2)
