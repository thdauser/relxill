cmake_minimum_required(VERSION 3.15)

project(relxill VERSION 1.4.20210104)


set(PROJECT_ROOT_DIR "-DRELXILL_SOURCE_DIR='\"${CMAKE_SOURCE_DIR}\"'")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic ${PROJECT_ROOT_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-std=c++17")

include_directories($ENV{HEADAS}/include/)
link_directories($ENV{HEADAS}/lib/)

configure_file(src/config.h.in config.h) # written into binary tree

add_subdirectory(src test)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(FILES "${PROJECT_BINARY_DIR}/config.h" DESTINATION model)
install(FILES "${PROJECT_BINARY_DIR}/lmodel_relxill.dat" DESTINATION model)


