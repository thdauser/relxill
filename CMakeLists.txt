cmake_minimum_required(VERSION 3.15)
project(relxill)

set(CMAKE_C_STANDARD 11)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-std=c++17")

include_directories($ENV{HEADAS}/include/)
link_directories($ENV{HEADAS}/lib/)

set(PROJECT_ROOT_DIR "-DRELXILL_SOURCE_DIR='\"${CMAKE_SOURCE_DIR}\"'")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic ${PROJECT_ROOT_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROJECT_ROOT_DIR}")

## testing setup
find_package(Catch2 REQUIRED)

set(SOURCE_FILES
        common.h
        donthcomp.c
        relbase.c
        relbase.h
        relcache.c
        relcache.h
        rellp.c
        rellp.h
        relmodels.c
        relmodels.h
        reltable.c
        reltable.h
        relutility.c
        relutility.h
        xilltable.c xilltable.h
        test_relxill.c test_relxill.h
        test_rellp.c test_rellp.h
        relphysics.c relphysics.h
        test_std_functions.c test_std_functions.h
         )

set(SOURCE_FILES_CPP
        ${SOURCE_FILES} cppmodels.cpp cppmodels.h cppparameters.h cppspectrum.h cppTypes.h cppModelDatabase.h cppparameters.cpp)

set(EXEC_FILES test_sta)

set(EXEC_FILES_CPP test_cpp)


foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${SOURCE_FILES})
    target_link_libraries(${execfile} m cfitsio)
endforeach (execfile ${EXEC_FILES})

foreach (execfile ${EXEC_FILES_CPP})
    add_executable(${execfile} ${execfile}.cpp ${SOURCE_FILES_CPP})
    target_link_libraries(${execfile} m cfitsio)
endforeach (execfile ${EXEC_FILES_CPP})


## testing setup
set(test_execname tests)
add_executable(${test_execname} test/unit/${test_execname}.cpp ${SOURCE_FILES_CPP} test/unit/tests-execmodel.cpp)
target_link_libraries(${test_execname} m cfitsio)
target_link_libraries(${test_execname} Catch2::Catch2)
