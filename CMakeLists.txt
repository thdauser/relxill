cmake_minimum_required(VERSION 3.15)
project(relxill C)

set(CMAKE_C_STANDARD 11)

include_directories($ENV{HEADAS}/include/)
link_directories($ENV{HEADAS}/lib/)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic ")


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
        xilltable.c
        xilltable.h
        test_relxill.c test_relxill.h)

set(EXEC_FILES
        test_sta
        test_xilltab
        )

foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${SOURCE_FILES})
    target_link_libraries(${execfile} m cfitsio)
endforeach (execfile ${EXEC_FILES})

target_link_libraries(test_sta m)
target_link_libraries(test_sta cfitsio)
