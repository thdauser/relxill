## testing setup (find_package(Catch2 REQUIRED) needs catch2 to be installed in gitlab)
set(CATCH2_FILES catch2/catch_amalgamated.cpp catch2/catch_amalgamated.cpp)

set(UNIT_TEST_FILES_CPP
        test-cppspectrum.cpp tests-execmodel.cpp
        )
      
set(UNIT_TEST_FILES
        test_rellp.c test_rellp.h
        test_relxill.c test_relxill.h
        test_std_functions.c test_std_functions.h
        test_xilltab.c test_xilltab.h
        )


set(EXEC_FILES run-tests)

set(EXEC_FILES_CPP tests)

foreach (execfile ${EXEC_FILES_CPP})
    add_executable(${execfile} ${execfile}.cpp ${UNIT_TEST_FILES_CPP} ${CATCH2_FILES})

    target_link_libraries(${execfile} PUBLIC ${EXTRA_LIBS})
    target_link_libraries(${execfile} PUBLIC m cfitsio)

    target_include_directories(${execfile} PUBLIC
            ${EXTRA_INCLUDES}
            ${PROJECT_BINARY_DIR}
            )

endforeach (execfile ${EXEC_FILES_CPP})


foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${UNIT_TEST_FILES})

    target_link_libraries(${execfile} PUBLIC ${EXTRA_LIBS})
    target_link_libraries(${execfile} PUBLIC m cfitsio)

    target_include_directories(${execfile} PUBLIC
            ${EXTRA_INCLUDES}
            ${PROJECT_BINARY_DIR}
            )

endforeach (execfile ${EXEC_FILES})


install(TARGETS ${EXEC_FILES} ${EXEC_FILES_CPP} DESTINATION bin)


#set(test_execname tests)
#add_executable(${test_execname} ${test_execname}.cpp ${UNIT_TEST_FILES} )
#target_link_libraries(${test_execname} m cfitsio)
#target_link_libraries(${test_execname} Catch2::Catch2)
