set(UNIT_TEST_FILES
        test_rellp.c  test_rellp.h
        test_relxill.c  test_relxill.h
        test_std_functions.c  test_std_functions.h
        test_xilltab.c  test_xilltab.h
        )


set(EXEC_FILES run-tests)

foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${UNIT_TEST_FILES} )

    target_link_libraries(${execfile} PUBLIC ${EXTRA_LIBS})
    target_link_libraries(${execfile} PUBLIC m cfitsio)

    target_include_directories(${execfile} PUBLIC
            ${EXTRA_INCLUDES}
            ${PROJECT_BINARY_DIR}
            )

endforeach (execfile ${EXEC_FILES})

install(TARGETS ${EXEC_FILES} DESTINATION bin)
