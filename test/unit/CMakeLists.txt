set(UNIT_TEST_FILES
        test_rellp.c  test_rellp.h
        test_relxill.c  test_relxill.h
        test_std_functions.c  test_std_functions.h
        test_xilltab.c  test_xilltab.h
        )



set(CONFIG_FILE ${PROJECT_BINARY_DIR}/config.h)

set(EXEC_FILES test-driver)

foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${CONFIG_FILE} ${UNIT_TEST_FILES} )
    target_link_libraries(${execfile} m cfitsio)

    target_link_libraries(${execfile} PUBLIC ${EXTRA_LIBS})
    target_include_directories(${execfile} PUBLIC
            "${EXTRA_INCLUDES}")

    target_include_directories(${execfile} PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
endforeach (execfile ${EXEC_FILES})

