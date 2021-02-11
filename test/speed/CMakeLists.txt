set(EXEC_FILES speed_test)

foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.cpp)

    target_link_libraries(${execfile} PUBLIC ${EXTRA_LIBS})
    target_link_libraries(${execfile} PUBLIC m cfitsio fftw3)

    target_include_directories(${execfile} PUBLIC
            ${EXTRA_INCLUDES}
            ${PROJECT_BINARY_DIR}
            )

endforeach (execfile ${EXEC_FILES})

install(TARGETS ${EXEC_FILES} DESTINATION bin)

