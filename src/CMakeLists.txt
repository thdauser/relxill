set(SOURCE_FILES
        common.h
        donthcomp.c
        relbase.c relbase.h
        relcache.c relcache.h
        rellp.c rellp.h
        relmodels.c relmodels.h
        reltable.c reltable.h
        relutility.c relutility.h
        xilltable.c xilltable.h
        relphysics.c relphysics.h
        cppModelDatabase.h
        cppmodels.cpp  cppmodels.h
        cppparameters.cpp  cppparameters.h
        cppspectrum.h
        cppTypes.h
        )

set(CONFIG_FILE ${PROJECT_BINARY_DIR}/config.h)

set(EXEC_FILES test_sta)

foreach (execfile ${EXEC_FILES})
    add_executable(${execfile} ${execfile}.c ${SOURCE_FILES} ${CONFIG_FILE} )
    target_link_libraries(${execfile} m cfitsio)
    target_include_directories(${execfile} PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
endforeach (execfile ${EXEC_FILES})

set(LIBNAME Relxill)
add_library(${LIBNAME} ${SOURCE_FILES})
target_include_directories(${LIBNAME} PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(${LIBNAME} m cfitsio)


# create the lmodel.dat file
set(LMODEL_NAME lmodel_relxill.dat)

if (DEFINED ENV{RELXILL_STABLE})
    add_custom_command(
            OUTPUT ${LMODEL_NAME}
            COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/modelfiles/lmodel_relxill_public.dat ${CMAKE_CURRENT_SOURCE_DIR}/${LMODEL_NAME}
            COMMAND echo " ****** Creating lmodel_relxill.dat ONLY for STABLE RELEASE MODELS ******"
            USES_TERMINAL
    )

else()
    add_custom_command(
            OUTPUT ${LMODEL_NAME}
            COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/modelfiles/lmodel_relxill_public.dat ${CMAKE_CURRENT_SOURCE_DIR}/${LMODEL_NAME}
            COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/modelfiles/lmodel_relxill_devel.dat >> ${CMAKE_CURRENT_SOURCE_DIR}/${LMODEL_NAME}
            COMMAND echo " ****** Creating lmodel_relxill.dat for ALL MODELS ******"
            USES_TERMINAL
    )
endif()

add_custom_target(lmodeldat ALL DEPENDS ${LMODEL_NAME})


add_custom_command(
        OUTPUT xspec_wrapper_lmodels.cpp xspec_wrapper_lmodels.h
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_wrapper_xspec.py
        ${CMAKE_CURRENT_SOURCE_DIR}/lmodel_relxill.dat ${CMAKE_CURRENT_SOURCE_DIR}/xspec_wrapper_lmodels.cpp
        DEPENDS lmodel_relxill.dat
        )

add_custom_target(xspec_wrapper ALL DEPENDS xspec_wrapper_lmodels.cpp xspec_wrapper_lmodels.h)



install(FILES ${SOURCE_FILES} DESTINATION model)
install(FILES
        modelfiles/compile_relxill.sh
        ${LMODEL_NAME}
        DESTINATION model)

install(TARGETS ${EXEC_FILES} DESTINATION bin)


