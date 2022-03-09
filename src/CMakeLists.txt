### SOURCE FILES: add all model files here ###
set(SOURCE_FILES
        common.h
        donthcomp.c
        Relbase.cpp Relbase.h
        Relcache.cpp Relcache.h
        Rellp.cpp Rellp.h
        Relmodels.cpp Relmodels.h
        reltable.c reltable.h
        relutility.c relutility.h
        xilltable.c xilltable.h
        relphysics.c relphysics.h
        writeOutfiles.c writeOutfiles.h
        Relprofile.cpp Relprofile.h
        Relxill.cpp Relxill.h
        ModelDatabase.h
        LocalModel.cpp LocalModel.h
        ModelParams.cpp ModelParams.h
        XspecSpectrum.h
        ModelInfo.h
        Relreturn_BlackBody.cpp Relreturn_BlackBody.h
        IonGradient.cpp IonGradient.h
        Relreturn_Corona.cpp Relreturn_Corona.h
        relreturn_datastruct.c relreturn_datastruct.h
        relreturn_table.c relreturn_table.h
        xspec_wrapper_lmodels.cpp xspec_wrapper_lmodels.h   #are created by the wrapper script
        Xillspec.cpp Xillspec.h)
############################################

set(CONFIG_FILE ${PROJECT_BINARY_DIR}/config.h)
set(SOURCE_FILES ${SOURCE_FILES} ${CONFIG_FILE})

set(EXEC_FILES_CPP test_sta)

foreach (execfile ${EXEC_FILES_CPP})
    add_executable(${execfile} ${execfile}.cpp ${SOURCE_FILES} ${CONFIG_FILE} )
    target_link_libraries(${execfile} cfitsio fftw3 m)
    target_include_directories(${execfile} PUBLIC "${PROJECT_BINARY_DIR}" )  # necessary to find config file
endforeach (execfile ${EXEC_FILES_CPP})

### RELXILL LIBRARY ###
set(LIBNAME Relxill)
add_library(${LIBNAME} ${SOURCE_FILES} ${CONFIG_FILE})
target_include_directories(${LIBNAME} PUBLIC "${PROJECT_BINARY_DIR}")
target_link_libraries(${LIBNAME} cfitsio fftw3 m)
########################



# create everything in the real source directly (and not as default the build)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

### Target: LMODELDAT ###
set(LMODEL_NAME lmodel_relxill.dat)
if (DEFINED ENV{RELXILL_STABLE})
    add_custom_command(
            OUTPUT ${SRC_DIR}/${LMODEL_NAME}
            COMMAND cp ${SRC_DIR}/modelfiles/lmodel_relxill_public.dat ${SRC_DIR}/${LMODEL_NAME}
            COMMAND echo " ****** Creating lmodel_relxill.dat ONLY for STABLE RELEASE MODELS ******"
            DEPENDS ${SRC_DIR}/modelfiles/lmodel_relxill_public.dat
            USES_TERMINAL
    )

else()
    add_custom_command(
            OUTPUT ${SRC_DIR}/${LMODEL_NAME}
            COMMAND cp ${SRC_DIR}/modelfiles/lmodel_relxill_public.dat ${SRC_DIR}/${LMODEL_NAME}
            COMMAND cat ${SRC_DIR}/modelfiles/lmodel_relxill_devel.dat >> ${SRC_DIR}/${LMODEL_NAME}
            COMMAND echo " ****** Creating lmodel_relxill.dat for ALL MODELS ******"
            DEPENDS ${SRC_DIR}/modelfiles/lmodel_relxill_public.dat
            DEPENDS ${SRC_DIR}/modelfiles/lmodel_relxill_devel.dat
            USES_TERMINAL
    )
endif()

add_custom_target(lmodeldat ALL DEPENDS ${LMODEL_NAME})
#############################



### Target: XSPEC_WRAPPER ###
add_custom_command(
        OUTPUT ${SRC_DIR}/xspec_wrapper_lmodels.cpp ${SRC_DIR}/xspec_wrapper_lmodels.h
        COMMAND ${SRC_DIR}/create_wrapper_xspec.py
        ${SRC_DIR}/${LMODEL_NAME} ${SRC_DIR}/xspec_wrapper_lmodels.cpp
        COMMAND echo `pwd`
        DEPENDS lmodeldat ${SRC_DIR}/create_wrapper_xspec.py ${SRC_DIR}/${LMODEL_NAME}
        )

add_custom_target(xspec_wrapper ALL DEPENDS xspec_wrapper_lmodels.cpp xspec_wrapper_lmodels.h)
#############################




###  INSTALL  ###
install(FILES ${SOURCE_FILES} DESTINATION model)
install(FILES modelfiles/compile_relxill.sh ${LMODEL_NAME} DESTINATION model)

install(TARGETS ${EXEC_FILES_CPP} DESTINATION bin)




