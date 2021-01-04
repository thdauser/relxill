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



add_custom_command(
        OUTPUT lmodel_relxill.dat
        COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/modelfiles/lmodel_relxill.dat ${PROJECT_BINARY_DIR}/
        USES_TERMINAL
)
add_custom_target(lmodeldat ALL DEPENDS lmodel_relxill.dat)

install(FILES ${SOURCE_FILES} DESTINATION model)
install(FILES
        modelfiles/compile_relxill.sh
        DESTINATION model)

install(TARGETS ${EXEC_FILES} DESTINATION bin)


