file(GLOB EXE_SOURCES *.f90)

foreach(EXE_SOURCE ${EXE_SOURCES})
    get_filename_component(EXE_NAME ${EXE_SOURCE} NAME_WE)
    add_executable(${EXE_NAME}.exe ${EXE_SOURCE})
    target_link_libraries(${EXE_NAME}.exe PRIVATE FastPoisson)
    # Set the module output directory for Fortran modules
    set_target_properties(${EXE_NAME}.exe PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

    # Install the executable to bin/testings
    install(TARGETS ${EXE_NAME}.exe
        RUNTIME DESTINATION bin
    )
endforeach()
