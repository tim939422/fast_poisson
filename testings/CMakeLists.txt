# Glob all test source files in the testings directory
file(GLOB TEST_SOURCES test_*.f90)

# Loop through all test source files
foreach(TEST_SOURCE ${TEST_SOURCES})
    # Get the executable name (remove directory and file extension)
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Define the executable target
    add_executable(${TEST_NAME}.exe ${TEST_SOURCE})

    # Link the FastPoisson library
    target_link_libraries(${TEST_NAME}.exe PRIVATE FastPoisson)

    # Set the module output directory for Fortran modules
    set_target_properties(${TEST_NAME}.exe PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

    # Install the executable to bin/testings
    install(TARGETS ${TEST_NAME}.exe
        RUNTIME DESTINATION bin/testings
    )
endforeach()
