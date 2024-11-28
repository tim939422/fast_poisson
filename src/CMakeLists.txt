# Glob all Fortran source files in the current directory
file(GLOB FORTRAN_SOURCES "*.f90" "*.F90" "*.f")

# Add the library using the globbed Fortran files
add_library(FastPoisson STATIC ${FORTRAN_SOURCES})

# Set the module output directory for Fortran modules
set_target_properties(FastPoisson PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Link FFTW3 to the FastPoisson library
target_link_libraries(FastPoisson PRIVATE ${FFTW_DOUBLE_LIB})

# Include FFTW3 headers if necessary
target_include_directories(FastPoisson PRIVATE ${FFTW_INCLUDE_DIRS})

# Install the library into ${CMAKE_INSTALL_PREFIX}/lib
install(TARGETS FastPoisson
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include
)
