# Set the source files
file(GLOB SRC_FILES *.f90)

# Add the executable target
add_executable(PoissonSolver ${SRC_FILES})

# Include directories for module files
target_include_directories(PoissonSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# Link the FFTW library
target_link_libraries(PoissonSolver PRIVATE ${FFTW_LIB})

# Set compiler flags for Fortran (adjust as needed)
target_compile_options(PoissonSolver PRIVATE -O3 -Wall -ffree-form)

# Ensure the FFTW include path is recognized
target_include_directories(PoissonSolver PRIVATE "${FFTW_ROOT}/include")

# Set the Fortran module output directory
set_target_properties(PoissonSolver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Specify the installation directory for the binary
install(TARGETS PoissonSolver DESTINATION ${CMAKE_SOURCE_DIR}/bin)