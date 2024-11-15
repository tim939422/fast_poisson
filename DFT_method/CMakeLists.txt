# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(PoissonSolver LANGUAGES Fortran)

# Set the output executable name
set(EXE_NAME poisson_1d.exe)

# Set source files
set(SRC_FILES
    main.f90
    fftw.f90
    m_poisson_solver.f90
    m_constants.f90
)

# Find FFTW library
set(FFTW_ROOT "/Users/duosifan/Codes/fftw/3.3.10")
find_library(FFTW_LIB fftw3 HINTS "${FFTW_ROOT}/lib")
find_path(FFTW_INCLUDE_DIR fftw3.h HINTS "${FFTW_ROOT}/include")

# Check if FFTW was found
if(NOT FFTW_LIB OR NOT FFTW_INCLUDE_DIR)
    message(FATAL_ERROR "FFTW library not found in ${FFTW_ROOT}")
endif()

# Add include directories for FFTW
include_directories(${FFTW_INCLUDE_DIR})

# Add the executable
add_executable(${EXE_NAME} ${SRC_FILES})

# Link FFTW library
target_link_libraries(${EXE_NAME} PRIVATE ${FFTW_LIB})

# Compiler flags
target_compile_options(${EXE_NAME} PRIVATE -O2)
