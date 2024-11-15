# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(PoissonSolver LANGUAGES Fortran)

# Set the output executable name
set(EXE_NAME poisson_1d.exe)

# Set source files
set(SRC_FILES
    main.f90
    m_poisson_solver.f90
    m_constants.f90
    m_tdma.f90
)

# Add the executable
add_executable(${EXE_NAME} ${SRC_FILES})

# Compiler flags
target_compile_options(${EXE_NAME} PRIVATE -O2)
