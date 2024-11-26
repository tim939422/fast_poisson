# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(PoissonSolver2D LANGUAGES Fortran)

# Find FFTW library
set(FFTW_ROOT "/Users/duosifan/Codes/fftw/3.3.10")
find_library(FFTW_LIB fftw3 HINTS "${FFTW_ROOT}/lib")

add_subdirectory(src)