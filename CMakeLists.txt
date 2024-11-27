cmake_minimum_required(VERSION 3.15)

# Project name and language
project(FastPoisson LANGUAGES Fortran)

# Set the FFTW root directory
set(FFTW_ROOT "/Users/duosifan/Codes/fftw/3.3.10")

# Add FFTW3 include and library paths
find_library(FFTW_LIB fftw3 HINTS "${FFTW_ROOT}/lib")
include_directories(${FFTW_INCLUDE_DIRS})
link_directories(${FFTW_LIBRARY_DIRS})

# Add the src directory as a subdirectory
add_subdirectory(src)

# Custom install directory inside the project
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})  # Root of the project