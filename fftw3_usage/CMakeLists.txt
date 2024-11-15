# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(FFTW3_Usage_Project Fortran)

# Find FFTW library
set(FFTW_ROOT "/Users/duosifan/Codes/fftw/3.3.10")
find_library(FFTW_LIB fftw3 HINTS "${FFTW_ROOT}/lib")

# Compiler flags
set(CMAKE_Fortran_FLAGS "-O2")

# Include FFTW3 directory
include_directories(${FFTW3_INCLUDE_DIR})

# Find all .f90 files
file(GLOB ALL_FORTRAN_FILES *.f90)

# Exclude fftw.f90 from the list
list(REMOVE_ITEM ALL_FORTRAN_FILES ${CMAKE_SOURCE_DIR}/fftw.f90)

# Compile fftw.f90 module
add_library(fftw_module OBJECT fftw.f90)
target_include_directories(fftw_module PRIVATE ${FFTW3_INCLUDE_DIR})

# Link FFTW3 library
link_directories(${FFTW3_LIB_DIR})

# Loop over each .f90 file to create an executable ending in .exe
foreach(SRC_FILE ${ALL_FORTRAN_FILES})
    # Get the filename without extension
    get_filename_component(EXEC_NAME ${SRC_FILE} NAME_WE)

    # Add .exe suffix
    set(EXEC_NAME_WITH_SUFFIX "${EXEC_NAME}.exe")

    # Add executable for each file
    add_executable(${EXEC_NAME_WITH_SUFFIX} ${SRC_FILE} fftw.f90)

    # Link FFTW3 library to each executable
    target_link_libraries(${EXEC_NAME_WITH_SUFFIX} PRIVATE ${FFTW_LIB})
endforeach()
