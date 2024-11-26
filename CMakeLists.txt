cmake_minimum_required(VERSION 3.15)

# Project name and language
project(FortranProject LANGUAGES Fortran)

# Add the src directory as a subdirectory
add_subdirectory(src)

# Custom install directory inside the project
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})  # Root of the project