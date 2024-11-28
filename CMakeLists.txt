cmake_minimum_required(VERSION 3.21)                        
cmake_policy(SET CMP0074 NEW)
project(FastPoisson LANGUAGES Fortran)


include(cmake/fft.cmake)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/.modules)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(application)
add_subdirectory(testings)
