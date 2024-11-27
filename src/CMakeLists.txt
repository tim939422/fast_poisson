# Add the executable target
add_executable(test_1d.exe
    m_grid.f90
    m_laplacian.f90
    m_numerics.f90
    m_tdma.f90
    test_1d.f90
)
target_link_libraries(test_1d.exe PRIVATE ${FFTW_LIB})

add_executable(test_2d.exe
    m_numerics.f90
    m_grid.f90
    m_poisson_2d.f90
    m_fft.f90
    m_laplacian.f90
    m_fftw3.f90
    m_tdma.f90
    test_2d.f90
)
target_link_libraries(test_2d.exe PRIVATE ${FFTW_LIB})


set(EXECUTABLES test_1d.exe test_2d.exe)
# Specify the install destination relative to the root
install(TARGETS ${EXECUTABLES}
    RUNTIME DESTINATION bin
)