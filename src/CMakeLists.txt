# Specify the executable name
set(EXECUTABLE_NAME test_1d.exe)

# Add the executable target
add_executable(${EXECUTABLE_NAME}
    m_grid.f90
    m_laplacian.f90
    m_numerics.f90
    m_tdma.f90
    test_1d.f90
)

# Specify the install destination relative to the root
install(TARGETS ${EXECUTABLE_NAME}
    RUNTIME DESTINATION bin
)