program test_1d
    use m_numerics
    use m_grid, only: twoside_stretch_grid, nonuniform_metric
    use m_laplacian
    use m_tdma
    implicit none
    real(rp), parameter :: L = two
    integer  :: nx
    real(rp) :: beta
    real(rp), allocatable, dimension(:) :: xf, xc, dxf, dxc ! grid and metric
    real(rp), allocatable, dimension(:) :: a, b, c ! tridiagonal coeff

    real(rp), allocatable, dimension(:) :: phi, dphidx, dphidx_a

    integer :: i, iunit
    character(len=20) :: fmt  ! Format string for printing
    real(rp) :: error

    nx = 32768
    beta = 1.2_rp
    ! allocate memory
    allocate(xf(0:nx+1), xc(0:nx+1), dxf(0:nx+1), dxc(0:nx+1))
    allocate(a(nx), b(nx), c(nx))
    allocate(phi(0:nx+1), dphidx(0:nx+1), dphidx_a(0:nx+1))

    ! create coordinate system
    ! face coordinate
    call twoside_stretch_grid(nx, beta, L, xf)
    xf(nx + 1) = two*xf(nx) - xf(nx - 1)
    ! cell coordinate
    do i = 1, nx + 1
        xc(i) = half*(xf(i - 1) + xf(i))
    end do
    xc(0) = two*xf(0) - xc(1)

    ! calculate metric
    call nonuniform_metric(nx, xf, xc, dxf, dxc)

    ! write grid and metric
    fmt = '(i5, 4(es23.15))'
    open(newunit=iunit, file="grid.dat", status="replace", action="write")
    write(iunit, '(a)') '# Index         xf                    xc                    dxf                   dxc'
    do i = 0, nx + 1
        write(iunit, fmt) i, xf(i), xc(i), dxf(i), dxc(i)
    end do
    close(iunit)

    ! build tridiagonal matrix for 1D laplacian
    call laplacian_1d_tridiagonal(nx, dxf, dxc, a, b, c)
    fmt = '(i5, 3(es23.15))'
    open(newunit=iunit, file="tridiagonal.dat", status="replace", action="write")
    write(iunit, '(a)') '# Index          a                      b                      c'
    do i = 1, nx
        write(iunit, fmt) i, a(i), b(i), c(i)
    end do
    close(iunit)

    do i = 1, nx
        phi(i) = -(half*PI)**2*cos(half*PI*xc(i))
    end do
    call tridag(nx, a, b, c, phi(1:nx))
    
    ! BC
    phi(0) = phi(1)
    phi(nx + 1) = phi(nx)

    ! calculate dphi/dx at face
    do i = 0, nx
        dphidx(i) = (phi(i + 1) - phi(i))/dxf(i)
        dphidx_a(i) = -half*PI*sin(half*PI*xf(i))
    end do

    error = norm2(dphidx(0:nx) - dphidx_a(0:nx))/norm2(dphidx_a(0:nx))
    print *, error

    ! Deallocate all allocated arrays
    if (allocated(xf)) deallocate(xf)
    if (allocated(xc)) deallocate(xc)
    if (allocated(dxf)) deallocate(dxf)
    if (allocated(dxc)) deallocate(dxc)
    if (allocated(a)) deallocate(a)
    if (allocated(b)) deallocate(b)
    if (allocated(c)) deallocate(c)
    if (allocated(phi)) deallocate(phi)
    if (allocated(dphidx)) deallocate(dphidx)
    if (allocated(dphidx_a)) deallocate(dphidx_a)

end program test_1d