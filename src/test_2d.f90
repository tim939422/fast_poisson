program test_2d
    use m_numerics
    use m_grid
    use m_poisson_2d
    implicit none
    
    real(rp), parameter :: L = one
    integer :: nx, ny

    real(rp) :: beta
    ! grid
    real(rp), allocatable, dimension(:) :: xf, xc, yf, yc
    ! metric
    real(rp) :: dx
    real(rp), allocatable, dimension(:) :: dyf, dyc

    real(rp), allocatable, dimension(:, :) :: phi, dphidx, dphidy, work

    integer :: i, j
    real(rp) :: error

    print *, "Enter the grid dimensions (nx, ny):"
    read(*, *) nx, ny

    allocate(xf(0:nx+1), xc(0:nx+1), yf(0:ny+1), yc(0:ny+1))
    allocate(dyf(0:ny+1), dyc(0:ny+1))
    allocate(phi(0:nx+1, 0:ny+1), dphidx(0:nx+1, 0:ny+1), dphidy(0:nx+1, 0:ny+1), work(0:nx+1, 0:ny+1))

    ! face coordinate
    call linspace(nx, L, xf)
    xf(nx + 1) = two*xf(nx) - xf(nx - 1)
    beta = 1.2_rp
    call twoside_stretch_grid(ny, beta, L, yf)
    yf(ny + 1) = two*yf(ny) - yf(ny - 1)

    ! cell coordinate
    do i = 1, nx + 1
        xc(i) = half*(xf(i - 1) + xf(i))
    end do
    xc(0) = two*xf(0) - xc(1)
    do j = 1, ny + 1
        yc(j) = half*(yf(j - 1) + yf(j))
    end do
    yc(0) = two*yf(0) - yc(1)

    ! calculate metric
    dx = L/real(nx, rp)
    call nonuniform_metric(ny, yf, yc, dyf, dyc)

    call init_poisson(nx, ny, dx, dyf, dyc)
    do j = 1, ny
        do i = 1, nx
            phi(i, j) = -5.0_rp*PI**2*sin(two*PI*xc(i))*cos(PI*yc(j))
        end do
    end do
    call solve_poisson(nx, ny, phi)
    ! test dphi/dx
    do j = 1, ny
        do i = 0, nx
            work(i, j)   = two*PI*cos(two*PI*xf(i))*cos(PI*yc(j))
            dphidx(i, j) = (phi(i + 1, j) - phi(i, j))/dx
        end do
    end do
    print *, norm2(work(0:nx, 1:ny) - dphidx(0:nx, 1:ny))/norm2(work(0:nx, 1:ny))


    call finalize_poisson
    
end program test_2d