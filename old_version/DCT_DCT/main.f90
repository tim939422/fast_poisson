program main
    use m_numerics, only: rp, third, half, one, two, PI
    use poisson_2d
    implicit none
    integer,  parameter  :: Nx = 8, Ny = 8
    real(rp), parameter  :: dx = one/real(Nx, rp), dy = one/real(Ny, rp)
    
    real(rp), dimension(0:Nx + 1) :: x
    real(rp), dimension(0:Ny + 1) :: y
    real(rp), dimension(0:Nx + 1, 0:Ny + 1) :: f, p

    integer :: i, j, iunit

    do i = 1, Nx
        x(i) = real(i - 1, rp)*dx
    end do
    x(:) = x(:) + half*dx

    do j = 1, Ny
        y(j) = real(j - 1, rp)*dy
    end do
    y(:) = y(:) + half*dx

    do j = 1, Ny
        do i = 1, Nx
            f(i, j) = (half*x(i)**2 - third*x(i)**3)*(half*y(j)**2 - third*y(j)**3)
        end do
    end do

    p(:, :) = f(:, :)
    call init_poisson_2d(nx, ny, dx, dy)
    call solve_poisson_2d(nx, ny, p)
    call finalize_poisson_2d

    open(newunit=iunit, file="p.bin", access="stream")
    write(iunit) p(1:Nx, 1:Ny)
    close(iunit)
end program main