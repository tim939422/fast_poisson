program demo_2d
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_rectilinear, only: t_rectilinear
    use m_poisson, only: t_poisson
    use m_gradients, only: t_gradient
    use m_io, only: array_write_3d

    implicit none
    real(rp), parameter :: Lx = 4.0_rp*PI, Ly = 2.0_rp
    real(rp), parameter :: beta = 1.2_rp
    integer :: nx, ny
    type(t_rectilinear) :: channel_grid
    type(t_poisson) :: potential_solver
    type(t_gradient) :: gradient
    real(rp), allocatable, dimension(:, :, :) :: phi, sol, ref
    integer :: i, j
    real(rp) :: relative_error

    print *, 'input N'
    read(*, *) nx
    ny = nx

    ! allocate memory for main program
    allocate(phi(0:nx + 1, 0:ny + 1, 0:1), sol(0:nx + 1, 0:ny + 1, 0:1), ref(0:nx + 1, 0:ny + 1, 0:1))

    ! initialize all objects
    call channel_grid%init([nx, ny], [Lx, Ly], beta)
    call potential_solver%init(channel_grid)
    call gradient%init(channel_grid)

    ! set source
    phi(:, :, :) = 0.0_rp
    associate(xc => channel_grid%xc, yc => channel_grid%yc)
        do j = 1, ny
            do i = 1, nx
                phi(i, j, 1) = -(1.0_rp + 0.25_rp*PI**2)*sin(xc(i))*cos(0.5_rp*PI*yc(j))
            end do
        end do
    end associate

    call potential_solver%solve_2d(phi)

    ! BC (periodic in x and Neumann in y)
    phi(0, 1:ny, 1) = phi(nx, 1:ny, 1)
    phi(nx + 1, 1:ny, 1) = phi(1, 1:ny, 1)
    phi(:, 0, 1) = phi(:, 1, 1)
    phi(:, ny + 1, 1) = phi(:, ny, 1)


    
    ! Verification
    ref(:, :, :) = 0.0_rp; sol(:, :, :) = 0.0_rp
    call gradient%gradpx(phi, sol)

    associate(xf => channel_grid%xf, yc => channel_grid%yc)
        do j = 1, ny
            do i = 0, nx
                ref(i, j, 1) = cos(xf(i))*cos(0.5*PI*yc(j))
            end do
        end do
    end associate
    relative_error = norm2(ref(0:nx, 1:ny, 1) - sol(0:nx, 1:ny, 1))/norm2(ref(0:nx, 1:ny, 1))
    write(*, '("Relative error in dphi/dx ", es23.15)') relative_error

    call gradient%gradpy_2d(phi, sol)
    associate(xc => channel_grid%xc, yf => channel_grid%yf)
        do j = 0, ny
            do i = 1, nx
                ref(i, j, 1) = -0.5_rp*PI*sin(xc(i))*sin(0.5_rp*PI*yf(j))
            end do
        end do
    end associate
    relative_error = norm2(ref(1:nx, 0:ny, 1) - sol(1:nx, 0:ny, 1))/norm2(ref(1:nx, 0:ny, 1))
    write(*, '("Relative error in dphi/dy ", es23.15)') relative_error

    
    call potential_solver%finalize
    call channel_grid%finalize()

    deallocate(phi, sol, ref)
end program demo_2d