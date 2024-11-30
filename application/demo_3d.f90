program demo_2d
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_rectilinear, only: t_rectilinear
    use m_poisson, only: t_poisson
    use m_gradients, only: t_gradient

    implicit none
    real(rp), parameter :: Lx = 4.0_rp*PI, Ly = 2.0_rp*PI, Lz = 2.0_rp 
    real(rp), parameter :: beta = 1.2_rp
    integer :: nx, ny, nz
    type(t_rectilinear) :: channel_grid
    type(t_poisson) :: potential_solver
    type(t_gradient) :: gradient
    real(rp), allocatable, dimension(:, :, :) :: phi, sol, ref
    integer :: i, j, k
    real(rp) :: relative_error

    print *, 'input N'
    read(*, *) nx
    ny = nx
    nz = nx

    ! allocate memory for main program
    allocate(phi(0:nx + 1, 0:ny + 1, 0:nz + 1), sol(0:nx + 1, 0:ny + 1, 0:nz + 1 ), ref(0:nx + 1, 0:ny + 1, 0:nz +1))

    ! initialize all objects
    call channel_grid%init([nx, ny, nz], [Lx, Ly, Lz], beta)
    

    call potential_solver%init(channel_grid)
    call gradient%init(channel_grid)

    ! set source
    phi(:, :, :) = 0.0_rp
    associate(xc => channel_grid%xc, yc => channel_grid%yc, zc => channel_grid%zc)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    phi(i, j, k) = -(2.0_rp + 0.25_rp*PI**2)*sin(xc(i))*sin(yc(j))*cos(0.5_rp*PI*zc(k))
                end do
            end do
        end do
    end associate


    call potential_solver%solve_3d(phi)

    ! BC
    ! periodic in x
    phi(0, 1:ny, 1:nz) = phi(nx, 1:ny, 1:nz)
    phi(nx + 1, 1:ny, 1:nz) = phi(1, 1:ny, 1:nz)
    ! periodic in y
    phi(:, 0, 1:nz) = phi(:, ny, 1:nz)
    phi(:, ny + 1, 1:nz) = phi(:, 1, 1:nz)
    ! Neumann in z
    phi(:, :, 0) = phi(:, :, 1)
    phi(:, :, nz + 1) = phi(:, :, nz)
    
    ! Verification
    call gradient%gradpx(phi, sol)
    associate(xf => channel_grid%xf, yc => channel_grid%yc, zc => channel_grid%zc)
        do k = 1, nz
            do j = 1, ny
                do i = 0, nx
                    ref(i, j, k) = cos(xf(i))*sin(yc(j))*cos(0.5_rp*PI*zc(k))
                end do
            end do 
        end do
    end associate
    relative_error = norm2(ref(0:nx, 1:ny, 1:nz) - sol(0:nx, 1:ny, 1:nz))/norm2(ref(0:nx, 1:ny, 1:nz))
    write(*, '("Relative error in dphi/dx ", es23.15)') relative_error

    call gradient%gradpy_3d(phi, sol)
    associate(xc => channel_grid%xc, yf => channel_grid%yf, zc => channel_grid%zc)
        do k = 1, nz
            do j = 0, ny
                do i = 1, nx
                ref(i, j, k) = sin(xc(i))*cos(yf(j))*cos(0.5_rp*PI*zc(k))
                end do
            end do
        end do
    end associate
    relative_error = norm2(ref(1:nx, 0:ny, 1:nz) - sol(1:nx, 0:ny, 1:nz))/norm2(ref(1:nx, 0:ny, 1:nz))
    write(*, '("Relative error in dphi/dy ", es23.15)') relative_error


    call gradient%gradpz(phi, sol)
    associate(xc => channel_grid%xc, yc => channel_grid%yc, zf => channel_grid%zf)
        do k = 0, nz
            do j = 1, ny
                do i = 1, nx
                ref(i, j, k) = -0.5_rp*PI*sin(xc(i))*sin(yc(j))*sin(0.5_rp*PI*zf(k))
                end do
            end do
        end do
    end associate
    relative_error = norm2(ref(1:nx, 1:ny, 0:nz) - sol(1:nx, 1:ny, 0:nz))/norm2(ref(1:nx, 1:ny, 0:nz))
    write(*, '("Relative error in dphi/dz ", es23.15)') relative_error



    
    call potential_solver%finalize()
    call channel_grid%finalize()
    call gradient%finalize()

    deallocate(phi, sol, ref)
end program demo_2d