module poisson_2d
    use, intrinsic :: iso_c_binding
    use m_numerics, only: rp, zero, quarter, one, two, PI
    use m_fftw3
    implicit none
    private

    ! memory for FFTW3
    real(rp), pointer, dimension(:, :) :: work
    type(c_ptr) :: work_ptr

    ! normalize factor
    real(rp) :: norm

    ! plan
    type(c_ptr) :: forward_plan_x, backward_plan_x, forward_plan_y, backward_plan_y

    ! laplacian operator
    real(rp), dimension(:, :), allocatable :: inv_laplacian

    public :: init_poisson_2d, finalize_poisson_2d, solve_poisson_2d
contains
    !> Initialize a 2D Poisson solver with FD2 method (accelerate by DCT)
    !! 
    !! @param nx number of cells in x direction
    !! @param ny number of cells in y direction
    !! @param dx cell width in x direction
    !! @param dy cell width in y direction
    subroutine init_poisson_2d(nx, ny, dx, dy)
        integer, intent(in) :: nx, ny
        real(rp), intent(in) :: dx, dy
        
        ! local
        integer :: i, j
        real(rp), dimension(nx) :: lambda_x
        real(rp), dimension(ny) :: lamdda_y
        
        ! allocate FFTW3 memory
        work_ptr = fftw_alloc_real(int(nx*ny, c_size_t))
        call c_f_pointer(work_ptr, work, [nx, ny])

        ! normalize factor
        norm = quarter/real(nx*ny, rp)

        ! plan DCT
        forward_plan_x = create_r2r_2d(nx, ny, FFTW_REDFT10, 'x')
        forward_plan_y = create_r2r_2d(nx, ny, FFTW_REDFT10, 'y')
        backward_plan_x = create_r2r_2d(nx, ny, FFTW_REDFT01, 'x')
        backward_plan_y = create_r2r_2d(nx, ny, FFTW_REDFT01, 'y')

        allocate(inv_laplacian(nx, ny))
        do i = 1, nx
            lambda_x(i) = two*(cos(real(i - 1)*PI/real(nx, rp)) - one)
        end do
        do j = 1, ny
            lamdda_y(j) = two*(cos(real(j - 1)*PI/real(ny, rp)) - one)
        end do
        lambda_x(:) = lambda_x(:)/dx**2
        lamdda_y(:) = lamdda_y(:)/dy**2

        do j = 1, ny
            do i = 1, nx
                inv_laplacian(i, j) = lambda_x(i) + lamdda_y(j)
            end do
        end do
        inv_laplacian(1, 1) = one
        inv_laplacian(:, :) = one/inv_laplacian(:, :)

    end subroutine init_poisson_2d

    subroutine solve_poisson_2d(nx, ny, p)
        integer, intent(in) :: nx, ny
        real(c_double), intent(inout), dimension(0:, 0:) :: p

        work(:, :) = p(1:nx, 1:ny)
        call fftw_execute_r2r(forward_plan_x, work, work)
        call fftw_execute_r2r(forward_plan_y, work, work)

        work(1, 1) = zero
        work(:, :) = work(:, :)*inv_laplacian(:, :)
        call fftw_execute_r2r(backward_plan_y, work, work)
        call fftw_execute_r2r(backward_plan_x, work, work)

        p(1:ny, 1:ny) = norm*work(:, :)
    end subroutine solve_poisson_2d

    subroutine finalize_poisson_2d()
        call fftw_free(work_ptr)
        call fftw_destroy_plan(forward_plan_x)
        call fftw_destroy_plan(backward_plan_x)
        call fftw_destroy_plan(forward_plan_y)
        call fftw_destroy_plan(backward_plan_y)
        deallocate(inv_laplacian)
    end subroutine finalize_poisson_2d
    

    ! create plan of 1D kind type real-to-real DFT along dir of a 2D array
    type(c_ptr) function create_r2r_2d(nx, ny, kind, dir) result(plan)
        ! parameters
        integer(c_int), parameter :: rank = 1, howmany_rank = 1

        ! interface
        integer, intent(in) :: nx, ny
        integer(c_int), intent(in) :: kind
        character(len=1), intent(in) :: dir ! only 'x', 'y'

        ! local variables
        integer(C_FFTW_R2R_KIND) :: mykind(1)
        real(c_double), dimension(nx, ny) :: work
        type(fftw_iodim) :: dims(rank), howmany_dims(howmany_rank)

        ! determine the access pattern
        if (dir == 'x') then
            dims(1)%n  = nx
            ! element of work(:, j) is continuous, a.k.a., stride = 1
            dims(1)%is = 1
            dims(1)%os = 1

            ! do j = 1, Ny
            howmany_dims(1)%n = ny
            ! stride of j increment is nx
            howmany_dims(1)%is = nx
            howmany_dims(1)%os = nx
        else if (dir == 'y') then
            ! work(i, :) has a stride = nx
            dims(1)%n  = ny
            dims(1)%is = nx
            dims(1)%os = nx

            ! do i = 1, Nx
            howmany_dims(1)%n = nx
            ! stride of i increment is 1 (continuous)
            howmany_dims(1)%is = 1
            howmany_dims(1)%os = 1
        end if

        mykind(1) = kind
        plan = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, work, work, mykind, FFTW_ESTIMATE)

    end function create_r2r_2d
end module poisson_2d