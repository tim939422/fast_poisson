module m_poisson_solver
    use, intrinsic :: iso_c_binding
    use fftw3
    use m_fft
    use m_constants
    use m_tdma
    use m_boundary, only: DIRICHLET, NEUMANN, PERIODIC
    implicit none
    private
    real(c_double), pointer, dimension(:, :) :: work
    type(c_ptr) :: work_ptr

    ! tranform stuff
    ! normalize factor for backwards
    real(c_double) :: norm_factor
    ! laplacian operator d^2/dx^2
    real(c_double), dimension(:), allocatable :: laplacian_x
    ! plan for forward and backward transform
    type(c_ptr) :: fwd_x_plan, bwd_x_plan

    ! lower, diagonal, upper of d^2/dy^2
    real(c_double), dimension(:), allocatable :: a, b, c
    real(c_double), dimension(:), allocatable :: bb ! plus the DFT laplacian
    logical :: is_periodic_y

    

    ! Declare public procedures
    public :: init_poisson, solve_poisson, finalize_poisson
contains
    subroutine init_poisson(nx, ny, dx, dy, BC)
        integer, intent(in) :: nx, ny
        real(c_double), intent(in) :: dx, dy
        ! BC(1, 1) : left
        ! BC(2, 1) : right
        ! BC(1, 2) : bottom
        ! BC(2, 2) : top
        integer, intent(in), dimension(:, :) :: BC

        ! local variables
        integer(C_FFTW_R2R_KIND) :: fwd_x_kind, bwd_x_kind
        real(c_double) :: norm_factor_x

        ! prepare DFT in x direction
        ! operator
        allocate(laplacian_x(nx))
        call create_laplacian_dft(nx, dx, BC(:, 1), laplacian_x)
        ! FFT
        work_ptr = fftw_alloc_real(int(nx*ny, c_size_t))
        call c_f_pointer(work_ptr, work, [nx, ny])
        call search_fft(nx, BC(:, 1), fwd_x_kind, bwd_x_kind, norm_factor_x)
        norm_factor = norm_factor_x
        fwd_x_plan = create_r2r_1d(nx, ny, fwd_x_kind, dir='x')
        bwd_x_plan = create_r2r_1d(nx, ny, bwd_x_kind, dir='x')



        ! prepare matrix method in y direction
        if (all(BC(:, 2) .eq. [PERIODIC, PERIODIC])) then 
            is_periodic_y = .true.
        end if
        allocate(a(ny), b(ny), c(ny))
        allocate(bb(ny))
        call create_laplacian_matrix(ny, dy, BC(:, 2), a, b, c)
        
    end subroutine init_poisson

    subroutine solve_poisson(nx, ny, p)
        integer, intent(in) :: nx, ny
        real(c_double), intent(inout), dimension(0:, 0:) :: p

        ! local variables
        integer :: i
        ! copy to work array
        work(:, :) = p(1:nx, 1:ny)

        ! perform forward FFT
        call fftw_execute_r2r(fwd_x_plan, work, work)

        if (is_periodic_y) then
            do i = 1, nx
                bb(:) = b(:) + laplacian_x(i)
                call TDMA_p(ny, a, bb, c, work(i, :))
            end do
        else
            do i = 1, nx
                bb(:) = b(:) + laplacian_x(i)
                call TDMA(ny, a, bb, c, work(i, :))
            end do
        end if

        ! perform backward FFT
        call fftw_execute_r2r(bwd_x_plan, work, work)
        p(1:nx, 1:ny) = norm_factor*work(:, :)
    end subroutine solve_poisson

    subroutine finalize_poisson()
        deallocate(laplacian_x)
        call fftw_free(work_ptr)
        call fftw_destroy_plan(fwd_x_plan)
        call fftw_destroy_plan(bwd_x_plan)
        
        deallocate(a, b, c)
    end subroutine finalize_poisson

    

    ! Private stuff
    subroutine search_fft(n, BC, fwd_kind, bwd_kind, factor)
        integer, intent(in) :: n
        integer, intent(in) :: BC(:)
        integer(C_FFTW_R2R_KIND), intent(out) :: fwd_kind, bwd_kind
        real(c_double), intent(out) :: factor ! normalize factor

        ! local variables
        real(c_double) :: invn
        
        invn = one/real(n, c_double)

        if (all(BC .eq. [PERIODIC, PERIODIC])) then
            fwd_kind = FFTW_R2HC; bwd_kind = FFTW_HC2R
        else if (all(BC .eq. [NEUMANN, NEUMANN])) then
            fwd_kind = FFTW_REDFT10; bwd_kind = FFTW_REDFT01
        else if (all(BC .eq. [NEUMANN, DIRICHLET]) .or. all(BC .eq. [DIRICHLET, NEUMANN])) then
            fwd_kind = FFTW_REDFT11; bwd_kind = FFTW_REDFT11
        else if (all(BC .eq. [DIRICHLET, DIRICHLET])) then
            fwd_kind = FFTW_RODFT10; bwd_kind = FFTW_RODFT01
        end if

        if (all(BC .eq. [PERIODIC, PERIODIC])) then
            factor = invn
        else
            factor = half*invn
        end if

    end subroutine search_fft

    ! No global data
    ! laplacian operator d^2/dx^2 with FD2 eigenvalues
    ! FD2 eigenvalues of d^2/dx^2 operator with different BC
    subroutine create_laplacian_dft(n, dx, BC, laplacian)
        integer, intent(in) :: n
        real(c_double), intent(in) :: dx
        integer, intent(in) :: BC(:)
        real(c_double), intent(out), dimension(n) :: laplacian

        ! local variables
        integer :: i
        real(c_double) :: invn, invdx2

        invn   = one/real(n, c_double)


        if (all(BC .eq. [PERIODIC, PERIODIC])) then
            do i = 1, n
                laplacian(i) = -two*(one - cos(two*PI*real(i - 1, c_double)*invn))
            end do
        else if (all(BC .eq. [NEUMANN, NEUMANN])) then
            do i = 1, n
                laplacian(i) = -two*(one - cos(PI*real(i - 1, c_double)*invn))
            end do
        else if (all(BC .eq. [NEUMANN, DIRICHLET]) .or. all(BC .eq. [DIRICHLET, NEUMANN])) then
            do i = 1, n
                laplacian(i) = -two*(one - cos(half*PI*real(2*i - 1, c_double)*invn))
            end do
        else if (all(BC .eq. [DIRICHLET, DIRICHLET])) then
            do i = 1, n
                laplacian(i) = -two*(one - cos(PI*real(i, c_double)*invn))
            end do
        end if

        ! scale by 1/dx^2
        invdx2 = one/dx**2
        do i = 1, n
            laplacian(i) = laplacian(i)*invdx2
        end do
    end subroutine create_laplacian_dft

    ! tridiagonal matrix for FD2 d^2/dx^2 operator matrix with different BC
    ! uniform grid only
    subroutine create_laplacian_matrix(n, dx, BC, lower, diag, upper)
        integer, intent(in) :: n
        real(c_double), intent(in) :: dx
        integer, intent(in) :: BC(:)
        real(C_DOUBLE), intent(out), dimension(n) :: lower, diag, upper
        
        ! local variables
        integer :: i
        real(c_double) :: invdx2

        ! interior points
        lower(:) = one
        upper(:) = one
        diag(:)  = -(lower(:) + upper(:))

        ! Boundary conditions
        ! 
        select case (BC(1))
        case (DIRICHLET)
            diag(1) = diag(1) - lower(1)
        case (NEUMANN)
            diag(1) = diag(1) + lower(1)
        end select

        select case (BC(2))
        case (DIRICHLET)
            diag(n) = diag(n) - upper(n)
        case (NEUMANN)
            diag(n) = diag(n) + upper(n)
        end select
        

        ! scale by grid spacing
        invdx2 = one/dx**2
        do i = 1, n
            lower(i) = lower(i)*invdx2
            upper(i) = upper(i)*invdx2
            diag(i) = diag(i)*invdx2
        end do
    end subroutine create_laplacian_matrix

end module m_poisson_solver