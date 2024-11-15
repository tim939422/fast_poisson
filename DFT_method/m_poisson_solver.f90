module m_poisson_solver
    use, intrinsic :: iso_c_binding
    use :: fftw3
    use :: m_constants
    implicit none
    private
    real(C_DOUBLE), pointer, dimension(:) :: work
    type(C_PTR) :: work_ptr
    real(C_DOUBLE) :: factor
    type(C_PTR)    :: forward, backward

    ! eigenvalues for Laplace
    real(C_DOUBLE), dimension(:), allocatable :: lambda_x


    ! Declare public procedures
    public :: init_poisson, solve_poisson, finalize_poisson

contains
    subroutine init_poisson(n, dx, BC)
        integer, intent(in) :: n ! Number of cells
        real(C_DOUBLE), intent(in) :: dx
        character(len=2), intent(in) :: BC ! Boundary condition string (e.g., 'PP', 'NN', 'ND')

        ! local variables
        integer(C_INT) :: fwd_kind, bwd_kind
        real(C_DOUBLE) :: invdx2 ! (1/dx^2)
        integer :: i

        ! allocate memory
        work_ptr = fftw_alloc_real(int(n, C_SIZE_T))
        call c_f_pointer(work_ptr, work, [n])
        allocate(lambda_x(n))

        ! select real to real DFT type
        select case (BC)
        case ('PP')
            fwd_kind = FFTW_R2HC
            bwd_kind = FFTW_HC2R
            factor = 1.0/real(n, C_DOUBLE)
        case ('NN')
            fwd_kind = FFTW_REDFT10
            bwd_kind = FFTW_REDFT01
            factor = 1.0/real(2*n, C_DOUBLE)
        case ('ND')
            fwd_kind = FFTW_REDFT11
            bwd_kind = FFTW_REDFT11
            factor = 1.0/real(2*n, C_DOUBLE)
        case default
            print *, "Error: Invalid BC: ", BC
        end select

        ! plan real to real DFT
        forward  = fftw_plan_r2r_1d(N, work, work, fwd_kind, FFTW_ESTIMATE)
        backward = fftw_plan_r2r_1d(N, work, work, bwd_kind, FFTW_ESTIMATE)

        ! eigenvalue for Laplace (FD2)
        call create_laplacian_dft(n, lambda_x, BC)
        invdx2 = one/dx**2
        do i = 1, N
            lambda_x(i) = lambda_x(i)*invdx2
        end do

    end subroutine

    ! in-place solve (cell center point)
    subroutine solve_poisson(n, p)
        integer, intent(in) :: n ! Number of cells
        real(C_DOUBLE), intent(inout), dimension(0:) :: p

        ! local
        integer :: i

        ! copy to work array
        work(:) = p(1:n)

        call fftw_execute_r2r(forward, work, work)
        work(1) = 0.0_C_DOUBLE
        do i = 2, N
            work(i) = work(i)/lambda_x(i)
        end do
        call fftw_execute_r2r(backward, work, work)

        ! copy it back
        p(1:N) = work(:)*factor
    end subroutine

    subroutine finalize_poisson()
        call fftw_free(work_ptr)
        deallocate(lambda_x)
        call fftw_destroy_plan(forward)
        call fftw_destroy_plan(backward)
    end subroutine finalize_poisson


    ! Private procedures
    ! FD2
    subroutine create_laplacian_dft(n, lambda, BC)
        integer, intent(in) :: n
        real(C_DOUBLE), intent(out), dimension(n) :: lambda
        character(len=2), intent(in) :: BC

        ! local variables
        integer :: i

        select case (BC)
        case ('PP')
            do i = 1, n
                lambda(i) = -two*(one - cos(two*PI*real(i - 1, C_DOUBLE)/real(N, C_DOUBLE)))
            end do
        case ('NN')
            do i = 1, n
                lambda(i) = -two*(one - cos(PI*real(i - 1, C_DOUBLE)/real(N, C_DOUBLE)))
            end do
        case ('ND')
            do i = 1, n
                lambda(i) = -two*(one - cos(PI*real(2*i - 1, C_DOUBLE)/real(2*N, C_DOUBLE)))
            end do
        case default
            print *, "Error: Invalid BC: ", BC
        end select

    end subroutine create_laplacian_dft
end module m_poisson_solver