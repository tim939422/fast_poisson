module m_poisson_solver
    use, intrinsic :: iso_c_binding
    use :: m_constants
    use :: m_tdma
    implicit none
    private
    real(C_DOUBLE), pointer, dimension(:) :: work

    ! eigenvalues for Laplace
    real(C_DOUBLE), dimension(:), allocatable :: a, b, c

    logical :: is_periodic

    ! Declare public procedures
    public :: init_poisson, solve_poisson, finalize_poisson

contains
    subroutine init_poisson(n, dx, BC)
        integer, intent(in) :: n ! Number of cells
        real(C_DOUBLE), intent(in) :: dx
        character(len=2), intent(in) :: BC ! Boundary condition string (e.g., 'PP', 'NN', 'ND')
        ! local
        real(C_DOUBLE) :: invdx2

        ! allocate memory
        allocate(work(n))
        allocate(a(n), b(n), c(n))

        is_periodic = .false.
        if ( BC == 'PP' ) then
            is_periodic = .true.
        end if

        call construct_tridiag_matrix(n, BC, a, b, c)
        invdx2 = one/dx**2
        a(:) = a(:)*invdx2
        b(:) = b(:)*invdx2
        c(:) = c(:)*invdx2

    end subroutine

    ! in-place solve (cell center point)
    subroutine solve_poisson(n, p)
        integer, intent(in) :: n ! Number of cells
        real(C_DOUBLE), intent(inout), dimension(0:) :: p

        ! copy to work array
        work(:) = p(1:n)
        
        if (is_periodic) then
            call TDMA_p(n, a, b, c, work)
        else
            call TDMA(n, a, b, c, work)
        end if

        p(1:n) = work(:)
    end subroutine

    subroutine finalize_poisson()
        deallocate(work, a, b, c)
    end subroutine finalize_poisson


    ! Private procedures
    ! FD2
    subroutine construct_tridiag_matrix(n, BC, lower, diag, upper)
        integer, intent(in) :: n
        character(len=2), intent(in) :: BC

        real(C_DOUBLE), intent(out), dimension(n) :: lower, diag, upper
        

        ! local variables
        lower(:) = one
        diag(:)  = -two
        upper(:) = one

        select case (BC)
        case ('PP')
            ! unchanged
            diag(1) = -two
            diag(n) = -two
        case ('NN')
            diag(1) = -one
            diag(n) = -one
        case ('ND')
            diag(1) = -one
            diag(n) = -three
        case default
            print *, "Error: Invalid BC: ", BC
        end select

    end subroutine construct_tridiag_matrix
end module m_poisson_solver