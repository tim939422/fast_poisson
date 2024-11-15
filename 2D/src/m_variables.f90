module m_variables
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=100) :: fname    

    ! grid
    real(c_double) :: dx, dy, dz
    real(c_double), allocatable, dimension(:) :: xc, yc ! cell center
    real(c_double), allocatable, dimension(:) :: xn, yn ! node center

    ! cell center variables
    real(c_double), allocatable, dimension(:, :) :: p

    ! edge center variables
    real(c_double), allocatable, dimension(:, :) :: dpdx, dpdy

    ! BC
    integer :: BC(2, 2)

contains
    subroutine init_variables(nx, ny)
        integer, intent(in) :: nx, ny

        write(*, '(A)') "Allocate memory for variables"
        allocate(xc(0:nx + 1), yc(0:ny + 1))
        allocate(xn(0:nx + 1), yn(0:ny + 1))

        allocate(p(0:nx + 1, 0:ny + 1))

        allocate(dpdx(0:nx + 1, 0:ny + 1))
        allocate(dpdy(0:nx + 1, 0:ny + 1))
        
    end subroutine init_variables

    subroutine finalize_variables()
        deallocate(xc, yc)
        deallocate(xn, yn)

        deallocate(p)

        deallocate(dpdx, dpdy)
    end subroutine finalize_variables
end module m_variables