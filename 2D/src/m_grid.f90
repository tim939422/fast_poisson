module m_grid
    use, intrinsic :: iso_c_binding
    use m_constants
    implicit none
contains
    ! domain start at (0, 0)
    subroutine create_regular_grid(nx, ny, Lx, Ly, dx, dy, xc, yc, xn, yn)
        integer, intent(in) :: nx, ny
        real(c_double), intent(in) :: Lx, Ly
        real(c_double), intent(out) :: dx, dy
        real(c_double), intent(out) :: xc(0:), yc(0:)
        real(c_double), intent(out) :: xn(0:), yn(0:)

        ! local
        integer :: i, j

        dx = Lx/real(nx, c_double); dy = Ly/real(ny, c_double)
        write(*, '("Creating regular grid")')
        write(*, '("     dx = ", e15.7, " dy = ", e15.7)') dx, dy
        do i = 0, nx + 1
            xn(i) = real(i, c_double)*dx
            xc(i) = xn(i) - half*dx
        end do

        do j = 0, ny + 1
            yn(j) = real(j, c_double)*dy
            yc(j) = yn(j) - half*dy
        end do

    end subroutine create_regular_grid
end module m_grid