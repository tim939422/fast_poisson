module m_RHS
    use, intrinsic :: iso_c_binding
    use m_constants
    use m_case, only: TAYLOR_GREEN, STEADY_CHANNEL, CAVITY, DEVELOPING_CHANNEL
    implicit none
    
contains
    ! some analytical RHS of Poisson equation
    subroutine fill_pressure_poisson_RHS_a(icase, nx, ny, xc, yc, RHS)
        integer, intent(in) :: icase, nx, ny
        real(c_double), intent(in) :: xc(0:), yc(0:)
        real(c_double), intent(out) :: RHS(0:, 0:)

        integer :: i, j
        real(c_double) :: factor

        select case (icase)
        case (TAYLOR_GREEN)
            factor = -eight*PI**2
            do j = 1, ny
                do i = 1, nx
                    RHS(i, j) = factor*sin(two*PI*xc(i))*sin(two*PI*yc(j))
                end do
            end do
        case (STEADY_CHANNEL)
            factor = -five*PI**2
            do j = 1, ny
                do i = 1, nx
                    RHS(i, j) = factor*sin(two*PI*xc(i))*cos(PI*yc(j))
                end do
            end do
        case (CAVITY)
            factor = -two*PI**2
            do j = 1, ny
                do i = 1, nx
                    RHS(i, j) = factor*cos(PI*xc(i))*cos(PI*yc(j))
                end do
            end do
        case (DEVELOPING_CHANNEL)
            factor = -13.0_c_double/four*PI**2
            do j = 1, ny
                do i = 1, nx
                    RHS(i, j) = factor*cos(three*PI/two*xc(i))*cos(PI*yc(j))
                end do
            end do
        end select
    end subroutine fill_pressure_poisson_RHS_a
end module m_RHS