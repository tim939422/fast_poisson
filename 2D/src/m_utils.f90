module m_utils
    use, intrinsic :: iso_c_binding
    use m_case, only: TAYLOR_GREEN, STEADY_CHANNEL, CAVITY, DEVELOPING_CHANNEL
    use m_constants
    implicit none
contains
    subroutine compare_pressure_gradient(icase, nx, ny, xc, yc, xn, yn, p, dpdx, dpdy)
        integer, intent(in) :: icase, nx, ny
        real(c_double), intent(in) :: xc(0:), yc(0:)
        real(c_double), intent(in) :: xn(0:), yn(0:)
        real(c_double), intent(in) :: p(0:, 0:), dpdx(0:, 0:), dpdy(0:, 0:)

        real(c_double) :: work(0:nx + 1, 0:ny + 1)
        integer :: i, j
        real(c_double) :: factor
        real(c_double) :: relative_error
        real(c_double) :: mean, mean_ref

        select case (icase)
        case (TAYLOR_GREEN)
            do j = 1, ny
                do i = 1, nx
                    work(i, j) = sin(two*PI*xc(i))*sin(two*PI*yc(j))
                end do
            end do
        case (STEADY_CHANNEL)
            do j = 1, ny
                do i = 1, nx
                    work(i, j) = sin(two*PI*xc(i))*cos(PI*yc(j))
                end do
            end do
        case (CAVITY)
            do j = 1, ny
                do i = 1, nx
                    work(i, j) = cos(PI*xc(i))*cos(PI*yc(j))
                end do
            end do
        case (DEVELOPING_CHANNEL)
            do j = 1, ny
                do i = 1, nx
                    work(i, j) = cos(three*PI/two*xc(i))*cos(PI*yc(j))
                end do
            end do
        end select

        ! extract mean
        mean = zero; mean_ref = zero
        do j = 1, ny
            do i = 1, nx
                mean = mean + p(i, j)
                mean_ref = mean_ref + work(i, j)
            end do
        end do
        mean = mean/real(nx*ny, c_double)
        mean_ref = mean_ref/real(nx*ny, c_double)
        
        relative_error = norm2(work(1:nx, 1:ny) - mean_ref + mean - p(1:nx, 1:ny))/norm2(work(1:nx, 1:ny) - mean_ref)
        write(*, '("relative error in p", e20.12)') relative_error


        ! compare dpdx
        select case (icase)
        case (TAYLOR_GREEN)
            factor = two*PI
            do j = 1, ny
                do i = 0, nx
                    work(i, j) = factor*cos(two*PI*xn(i))*sin(two*PI*yc(j))
                end do
            end do
        case (STEADY_CHANNEL)
            factor = two*PI
            do j = 1, ny
                do i = 0, nx
                    work(i, j) = factor*cos(two*PI*xn(i))*cos(PI*yc(j))
                end do
            end do
        case (CAVITY)
            factor = -PI
            do j = 1, ny
                do i = 0, nx
                    work(i, j) = factor*sin(PI*xn(i))*cos(PI*yc(j))
                end do
            end do
        case (DEVELOPING_CHANNEL)
            factor = -three*PI/two
            do j = 1, ny
                do i = 0, nx
                    work(i, j) = factor*sin(three*PI/two*xn(i))*cos(PI*yc(j))
                end do
            end do
        end select

        relative_error = norm2(work(0:nx, 1:ny) - dpdx(0:nx, 1:ny))/norm2(work(0:nx, 1:ny))
        write(*, '("relative error in dp/dx", e20.12)') relative_error


        ! compare dpdy
        select case (icase)
        case (TAYLOR_GREEN)
            factor = two*PI
            do j = 0, ny
                do i = 1, nx
                    work(i, j) = factor*sin(two*PI*xc(i))*cos(two*PI*yn(j))
                end do
            end do
        case (STEADY_CHANNEL)
            factor = -PI
            do j = 0, ny
                do i = 1, nx
                    work(i, j) = factor*sin(two*PI*xc(i))*sin(PI*yn(j))
                end do
            end do
        case (CAVITY)
            factor = -PI
            do j = 0, ny
                do i = 1, nx
                    work(i, j) = factor*cos(PI*xc(i))*sin(PI*yn(j))
                end do
            end do
        case (DEVELOPING_CHANNEL)
            factor = -PI
            do j = 0, ny
                do i = 1, nx
                    work(i, j) = factor*cos(three*PI/two*xc(i))*sin(PI*yn(j))
                end do
            end do
        end select

        relative_error = norm2(work(1:nx, 0:ny) - dpdy(1:nx, 0:ny))/norm2(work(1:nx, 0:ny))
        write(*, '("relative error in dp/dy", e20.12)') relative_error

    end subroutine compare_pressure_gradient
end module m_utils