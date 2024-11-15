module m_case
    use, intrinsic :: iso_c_binding
    use m_constants
    implicit none
    integer, parameter :: TAYLOR_GREEN   = 0
    integer, parameter :: STEADY_CHANNEL = 1
    integer, parameter :: CAVITY         = 2
    integer, parameter :: DEVELOPING_CHANNEL = 3

    character(len=*), dimension(0:3), parameter :: case_string = ['Taylor Green      ', &
                                                                   'Steady Channel    ', &
                                                                   'Cavity            ', &
                                                                   'Developing Channel']
contains
    subroutine init_case(icase, nx, ny, p)
        integer, intent(in) :: icase
        integer, intent(in) :: nx
        integer, intent(in) :: ny

        real(c_double), intent(out) :: p(0:, 0:)

        integer :: i, j
        i = nx; j = ny
        select case (icase)
        case (TAYLOR_GREEN)
            p(:, :) = zero
        case (STEADY_CHANNEL)
            p(:, :) = zero
        case (CAVITY)
            p(:, :) = zero
        case (DEVELOPING_CHANNEL)
            p(:, :) = zero
        end select
    end subroutine init_case
end module m_case