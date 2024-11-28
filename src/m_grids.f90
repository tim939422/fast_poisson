module m_grids

    !> module of Cartesian grid generator and metric calculation
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28

    use m_kinds, only: rp

    implicit none

    private

    public :: staggered_twoside_stretched, staggered_uniform
contains

    subroutine staggered_twoside_stretched(n, L, beta, xf, xc)
        !> staggered hyperbolic tangent twoside stretched grid on [0, L]
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> domain length
        real(rp), intent(in) :: L
        !> stretching factor
        real(rp), intent(in) :: beta
        !> grid face and cell coordinate
        real(rp), intent(out), dimension(0:) :: xf, xc

        ! local

        ! work
        ! face coordinate
        call twoside_stretched(n, L, beta, xf)
        xf(n + 1) = 2.0_rp*xf(n) - xf(n - 1)

        ! cell coordinate
        call face2cell(n, xf, xc)

    end subroutine staggered_twoside_stretched

    subroutine staggered_uniform(n, L, xf, xc)
        !> staggered uniform grid on [0, L]
        !>
        !> note -
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> domain length
        real(rp), intent(in) :: L
        !> grid face and cell coordinate
        real(rp), intent(out), dimension(0:) :: xf, xc

        ! local

        ! work
        ! face coordinate
        call uniform(n, L, xf)
        xf(n + 1) = 2.0_rp*xf(n) - xf(n - 1)

        call face2cell(n, xf, xc)
    end subroutine staggered_uniform

    ! Private

    subroutine twoside_stretched(n, L, beta, x)
        !> create hyperbolic tangent twoside stretched grid on [0, L]
        !>
        !> note - include endpoints
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> domain length
        real(rp), intent(in) :: L
        !> stretching factor
        real(rp), intent(in) :: beta
        !> grid coordinate (at least of length n + 1)
        real(rp), intent(out) :: x(0:)

        ! local
        integer :: i

        ! work
        do i = 0, n
            x(i) = 0.5_rp*L*(1.0_rp - tanh(beta*(1.0_rp - 2.0_rp*real(i, rp)/real(n, rp)))/tanh(beta))
        end do

    end subroutine twoside_stretched

    subroutine uniform(n, L, x)
        !> create uniform grid on [0, L]
        !>
        !> note - include endpoints
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> domain length
        real(rp), intent(in) :: L
        !> grid coordinate (at least of length n + 1)
        real(rp), intent(out) :: x(0:)

        ! local
        integer :: i
        real(rp) :: dx

        ! work
        dx = L/real(n, rp)

        do i = 0, n
            x(i) = real(i, rp)*dx
        end do

    end subroutine uniform

    subroutine face2cell(n, xf, xc)
        !> convert face coordinate to cell coordinate
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> face coordinate
        real(rp), intent(in) :: xf(0:)
        !> cell coordinate
        real(rp), intent(out) :: xc(0:)

        ! local
        integer :: i

        ! work
        do i = 1, n + 1
            xc(i) = 0.5_rp*(xf(i - 1) + xf(i))
        end do
        xc(0) = 2.0_rp*xf(0) - xc(1)

    end subroutine face2cell
end module m_grids