module m_metrics
    !> Module of 1D grid metrics evaluation
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28
    use m_kinds
    implicit none

    private

    public :: staggered_metric
    
contains
    subroutine staggered_metric(n, zf, zc, dzf, dzc)
        !> calculate grid metric dx/dxi by FD2
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        !> number of segments
        integer, intent(in)  :: n
        !> grid face and cell coordinate
        real(rp), intent(in), dimension(0:) :: zf, zc
        !> grid metric at face and cell coordinate
        real(rp), intent(out), dimension(0:) :: dzf, dzc

        ! local
        integer :: i

        ! work

        ! face metric
        do i = 0, n
            dzf(i) = zc(i + 1) - zc(i)
        end do
        dzf(n + 1) = dzf(n)

        ! cell metric
        do i = 1, n + 1
            dzc(i) = zf(i) - zf(i - 1)
        end do
        dzc(0) = dzc(1)

    end subroutine staggered_metric
end module m_metrics
