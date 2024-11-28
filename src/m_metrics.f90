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
    subroutine staggered_metric(n, xf, xc, dxf, dxc)
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
        real(rp), intent(in), dimension(0:) :: xf, xc
        !> grid metric at face and cell coordinate
        real(rp), intent(out), dimension(0:) :: dxf, dxc

        ! local
        integer :: i

        ! work

        ! face metric
        do i = 0, n
            dxf(i) = xc(i + 1) - xc(i)
        end do
        dxf(n + 1) = dxf(n)

        ! cell metric
        do i = 1, n + 1
            dxc(i) = xf(i) - xf(i - 1)
        end do
        dxc(0) = dxc(1)

    end subroutine staggered_metric
end module m_metrics
