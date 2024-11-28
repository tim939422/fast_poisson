module m_kinds
    !> Module setting the floating point precision (default double precision)
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan Nov 2024

    implicit none

    private
    
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))

    integer, parameter, public :: rp = dp
end module m_kinds
