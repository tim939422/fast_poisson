module m_constants
    !> Module setting some mathematics and physics constants
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28

    use m_kinds, only: rp

    implicit none
    private

    real(rp), parameter, public :: PI  = 4.0_rp*atan(1.0_rp)
    real(rp), parameter, public :: eps = epsilon(1.0_rp)
end module m_constants
