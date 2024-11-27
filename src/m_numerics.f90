module m_numerics
    implicit none
    ! single, double precision reals
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))

    ! I will only use double precision
    integer, parameter :: rp = dp

    ! numerical constants with selected precision
    real(rp), parameter :: zero = 0.0_rp
    real(rp), parameter :: half = 0.5_rp
    real(rp), parameter :: one = 1.0_rp
    real(rp), parameter :: two = 2.0_rp
    real(rp), parameter :: four = 4.0_rp

    real(rp), parameter :: eps = epsilon(one)
    real(rp), parameter :: PI = four*atan(one)
end module m_numerics