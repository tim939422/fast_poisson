module m_numerics
    use, intrinsic :: iso_c_binding
    implicit none

    integer, parameter :: rp = c_double
    real(rp), parameter :: zero = 0.0_rp
    real(rp), parameter :: quarter = 0.25_rp
    real(rp), parameter :: third = 1.0_rp/3.0_rp
    real(rp), parameter :: half = 0.5_rp
    real(rp), parameter :: one = 1.0_rp
    real(rp), parameter :: two = 2.0_rp
    real(rp), parameter :: three = 3.0_rp
    real(rp), parameter :: four = 4.0_rp
    real(rp), parameter :: five = 5.0_rp
    real(rp), parameter :: eight = 8.0_rp

    real(rp), parameter :: PI = four*atan(one)

    real(rp), parameter :: eps = epsilon(one)
    real(rp), parameter :: large  = huge(one) - eight
    real(rp), parameter :: mlarge = -large

end module m_numerics