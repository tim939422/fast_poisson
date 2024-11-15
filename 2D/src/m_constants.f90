module m_constants
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), parameter :: zero = 0.0_c_double
    real(c_double), parameter :: half = 0.5_c_double
    real(c_double), parameter :: one = 1.0_c_double
    real(c_double), parameter :: two = 2.0_c_double
    real(c_double), parameter :: three = 3.0_c_double
    real(c_double), parameter :: four = 4.0_c_double
    real(c_double), parameter :: five = 5.0_c_double
    real(c_double), parameter :: eight = 8.0_c_double
    real(c_double), parameter :: PI = four*atan(one)

    real(C_DOUBLE), parameter :: eps = epsilon(one)
end module m_constants