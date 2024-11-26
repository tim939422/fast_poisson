module m_constants
    use, intrinsic :: iso_c_binding
    implicit none
    real(C_DOUBLE), parameter :: one = 1.0_C_DOUBLE
    real(C_DOUBLE), parameter :: two = 2.0_C_DOUBLE
    real(C_DOUBLE), parameter :: three = 3.0_C_DOUBLE
    real(C_DOUBLE), parameter :: four = 4.0_C_DOUBLE
    real(C_DOUBLE), parameter :: PI = four*atan(one)
end module m_constants