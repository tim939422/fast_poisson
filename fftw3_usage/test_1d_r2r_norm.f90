program test_1d_r2r_norm
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none
    integer, parameter :: NX = 8, NY = 16
    real(c_double), dimension(NX) :: in_x
    real(c_double), pointer, dimension(:, :) :: out
    type(c_ptr) :: out_ptr

    

end program test_1d_r2r_norm