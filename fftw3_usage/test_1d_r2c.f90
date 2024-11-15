program test_1d_r2c
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none
    
    integer, parameter :: N = 8
    integer, parameter :: NC = N/2 + 1
    type(C_PTR)    :: plan
    real(C_DOUBLE), dimension(N) :: in
    complex(C_DOUBLE_COMPLEX), dimension(NC) :: out
    integer :: i, iunit

    plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE)

    in = [1.0_C_DOUBLE,  &
          2.3_C_DOUBLE,  &
          1.4_C_DOUBLE,  &
          4.0_C_DOUBLE,  &
          1.32_C_DOUBLE, &
          3.0_C_DOUBLE,  &
          1.0_C_DOUBLE,  &
          4.2_C_DOUBLE]
    call fftw_execute_dft_r2c(plan, in, out)


    ! Print the output

    open(newunit=iunit, file="test_1d_r2c.out")
    write(iunit, '(A)')  "# Output of FFT:"
    write(iunit, '(A)')  "#   i      out(i)%Re            out(i)%Im"
    do i = 1, NC
        write(iunit, '(i5, 2E20.12)') i, out(i)%re, out(i)%im
    end do
    close(iunit)

    call fftw_destroy_plan(plan)

end program test_1d_r2c