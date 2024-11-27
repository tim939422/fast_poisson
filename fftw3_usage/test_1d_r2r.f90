program test_1d_r2r
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none
    
    integer, parameter :: N = 8
    type(C_PTR)    :: plan
    real(C_DOUBLE), dimension(N) :: in
    real(C_DOUBLE), dimension(N) :: out
    integer :: i, iunit

    plan = fftw_plan_r2r_1d(N, in, out, FFTW_R2HC, FFTW_ESTIMATE)

    in = [1.0_C_DOUBLE,  &
          2.3_C_DOUBLE,  &
          1.4_C_DOUBLE,  &
          4.0_C_DOUBLE,  &
          1.32_C_DOUBLE, &
          3.0_C_DOUBLE,  &
          1.0_C_DOUBLE,  &
          4.2_C_DOUBLE]
    call fftw_execute_r2r(plan, in, out)


    ! Print the output

    open(newunit=iunit, file="test_1d_r2r.out")
    write(iunit, '(A)')  "# Output of FFT:"
    write(iunit, '(A)')  "#   k         out(i)"
    do i = 1, N
        write(iunit, '(i5, es20.12)') i, out(i)
    end do
    close(iunit)

    call fftw_destroy_plan(plan)
end program test_1d_r2r