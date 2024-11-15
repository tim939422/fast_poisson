! Solve 3D Poisson equation with Pseudo Spectral Method
program test_3d_r2c
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none
    real(c_double), parameter :: PI = 4.0_c_double*atan(1.0_c_double)
    real(c_double), parameter :: PI2 = 2.0_c_double*PI
    integer, parameter :: N = 128
    real(c_double), parameter :: dx = PI2/real(N, c_double)
    real(c_double), parameter :: factor = 1.0_c_double/real(N*N*N, c_double)
    type(c_ptr) :: fwd_plan, bwd_plan
    integer :: iunit, i, j, k

    ! reference data
    real(c_double), dimension(N) :: x, y, z
    real(c_double), dimension(N, N, N) :: p_ref
    real(c_double), dimension(N) :: akx, aky, akz
    real(c_double), dimension(N/2 + 1, N, N) :: laplace

    ! fftw stuff
    real(c_double), pointer :: p(:, :, :)
    complex(c_double_complex), pointer :: p_hat(:, :, :)
    type(c_ptr) :: work_r, work_c
    
    ! allocate meomery
    work_r = fftw_alloc_real(int(N*N*N, c_size_t))
    call c_f_pointer(work_r, p, [N, N, N])
    work_c = fftw_alloc_complex(int((N/2 + 1)*N*N, c_size_t))
    call c_f_pointer(work_c, p_hat, [N/2 + 1, N, N])

    ! plan DFT
    fwd_plan = fftw_plan_dft_r2c_3d(N, N, N, p, p_hat, FFTW_ESTIMATE)
    bwd_plan = fftw_plan_dft_c2r_3d(N, N, N, p_hat, p, FFTW_ESTIMATE)


    ! geometry
    do i = 1, N
        x(i) = real(i - 1, c_double)*dx; y(i) = x(i); z(i) = x(i)
    end do
    do k = 1, N
        do j = 1, N
            do i = 1, N
                p_ref(i, j, k) = (cos(x(i)) + cos(y(j)))*cos(z(k))
                ! RHS of nabla^2 p = f
                p(i, j, k) = -2.0_c_double*p_ref(i, j, k)
            end do
        end do
    end do

    call fftw_execute_dft_r2c(fwd_plan, p, p_hat)

    call fftfreq(n, akx); call fftfreq(n, aky); call fftfreq(n, akz)
    do k = 1, N
        do j = 1, N
            do i = 1, N/2 + 1
                if (i == 1 .and. j == 1 .and. k == 1) then
                    laplace(i, j, k) = 0.0_c_double
                else
                laplace(i, j, k) = -1.0_c_double/(akx(i)**2 + aky(j)**2 + akz(k)**2)
                end if
            end do
        end do
    end do

    do k = 1, N
        do j = 1, N
            do i = 1, N/2 + 1
                p_hat(i, j, k) = laplace(i, j, k)*p_hat(i, j, k)
            end do
        end do
    end do

    call fftw_execute_dft_c2r(bwd_plan, p_hat, p)
    p(:, :, :) = factor*p(:, :, :)

    open(newunit=iunit, file="test_3d_r2c.out")
    write(iunit, *) norm2(p - p_ref)/norm2(p_ref)
    close(iunit)

    call fftw_destroy_plan(fwd_plan)
    call fftw_destroy_plan(bwd_plan)
    call fftw_free(work_r)
    call fftw_free(work_c)
contains
    subroutine fftfreq(n, ak)
        integer, intent(in) :: n
        real(c_double), intent(out), dimension(:) :: ak

        integer :: l
        do l = 1, n/2
            ak(l) = real(l - 1, c_double)
        end do

        do l = n/2 + 1, n
            ak(l) = real(-(n - l + 1), c_double)
        end do
    end subroutine fftfreq
end program test_3d_r2c