module m_poisson_2d
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t, c_associated
    use m_fftw3
    use m_numerics
    use m_fft, only: create_r2r_2d
    use m_laplacian
    use m_tdma
    implicit none
    private

    ! memory for FFTW3
    real(rp), pointer, dimension(:, :) :: work
    type(c_ptr) :: work_ptr

    ! normalize factor
    real(rp) :: norm

    ! plan
    type(c_ptr) :: forward_plan_x, backward_plan_x

    ! laplacian operator
    real(rp), dimension(:), allocatable :: laplacian_x
    ! tridiagonal for laplacian in y
    real(rp), dimension(:), allocatable :: a, b, c, bb

    public :: init_poisson, solve_poisson, finalize_poisson
contains
    subroutine init_poisson(nx, ny, dx, dyf, dyc)
        integer, intent(in) :: nx, ny
        real(rp), intent(in) :: dx
        real(rp), intent(in), dimension(0:) :: dyf, dyc
        
        ! allocate memory
        work_ptr = fftw_alloc_real(int(nx*ny, c_size_t))
        call c_f_pointer(work_ptr, work, [nx, ny])
        allocate(laplacian_x(nx))
        allocate(a(ny), b(ny), c(ny), bb(ny))

        ! plan DFT
        forward_plan_x  = create_r2r_2d(nx, ny, FFTW_R2HC, 'x')
        backward_plan_x = create_r2r_2d(nx, ny, FFTW_HC2R, 'x')
        norm = one/(real(nx, rp))

        ! setup Laplacian operator
        call laplacian_dft(nx, dx, laplacian_x)
        call laplacian_1d_tridiagonal(ny, dyf, dyc, a, b, c)

    end subroutine init_poisson

    subroutine solve_poisson(nx, ny, phi)
        integer, intent(in) :: nx, ny
        real(rp), intent(inout), dimension(0:, 0:) :: phi

        ! local
        integer :: i

        work(:, :) = phi(1:nx, 1:ny)
        call fftw_execute_r2r(forward_plan_x, work, work)
        do i = 1, nx
            bb(:) = b(:) + laplacian_x(i)
            call tridag(ny, a, bb, c, work(i, :))
        end do
        call fftw_execute_r2r(backward_plan_x, work, work)

        phi(1:nx, 1:ny) = norm*work(:, :)

    end subroutine solve_poisson

    subroutine finalize_poisson
        ! destroy all c pointer if exists
        if (c_associated(work_ptr)) call fftw_free(work_ptr)
        if (c_associated(forward_plan_x)) call fftw_destroy_plan(forward_plan_x)
        if (c_associated(backward_plan_x)) call fftw_destroy_plan(backward_plan_x)

        ! destroy fortran array
        if (allocated(laplacian_x)) deallocate(laplacian_x)
        if (allocated(a)) deallocate(a)
        if (allocated(b)) deallocate(b)
        if (allocated(bb)) deallocate(bb)
        if (allocated(c)) deallocate(c)
    end subroutine finalize_poisson
end module m_poisson_2d