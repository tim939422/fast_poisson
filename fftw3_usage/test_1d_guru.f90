program test_1d_guru
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none
    integer, parameter :: NX = 8, NY = 16, NZ = 4
    real(c_double), parameter :: eps = epsilon(1.0_C_DOUBLE)
    real(c_double), dimension(NX) :: in_x, out_x
    real(c_double), dimension(NY) :: in_y, out_y
    type(c_ptr) :: plan_x_1d, plan_y_1d
    integer :: iunit, i, j, k

    real(c_double), pointer, dimension(:, :, :) :: out
    type(c_ptr) :: out_ptr

    ! guru information
    integer(C_INT) :: rank, howmany_rank
    type(fftw_iodim), allocatable, dimension(:) :: dims, howmany_dims
    integer(C_FFTW_R2R_KIND), allocatable, dimension(:) :: kind_x, kind_y
    type(c_ptr) :: plan_x, plan_y

    ! check equivalence
    logical :: is_same
    



    ! Perform 1D transform to the 1D data (baseline)
    plan_x_1d = fftw_plan_r2r_1d(NX, in_x, out_x, FFTW_R2HC, FFTW_ESTIMATE)
    plan_y_1d = fftw_plan_r2r_1d(NY, in_y, out_y, FFTW_R2HC, FFTW_ESTIMATE)

    ! I put some number in out_x and out_y
    in_x = [1.0_C_DOUBLE, 2.3_C_DOUBLE, 1.4_C_DOUBLE, 4.0_C_DOUBLE,  &
            1.32_C_DOUBLE, 3.0_C_DOUBLE, 1.0_C_DOUBLE,  4.2_C_DOUBLE]
    in_y = [1.0_C_DOUBLE, 2.3_C_DOUBLE, 1.4_C_DOUBLE, 4.0_C_DOUBLE,  &
            1.32_C_DOUBLE, 3.0_C_DOUBLE, 1.0_C_DOUBLE,  4.2_C_DOUBLE, &
            2.45_C_DOUBLE, 0.32_C_DOUBLE, 3.4_C_DOUBLE, 0.22_C_DOUBLE, &
            1.45_C_DOUBLE, 0.98_C_DOUBLE, 2.23_C_DOUBLE, 1.02_C_DOUBLE]

    call fftw_execute_r2r(plan_x_1d, in_x, out_x)
    call fftw_execute_r2r(plan_y_1d, in_y, out_y)

    call fftw_destroy_plan(plan_x_1d)
    call fftw_destroy_plan(plan_y_1d)

    ! memory for out
    out_ptr = fftw_alloc_real(int(NX*NY*NZ, c_size_t))
    call c_f_pointer(out_ptr, out, [NX, NY, NZ])

    rank = 1
    howmany_rank = 2
    allocate(dims(rank), howmany_dims(howmany_rank))
    allocate(kind_x(rank), kind_y(rank))

    ! We first do the 1D transform along x
    ! The logical is as
    ! do k = 1, Nz
    !   do j = 1, Ny
    !     work tranform on u(:, j, k)
    !   end do
    ! end do
    dims(1)%n = NX
    dims(1)%is = 1 ! element of u(:, j, k) is continuous, a.k.a, stride = 1
    dims(1)%os = 1
    
    ! do j = 1, Ny
    howmany_dims(1)%n  = NY
    ! increment j by 1, the stride is NX
    howmany_dims(1)%is = NX 
    howmany_dims(1)%os = NX

    ! do k = 1, Nz
    howmany_dims(2)%n = NZ
    ! increment k by 1, the stride is NX*NY
    howmany_dims(2)%is = NX*NY
    howmany_dims(2)%os = NX*NY

    kind_x(1) = FFTW_R2HC
    plan_x = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, out, out, kind_x, FFTW_ESTIMATE)
    do k = 1, NZ
        do j = 1, NY
            out(:, j, k) = in_x(:)
        end do
    end do
    call fftw_execute_r2r(plan_x, out, out)
    open(newunit=iunit, file="test_1d_guru_x.txt")
    write(iunit, '(8E20.12)') out_x
    is_same = .true.
    do k = 1, NZ
        do j = 1, NY
            write(iunit, '(8E20.12)') out(:, j, k)
            if (norm2(out_x - out(:, j, k)) > eps) then
                is_same = .false.
                exit
            end if
        end do
    end do
    if (is_same) then
        print *, "1D in x-direction: PASS"
    else
        print *, "1D in x-direction: FAIL"
    end if
    close(iunit)
    call fftw_destroy_plan(plan_x)

    ! Now, we do the 1D transform along y
    ! The logical is as
    ! do k = 1, Nz
    !   do i = 1, Nx
    !     work transform on u(i, :, k)
    !   end do
    ! end do
    ! for u(i, :, k)
    ! u(i, j, k) and u(i, j + 1, k) is separated by NX
    dims(1)%n  = NY
    dims(1)%is = NX
    dims(1)%os = NX

    ! do i = 1, Nx
    howmany_dims(1)%n = NX
    ! increment i by 1, the stride is 1
    howmany_dims(1)%is = 1
    howmany_dims(1)%os = 1

    ! do k = 1, Nz
    howmany_dims(2)%n = NZ
    ! increment k by 1, the stride is NX*NY
    howmany_dims(2)%is = NX*NY
    howmany_dims(2)%os = NX*NY

    kind_y(1) = FFTW_R2HC
    plan_y = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, out, out, kind_y, FFTW_ESTIMATE)
    do k = 1, NZ
        do i = 1, NX
            out(i, :, k) = in_y(:)
        end do
    end do
    call fftw_execute_r2r(plan_y, out, out)
    open(newunit=iunit, file="test_1d_guru_y.txt")
    write(iunit, '(16E20.12)') out_y
    is_same = .true.
    do k = 1, NZ
        do i = 1, NX
            write(iunit, '(8E20.12)') out(i, :, k)
            if (norm2(out_y - out(i, :, k)) > eps) then
                is_same = .false.
                exit
            end if
        end do
    end do
    if (is_same) then
        print *, "1D in y-direction: PASS"
    else
        print *, "1D in y-direction: FAIL"
    end if
    close(iunit)
    call fftw_destroy_plan(plan_y)
    call fftw_free(out_ptr)
    deallocate(dims, howmany_dims)
    deallocate(kind_x, kind_y)
end program test_1d_guru