program test_m_fft_2d
    use, intrinsic :: iso_c_binding, only: c_ptr
    use m_kinds, only: rp
    use m_fft, only: DFT, IDFT, create_r2r_2d, execute_fft_2d, destroy_plan
    implicit none
    integer :: nx, ny
    real(rp), dimension(:, :), allocatable :: in, work

    type(c_ptr) :: plans(4)
    integer :: itypes(4), dirs(4)
    character(len=12) :: fnames(4)
    integer :: i, iunit
    real(rp) :: factor ! FFT normalize factor

    nx = 128
    ny = 64
    allocate(in(nx, ny), work(nx, ny))

    open(newunit=iunit, file='data.bin', access='stream', form='unformatted')
    read(iunit) in
    close(iunit)

    ! Test the transform result with Python Implementation
    fnames = ['sol/fwdx.bin', 'sol/bwdx.bin', 'sol/fwdy.bin', 'sol/bwdy.bin']
    itypes = [DFT, IDFT, DFT, IDFT]
    dirs = [0, 0, 1, 1]
    do i = 1, 4
        plans(i) = create_r2r_2d(nx, ny, itypes(i), dirs(i))
        work(:, :) = in(:, :)
        call execute_fft_2d(plans(i), work)

        open(newunit=iunit, file=fnames(i), access='stream', form='unformatted')
        write(iunit) work
        close(iunit)
    end do

    ! Test transform in x
    factor = 1.0_rp/real(nx, rp)
    work(:, :) = in(:, :)
    call execute_fft_2d(plans(1), work)
    call execute_fft_2d(plans(2), work)
    work(:, :) = factor*work(:, :)
    print *, norm2(in - work)

    ! Test transform in y
    factor = 1.0_rp/real(ny, rp)
    work(:, :) = in(:, :)
    call execute_fft_2d(plans(3), work)
    call execute_fft_2d(plans(4), work)
    work(:, :) = factor*work(:, :)
    print *, norm2(in - work)

    ! Destroy all plan
    do i = 1, 4
        call destroy_plan(plans(i))
    end do

    deallocate(in, work)
end program test_m_fft_2d