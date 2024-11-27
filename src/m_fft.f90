module m_fft
    use, intrinsic :: iso_c_binding
    use m_fftw3
    implicit none
    
contains
    ! create plan of 1D kind type real-to-real DFT along dir of a 2D array
    type(c_ptr) function create_r2r_2d(nx, ny, kind, dir) result(plan)
        ! parameters
        integer(c_int), parameter :: rank = 1, howmany_rank = 1

        ! interface
        integer, intent(in) :: nx, ny
        integer(c_int), intent(in) :: kind
        character(len=1), intent(in) :: dir ! only 'x', 'y'

        ! local variables
        real(c_double), dimension(nx, ny) :: in, out
        type(fftw_iodim) :: dims(rank), howmany_dims(howmany_rank)

        ! determine the access pattern
        if (dir == 'x') then
            dims(1)%n  = nx
            ! element of work(:, j) is continuous, a.k.a., stride = 1
            dims(1)%is = 1
            dims(1)%os = 1

            ! do j = 1, Ny
            howmany_dims(1)%n = ny
            ! stride of j increment is nx
            howmany_dims(1)%is = nx
            howmany_dims(1)%os = nx
        else if (dir == 'y') then
            ! work(i, :) has a stride = nx
            dims(1)%n  = ny
            dims(1)%is = nx
            dims(1)%os = nx

            ! do i = 1, Nx
            howmany_dims(1)%n = nx
            ! stride of i increment is 1 (continuous)
            howmany_dims(1)%is = 1
            howmany_dims(1)%os = 1
        end if

        plan = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, &
                                  in, out, &
                                  [kind], FFTW_ESTIMATE)

    end function create_r2r_2d

    
end module m_fft