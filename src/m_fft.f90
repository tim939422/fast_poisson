module m_fft
    !> module of wrapper procedures of FFTW 3 (isolate the interface AS MUCH AS POSSIBLE)
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28
    use, intrinsic :: iso_c_binding, only: c_ptr
    use m_kinds, only: rp
    use m_fftw3, only: FFTW_R2HC, FFTW_HC2R, FFTW_ESTIMATE
    use m_fftw3, only: fftw_plan_guru_r2r, fftw_execute_r2r, fftw_destroy_plan
    use m_fftw3, only: fftw_iodim
    implicit none

    private
    
    !> public available alias of FFTW 3 real to real kind (I use itype)
    integer, parameter, public :: DFT  = FFTW_R2HC
    integer, parameter, public :: IDFT = FFTW_HC2R

    !> private macros for FFT direction
    integer, parameter :: X = 0
    integer, parameter :: Y = 1
    integer, parameter :: Z = 2
    
    !> Declare public interface
    public :: create_r2r_2d, create_r2r, execute_fft_2d, execute_fft, destroy_plan
contains
    function create_r2r_2d(nx, ny, itype, dir) result(plan)
        !> function to create FFTW 3 plan of 1D real to real FFT of 2D array with guru interface
        !>
        !> note - everything is unnormalized
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! return type
        type(c_ptr) :: plan

        ! interface
        !> dimension in x
        integer, intent(in) :: nx
        !> dimension in y
        integer, intent(in) :: ny
        !> FFT type available: DFT, IDFT
        integer, intent(in) :: itype
        !> direction of FFT available: X-0, Y-1
        integer, intent(in) :: dir

        ! local
        real(rp) :: in(nx, ny)
        integer :: rank, howmany_rank
        type(fftw_iodim), allocatable :: dims(:), howmany_dims(:)

        ! work
        rank = 1
        howmany_rank = 1
        allocate(dims(rank), howmany_dims(howmany_rank))

        if (dir == X) then
            ! 1D FFT of size nx
            dims(1)%n  = nx
            ! stride of i increment: 1
            dims(1)%is = 1
            dims(1)%os = 1

            ! Perform 1D FFT ny times
            howmany_dims(1)%n = ny
            ! stride of j increment: nx
            howmany_dims(1)%is = nx
            howmany_dims(1)%os = nx

        else if (dir == Y) then
            ! 1D FFT of size ny
            dims(1)%n  = ny
            ! stride of j increment: nx
            dims(1)%is = nx
            dims(1)%os = nx

            ! Perform 1D FFT nx times
            howmany_dims(1)%n = nx
            ! stride of i increment: 1
            howmany_dims(1)%is = 1
            howmany_dims(1)%os = 1
        end if

        plan = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, in, in, [itype], FFTW_ESTIMATE)

        deallocate(dims, howmany_dims)
    end function create_r2r_2d

    function create_r2r(nx, ny, nz, itype, dir) result(plan)
        !> function to create FFTW 3 plan of 1D real to real FFT of 3D array with guru interface
        !>
        !> note - everything is unnormalized
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! return type
        type(c_ptr) :: plan

        ! interface
        !> dimension in x
        integer, intent(in) :: nx
        !> dimension in y
        integer, intent(in) :: ny
        !> dimension in z
        integer, intent(in) :: nz
        !> FFT type available: DFT, IDFT
        integer, intent(in) :: itype
        !> direction of FFT available: X-0, Y-1, Z-2
        integer, intent(in) :: dir

        ! local
        real(rp) :: in(nx, ny, nz)
        integer :: rank, howmany_rank
        type(fftw_iodim), allocatable :: dims(:), howmany_dims(:)

        ! work
        rank = 1
        howmany_rank = 2
        allocate(dims(rank), howmany_dims(howmany_rank))

        if (dir == X) then
            ! 1D FFT of size nx
            dims(1)%n  = nx
            ! stride of i increment: 1
            dims(1)%is = 1
            dims(1)%os = 1

            ! Perform 1D FFT ny*nz times
            howmany_dims(1)%n = ny
            ! stride of j increment: nx
            howmany_dims(1)%is = nx
            howmany_dims(1)%os = nx

            howmany_dims(2)%n = nz
            ! stride of k increment: nx*ny
            howmany_dims(2)%is = nx*ny
            howmany_dims(2)%os = nx*ny
        else if (dir == Y) then
            ! 1D FFT of size ny
            dims(1)%n  = ny
            ! stride of j increment: nx
            dims(1)%is = nx
            dims(1)%os = nx

            ! Perform 1D FFT nx*nz times
            howmany_dims(1)%n = nx
            ! stride of i increment: 1
            howmany_dims(1)%is = 1
            howmany_dims(1)%os = 1

            howmany_dims(2)%n = nz
            ! stride of k increment: nx*ny
            howmany_dims(2)%is = nx*ny
            howmany_dims(2)%os = nx*ny
        end if


        plan = fftw_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, in, in, [itype], FFTW_ESTIMATE)

        deallocate(dims, howmany_dims)
    end function create_r2r

    subroutine execute_fft_2d(plan, work)
        !> Wrapper for fftw_execute_r2r
        !>
        !> note - work should be contiguous and aligned
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        !> FFTW 3 plan
        type(c_ptr), intent(in) :: plan
        !> work array (in-place transform)
        real(rp), intent(inout) :: work(:, :)

        call fftw_execute_r2r(plan, work, work)

    end subroutine execute_fft_2d

    subroutine execute_fft(plan, work)
        !> Wrapper for fftw_execute_r2r
        !>
        !> note - work should be contiguous and aligned
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        !> FFTW 3 plan
        type(c_ptr), intent(in) :: plan
        !> work array (in-place transform)
        real(rp), intent(inout) :: work(:, :, :)

        call fftw_execute_r2r(plan, work, work)

    end subroutine execute_fft

    subroutine destroy_plan(plan)
        !> Wrapper for fftw_destroy_plan
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        type(c_ptr), intent(in) :: plan

        call fftw_destroy_plan(plan)

    end subroutine destroy_plan
end module m_fft
