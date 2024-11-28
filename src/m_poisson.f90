module m_poisson
    !> module of FFT acceralated Poisson Solver Class
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28

    use m_kinds, only: rp
    use m_rectilinear, only: t_rectilinear_2d
    use m_laplacian, only: fft_laplacian, matrix_laplacian
    use m_fft, only: DFT, IDFT
    use m_fft, only: create_r2r_2d, execute_fft_2d
    use m_fft, only: create_r2r_3d, execute_fft_3d
    use m_fft, only: destroy_plan
    use, intrinsic :: iso_c_binding, only: c_ptr
    use m_tdma, only: tridag

    implicit none

    private

    type, public :: t_poisson_2d
        private
        !> dimension in x (# of cells)
        integer, public :: nx
        !> dimension in y
        integer, public :: ny
        !> laplacian operator in x
        real(rp), allocatable :: laplacian_x(:)
        !> tridiagonal matrix in y (bb is the b + laplacian_x(i))
        real(rp), allocatable :: a(:), b(:), c(:), bb(:)
        !> work array for FFT
        real(rp), allocatable :: work(:, :)
        !> forward FFT plan
        type(c_ptr) :: forward
        !> backward FFT plan
        type(c_ptr) :: backward
        !> FFT normalization factor
        real(rp) :: factor
    contains
        private
        procedure, public :: init => init_poisson_2d
        procedure, public :: solve => solve_poisson_2d
        procedure :: alloc => allocate_poisson_2d
        procedure, public :: finalize => finalize_poisson_2d
    end type t_poisson_2d

    

contains

    subroutine init_poisson_2d(self, grid)
        ! interface
        class(t_poisson_2d) :: self
        type(t_rectilinear_2d) :: grid

        ! local
        integer  :: nx, ny
        real(rp) :: dx

        ! work
        nx = grid%nx; ny = grid%ny; dx = grid%dx
        self%nx = nx; self%ny = ny
        call self%alloc

        ! setup operator
        call fft_laplacian(nx, dx, self%laplacian_x)
        call matrix_laplacian(ny, grid%dyf, grid%dyc, self%a, self%b, self%c)

        ! plan FFT
        self%forward = create_r2r_2d(nx, ny, DFT, 0)
        self%backward = create_r2r_2d(nx, ny, IDFT, 0)
        self%factor = 1.0_rp/real(nx, rp)

    end subroutine init_poisson_2d

    subroutine solve_poisson_2d(self, phi)
        ! interface
        class(t_poisson_2d) :: self
        real(rp), dimension(0:, 0:) :: phi

        ! local
        integer :: i

        ! work
        associate(nx => self%nx, ny=>self%ny, work=>self%work, &
                  forward=>self%forward, backward=>self%backward, &
                  a=>self%a, b=>self%b, c=>self%c, bb=>self%bb, &
                  laplacian_x=>self%laplacian_x)
            
            ! Copy to work array
            work(:, :) = phi(1:nx, 1:ny)

            call execute_fft_2d(forward, work)

            ! now we have
            ! (\lambda_x)_i/dx^2 \tilde{\phi}_{i, j} + a_j* ..

            do i = 1, nx
                bb(:) = b(:) + laplacian_x(i)
                call tridag(a, bb, c, work(i, :), ny)
            end do

            call execute_fft_2d(backward, work)

            ! Copy back to original array
            phi(1:nx, 1:ny) = self%factor*work(:, :)

        end associate

    end subroutine solve_poisson_2d

    subroutine allocate_poisson_2d(self)
        ! interface
        class(t_poisson_2d) :: self

        associate(nx => self%nx, ny=>self%ny)
            allocate(self%laplacian_x(nx))
            allocate(self%a(ny), self%b(ny), self%c(ny), self%bb(ny))
            allocate(self%work(nx, ny))
        end associate
        
        print *, "t_poisson_2d resource allocated"

    end subroutine allocate_poisson_2d

    subroutine finalize_poisson_2d(self)
        ! interface
        class(t_poisson_2d) :: self

        ! Deallocate laplacian_x
        if (allocated(self%laplacian_x)) then
            deallocate(self%laplacian_x)
        endif

        ! Deallocate a, b, c
        if (allocated(self%a)) then
            deallocate(self%a)
        endif
        if (allocated(self%b)) then
            deallocate(self%b)
        endif
        if (allocated(self%c)) then
            deallocate(self%c)
        endif
        if (allocated(self%bb)) deallocate(self%bb)

        ! Deallocate work
        if (allocated(self%work)) then
            deallocate(self%work)
        endif


        call destroy_plan(self%forward)
        call destroy_plan(self%backward)

    end subroutine finalize_poisson_2d
end module m_poisson