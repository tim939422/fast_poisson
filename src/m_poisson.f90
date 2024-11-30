module m_poisson
    !> module of FFT acceralated Poisson Solver Class
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28

    use m_kinds, only: rp
    use m_rectilinear, only: t_rectilinear
    use m_laplacian, only: fft_laplacian, matrix_laplacian
    use m_fft, only: DFT, IDFT
    use m_fft, only: create_r2r, execute_fft, destroy_plan
    use, intrinsic :: iso_c_binding, only: c_ptr
    use m_tdma, only: tridag
    use m_io, only: array_write_3d
    

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
        real(rp), allocatable :: work(:, :, :)
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

    type, public :: t_poisson
        private
        !> dimension in x (# of cells)
        integer, public :: nx
        !> dimension in y
        integer, public :: ny
        !> dimension in y
        integer, public :: nz
        !> laplacian operator in x
        real(rp), allocatable :: laplacian_x(:)
        !> laplacian operator in y
        real(rp), allocatable :: laplacian_y(:)
        !> pre-calculate laplacian_x(i) + laplacian_x(j)
        real(rp), allocatable :: laplacian_xy(:, :)
        !> tridiagonal matrix in z (bb is the b + laplacian_x(i) + laplacian_x(j))
        real(rp), allocatable :: a(:), b(:), c(:), bb(:)
        !> work array for FFT
        real(rp), allocatable :: work(:, :, :)
        !> forward FFT plan x -> y
        type(c_ptr) :: forward(2)
        !> backward FFT plan y -> x
        type(c_ptr) :: backward(2)
        !> FFT normalization factor
        real(rp) :: factor
    contains
        private
        procedure, public :: init => init_poisson
        procedure, public :: solve => solve_poisson_3d
        procedure :: alloc => allocate_poisson
        procedure, public :: finalize => finalize_poisson
    end type t_poisson

contains

    subroutine init_poisson_2d(self, grid)
        ! interface
        class(t_poisson_2d) :: self
        type(t_rectilinear) :: grid

        ! local
        integer  :: nx, ny
        real(rp) :: dx

        ! work
        nx = grid%nx; ny = grid%ny
        dx = grid%dx
        self%nx = nx; self%ny = ny
        call self%alloc

        ! setup operator
        call fft_laplacian(nx, dx, self%laplacian_x)
        associate(dyf => grid%dsf, dyc => grid%dsc)
            call matrix_laplacian(ny, dyf, dyc, self%a, self%b, self%c)
        end associate

        ! plan FFT
        self%forward = create_r2r(nx, ny, 1, DFT, 0)
        self%backward = create_r2r(nx, ny, 1, IDFT, 0)
        self%factor = 1.0_rp/real(nx, rp)

    end subroutine init_poisson_2d

    subroutine solve_poisson_2d(self, phi)
        ! interface
        class(t_poisson_2d) :: self
        real(rp), dimension(0:, 0:, 0:) :: phi

        ! local
        integer :: i

        ! work
        associate(nx => self%nx, ny=>self%ny, work=>self%work, &
                  forward=>self%forward, backward=>self%backward, &
                  a=>self%a, b=>self%b, c=>self%c, bb=>self%bb, &
                  laplacian_x=>self%laplacian_x)
            
            ! Copy to work array
            work(:, :, :) = phi(1:nx, 1:ny, 1:1)

            call execute_fft(forward, work)

            do i = 1, nx
                bb(:) = b(:) + laplacian_x(i)
                call tridag(a, bb, c, work(i, :, 1), ny)
            end do

            call execute_fft(backward, work)

            ! Copy back to original array
            phi(1:nx, 1:ny, 1:1) = self%factor*work(:, :, :)

        end associate

    end subroutine solve_poisson_2d

    subroutine allocate_poisson_2d(self)
        ! interface
        class(t_poisson_2d) :: self

        associate(nx => self%nx, ny=>self%ny)
            allocate(self%laplacian_x(nx))
            allocate(self%a(ny), self%b(ny), self%c(ny), self%bb(ny))
            allocate(self%work(nx, ny, 1))
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



    subroutine init_poisson(self, grid)
        ! interface
        class(t_poisson) :: self
        type(t_rectilinear) :: grid

        ! local
        integer  :: nx, ny, nz
        real(rp) :: dx, dy
        integer :: i, j


        ! work
        nx = grid%nx; ny = grid%ny; nz = grid%nz
        dx = grid%dx; dy = grid%dy
        self%nx = nx; self%ny = ny; self%nz = nz
        call self%alloc

        ! setup operator
        call fft_laplacian(nx, dx, self%laplacian_x)
        call fft_laplacian(ny, dy, self%laplacian_y)
        ! pre-calculate
        do j = 1, ny
            do i = 1, nx
                self%laplacian_xy(i, j) = self%laplacian_x(i) + self%laplacian_y(j)
            end do
        end do

        call matrix_laplacian(nz, grid%dsf, grid%dsc, self%a, self%b, self%c)


        ! plan FFT
        self%forward(1) = create_r2r(nx, ny, nz, DFT, 0) ! X
        self%forward(2) = create_r2r(nx, ny, nz, DFT, 1) ! Y
        self%backward(1) = create_r2r(nx, ny, nz, IDFT, 1) ! Y
        self%backward(2) = create_r2r(nx, ny, nz, IDFT, 0) ! X
        self%factor = 1.0_rp/real(nx*ny, rp)

    end subroutine init_poisson

    subroutine solve_poisson_3d(self, phi)
        ! interface
        class(t_poisson) :: self
        real(rp), dimension(0:, 0:, 0:) :: phi

        ! local
        integer :: i, j

        ! work
        associate(nx => self%nx, ny => self%ny, nz => self%nz, &
                  work=>self%work, &
                  forward=>self%forward, backward=>self%backward, &
                  a=>self%a, b=>self%b, c=>self%c, bb=>self%bb, &
                  laplacian_xy=>self%laplacian_xy)
            
            ! Copy to work array
            work(:, :, :) = phi(1:nx, 1:ny, 1:nz)

            ! forward transform X -> Y
            call execute_fft(forward(1), work)
            call execute_fft(forward(2), work)
            
            do j = 1, ny
                do i = 1, nx
                    bb(:) = b(:) + laplacian_xy(i, j)
                    call tridag(a, bb, c, work(i, j, :), nz)
                end do
            end do

            ! backward transform Y -> X
            call execute_fft(backward(1), work)
            call execute_fft(backward(2), work)

            ! Copy back to original array
            phi(1:nx, 1:ny, 1:nz) = self%factor*work(:, :, :)
        end associate

    end subroutine solve_poisson_3d

    subroutine allocate_poisson(self)
        ! interface
        class(t_poisson) :: self

        associate(nx => self%nx, ny => self%ny, nz => self%nz)
            allocate(self%laplacian_x(nx))
            allocate(self%laplacian_y(ny))
            allocate(self%laplacian_xy(nx, ny))
            allocate(self%a(nz), self%b(nz), self%c(nz), self%bb(nz))
            allocate(self%work(nx, ny, nz))
        end associate

        print *, "t_poisson_3d resource allocated"

    end subroutine allocate_poisson

    subroutine finalize_poisson(self)
        ! interface
        class(t_poisson) :: self

        ! Deallocate laplacian_x and laplacian y
        if (allocated(self%laplacian_x)) then
            deallocate(self%laplacian_x)
        endif
        if (allocated(self%laplacian_y)) then
            deallocate(self%laplacian_y)
        endif

        if (allocated(self%laplacian_xy)) deallocate(self%laplacian_xy)

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


        call destroy_plan(self%forward(1)); call destroy_plan(self%forward(2))
        call destroy_plan(self%backward(1)); call destroy_plan(self%backward(2))
    end subroutine finalize_poisson
end module m_poisson