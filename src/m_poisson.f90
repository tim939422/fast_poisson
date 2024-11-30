module m_poisson
    !> module of FFT acceralated Poisson Solver Class
    !>
    !> - see docs for details
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
        procedure, public :: solve_3d => solve_poisson_3d
        procedure, public :: solve_2d => solve_poisson_2d
        procedure :: alloc => allocate_poisson
        procedure, public :: finalize => finalize_poisson
    end type t_poisson

contains


    subroutine init_poisson(self, grid)
        !> constructor for type t_poisson
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_poisson) :: self
        !> rectilinear grid with metric information
        type(t_rectilinear) :: grid

        ! work
        copy_parameter: associate(nx => self%nx, ny => self%ny, nz => self%nz)
            nx = grid%nx; ny = grid%ny; nz = grid%nz
        end associate copy_parameter
        
        call self%alloc

        setup_operator: block

            integer :: i, j

            associate(nx => self%nx, ny => self%ny, nz => self%nz, &
                laplacian_x => self%laplacian_x, laplacian_y => self%laplacian_y, &
                laplacian_xy => self%laplacian_xy)
                
                call fft_laplacian(nx, grid%dx, laplacian_x)
                if (nz > 1) then

                    call fft_laplacian(nx, grid%dy, laplacian_y)
                    call matrix_laplacian(nz, grid%dsf, grid%dsc, self%a, self%b, self%c)

                    ! pre-calculate
                    do j = 1, ny
                        do i = 1, nx
                            laplacian_xy(i, j) = laplacian_x(i) + laplacian_y(j)
                        end do
                    end do
                    call matrix_laplacian(nz, grid%dsf, grid%dsc, self%a, self%b, self%c)
                else
                    call matrix_laplacian(ny, grid%dsf, grid%dsc, self%a, self%b, self%c)
                end if
            
            end associate

        end block setup_operator


        ! plan FFT
        plan_fft: associate(nx => self%nx, ny => self%ny, nz => self%nz)

            self%forward(1) = create_r2r(nx, ny, nz, DFT, 0) ! X
            if (nz > 1) then
                self%forward(2) = create_r2r(nx, ny, nz, DFT, 1) ! Y
                self%backward(1) = create_r2r(nx, ny, nz, IDFT, 1) ! Y
            end if
            self%backward(2) = create_r2r(nx, ny, nz, IDFT, 0) ! X

            if (nz > 1) then
                ! 3D
                self%factor = 1.0_rp/real(nx*ny, rp)
            else
                ! 2D
                self%factor = 1.0_rp/real(nx, rp)
            end if

        end associate plan_fft

    end subroutine init_poisson

    subroutine solve_poisson_3d(self, phi)
        !> solver 3D Poisson equation on a rectilinear grid with possible
        !> stretched grid in z
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_poisson) :: self
        ! source term S in nabla^2 phi = S
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
            
            ! solve linear system
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

    subroutine solve_poisson_2d(self, phi)
        !> solver 2D Poisson equation on a rectilinear grid with possible
        !> stretched grid in y
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30

        ! interface
        class(t_poisson) :: self
        real(rp), dimension(0:, 0:, 0:) :: phi

        ! local
        integer :: i, j

        ! work
        associate(nx => self%nx, ny => self%ny, &
                  work=>self%work, &
                  forward=>self%forward, backward=>self%backward, &
                  a=>self%a, b=>self%b, c=>self%c, bb=>self%bb, &
                  laplacian_x=>self%laplacian_x)
            
            ! Copy to work array
            work(:, :, :) = phi(1:nx, 1:ny, 1:1)

            ! forward transform X
            call execute_fft(forward(1), work)
            
            ! solve linear system
            do i = 1, nx
                bb(:) = b(:) + laplacian_x(i)
                call tridag(a, bb, c, work(i, :, 1), ny)
            end do

            ! backward transform X
            call execute_fft(backward(2), work)

            ! Copy back to original array
            phi(1:nx, 1:ny, 1:1) = self%factor*work(:, :, :)
            
        end associate

    end subroutine solve_poisson_2d

    subroutine allocate_poisson(self)
        ! interface
        class(t_poisson) :: self

        associate(nx => self%nx, ny => self%ny, nz => self%nz)

            allocate(self%laplacian_x(nx))
            allocate(self%work(nx, ny, nz))

            if (nz > 1) then
                allocate(self%laplacian_y(ny))
                allocate(self%laplacian_xy(nx, ny))
                allocate(self%a(nz), self%b(nz), self%c(nz), self%bb(nz))
            else
                allocate(self%a(ny), self%b(ny), self%c(ny), self%bb(ny))
            end if

        end associate

        print *, "t_poisson_3d resource allocated"

    end subroutine allocate_poisson

    subroutine finalize_poisson(self)
        ! interface
        class(t_poisson) :: self

        if (allocated(self%laplacian_x)) deallocate(self%laplacian_x)
        if (allocated(self%laplacian_y)) deallocate(self%laplacian_y)
        if (allocated(self%laplacian_xy)) deallocate(self%laplacian_xy)

        if (allocated(self%a)) deallocate(self%a)
        if (allocated(self%b)) deallocate(self%b)
        if (allocated(self%c)) deallocate(self%c)
        if (allocated(self%bb)) deallocate(self%bb)

        if (allocated(self%work)) then
            deallocate(self%work)
        endif

        associate(nz => self%nz)

            call destroy_plan(self%forward(1))
            if (nz > 1) then
                call destroy_plan(self%forward(2))
                call destroy_plan(self%backward(1))
            end if
            call destroy_plan(self%backward(2))

        end associate

    end subroutine finalize_poisson

end module m_poisson