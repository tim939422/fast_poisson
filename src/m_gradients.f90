module m_gradients
    !> module of gradient evaluation of field stored in cell coordinate
    !> like pressure, passive scalar and conformation tensor
    !>
    !> note - resulted gradient is stored at face coordinate (see below)
    !>
    !>        gradpx - xf (i-face) for i: 0 -> nx, j: 1 -> ny, k: 1 -> nz
    !>        gradpy - yf (j-face) for i: 1 -> nx, j: 0 -> ny, k: 1 -> nz
    !>        gradpz - zf (k-face) for i: 1 -> nx, j: 1 -> ny, k: 0 -> nz
    !>
    !>        Example of (dphi/dx)(i) = (phi(i + 1) - phi(i))/dxf(i)
    !>
    !>            phi(i)    phi(i + 1)
    !>        |-----o-----|-----o-----|
    !>              |           |
    !>              |           |
    !>              |<  dxf(i) >|
    !>              
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-30

    use m_kinds, only: rp
    use m_rectilinear, only: t_rectilinear

    implicit none
    
    private

    type, public :: t_gradient
        private
        !> inverse of dx (uniform)
        real(rp) :: invdx
        !> inverse of dy (uniform)
        real(rp) :: invdy
        !> inverse of dyf (2D) or dzf (3D)
        real(rp), allocatable :: invdsf(:)
        !> number of interior cells - x
        integer :: nx
        !> number of interior cells - y
        integer :: ny
        !> number of interior cells - z
        integer :: nz
    contains
        private
        procedure, public :: init => init_gradient
        procedure, public :: gradpx
        procedure, public :: gradpy_2d
        procedure, public :: gradpy_3d
        procedure, public :: gradpz
        procedure, public :: finalize => finalize_gradient
    end type t_gradient

contains
    subroutine init_gradient(self, grid)
        !> constructor for type t_gradient
        !>
        !> note -
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self
        !> rectilinear grid with metric information
        type(t_rectilinear) :: grid

        ! work
        copy_parameter: associate(nx => self%nx, ny => self%ny, nz => self%nz)
            nx = grid%nx; ny = grid%ny; nz = grid%nz
        end associate copy_parameter
        
        allocate_array: associate(nx => self%nx, ny => self%ny, nz => self%nz)
            if (nz > 1) then
                ! 3D
                allocate(self%invdsf(0:nz + 1))
            else
                ! 2D
                allocate(self%invdsf(0:ny + 1))
            end if
        end associate allocate_array
        
        fd_coefficient: associate(nz => self%nz, &
            invdx => self%invdx, invdy => self%invdy, invdsf => self%invdsf)

            invdx = 1.0_rp/grid%dx
            if (nz > 1) then
                invdy = 1.0_rp/grid%dy
            end if
            invdsf(:) = 1.0_rp/grid%dsf(:)
        end associate fd_coefficient

    end subroutine init_gradient


    subroutine gradpx(self, phi, work)
        !> evaluate dp/dx at i-face
        !>
        !> note - for 2D, phi should be of shape (0:nx + 1, 0:ny + 1, 0:1)
        !>        where only phi(0:nx +1, 0:ny + 1, 1) is meaningful
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        ! local
        integer :: i, j, k

        associate(nx => self%nx, ny=>self%ny, nz => self%nz, invdx=>self%invdx)
            do k = 1, nz
                do j = 1, ny
                    do i = 0, nx
                        work(i, j, k) = (phi(i + 1, j, k) - phi(i,  j, k))*invdx
                    end do
                end do
            end do
        end associate

    end subroutine gradpx

    subroutine gradpy_2d(self, phi, work)
        !> evaluate dp/dy at j-face
        !>
        !> note - same as gradpx
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        ! local
        integer :: i, j, k

        ! work
        associate(nx => self%nx, ny=>self%ny, nz => self%nz, invdyf=>self%invdsf)
            do k = 1, nz
                do j = 0, ny
                    do i = 1, nx
                        work(i, j, k) = (phi(i, j + 1, k) - phi(i,  j, k))*invdyf(j)
                    end do
                end do
            end do
        end associate

    end subroutine gradpy_2d   

    subroutine gradpy_3d(self, phi, work)
        !> evaluate dp/y at j-face
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        ! local
        integer :: i, j, k

        ! work
        associate(nx => self%nx, ny=>self%ny, nz => self%nz, invdy=>self%invdy)
            do k = 1, nz
                do j = 0, ny
                    do i = 1, nx
                        work(i, j, k) = (phi(i, j + 1, k) - phi(i,  j, k))*invdy
                    end do
                end do
            end do
        end associate

    end subroutine gradpy_3d

    subroutine gradpz(self, phi, work)
        !> evaluate dp/dz at k-face
        !>
        !> note -
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        ! local
        integer :: i, j, k

        ! work
        associate(nx => self%nx, ny=>self%ny, nz => self%nz, invdzf=>self%invdsf)
            do k = 0, nz
                do j = 1, ny
                    do i = 1, nx
                        work(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k))*invdzf(k)
                    end do
                end do
            end do
        end associate

    end subroutine gradpz

    subroutine finalize_gradient(self)
        !> free resource for type t_gradient
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        class(t_gradient) :: self

        ! work
        if (allocated(self%invdsf)) deallocate(self%invdsf)

    end subroutine finalize_gradient
end module m_gradients