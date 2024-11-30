module m_gradients

    use m_kinds, only: rp
    use m_rectilinear, only: t_rectilinear

    implicit none
    
    private

    type, public :: t_gradient
        private
        real(rp) :: invdx
        real(rp) :: invdy
        !> inverse of dyf (2D) or dzf (3D)
        real(rp), allocatable :: invdsf(:)
        integer :: nx
        integer :: ny
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
        ! interface
        class(t_gradient) :: self
        type(t_rectilinear) :: grid

        ! local
        integer  :: nx, ny, nz

        ! work
        nx = grid%nx; ny = grid%ny; nz = grid%nz
        self%nx = nx; self%ny = ny; self%nz = grid%nz

        self%invdx = 1.0_rp/grid%dx

        if (grid%is_3d) then
            self%invdy = 1.0_rp/grid%dy
            allocate(self%invdsf(0:nz + 1))
        else
            allocate(self%invdsf(0:ny + 1))
        end if
        self%invdsf(:) = 1.0_rp/grid%dsf(:)

    end subroutine init_gradient


    subroutine gradpx(self, phi, work)
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
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        integer :: i, j, k

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
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        integer :: i, j, k

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
        ! interface
        class(t_gradient) :: self
        real(rp), intent(in) :: phi(0:, 0:, 0:)
        real(rp), intent(out) :: work(0:, 0:, 0:)

        
        integer :: i, j, k

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
        ! interface
        class(t_gradient) :: self

        if (allocated(self%invdsf)) deallocate(self%invdsf)

    end subroutine finalize_gradient
end module m_gradients