module m_gradients

    use m_kinds, only: rp
    use m_rectilinear, only: t_rectilinear_2d

    implicit none
    
    private

    type, public :: t_gradient_2d
        private
        real(rp) :: invdx
        real(rp), allocatable :: invdyf(:)
        integer :: nx
        integer :: ny
    contains
        private
        procedure, public :: init => init_gradient_2d
        procedure, public :: dpdx => gradpx_2d
        procedure, public :: dpdy => gradpy_2d
        procedure, public :: finalize => finalize_gradient_2d
    end type t_gradient_2d

contains

    subroutine init_gradient_2d(self, grid)
        ! interface
        class(t_gradient_2d) :: self
        type(t_rectilinear_2d) :: grid

        ! local
        integer  :: nx, ny

        ! work
        nx = grid%nx; ny = grid%ny
        self%nx = nx; self%ny = ny
        self%invdx = 1.0_rp/grid%dx
        allocate(self%invdyf(0:ny + 1))
        self%invdyf(:) = 1.0_rp/grid%dyf(:)

    end subroutine init_gradient_2d

    subroutine gradpx_2d(self, phi, work)
        ! interface
        class(t_gradient_2d) :: self
        real(rp), intent(in) :: phi(0:, 0:)
        real(rp), intent(out) :: work(0:, 0:)

        ! local
        integer :: i, j

        associate(nx => self%nx, ny=>self%ny, invdx=>self%invdx)
            do j = 1, ny
                do i = 0, nx
                    work(i, j) = (phi(i + 1, j) - phi(i,  j))*invdx
                end do
            end do
        end associate
        
    end subroutine gradpx_2d

    subroutine gradpy_2d(self, phi, work)
        ! interface
        class(t_gradient_2d) :: self
        real(rp), intent(in) :: phi(0:, 0:)
        real(rp), intent(out) :: work(0:, 0:)
        
        ! local
        integer :: i, j

        associate(nx => self%nx, ny=>self%ny, invdyf=>self%invdyf)
            do j = 0, ny
                do i = 1, nx
                    work(i, j) = (phi(i, j + 1) - phi(i,  j))*invdyf(j)
                end do
            end do
        end associate

    end subroutine gradpy_2d

    subroutine finalize_gradient_2d(self)
        ! interface
        class(t_gradient_2d) :: self

        if (allocated(self%invdyf)) deallocate(self%invdyf)

    end subroutine finalize_gradient_2d
    
end module m_gradients