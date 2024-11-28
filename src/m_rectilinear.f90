module m_rectilinear
    !> module of rectilinear grid type 2D/3D
    !>
    !> - note only allow non-uniform in the last axis
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28
    
    use m_kinds, only: rp
    use m_connectors, only: staggered_twoside_stretched, staggered_uniform
    use m_metrics, only: staggered_metric

    implicit none

    private

    type, public :: t_rectilinear_2d
        private
        !> dimension in x (# of cells)
        integer, public :: nx
        !> dimension in y
        integer, public :: ny
        !> face coordinate in x
        real(rp), allocatable, public :: xf(:)
        !> face coordinate in y
        real(rp), allocatable, public :: yf(:)
        !> cell coordinate in x
        real(rp), allocatable, public :: xc(:)
        !> cell coordinate in y
        real(rp), allocatable, public :: yc(:)
        !> metric dx/dxi (uniform)
        real(rp), public :: dx
        !> metric dy/deta (face)
        real(rp), allocatable, public :: dyf(:)
        !> metric dy/deta (cell)
        real(rp), allocatable, public :: dyc(:)
    contains
        private
        procedure, public :: init => init_rectilinear_2d
        procedure :: alloc => allocate_rectilinear_2d
        procedure, public :: finalize => finalize_rectilinear_2d
    end type t_rectilinear_2d

    type, public :: t_rectilinear_3d
        private
        !> dimension in x (# of cells)
        integer, public :: nx
        !> dimension in y
        integer, public :: ny
        !> dimension in z
        integer, public :: nz
        !> face coordinate in x
        real(rp), allocatable, public :: xf(:)
        !> face coordinate in y
        real(rp), allocatable, public :: yf(:)
        !> face coordinate in z
        real(rp), allocatable, public :: zf(:)
        !> cell coordinate in x
        real(rp), allocatable, public :: xc(:)
        !> cell coordinate in y
        real(rp), allocatable, public :: yc(:)
        !> cell coordinate in z
        real(rp), allocatable, public :: zc(:)

        !> metric dx/dxi (uniform)
        real(rp), public :: dx
        !> metric dy/deta (uniform)
        real(rp), public :: dy
        !> metric dz/dzeta (face)
        real(rp), allocatable, public :: dzf(:)
        !> metric dz/dzeta (cell)
        real(rp), allocatable, public :: dzc(:)
    contains
        private
        procedure, public :: init => init_rectilinear_3d
        procedure :: alloc => allocate_rectilinear_3d
        procedure, public :: finalize => finalize_rectilinear_3d
    end type t_rectilinear_3d
contains
    
    subroutine init_rectilinear_2d(self, nx, ny, Lx, Ly, beta)
        !> constructor for type t_rectlinear_2d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear_2d) :: self
        !> dimension in x
        integer, intent(in) :: nx
        !> dimension in y
        integer, intent(in) :: ny
        !> domain length in x
        real(rp), intent(in) :: Lx
        !> domain length in y
        real(rp), intent(in) :: Ly
        !> stretching factor in y
        real(rp), intent(in), optional :: beta

        ! local

        ! work
        self%nx = nx; self%ny = ny
        call self%alloc

        ! x grid and metric
        call staggered_uniform(nx, Lx, self%xf, self%xc)
        self%dx = Lx/real(nx, rp)

        ! y grid and metric
        if (present(beta)) then
            call staggered_twoside_stretched(ny, Ly, beta, self%yf, self%yc)
        else
            call staggered_uniform(ny, Ly, self%yf, self%yc)
        end if
        call staggered_metric(ny, self%yf, self%yc, self%dyf, self%dyc)
        
    end subroutine init_rectilinear_2d

    subroutine allocate_rectilinear_2d(self)
        !> allocate resource for type t_rectilinear_22d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear_2d) :: self

        associate(nx => self%nx, ny=>self%ny)
            allocate(self%xf(0:nx + 1), self%xc(0:nx + 1))
            allocate(self%yf(0:ny + 1), self%yc(0:ny + 1), self%dyf(0:ny + 1), self%dyc(0:ny + 1))
        end associate

        print *, "t_rectilinear_2d resource allocated"
    end subroutine allocate_rectilinear_2d

    subroutine finalize_rectilinear_2d(self)
        !> free resource for type t_rectilinear_22d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        class(t_rectilinear_2d) :: self
        
        if (allocated(self%xf)) deallocate(self%xf)
        if (allocated(self%yf)) deallocate(self%yf)
        if (allocated(self%xc)) deallocate(self%xc)
        if (allocated(self%yc)) deallocate(self%yc)
        if (allocated(self%dyf)) deallocate(self%dyf)
        if (allocated(self%dyc)) deallocate(self%dyc)

        print *, "t_rectilinear_2d resource freed"
        
    end subroutine finalize_rectilinear_2d


    subroutine init_rectilinear_3d(self, nx, ny, nz, Lx, Ly, Lz, beta)
        !> constructor for type t_rectlinear_2d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear_3d) :: self
        !> dimension in x
        integer, intent(in) :: nx
        !> dimension in y
        integer, intent(in) :: ny
        !> dimension in z
        integer, intent(in) :: nz
        !> domain length in x
        real(rp), intent(in) :: Lx
        !> domain length in y
        real(rp), intent(in) :: Ly
        !> domain length in z
        real(rp), intent(in) :: Lz
        !> stretching factor in z
        real(rp), intent(in), optional :: beta

        ! local

        ! work
        self%nx = nx; self%ny = ny; self%nz = nz
        call self%alloc

        ! x grid and metric
        call staggered_uniform(nx, Lx, self%xf, self%xc)
        self%dx = Lx/real(nx, rp)

        ! y grid and metric
        call staggered_uniform(ny, Ly, self%yf, self%yc)
        self%dy = Ly/real(ny, rp)

        ! y grid and metric
        if (present(beta)) then
            call staggered_twoside_stretched(nz, Lz, beta, self%zf, self%zc)
        else
            call staggered_uniform(nz, Lz, self%zf, self%zc)
        end if
        call staggered_metric(nz, self%zf, self%zc, self%dzf, self%dzc)
        
    end subroutine init_rectilinear_3d

    subroutine allocate_rectilinear_3d(self)
        !> allocate resource for type t_rectilinear_22d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear_3d) :: self

        associate(nx => self%nx, ny=>self%ny, nz => self%nz)
            allocate(self%xf(0:nx + 1), self%xc(0:nx + 1))
            allocate(self%yf(0:ny + 1), self%yc(0:ny + 1))
            allocate(self%zf(0:nz + 1), self%zc(0:nz + 1), self%dzf(0:nz + 1), self%dzc(0:nz + 1))
        end associate

        print *, "t_rectilinear_3d resource allocated"
    end subroutine allocate_rectilinear_3d

    subroutine finalize_rectilinear_3d(self)
        !> free resource for type t_rectilinear_22d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        class(t_rectilinear_3d) :: self
        
        if (allocated(self%xf)) deallocate(self%xf)
        if (allocated(self%yf)) deallocate(self%yf)
        if (allocated(self%xc)) deallocate(self%xc)
        if (allocated(self%yc)) deallocate(self%yc)
        if (allocated(self%dzf)) deallocate(self%dzf)
        if (allocated(self%dzc)) deallocate(self%dzc)

        print *, "t_rectilinear_3d resource freed"
        
    end subroutine finalize_rectilinear_3d



end module m_rectilinear