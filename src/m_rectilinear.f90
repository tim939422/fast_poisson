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

    type, public :: t_rectilinear
        private
        !> grid dimension flag
        logical, public :: is_3d
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
        !> face metric - dy/deta (2D): dyf or dz/dzeta (3D): dzf
        real(rp), allocatable, public :: dsf(:)
        !> cell metric - dy/deta (2D): dyc or dz/dzeta (3D): dzc
        real(rp), allocatable, public :: dsc(:)

    contains
        private
        procedure, public :: init => init_rectilinear
        procedure :: alloc => allocate_rectilinear
        procedure, public :: finalize => finalize_rectilinear
    end type t_rectilinear
contains

    subroutine init_rectilinear(self, dims, domain, beta)
        !> constructor for type t_rectlinear_2d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear) :: self
        !> dimension in cell: dims(1) - nx, dims(2) - ny, dims(3) - nz
        integer, intent(in) :: dims(:)
        !> domain length: domain(1) - Lx, domain(2) - Ly, domain(3) - Lz
        real(rp), intent(in) :: domain(:)
        !> stretching factor in z
        real(rp), intent(in), optional :: beta

        ! local
        integer :: nx, ny, nz
        real(rp) :: Lx, Ly, Lz
        logical :: is_3d

        ! work
        is_3d = size(dims) > 2; self%is_3d = is_3d
        nx = dims(1); self%nx = nx
        ny = dims(2); self%ny = ny
        if (is_3d) then
            nz = dims(3); self%nz = nz
        end if
        call self%alloc

        Lx = domain(1); Ly = domain(2)
        ! x grid and metric
        call staggered_uniform(nx, Lx, self%xf, self%xc)
        self%dx = Lx/real(nx, rp)

        if (is_3d) then
            ! y grid and metric
            call staggered_uniform(ny, Ly, self%yf, self%yc)
            self%dy = Ly/real(ny, rp)

            ! z grid and metric
            Lz = domain(3)
            if (present(beta)) then
                call staggered_twoside_stretched(nz, Lz, beta, self%zf, self%zc)
            else
                call staggered_uniform(nz, Lz, self%zf, self%zc)
            end if
            call staggered_metric(nz, self%zf, self%zc, self%dsf, self%dsc)
        else
            if (present(beta)) then
                call staggered_twoside_stretched(ny, Ly, beta, self%yf, self%yc)
            else
                call staggered_uniform(ny, Ly, self%yf, self%yc)
            end if
            call staggered_metric(ny, self%yf, self%yc, self%dsf, self%dsc)
        end if

    end subroutine init_rectilinear

    subroutine allocate_rectilinear(self)
        !> allocate resource for type t_rectilinear_22d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        class(t_rectilinear) :: self

        associate(nx => self%nx, ny=>self%ny, nz => self%nz, is_3d => self%is_3d)
            allocate(self%xf(0:nx + 1), self%xc(0:nx + 1))
            allocate(self%yf(0:ny + 1), self%yc(0:ny + 1))
            
            if (is_3d) then
                allocate(self%zf(0:nz + 1), self%zc(0:nz + 1))
                allocate(self%dsf(0:nz + 1), self%dsc(0:nz + 1))
            else
                allocate(self%dsf(0:ny + 1), self%dsc(0:ny + 1))
            end if

        end associate

        print *, "t_rectilinear_3d resource allocated"
    end subroutine allocate_rectilinear

    subroutine finalize_rectilinear(self)
        !> free resource for type t_rectilinear_2d
        !>
        !> note - 
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28

        ! interface
        class(t_rectilinear) :: self
        
        if (allocated(self%xf)) deallocate(self%xf)
        if (allocated(self%yf)) deallocate(self%yf)
        if (allocated(self%xc)) deallocate(self%xc)
        if (allocated(self%yc)) deallocate(self%yc)
        if (allocated(self%dsf)) deallocate(self%dsf)
        if (allocated(self%dsc)) deallocate(self%dsc)

        print *, "t_rectilinear_3d resource freed"
        
    end subroutine finalize_rectilinear



end module m_rectilinear