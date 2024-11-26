module m_parameters
    use, intrinsic :: iso_c_binding
    implicit none
    ! case
    integer :: icase

    ! grid size
    integer :: nx, ny

    ! geometry
    real(c_double) :: Lx, Ly


    private :: print_info
contains
    subroutine parse_input(fname)
        character(len=*), intent(in) :: fname

        ! local
        integer :: iunit

        namelist /case/ icase
        namelist /gridsize/ nx, ny
        namelist /geometry/ Lx, Ly

        open(newunit=iunit, file=trim(fname), status="old")
        read(iunit, nml=case)
        read(iunit, nml=gridsize)
        read(iunit, nml=geometry)
        close(iunit)

        call print_info
    end subroutine parse_input

    subroutine print_info()
        use m_case, only: case_string
        write(*, '("Case     : ", A)') trim(case_string(icase))
        write(*, '("(Nx, Ny) : ", 2i5)') nx, ny
        write(*, '("Domain   : ", "x 0 -> ", f15.7, " y 0 ->", f15.7)') Lx, Ly
    end subroutine print_info
end module m_parameters