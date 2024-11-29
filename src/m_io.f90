module m_io

    use m_kinds, only: rp

    implicit none
    private

    public :: array_write_1d, array_write_2d, array_write_3d
contains
    subroutine array_write_1d(fname, data)
        character(len=*), intent(in) :: fname
        real(rp), dimension(:) :: data

        integer :: iunit
        open(newunit=iunit, file=trim(fname), status='replace', access='stream', form='unformatted')
        write(iunit) data
        close(iunit)
    end subroutine array_write_1d

    subroutine array_write_2d(fname, data)
        character(len=*), intent(in) :: fname
        real(rp), dimension(:, :) :: data

        integer :: iunit
        open(newunit=iunit, file=trim(fname), status='replace', access='stream', form='unformatted')
        write(iunit) data
        close(iunit)
    end subroutine array_write_2d

    subroutine array_write_3d(fname, data)
        character(len=*), intent(in) :: fname
        real(rp), dimension(:, :, :) :: data

        integer :: iunit
        open(newunit=iunit, file=trim(fname), status='replace', access='stream', form='unformatted')
        write(iunit) data
        close(iunit)
    end subroutine array_write_3d
end module m_io