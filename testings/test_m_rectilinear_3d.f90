program test_m_rectilinear_3d
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_rectilinear, only: t_rectilinear
    implicit none
    
    integer, parameter :: nx = 128, ny = 128, nz = 128
    real(rp), parameter :: Lx = 4.0_rp*PI, Ly = 2.0_rp*PI, Lz = 2.0_rp
    real(rp), parameter :: beta = 1.2_rp

    type(t_rectilinear) :: channel_grid
    integer :: iunit
    
    call channel_grid%init([nx, ny, nz], [Lx, Ly, Lz], beta=beta)

    ! Write xf
    open(newunit=iunit, file='sol/xf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%xf
    close(iunit)

    ! Write xc
    open(newunit=iunit, file='sol/xc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%xc
    close(iunit)

    ! Write dx
    open(newunit=iunit, file='sol/dx.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dx
    close(iunit)

    ! Write yf
    open(newunit=iunit, file='sol/yf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%yf
    close(iunit)

    ! Write yc
    open(newunit=iunit, file='sol/yc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%yc
    close(iunit)

    ! Write dy
    open(newunit=iunit, file='sol/dy.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dy
    close(iunit)

    ! Write zf
    open(newunit=iunit, file='sol/zf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%zf
    close(iunit)

    ! Write zc
    open(newunit=iunit, file='sol/zc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%zc
    close(iunit)

    ! Write dzf
    open(newunit=iunit, file='sol/dzf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dsf
    close(iunit)

    ! Write dzc
    open(newunit=iunit, file='sol/dzc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dsc
    close(iunit)    


    call channel_grid%finalize
end program test_m_rectilinear_3d