program test_m_rectilinear_2d
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_rectilinear, only: t_rectilinear
    implicit none
    
    integer, parameter :: nx = 128, ny = 128
    real(rp), parameter :: Lx = 4.0_rp*PI, Ly = 2.0_rp
    real(rp), parameter :: beta = 1.2_rp

    type(t_rectilinear) :: channel_grid
    integer :: iunit
    
    call channel_grid%init([nx, ny], [Lx, Ly], beta=beta)

    ! Write xf
    open(newunit=iunit, file='sol/xf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%xf
    close(iunit)

    ! Write yf
    open(newunit=iunit, file='sol/yf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%yf
    close(iunit)

    ! Write xc
    open(newunit=iunit, file='sol/xc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%xc
    close(iunit)

    ! Write yc
    open(newunit=iunit, file='sol/yc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%yc
    close(iunit)

    ! Write dx
    open(newunit=iunit, file='sol/dx.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dx
    close(iunit)

    ! Write dyf
    open(newunit=iunit, file='sol/dyf.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dsf
    close(iunit)

    ! Write dyc
    open(newunit=iunit, file='sol/dyc.bin', access='stream', form='unformatted')
    write(iunit) channel_grid%dsc
    close(iunit)


    call channel_grid%finalize
end program test_m_rectilinear_2d