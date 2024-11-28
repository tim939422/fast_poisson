program test_m_laplacian
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_laplacian, only: fft_laplacian, matrix_laplacian
    implicit none
    
    integer, parameter :: nx = 128, ny = 128
    real(rp), parameter :: Lx = 4.0_rp*PI, Ly = 2.0_rp
    real(rp), dimension(nx) :: laplacian
    real(rp), dimension(ny) :: a, b, c
    real(rp), dimension(0:ny + 1) :: dyf, dyc
    real(rp) :: dx
    integer :: iunit

    dx = Lx/real(nx, rp)
    call fft_laplacian(nx, dx, laplacian)
    open(newunit=iunit, file='sol/laplacian.bin', access='stream', form='unformatted')
    write(iunit) laplacian
    close(iunit)


    open(newunit=iunit, file='dyf.bin', access='stream', form='unformatted')
    read(iunit) dyf
    close(iunit)

    open(newunit=iunit, file='dyc.bin', access='stream', form='unformatted')
    read(iunit) dyc
    close(iunit)

    call matrix_laplacian(ny, dyf, dyc, a, b, c)

        ! Write a
    open(newunit=iunit, file='sol/a.bin', access='stream', form='unformatted')
    write(iunit) a
    close(iunit)

    ! Write b
    open(newunit=iunit, file='sol/b.bin', access='stream', form='unformatted')
    write(iunit) b
    close(iunit)

    ! Write c
    open(newunit=iunit, file='sol/c.bin', access='stream', form='unformatted')
    write(iunit) c
    close(iunit)

    
end program test_m_laplacian