program test_m_metrics
    use m_kinds, only: rp
    use m_metrics, only: staggered_metric
    implicit none
    
    integer :: n
    real(rp), allocatable, dimension(:) :: xf, xc, dxf, dxc
    integer :: iunit

    n = 128
    allocate(xf(0:n + 1), xc(0:n + 1))
    allocate(dxf(0:n + 1), dxc(0:n + 1))

    open(newunit=iunit, file='xf.bin', access='stream', form='unformatted')
    read(iunit) xf
    close(iunit)

    open(newunit=iunit, file='xc.bin', access='stream', form='unformatted')
    read(iunit) xc
    close(iunit)

    call staggered_metric(n, xf, xc, dxf, dxc)

    open(newunit=iunit, file='sol/dxf.bin', access='stream', form='unformatted')
    write(iunit) dxf
    close(iunit)
    open(newunit=iunit, file='sol/dxc.bin', access='stream', form='unformatted')
    write(iunit) dxc
    close(iunit)

    deallocate(xf, xc, dxf, dxc)
end program test_m_metrics