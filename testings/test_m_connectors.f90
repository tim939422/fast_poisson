program test_m_connectors
    use m_kinds, only: rp
    use m_constants, only: PI
    use m_connectors, only: staggered_twoside_stretched, staggered_uniform
    implicit none
    
    integer :: n
    real(rp) :: L, beta
    real(rp), allocatable, dimension(:) :: xf, xc
    integer :: iunit
    
    ! test uniform grid
    n = 128
    L = 4.0_rp*PI
    allocate(xf(0:n + 1), xc(0:n + 1))
    call staggered_uniform(n, L, xf, xc)
    open(newunit=iunit, file='sol/staggered_uniform_grid_xf.bin', access='stream', form='unformatted')
    write(iunit) xf
    close(iunit)
    open(newunit=iunit, file='sol/staggered_uniform_grid_xc.bin', access='stream', form='unformatted')
    write(iunit) xc
    close(iunit)
    deallocate(xf, xc)

    ! test twoside stretched grid
    n = 128
    L = 2.0_rp
    beta = 1.2_rp
    allocate(xf(0:n + 1), xc(0:n + 1))
    call staggered_twoside_stretched(n, L, beta, xf, xc)
    open(newunit=iunit, file='sol/staggered_twoside_stretched_grid_xf.bin', access='stream', form='unformatted')
    write(iunit) xf
    close(iunit)
    open(newunit=iunit, file='sol/staggered_twoside_stretched_grid_xc.bin', access='stream', form='unformatted')
    write(iunit) xc
    close(iunit)
    deallocate(xf, xc)

end program test_m_connectors