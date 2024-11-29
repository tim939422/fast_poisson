module m_laplacian
    !> module of eigenvalues and difference matrix of Laplacian operator subroutines
    !>
    !> copyright - Yang Group, BUAA
    !>
    !> author - D. Fan, 2024-11-28

    use m_kinds, only: rp
    use m_constants, only: PI
    implicit none

    private

    public :: fft_laplacian, matrix_laplacian

contains

    subroutine fft_laplacian(n, dx, laplacian)
        !> create Laplacian operator with FFT method
        !>
        !> note - only for FD2 on a staggered grid and periodic BC now
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        !> number of cells
        integer, intent(in) :: n
        !> uniform cell width
        real(rp), intent(in) :: dx
        !> eigenvalue divided by squared spacing: lambda/dx^2 
        real(rp), intent(out) :: laplacian(:)

        ! local
        integer :: i
        real(rp) :: invdx2

        do i = 1, n
            laplacian(i) = 2.0_rp*(cos(2.0_rp*PI*real(i - 1, rp)/real(n, rp)) - 1)
        end do

        ! scale by dx^2
        invdx2 = 1.0_rp/dx**2
        laplacian(:) = invdx2*laplacian(:)

    end subroutine fft_laplacian

    subroutine matrix_laplacian(n, dzf, dzc, a, b, c)
        !> create tridiagonal matrix for 1D Laplacian operator on a non-uniform grid
        !>
        !> note - only work for Neumann-Neumann BC now
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-28
        
        ! interface
        !> number of cells
        integer, intent(in) :: n
        !> grid metric at face coordinate
        real(rp), intent(in) :: dzf(0:)
        !> grid metric at cell coordinate
        real(rp), intent(in) :: dzc(0:)
        !> lower diagonal
        real(rp), intent(out) :: a(:)
        !> main diagonal
        real(rp), intent(out) :: b(:)
        !> upper diagonal
        real(rp), intent(out) :: c(:)

        ! local
        integer :: i

        ! work
        do i = 1, n
            a(i) = 1.0_rp/(dzf(i - 1)*dzc(i))
            c(i) = 1.0_rp/(dzf(i)*dzc(i))
            b(i) = -(a(i) + c(i))
        end do

        ! Neumann BC
        b(1) = b(1) + a(1)
        b(n) = b(n) + c(n)

    end subroutine matrix_laplacian

    
end module m_laplacian