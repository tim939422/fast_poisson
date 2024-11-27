module m_laplacian
    use m_numerics
    implicit none
    private
    public :: laplacian_1d_tridiagonal, laplacian_dft
contains
    subroutine laplacian_1d_tridiagonal(n, dxf, dxc, a, b, c)
        integer, intent(in) :: n
        real(rp), intent(in), dimension(0:) :: dxf, dxc
        real(rp), intent(out), dimension(:) :: a, b, c

        integer :: i

        do i = 1, n
            a(i) = 1/(dxf(i - 1)*dxc(i))
            c(i) = 1/(dxf(i)*dxc(i))
            b(i) = -(a(i) + c(i))
        end do

        call tridiagonal_bc(n, a, b, c)
    end subroutine laplacian_1d_tridiagonal

    subroutine laplacian_dft(n, dx, laplacian)
        integer, intent(in)   :: n
        real(rp), intent(in)  :: dx
        real(rp), intent(out), dimension(:) :: laplacian

        integer :: i
        real(rp) :: invdx2
        ! eigenvalues
        ! P-P
        do i = 1, n
            laplacian(i) = two*(cos(two*PI*real(i - 1, rp)/real(n, rp)) - one)
        end do

        invdx2 = one/dx**2
        do i = 1, n
            laplacian(i) = laplacian(i)*invdx2
        end do

    end subroutine laplacian_dft



    ! Private subroutine
    subroutine tridiagonal_bc(n, a, b, c)
        integer, intent(in) :: n
        real(rp), intent(inout), dimension(:) :: a, b, c

        ! homogeneous Neumann BC
        b(1) = b(1) + a(1)
        b(n) = b(n) + c(n)

    end subroutine tridiagonal_bc
end module m_laplacian