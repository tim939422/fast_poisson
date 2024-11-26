module m_laplacian
    use m_numerics
    implicit none
    private
    public :: laplacian_1d_tridiagonal
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



    ! Private subroutine
    subroutine tridiagonal_bc(n, a, b, c)
        integer, intent(in) :: n
        real(rp), intent(inout), dimension(:) :: a, b, c

        ! homogeneous Neumann BC
        b(1) = b(1) + a(1)
        b(n) = b(n) + c(n)

    end subroutine tridiagonal_bc
end module m_laplacian