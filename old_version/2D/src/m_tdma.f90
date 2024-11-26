module m_tdma
    use, intrinsic :: iso_c_binding
    use m_constants, only: eps, one
    implicit none
    ! Declare public procedures
    public :: TDMA, TDMA_p
contains
    subroutine TDMA(n, a, b, c, p)
        use, intrinsic :: iso_c_binding
        use m_constants, only: eps, one
        implicit none
        integer, intent(in) :: n
        real(c_double), intent(in), dimension(:) :: a, b, c 
        real(c_double), intent(inout), dimension(:) :: p

        ! local
        real(c_double), dimension(n) :: d
        real(c_double) :: z
        integer :: l

        ! Gauss elimination
        z    = one/b(1)
        d(1) = c(1)*z
        p(1) = p(1)*z
        do l = 2, n
            z = one/(b(l) - a(l)*d(l - 1) + eps)
            d(l) = c(l)*z
            p(l) = (p(l) - a(l)*p(l - 1))*z
        end do

        ! Back substitution
        do l = n - 1, 1, -1
            p(l) = p(l) - d(l)*p(l + 1)
        end do

    end subroutine TDMA

    subroutine TDMA_p(n, a, b, c, p)
        use, intrinsic :: iso_c_binding
        use m_constants, only: eps, zero, one
        implicit none
        integer, intent(in) :: n
        real(c_double), intent(in), dimension(:) :: a, b, c 
        real(c_double), intent(inout), dimension(:) :: p
    
        ! local
        real(c_double), dimension(n) :: p1, p2
    
        ! Solve for interior points
        p1(1:n-1) = p(1:n-1)
        call TDMA(n - 1, a, b, c, p1)

        ! Solve the homogeneous system for correction
        p2(:) = zero
        p2(1) = -a(1)
        p2(n-1) = -c(n-1)
        call TDMA(n - 1, a, b, c, p2)

        ! Update
        p(n) = (p(n) - c(n)*p1(1) - a(n)*p1(n-1))/(b(n) + c(n)*p2(1) + a(n)*p2(n - 1) + eps)
        p(1:n-1) = p1(1:n-1) + p2(1:n-1)*p(n)

    end subroutine TDMA_p
end module m_tdma