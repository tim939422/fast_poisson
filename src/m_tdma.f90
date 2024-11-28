module m_tdma
    use m_kinds, only: rp
    use m_constants, only: eps
    implicit none
    private

    public :: tridag
contains
    
    subroutine tridag(a, b, c, r, n)
        real(rp), intent(in), dimension(:) :: a, b, c
        real(rp), intent(inout), dimension(:) :: r
        integer, intent(in) :: n

        ! local
        real(rp) :: gam(n), bet
        integer  :: j

        ! Decomposition and forward substitution
        bet  = b(1)
        r(1) = r(1)/bet
        do j = 2, n
            gam(j) = c(j - 1)/bet
            bet = b(j) - a(j)*gam(j) + eps
            r(j) = (r(j) - a(j)*r(j - 1))/bet
        end do

        ! Backsubstitution
        do j = n - 1, 1, -1
            r(j) = r(j) - gam(j + 1)*r(j + 1)
        end do
    
    end subroutine tridag

end module m_tdma