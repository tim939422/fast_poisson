module m_tdma
    use m_kinds, only: rp
    use m_constants, only: eps
    implicit none
    private

    public :: tridag
contains
    
    subroutine tridag(a, b, c, r, n)
        !> TDMA algorithm of a tridiagonal system (adopted from numerical recipes
        !> pp. 43)
        !>
        !> note - we solve the following system
        !>
        !> [ b1    c1     0     ...       0     ] [ u1   ]   = [ r1   ]
        !> [ a2    b2    c2     ...       0     ] [ u2   ]     [ r2   ]
        !> [  0    a3    b3    c3     ...       ] [ u3   ]     [ r3   ]
        !> [ ...   ...   ...   ...     ...      ] [ ...  ]     [ ...  ]
        !> [  0     0   aN-1  bN-1   cN-1       ] [ uN-1 ]     [ rN-1 ]
        !> [  0     0     0    aN     bN        ] [ uN   ]     [ rN   ]
        !>
        !> singularity treatment - use regularization of machine zero, which
        !>                         is useful in Neumann-Neumann BC. The
        !>                         algorithm can solve for u up to a constant
        !>                         vector. This will suffice if we only care
        !>                         about the derivatives of u
        !>
        !>
        !> copyright - Yang Group, BUAA
        !>
        !> author - D. Fan, 2024-11-30
        
        ! interface
        real(rp), intent(in), dimension(:) :: a, b, c
        real(rp), intent(inout), dimension(:) :: r
        integer, intent(in) :: n

        ! local
        real(rp) :: gam(n), bet
        integer  :: j

        ! work
        ! decomposition and forward substitution
        bet  = b(1) + eps
        r(1) = r(1)/bet
        do j = 2, n
            gam(j) = c(j - 1)/bet
            bet = b(j) - a(j)*gam(j) + eps
            r(j) = (r(j) - a(j)*r(j - 1))/bet
        end do

        ! back substitution
        do j = n - 1, 1, -1
            r(j) = r(j) - gam(j + 1)*r(j + 1)
        end do
    
    end subroutine tridag
    
end module m_tdma
