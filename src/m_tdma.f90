module m_tdma
    use m_numerics
    implicit none
contains
    !> Solve a tridiagonal linear system with TDMA algorithm
    !>
    !> @param[in]    n   dimension of matrix (n x n)
    !> @param[in]    a   lower diagonal
    !> @param[in]    b   main diagonal
    !> @param[in]    c   upper diagonal
    !> @param[inout] r   rhs (in), solution (out)
    !>
    !> This is a modified version of tridag in Numeical Recipes in Fortran 77
    !> (see p. 43). A small regularization eps is added to deal with singular
    !> A in the Neumann-Neumann boundary condition
    !>
    !>                A                      *   u       =    r
    !> [ b1    c1     0     ...       0     ] [ u1   ]     [ r1   ]
    !> [ a2    b2    c2     ...       0     ] [ u2   ]     [ r2   ]
    !> [  0    a3    b3    c3     ...       ] [ u3   ]   = [ r3   ]
    !> [ ...   ...   ...   ...     ...      ] [ ...  ]     [ ...  ]
    !> [  0     0   aN-1  bN-1   cN-1       ] [ uN-1 ]     [ rN-1 ]
    !> [  0     0     0    aN     bN        ] [ uN   ]     [ rN   ]
    subroutine tridag(n, a, b, c, r)
        integer, intent(in) :: n
        real(rp), intent(in), dimension(:) :: a, b, c
        real(rp), intent(inout), dimension(:) :: r

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