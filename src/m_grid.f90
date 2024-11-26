module m_grid
    use m_numerics
    implicit none
contains

    !> Generate a two-sided stretched grid using a hyperbolic tangent distribution.
    !>
    !> @param[in] n     number of cells (n + 1 faces)
    !> @param[in] beta  stretching factor controlling the clustering of points
    !> @param[in] L     domain length
    !> @param[out] x    face coordinate 0:n
    subroutine twoside_stretch_grid(n, beta, L, x)
        integer, intent(in) :: n
        real(rp), intent(in) :: beta, L
        real(rp), intent(out) :: x(0:)

        integer :: i

        do i = 0, n
            x(i) = half*L*(one - tanh(beta*(one - two*real(i, rp)/real(N, rp)))/tanh(beta))
        end do

    end subroutine twoside_stretch_grid

    
    !> Evaluate grid metric for 1D non-uniform grid with 2nd order central difference
    !>
    !> @param[in] n     number of interior cells (n + 1 faces)
    !> @param[in] xf    face coordinate 0:n+1
    !> @param[in] xc    cell coorindate 0:n+1
    !> @param[out] dxf  metric at face
    !> @param[out] dxc  metric at cell
    subroutine nonuniform_metric(n, xf, xc, dxf, dxc)
        integer, intent(in) :: n
        real(rp), intent(in), dimension(0:) :: xf, xc
        real(rp), intent(out), dimension(0:) :: dxf, dxc

        integer :: i

        ! \(x_{\xi}\) at cell center
        do i = 1, N + 1
            dxc(i) = xf(i) - xf(i - 1)
        end do
        dxc(0) = dxc(1)

        ! \(x_{\xi}\) at face center
        do i = 0, N
            dxf(i) = xc(i + 1) - xc(i)
        end do
        dxf(N + 1) = dxf(N)

    end subroutine nonuniform_metric
end module m_grid