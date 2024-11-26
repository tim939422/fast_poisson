program main
    use m_parameters
    use m_variables
    use m_grid
    use m_case
    use m_boundary
    use m_poisson_solver
    use m_RHS
    use m_utils, only: compare_pressure_gradient
    implicit none
    integer :: i, j

    call get_command_argument(1, fname)
    if (len_trim(fname) == 0) then
        print *, "Error: No input file specified. Usage: PoissonSolver <input>"
        stop
    end if

    call parse_input(fname)
    call init_variables(nx, ny)
    call create_regular_grid(nx, ny, Lx, Ly, dx, dy, xc, yc, xn, yn)

    ! initialize case
    call init_case(icase, nx, ny, p)

    ! initialize boundary condition
    call init_boundary(icase, BC)

    ! initialize Poisson equation
    call init_poisson(nx, ny, dx, dy, BC)

    ! fill RHS of pressure
    call fill_pressure_poisson_RHS_a(icase, nx, ny, xc, yc, p)

    call solve_poisson(nx, ny, p)

    ! set pressure boundary
    call set_pressure_boundary(nx, ny, BC, p)

    ! compute pressure gradient
    ! dp/dx
    do j = 0, ny + 1
        do i = 0, nx
            dpdx(i, j) = (p(i + 1, j) - p(i, j))/dx
        end do
    end do

    ! dp/dy
    do j = 0, ny
        do i = 0, nx + 1
            dpdy(i, j) = (p(i, j + 1) - p(i, j))/dy
        end do
    end do
    
    call compare_pressure_gradient(icase, nx, ny, xc, yc, xn, yn, p, dpdx, dpdy)

    call finalize_poisson
    call finalize_variables
end program main