program main
    use m_constants
    use m_poisson_solver
    implicit none
    real(C_DOUBLE), parameter :: L = 1.0_C_DOUBLE

    

    real(C_DOUBLE), dimension(:), allocatable :: xc, xe
    real(C_DOUBLE), dimension(:), allocatable :: p, dpdx, dpdx_a
    real(C_DOUBLE) :: dx, invdx
    integer :: i, iunit

    ! configuration
    integer            :: N
    character(len=2)   :: BC
    character(len=100) :: fname
    namelist /config/ N, BC

    call get_command_argument(1, fname)
    if (len_trim(fname) == 0) then
        print *, "Error: No input file specified. Usage: ./poisson_1d <input>"
        stop
    end if

    ! read input
    ! default value
    N = 256
    BC = 'PP'
    open(newunit=iunit, file=trim(fname), status="old")
    read(iunit, nml=config)
    close(iunit)

    
    ! Allocate memeory for Poisson
    allocate(xc(0:N + 1), xe(0:N + 1))
    allocate(p(0:N + 1), dpdx(0:N + 1), dpdx_a(0:N + 1))

    ! Geometry
    dx = L/real(N, C_DOUBLE)
    invdx = 1.0_C_DOUBLE/dx
    do i = 1, N
        xc(i) = real(i - 1, C_DOUBLE)*dx + 0.5_C_DOUBLE*dx
    end do

    call init_poisson(N, dx, BC)

    ! set RHS
    select case (BC)
    case ('PP')
        do i = 1, N
            p(i) = -(two*PI)**2*sin(two*PI*xc(i))
        end do
    case ('NN')
        do i = 1, N
            p(i) = -PI**2*cos(PI*xc(i))
        end do
    case ('ND')
        do i = 1, N
            p(i) = -(three*PI/two)**2*cos(three*PI/two*xc(i))
        end do
    case default
        print *, "Error: Invalid BC: ", BC
    end select

    call solve_poisson(N, p)

    ! Set BC
    select case (BC)
    case ('PP')
        p(0)     = p(N)
        p(N + 1) = p(1)
    case ('NN')
        p(0)     = p(1)
        p(N + 1) = p(N)
    case ('ND')
        p(0)     = p(1)
        p(N + 1) = -p(N)
    case default
        print *, "Error: Invalid BC: ", BC
    end select 

    ! dp/dx by FD2
    do i = 0, N
        dpdx(i) = (p(i + 1) - p(i))*invdx
    end do

    ! edge coordinate and analytical solution
    do i = 0, N
        xe(i)     = real(i, C_DOUBLE)*dx
    end do

    ! Analytical solution
    select case (BC)
    case ('PP')
        do i = 0, N
            dpdx_a(i) = two*PI*cos(two*PI*xe(i))
        end do
    case ('NN')
        do i = 0, N
            dpdx_a(i) = -PI*sin(PI*xe(i))
        end do
    case ('ND')
        do i = 0, N
            dpdx_a(i) = -three*PI/two*sin(three*PI/two*xe(i))
        end do
    case default
        print *, "Error: Invalid BC: ", BC
    end select

    ! evaluate error
    print *, norm2(dpdx(0:N) - dpdx_a(0:N))/norm2(dpdx_a(0:N))

    ! Clean up
    deallocate(xc, xe)
    deallocate(p, dpdx, dpdx_a)
    call finalize_poisson()

end program main