module m_boundary
    use, intrinsic :: iso_c_binding
    use m_case, only: TAYLOR_GREEN, STEADY_CHANNEL, CAVITY, DEVELOPING_CHANNEL
    implicit none
    integer, parameter :: DIRICHLET = 0, NEUMANN = 1, PERIODIC = 2
    character(len=*), private, dimension(0:2), parameter :: pbc_string = ["DIRICHLET", &
                                                                          "NEUMANN  ", &
                                                                          "PERIODIC "]

contains
    ! setup boundary condition for a case
    subroutine init_boundary(icase, BC)
        integer, intent(in)  :: icase
        ! BC(1, 1) : left
        ! BC(2, 1) : right
        ! BC(1, 2) : bottom
        ! BC(2, 2) : top
        integer, intent(out) :: BC(:, :)
        integer :: i, j
        select case (icase)
        case (TAYLOR_GREEN)
            BC = reshape(source = [PERIODIC, PERIODIC, PERIODIC, PERIODIC], &
                         shape  = [2, 2])
        case (STEADY_CHANNEL)
            BC = reshape(source = [PERIODIC, PERIODIC, NEUMANN, NEUMANN], &
                         shape  = [2, 2])
        case (CAVITY)
            BC = reshape(source = [NEUMANN, NEUMANN, NEUMANN, NEUMANN], &
                         shape  = [2, 2])
        case (DEVELOPING_CHANNEL)
            BC = reshape(source = [NEUMANN, DIRICHLET, NEUMANN, NEUMANN], &
                         shape  = [2, 2])
        end select

        write(*, '(A)') "Boundary condition for left, right, bottom and top:"
        do j = 1, 2
            do i = 1, 2
                write(*, '(A)') trim(pbc_string(BC(i, j)))
            end do
        end do

    end subroutine init_boundary
    
    subroutine set_pressure_boundary(nx, ny, BC, p)
        integer, intent(in) :: nx, ny
        ! BC(1, 1) : left
        ! BC(2, 1) : right
        ! BC(1, 2) : bottom
        ! BC(2, 2) : top
        integer, intent(in), dimension(:, :) :: BC
        real(c_double), intent(inout) :: p(0:, 0:)
        
        ! left and right
        select case(BC(1, 1))
        case (PERIODIC)
            p(0, 1:ny) = p(nx, 1:ny)
        case (DIRICHLET)
            p(0, 1:ny) = -p(1, 1:ny)
        case (NEUMANN)
            p(0, 1:ny) = p(1, 1:ny)
        end select

        select case(BC(2, 1))
        case (PERIODIC)
            p(nx + 1, 1:ny) = p(1, 1:ny)
        case (DIRICHLET)
            p(nx + 1, 1:ny) = -p(nx, 1:ny)
        case (NEUMANN)
            p(nx + 1, 1:ny) = p(nx, 1:ny)  
        end select

        ! bottom and top
        select case(BC(1, 2))
        case (PERIODIC)
            p(:, 0) = p(:, ny)
        case (DIRICHLET)
            p(:, 0) = -p(:, 1)
        case (NEUMANN)
            p(:, 0) = p(:, 1) 
        end select

        select case(BC(2, 2))
        case (PERIODIC)
            p(:, ny + 1) = p(:, 1)
        case (DIRICHLET)
            p(:, ny + 1) = -p(:, ny)
        case (NEUMANN)
            p(:, ny + 1) = p(:, ny) 
        end select


    end subroutine set_pressure_boundary
    
end module m_boundary