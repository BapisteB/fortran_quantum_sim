module m_gate
    implicit none

    type :: t_base_gate
        character(len=10) :: name
        complex, dimension(:, :), allocatable :: matrix
    end type

    type :: t_gate
        integer, dimension(:), allocatable :: tgt
        type(t_base_gate) :: base_gate
    end type 

    real, parameter :: inv_sqrt2 = 1.0 / sqrt(2.0)

    contains

        pure function make_gate(name, mat) result(g)
            character(len=*), intent(in) :: name
            complex, intent(in) :: mat(:,:)
            type(t_base_gate) :: g

            g%name = name
            allocate(g%matrix(size(mat,1), size(mat,2)))
            g%matrix = mat
        end function make_gate

end module m_gate
