module m_gate
    implicit none

    type :: t_base_gate
        character(len=10) :: name
        complex, dimension(2, 2) :: matrix
    end type


    type :: t_gate
        integer :: tgt
        type(t_base_gate) :: base_gate
    end type 

    type(t_base_gate), parameter :: XBaseGate = t_base_gate("X", &
        reshape([cmplx(0, 0), cmplx(1, 0), &
        cmplx(1, 0), cmplx(0, 0)], [2, 2]))

    type(t_base_gate), parameter :: YBaseGate = t_base_gate("Y", &
        reshape([cmplx(0, 0), cmplx(0, -1), &
        cmplx(0, 1), cmplx(0, 0)], [2, 2]))
end module m_gate
