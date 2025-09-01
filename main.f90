program test_hashmap
    use m_statevec
    use m_gate
    !implicit none

    complex, dimension(:), allocatable :: amps
    type(t_gate), dimension(:), allocatable :: gates
    complex, dimension(2, 2) :: mat
    integer :: s
    type(t_statevec) :: stv

    type(t_gate) :: XGate = t_gate(0, XBaseGate)
    type(t_gate) :: YGate = t_gate(0, YBaseGate)
    integer :: i
    integer, dimension(2) :: y

    stv%amplitudes = amps
    allocate(stv%amplitudes(2 ** 2))
    
    call add_gate(stv, XGate)
    call add_gate(stv, XGate)
    print *, size(stv%gates)

    do i = 1, size(stv%gates)
        print *, stv%gates(i)%base_gate%name
    end do

    do i = 1, size(stv%amplitudes)
        print *, i, stv%amplitudes(i)
    end do

    stv%amplitudes(1) = cmplx(1, 0)
    call compute_stv(stv)

    do i = 1, size(stv%amplitudes)
        print *, i, stv%amplitudes(i)
    end do


end program test_hashmap

