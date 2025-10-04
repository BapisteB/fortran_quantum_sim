program test_hashmap
    use m_statevec
    use m_gate
    !implicit none

    integer :: i
    type(t_statevec) :: stv
    integer, dimension(:), allocatable :: tgt_1, tgt_2
    type(t_base_gate) :: XBaseGate, YBaseGate, HBaseGate, CNOTBaseGate

    logical :: compiled_with_openmp = .false.

    !$ compiled_with_openmp = .true.

    if (compiled_with_openmp) then 
        write(*,*) 'OpenMP used...'
    else
        write(*,*) 'Openmp not used...'
    end if

    XBaseGate = make_gate("X", &
        reshape([cmplx(0, 0), cmplx(1, 0), &
        cmplx(1, 0), cmplx(0, 0)], [2, 2]))

    YBaseGate = make_gate("Y", &
        reshape([cmplx(0, 0), cmplx(0, -1), &
        cmplx(0, 1), cmplx(0, 0)], [2, 2]))

    HBaseGate = make_gate("H", &
        reshape([cmplx(inv_sqrt2, 0), cmplx(inv_sqrt2, 0), &
        cmplx(inv_sqrt2, 0), cmplx(-inv_sqrt2, 0)], [2,2]))

    CNOTBaseGate = make_gate("CNOT", &
        reshape([cmplx(1, 0), cmplx(0, 0), cmplx(0, 0), cmplx(0, 0), &
        cmplx(0, 0), cmplx(1, 0), cmplx(0, 0), cmplx(0, 0), &
        cmplx(0, 0), cmplx(0, 0), cmplx(0, 0), cmplx(1, 0), &
        cmplx(0, 0), cmplx(0, 0), cmplx(1, 0), cmplx(0, 0)], &
        [4, 4]))

    allocate(stv%amplitudes(2 ** 2))

    allocate(tgt_1(1))
    allocate(tgt_2(2))
    tgt_1(1) = 0
    tgt_2(1) = 0
    tgt_2(2) = 1

    do i = 1, 10000000
        call add_gate(stv, t_gate(tgt_1, HBaseGate))
    end do
    !call add_gate(stv, t_gate(tgt_2, CNOTBaseGate))
    print *, stv%last_gate_placed

    stv%amplitudes(1) = cmplx(1, 0)
    do i = 1, size(stv%amplitudes)
        print *, i - 1, stv%amplitudes(i)
    end do

    call compute_stv(stv)

    do i = 1, size(stv%amplitudes)
        print *, i - 1, stv%amplitudes(i)
    end do


end program test_hashmap

