module m_statevec
    use m_gate
    implicit none

    type :: t_statevec
        ! Basically a circuit
        type(t_gate), allocatable :: gates(:)
        ! The statevector
        complex, allocatable :: amplitudes(:)

        integer :: last_gate_placed = 1
    end type

    contains
        
        ! probably better ways to do this but it's not the focus
        subroutine add_gate(stv, gate)
            type(t_statevec) :: stv
            type(t_gate), intent(in) :: gate
            type(t_gate), dimension(:), allocatable :: tmp
            integer :: capacity

            if (.not. allocated(stv%gates)) then
                allocate(stv%gates(1))
                stv%gates(stv%last_gate_placed) = gate
                return
            endif 

            capacity = size(stv%gates)

            if (stv%last_gate_placed == capacity) then
                allocate(tmp(capacity * 4))
                tmp(:capacity) = stv%gates
                deallocate(stv%gates)
                call move_alloc(tmp, stv%gates)
            endif 

            stv%last_gate_placed = stv%last_gate_placed + 1
            stv%gates(stv%last_gate_placed) = gate

        end subroutine

        subroutine compute_stv(stv)
            type(t_statevec) :: stv
            !complex, dimension(:), allocatable, intent(out) :: result
            type(t_gate) :: gate
            integer :: i, j, state, s_idxs
            integer, dimension(:), allocatable:: idxs
            complex, dimension(:), allocatable :: amps, new_amps
            complex, dimension(size(stv%amplitudes)) :: after_compute

            ! Foreach in gates
            do i = 1, stv%last_gate_placed
                after_compute(:) = 0
                gate = stv%gates(i)
                do state = 1, size(stv%amplitudes)
                    ! 1. Recover amplitudes
                    if (abs(real(stv%amplitudes(state))) < 1.0e-8 .and. abs(aimag(stv%amplitudes(state))) < 1.0e-8) then
                        cycle 
                    end if
                    ! 1.1 Get indices
                    ! Maybe this should return a mask list and find a way to do everything with only one operation afterwards ==>
                    ! list comprehension or something alike
                    idxs = get_indices(gate%tgt, state - 1)
                    ! 1.2 Initialize state to multiply
                    s_idxs = size(idxs)
                    allocate(amps(s_idxs))
                    allocate(new_amps(s_idxs))
                    !$OMP DO
                    do j = 1, s_idxs
                        amps(j) = stv%amplitudes(idxs(j) + 1)
                    end do
                    ! 1.3 Multiply
                    new_amps = matmul(gate%base_gate%matrix, amps)
                    ! 1.4 Put back
                    !$OMP DO
                    do j = 1, s_idxs
                        after_compute(idxs(j) + 1) = after_compute(idxs(j) + 1) + new_amps(j)
                    end do
                    deallocate(amps)
                    deallocate(new_amps)
                end do
                stv%amplitudes(:) = after_compute / sum(after_compute)
                ! Find better way to normalize
            end do

        end subroutine

        function get_indices(targeted_bits, base_state)
            integer, dimension(:) :: targeted_bits
            integer :: base_state, bit_idx, i, s, s_pow, current_state, mask
            integer, dimension(:), allocatable :: get_indices
            
            s = size(targeted_bits)
            s_pow = 2 ** s
            allocate(get_indices(s_pow))
            
            !$OMP DO
            do bit_idx = 0, s_pow - 1
                current_state = base_state
                !!$OMP DO
                do i = 1, s
                    mask = shiftl(1, targeted_bits(s - i + 1))
                    if (btest(bit_idx, i-1)) then
                        current_state = ior(current_state, mask)
                    else
                        current_state = iand(current_state, not(mask))
                    end if
                end do
                get_indices(bit_idx + 1) = current_state
            end do

        end function

end module m_statevec
