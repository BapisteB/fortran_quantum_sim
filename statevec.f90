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
                allocate(tmp(capacity * 2))
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
            integer :: i, j, state
            integer, dimension(2) :: idxs
            complex, dimension(2) :: amps, new_amps
            complex, dimension(size(stv%amplitudes)) :: after_compute
            
            ! Foreach in gates
            !!! TODO
            do i = 1, stv%last_gate_placed
                gate = stv%gates(i)
                do state = 1, size(stv%amplitudes)
                    ! 1. Recover amplitudes
                        ! 1.1 Get indices
                        idxs = get_indices(gate%tgt, state - 1)
                        ! 1.2 Initialize state to multiply
                        ! TODO: SIMD copy ?
                        do j = 1, size(idxs)
                            amps(j) = stv%amplitudes(idxs(j) + 1)
                        end do
                        ! 1.3 Multiply
                        new_amps = matmul(gate%base_gate%matrix, amps)
                        ! 1.4 Put back
                        do j = 1, size(idxs)
                            after_compute(idxs(j) + 1) = new_amps(j)
                        end do
                end do
                stv%amplitudes = after_compute
            end do

        end subroutine

        function get_indices(targeted_bits, base_state)
            integer :: targeted_bits
            integer :: base_state
            integer, dimension(2) :: get_indices

            get_indices(1) = ibclr(base_state, targeted_bits)
            get_indices(2) = ibset(base_state, targeted_bits)
        end function

end module m_statevec
