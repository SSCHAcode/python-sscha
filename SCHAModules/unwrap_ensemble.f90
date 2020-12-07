!
! This subroutine exploits symmetries to unwrap the ensemble.
! In this way the calculation with symmetries is much easier.
!
! It needs the crystal coordinates. 
subroutine unwrap_ensemble(coords, syms. irt, n_configs, new_coords, n_at, n_sym)
    implicit none
    double precision, dimension(n_configs, n_at * 3), intent(in) :: coords 
    integer, dimension(n_sym, 3, 3), intent(in) :: syms
    integer, dimension(n_sym, n_at), intent(in) :: irt 
    double precision, dimension(n_configs * n_sym, 3 * n_at), intent(out) :: new_coords
    integer :: n_configs, n_at, n_sym

    ! Apply the symmetries for each coords
    integer i_sym, i_conf, i_atm, x_atm
    integer k, h

    double precision, dimension(3) :: work
    
    do i_conf = 1, n_configs
        do i_sym = 1, n_sym
            do i_atm = 1, n_at
                x_atm = irt(i_sym, i_atm)

                ! Apply the symmetry and save the result to work
                work(:) = 0.0d0
                do k = 1,3 
                    do h = 1,3 
                        work(k) = work(k) + syms(i_sym, k, h) * coords(i_conf, start_id(i_atm) + h)
                    enddo
                enddo

                ! Now put the result on work in the correct atomic index
                new_coords(n_sym * (i_conf - 1) + i_sym, start_id(x_atm) : end_id(x_atm)) = work(:)
            enddo
        enddo
    enddo

end subroutine

function start_id(x) result (y)
    integer, intent(in) :: x
    integer, intent(out) :: y

    y = 3* (x - 1) + 1
end function 

function end_id(x) result (y)
    integer, intent(in) :: x
    integer :: y

    y = 3* (x - 1) + 3
end function 