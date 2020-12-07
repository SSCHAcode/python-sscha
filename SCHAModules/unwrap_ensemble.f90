!
! This subroutine exploits symmetries to unwrap the ensemble.
! In this way the calculation with symmetries is much easier.
!
! It needs the crystal coordinates. 

subroutine unwrap_ensemble(coords, syms, irt, n_configs, new_coords, n_at, n_sym)
    implicit none
    double precision, dimension(n_configs, n_at * 3), intent(in) :: coords 
    integer, dimension(n_sym, 3, 3), intent(in) :: syms
    integer, dimension(n_sym, n_at), intent(in) :: irt 
    double precision, dimension(n_configs * n_sym, 3 * n_at), intent(out) :: new_coords
    integer :: n_configs, n_at, n_sym

    ! Apply the symmetries for each coords
    integer i_sym, i_conf, i_atm, x_atm
    integer k, h, start_id, end_id

    double precision, dimension(3) :: work
    
    do i_conf = 1, n_configs
        do i_sym = 1, n_sym
            do i_atm = 1, n_at
                x_atm = irt(i_sym, i_atm)
                start_id = 3* (i_atm - 1) + 1
                end_id = 3* (i_atm - 1) + 3



                ! Apply the symmetry and save the result to work
                work(:) = 0.0d0
                do k = 1,3 
                    work(k) = sum(syms(i_sym, k, :) * coords(i_conf, start_id : end_id))
                enddo

                ! if (i_conf == 1 .and. i_sym == 0) then
                !     print *, "ATM:", i_atm , "into", x_atm
                !     print *, "Vectors change from"
                !     print *, coords(i_conf, start_id : end_id)
                !     print *, "TO:"
                !     print *, work(:)
                !     print *, "SYM"
                !     print *, syms(i_sym, 1, :)
                !     print *, syms(i_sym, 2, :)
                !     print *, syms(i_sym, 3, :)
                !     print *, "TRY first:"
                !     print *, sum(syms(i_sym, 1, :) * coords(i_conf, start_id : end_id))
                ! endif   


                start_id = 3* (x_atm - 1) + 1
                end_id =  3* (x_atm - 1) + 3

                ! Now put the result on work in the correct atomic index
                new_coords(n_sym * (i_conf - 1) + i_sym, start_id : end_id) = work(:)
            enddo
        enddo
    enddo

end subroutine
