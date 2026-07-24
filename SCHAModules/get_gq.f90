
subroutine get_gq (a, wr, transmode, T, g, da, n_mode, iq)

    use thermodynamic
    implicit none
  

    double precision, dimension(iq, n_mode), intent(in) :: a, wr
    logical, dimension(iq, n_mode), intent(in) :: transmode
    double precision, intent(in) :: T
    double precision, dimension(iq, iq, n_mode,n_mode), intent(out) :: g                 
  
    integer :: mu, nu, iq, qmu, qnu
    integer :: n_mode
    double precision, dimension(iq,n_mode), intent(out) :: da

    logical, parameter :: debug = .false.


    if (debug) then
      print *, "=== DEBUG GET_G ==="
      print *, "N_MODE:", n_mode
      call flush() 
    end if

    call w_to_daq(wr, T, da, iq, n_mode)

    ! do i=1, iq
    !     print*, da(iq,:)
    ! end do
  
    ! Calculate  the matrix g that will enter in the final equation
    do qmu = 1, iq
        do mu = 1, n_mode
            do qnu = 1, iq
                do nu = 1, n_mode
                    if (transmode(qmu, mu) .or. transmode(qnu, nu)) then
                        g(qmu, qnu, mu, nu) = 0.0d0

                    else
                        ! 2. Calculate the difference between the two frequencies
                        ! We use a relative difference (abs(w1-w2)/w2) to be scale-independent
                        if (abs(wr(qmu, mu) - wr(qnu, nu)) < 0.00001d0 * wr(qnu, nu)) then
                            
                            ! DEGENERACY CASE: Frequencies are so close that we must use the Limit
                            ! This handles both (qmu,mu)==(qnu,nu) AND different modes with same frequency
                            g(qmu, qnu, mu, nu) = da(qmu, mu) / (wr(qmu, mu) * a(qmu, mu)**3.0d0)

                        else
                            ! GENERAL CASE: Frequencies are different enough to divide safely
                            g(qmu, qnu, mu, nu) = (a(qmu, mu)**2.0d0 - a(qnu, nu)**2.0d0) / &
                                                (wr(qmu, mu)**2.0d0 - wr(qnu, nu)**2.0d0) / &
                                                (a(qnu, nu)**2.0d0 * a(qmu, mu)**2.0d0)
                        end if
                    end if
!                     if (transmode(qmu, mu) .or. transmode(qnu, nu)) then
!                         g(qmu, qnu, mu,nu) = 0.0d0
!                     else if ((mu .eq. nu) .and. (qmu .eq. qnu)) then
!                         g(qmu, qnu, mu,nu) = da(qmu, mu) / wr(qmu, mu) / a(qmu, mu)**3.0d0
!                     else if &
! (((mu .ne. nu ) .or. (qmu .ne. qnu)) .and. abs((wr(qmu, mu)-wr(qnu, nu))/wr(qnu, nu)) .lt. 0.00001d0) then
!                     !g(mu,nu) = 0.0d0
!                         g(qmu, qnu, mu,nu) = da(qmu, mu) / wr(qmu, mu) / a(qmu, mu)**3.0d0
!                     else
!                         g(qmu, qnu, mu,nu) = (a(qmu, mu)**2.0d0 - a(qnu, nu)**2.0d0) / &
!                         (wr(qmu, mu)**2.0d0 - wr(qnu, nu)**2.0d0)  / &
!                         ( a(qnu, nu)**2.0d0 * a(qmu, mu)**2.0d0 )
!                     end if
                end do
            end do
        end do
    end do
  
  end subroutine get_gq
  