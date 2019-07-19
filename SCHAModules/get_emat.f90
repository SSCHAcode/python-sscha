
! This subroutine calculates the matrix containing polarization vectors
! masses and the lengths

subroutine get_emat (er, a, amass, ityp_sc, v3_log, transmode, e, n_mode, nat_sc, ntyp)

    implicit none
    
      double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
      double precision, dimension(n_mode), intent(in) :: a
      double precision, dimension(ntyp), intent(in) :: amass
      integer, dimension(nat_sc), intent(in) :: ityp_sc
      logical, intent(in) :: v3_log
      logical, dimension(n_mode), intent(in) :: transmode
      double precision, dimension(n_mode,3*nat_sc), intent(out) :: e
    
      integer :: mu, ka, i, nat_sc, n_mode, alpha, ntyp

      logical, parameter :: debug = .true.
    
      !nat_sc = size(er(:,1,1))
      !n_mode = size(a)

      if (debug) then
        print *, "=== DEBUG GET_EMAT ==="
        print *, "NAT_SC:", nat_sc 
        print *, "N_MODE:", n_mode 
        call flush()
      endif
    
      do mu = 1, n_mode
        ka = 0
        do i = 1, nat_sc
          do alpha = 1, 3
            ka = ka + 1
            if (transmode(mu)) then
              e(mu,ka) = 0.0d0
            else if (v3_log) then
              e(mu,ka) = er(i,mu,alpha) * a(mu) / sqrt(amass(ityp_sc(i)))
            else
              e(mu,ka) = er(i,mu,alpha) * sqrt(amass(ityp_sc(i))) / a(mu)
            end if
          end do
        end do
      end do
    
    end subroutine get_emat
    