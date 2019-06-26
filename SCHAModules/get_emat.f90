
! This subroutine calculates the matrix containing polarization vectors
! masses and the lengths

subroutine get_emat (er, a, amass, ityp_sc, v3_log, transmode, e)

    implicit none
    
      double precision, dimension(:,:,:), intent(in) :: er
      double precision, dimension(:), intent(in) :: a
      double precision, dimension(:), intent(in) :: amass
      integer, dimension(:), intent(in) :: ityp_sc
      logical, intent(in) :: v3_log
      logical, dimension(:), intent(in) :: transmode
      double precision, dimension(:,:), intent(out) :: e
    
      integer :: mu, ka, i, nat_sc, n_mode, alpha
    
      nat_sc = size(er(:,1,1))
      n_mode = size(a)
    
      do mu = 1, n_mode
        ka = 0
        do i = 1, nat_sc
          do alpha = 1, 3
            ka = ka + 1
            if (v3_log) then
              e(mu,ka) = er(i,mu,alpha) * a(mu) / sqrt(amass(ityp_sc(i)))
            else if (.not. v3_log .and. transmode(mu)) then
              e(mu,ka) = 0.0d0
            else
              e(mu,ka) = er(i,mu,alpha) * sqrt(amass(ityp_sc(i))) / a(mu)
            end if
          end do
        end do
      end do
    
    end subroutine get_emat
    