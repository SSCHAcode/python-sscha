! This subroutine calculates the rotated displacement and
! upsilon matrices needed for the third order force constants

subroutine get_ur_upsilon_matrices (a, er, transmode, amass, ityp_sc, u, &
        ur, eprod, nat_sc, n_mode, n_random, ntyp)

    use omp_lib
    use stochastic
    implicit none

    double precision, dimension(n_mode), intent(in) :: a
    double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
    logical, dimension(n_mode), intent(in) :: transmode
    double precision, dimension(ntyp), intent(in) :: amass
    integer, dimension(nat_sc), intent(in) :: ityp_sc
    double precision, dimension(n_random,nat_sc,3), intent(in) :: u
    double precision, dimension(n_random,n_mode), intent(out) :: ur
    double precision, dimension(n_mode,n_mode), intent(out) :: eprod

    integer :: nat_sc, n_mode, n_random, ntyp
  
    double precision, dimension(:,:), allocatable :: e, u2
    
    integer :: i, j, mu, alpha
    integer :: ka
  
    allocate(e(n_mode,n_mode))
    allocate(u2(n_random,n_mode))

    ! Calculate e polarization vectors with lengths and masses

    do mu = 1, n_mode
        ka = 0
        do i = 1, nat_sc
          do alpha = 1, 3
            ka = ka + 1
            if (transmode(mu)) then
              e(mu,ka) = 0.0d0
            else
              e(mu,ka) = er(i,mu,alpha) * sqrt(amass(ityp_sc(i))) / a(mu)
            end if
          end do
        end do
    end do

    ! Get displacements in a rank two tensor

    ka = 0

    do i = 1, nat_sc
      do alpha = 1, 3
        ka = ka + 1
        u2(:,ka) = u(:,i,alpha)
      end do
    end do

    ! Calculate product between two e matrices
    ! eprod is the Upsilon matrix (see Raffaello's paper)
    call dgemm('T','N',n_mode,n_mode,n_mode,1.0d0,e,n_mode,&
               e,n_mode,0.0d0,eprod,n_mode)

    ! Rotate displacementes
    ! Here ur are Upsilon * u2
    call dgemm('N', 'T', n_random, n_mode, n_mode, 1.0d0, &
            u2, n_random, eprod, n_mode, 0.0d0, ur, n_random)
    
    deallocate(e, u2)

  end subroutine get_ur_upsilon_matrices
