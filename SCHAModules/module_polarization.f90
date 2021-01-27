module polarization

use lattice
use settings

public :: pol_vectors
public :: pol_vectors_2
public :: pol_vectors_check
public :: get_phase
public :: inequivalent_to_all
public :: inequivalent_to_all_NN
public :: write_dyn
public :: inequivalent_to_all_jnuq
public :: impose_degeneracy
public :: write_fc 
public :: read_fc_supercell
public :: diag_fc
public :: test_dyn_mat
public :: phiqtophistar
public :: ewqtoewr
public :: get_harmonic_energy
public :: get_harmonic_pressure

contains

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the symmetrized generators in all
! the q's with all the coefficients.

subroutine fc_supercell( ghrtot, q, dyn_coeff, nred, nqs, tau, tau_sc, &
                         ityp, itau, phi_sc) 

  implicit none

  double complex, dimension(:,:,:,:,:,:), intent(in) :: ghrtot
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: nred, nqs
  double precision, dimension(:,:), intent(in) :: dyn_coeff
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(:,:,:,:), intent(out) :: phi_sc

  integer :: nq, nqirr, nat, natsc
  integer :: i, j, alpha, beta, iq, lim1, lim2, qtot, qirr, sigma 
  double precision, dimension(3) :: latvec
  double complex :: im, one, complex_number
  double precision :: twopi
  double complex, allocatable :: expf(:)

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0 

  nq    = size(ghrtot(:,1,1,1,1,1))
  nqirr = size(nred(:))
  nat   = size(ityp)
  natsc = nq * nat

  ! profile.3 : use expf
  allocate(expf(nq))


  do j = 1, natsc
    do i = 1, natsc
      latvec(:) = tau_sc(:,i) - tau(:,itau(i)) - tau_sc(:,j) + tau(:,itau(j))
      
      ! precompute the exponential and the division by nq which are the most 
      ! expensive parts of this loop. Also the loops are now in a muche better
      ! order, but a decent compiler can reorder them anyway
      do qtot = 1,nq
        expf(qtot) = exp( im * twopi * dot_product(q(:,qtot),latvec)) / DBLE(nq)
      enddo
      
      do beta = 1, 3
        do alpha = 1, 3
          complex_number = (0.0d0,0.0d0)
          do qirr = 1, nqirr    
            if (qirr .eq. 1) then
              lim1 = 1
              lim2 = nqs(1)
            else
              lim1 = sum(nqs(1:qirr-1)) + 1
              lim2 = sum(nqs(1:qirr))
            end if
            do sigma = 1, nred(qirr)
              do qtot = lim1, lim2
                complex_number = complex_number + dyn_coeff(qirr,sigma) * &
                                 expf(qtot) * ghrtot(qtot,sigma,alpha,beta,itau(i),itau(j))
              end do 
            end do
          end do 
          if (abs(aimag(complex_number)) .gt. 1.0d-5) then
            print *, complex_number
            print *, ''
            print *, ' ERROR: There are force constants in the supercell that   '
            print *, '        are complex. This is not possible.              '
            print *, '        Stopping...                                     '
            print *, ' 1 '             
            stop
          end if
          phi_sc(alpha,beta,i,j) = DBLE(complex_number)
        end do
      end do
    end do
  end do

  deallocate(expf)

end subroutine fc_supercell

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the dynamical matrices at each q.

subroutine fc_supercell_from_dyn ( phitot, q, nqs, tau, tau_sc, &
                                   ityp, itau, phitot_sc) 

  implicit none

  double complex, dimension(:,:,:,:,:), intent(in) :: phitot
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: nqs
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(:,:,:,:), intent(out) :: phitot_sc

  integer :: nq, nat, natsc
  integer :: i, j, alpha, beta, qtot
  double precision, dimension(3) :: latvec
  double complex :: im, one, complex_number
  double precision :: twopi

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0 

  nq    = size(phitot(:,1,1,1,1))
  nat   = size(ityp)
  natsc = nq * nat

  do i = 1, natsc
    do j = 1, natsc
      latvec(:) = tau_sc(:,i) - tau(:,itau(i)) - tau_sc(:,j) + tau(:,itau(j))
      do alpha = 1, 3
        do beta = 1, 3
          complex_number = (0.0d0,0.0d0)
          do qtot = 1, nq
            complex_number = complex_number +  &
                             exp( im * twopi * dot_product(q(:,qtot),latvec)) * &
                             phitot(qtot,alpha,beta,itau(i),itau(j)) / &
                             dble(nq) 
          end do 
          if (abs(aimag(complex_number)) .gt. 1.0d-5) then
            print *, complex_number
            print *, ''
            print *, ' ERROR: There are force constants in the supercell that   '
            print *, '        are complex. This is not possible.              '
            print *, '        Stopping...                                     '
            print *, ' 2 '             
            stop
          end if
          phitot_sc(alpha,beta,i,j) = real(complex_number)
        end do
      end do
    end do
  end do

end subroutine fc_supercell_from_dyn

! This subroutine gets the translation vectors of the
! supercell

subroutine get_latvec ( tau_sc, tau, itau, latvec )

  implicit none

  double precision, dimension(:,:), intent(in) :: tau_sc, tau
  integer, dimension(:), intent(in) :: itau
  double precision, dimension(:,:), intent(out) :: latvec

  integer :: ka, i, nat_sc

  nat_sc = size(tau_sc(1,:))

  ! Prepare list of lattice vectors

  ka = 0

  do i = 1, nat_sc
    if (itau(i) .ne. 1) cycle
    ka = ka + 1
    latvec(ka,:) = tau_sc(:,i) - tau(:,1)
  end do

end subroutine get_latvec

! This subroutine imposes the force constant matrix
! to have translational symmetry, \phi_st (R) = \phi_st (-R).
! It also averages all the values for all the atoms in the unit
! cell that are equivalent.

subroutine impose_trans ( phi, tau, tau_sc, ityp, itau, at_sc, latvec, phi_r )

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: phi
  double precision, dimension(:,:), intent(in) :: tau, tau_sc
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(3,3), intent(in) :: at_sc
  double precision, dimension(:,:), intent(in) :: latvec
  double precision, dimension(:,:,:,:,:), intent(out) :: phi_r

  integer :: nat, nat_sc, nr
  double precision, dimension(3) :: cholat, vect, diff
  double precision, dimension(27,3) :: superlatvec
  double precision :: prec
  logical, dimension(:), allocatable :: assigned
  integer :: ka, i, j, k, l, r, is, js, la, r1, r2
  
  prec = 1.0d-6

  nat    = size(tau(1,:))
  nat_sc = size(tau_sc(1,:))

  nr = nat_sc / nat
  
  allocate( assigned(nr) )

  ! Create the supercell lattice vectors

  ka = 0

  do i = -1, 1
    do j = -1, 1
      do k = -1, 1
        ka = ka + 1
        superlatvec(ka,:) = dble(i) * at_sc(:,1) + dble(j) * at_sc(:,2) + dble(k) * at_sc(:,3)
!        print *, superlatvec(ka,:)
      end do
    end do
  end do

  ! Average over all the possible replicas in the supercell

  do i = 1, nat
    do j = 1, nat
      diff(:) = tau(:,i) - tau(:,j)
      print *, ''
      print *, ' i j :: ', i, j
      print *, ''
      do r = 1, nr
        ka = 0
        print *, ''
        print *, ' latvec :: ', r, latvec(r,:)
        print *, ''
        phi_r (r,:,:,i,j) = 0.0d0
        do is = 1, nat_sc
          if (itau(is) .eq. i) then
            do js = 1, nat_sc
              do la = 1, 27
                vect(:) = diff(:) - (tau_sc(:,is) - latvec(r,:) - tau_sc(:,js) + superlatvec(la,:))
                if (dot_product(vect,vect) .lt. prec) then
                  ka = ka + 1
                  print *, ka
                  print *, is, js
                  phi_r (r,:,:,i,j) = phi_r (r,:,:,i,j) + phi(:,:,is,js)
                end if
              end do
            end do
          end if   
        end do  
        if ( ka .ne. nr ) then
          print *, ''
          print *, ' ERROR: The number of points in the supercell lattice'
          print *, '        is different from the one expected.  '
          print *, '        Stopping...                                     '
          stop
        end if
        phi_r(r,:,:,i,j) = phi_r(r,:,:,i,j) / dble(nr)
      end do
    end do
  end do
    
  ! Impose that the phi_r (r,:,:,i,j) = phi_r (-r,:,:,j,i)

  do i = 1, nat
    do j = i, nat
      assigned(:) = .false.
      do r1 = 1, nr
        if (assigned(r1)) cycle
        print *, '  R ', latvec(r1,:) 
        do r2 = 1, nr
          do la = 1, 27
            vect(:) = latvec(r1,:) + latvec(r2,:) + superlatvec(la,:)
            if (dot_product(vect,vect) .lt. prec) then
              print *, ' -R ', latvec(r2,:) 
              phi_r(r1,:,:,i,j) = 0.5d0*(phi_r(r1,:,:,i,j) + phi_r(r2,:,:,j,i)) 
              phi_r(r2,:,:,j,i) = phi_r(r1,:,:,i,j)                             
              assigned(r2) = .true.
            end if
          end do
        end do
      end do
    end do
  end do

  
!  ! Apply the symmetrization
!
!  do i = 1, nat_sc
!    do j = 1, nat_sc
!       cholat(:) = tau_sc(:,j) - tau_sc(:,i)
!       ! Search atom with - R
!       do k = 1, nat_sc
!         do l = 1, 27
!            vect(:) = (tau_sc(:,i) - cholat(:) - superlatvec(l,:)) - tau_sc(:,k) 
!            if (dot_product(vect,vect) .lt. prec)  then
!              print *, ' i j -j :: ', i, j, k
!              phi(:,:,i,j) = 0.5d0*(phi(:,:,i,j) + phi(:,:,i,k))
!            end if
!         end do 
!       end do 
!    end do
!  end do
       
end subroutine impose_trans

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the dynamical matrices at each q.

subroutine dyn_from_fc_r ( phi, q, nqs, tau, tau_sc, &
                           ityp, itau, latvec, dyn)

  implicit none

  double precision, dimension(:,:,:,:,:), intent(in) :: phi
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: nqs
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(:,:), intent(in) :: latvec
  double complex, dimension(:,:,:,:,:), intent(out) :: dyn

  integer :: nq, nat, natsc
  integer :: i, j, k, alpha, beta, qtot, R
  integer :: ka
  double precision, dimension(3) :: vecaux
  double complex :: im, one, complex_number
  double precision :: twopi, prec

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0

  prec = 1.0d-6

  nq    = size(dyn(:,1,1,1,1))
  nat   = size(ityp)
  natsc = nq * nat

  do qtot = 1, nq
    do i = 1, nat
      do j = 1, nat
        do alpha = 1, 3
          do beta = 1, 3
            complex_number = (0.0d0,0.0d0)
            do R = 1, nq
              complex_number = complex_number +  &
                               exp( - im * twopi * dot_product(q(:,qtot),latvec(R,:))) * &
                               phi(R,alpha,beta,i,j)
            end do
            dyn(qtot,alpha,beta,i,j) = complex_number
          end do
        end do
      end do
    end do
  end do

end subroutine dyn_from_fc_r

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the dynamical matrices at each q.

subroutine dyn_from_fc ( phitot_sc, q, nqs, tau, tau_sc, &
                                   ityp, itau, dyn) 

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: phitot_sc
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: nqs
  integer, dimension(:), intent(in) :: ityp, itau
  double complex, dimension(:,:,:,:,:), intent(out) :: dyn 

  integer :: nq, nat, natsc
  integer :: i, j, k, alpha, beta, qtot, R
  integer :: ka
  double precision, dimension(:,:), allocatable :: latvec
  double precision, dimension(3) :: vecaux
  double complex :: im, one, complex_number
  double precision :: twopi, prec

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0 

  prec = 1.0d-6

  nq    = size(dyn(:,1,1,1,1))
  nat   = size(ityp)
  natsc = nq * nat

  allocate(latvec(nq,3))

  ! Prepare list of lattice vectors

  ka = 0

  do i = 1, natsc
    if (itau(i) .ne. 1) cycle
    ka = ka + 1
    latvec(ka,:) = tau_sc(:,i) - tau(:,1)
  end do 

  ! Print list of lattice vectors

  do i = 1, nq
    print *, latvec(i,:)
  end do 

  do qtot = 1, nq
    do i = 1, nat
      do j = 1, nat
        do alpha = 1, 3
          do beta = 1, 3
            complex_number = (0.0d0,0.0d0)
            do R = 1, nq
              ! Check what the atom in the supercell is
              do k = 1, natsc
                vecaux = tau(:,i) + latvec(R,:) - tau_sc(:,k)
                if ( sqrt(dot_product(vecaux,vecaux)) .lt. prec ) then
                  complex_number = complex_number +  &
                                  exp( - im * twopi * dot_product(q(:,qtot),latvec(R,:))) * &
                                  phitot_sc(alpha,beta,k,j)
                end if
              end do
            end do
            dyn(qtot,alpha,beta,i,j) = complex_number
          end do
        end do
      end do
    end do
  end do

end subroutine dyn_from_fc

! This subroutine writes the force constants matrix in a suitable
! format for the Fourier interpolation

subroutine get_frc( phi_sc, supercell_size, tau, tau_sc, at, itau, frc)

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: phi_sc
  integer, dimension(3), intent(in) :: supercell_size
  double precision, dimension(:,:), intent(in) :: tau, tau_sc
  double precision, dimension(3,3), intent(in) :: at
  integer, dimension(:), intent(in) :: itau  
  double precision, dimension(:,:,:,:,:,:,:), intent(out) :: frc

  integer :: alpha, beta, i, j, l, m, n, sup1, sup2, sup3
  integer :: natsc, nat
  double precision, dimension(3) :: vect
  double precision, dimension(:,:,:,:,:), allocatable :: phi_auxx

  natsc = size(tau_sc(1,:))
  nat   = size(tau(1,:))

  allocate(phi_auxx(nat,nat,supercell_size(1),supercell_size(2),supercell_size(3)))

  do alpha = 1, 3
    do beta = 1, 3
      do i = 1, natsc
        do j = 1, nat
          vect(:) = tau_sc(:,i)-tau(:,itau(i))
          call asign_supercell_index_new(vect,at,l,m,n)
          phi_auxx(itau(i),j,l,m,n) = phi_sc(alpha,beta,i,j)
        end do
      end do
      do i = 1, nat 
        do j = 1, nat
          do sup3 = 1, supercell_size(3) 
            do sup2 = 1, supercell_size(2) 
              do sup1 = 1, supercell_size(1)
                frc(sup1,sup2,sup3,alpha,beta,i,j) = phi_auxx(i,j,sup1,sup2,sup3) 
              end do
            end do
          end do 
        end do  
      end do 
    end do
  end do

  deallocate(phi_auxx)

end subroutine get_frc

!! This subroutine calculates the force constant matrices in the
!! real space in the supercell from the dynamical matrices at each q.
!
!subroutine dyn_from_fc ( frc, q, tau, tau_sc, at, &
!                         itau, phi) 
!
!  implicit none
!
!  double precision, dimension(:,:,:,:,:,:,:), intent(in) :: frc
!  double precision, dimension(3), intent(in) :: q
!  double precision, dimension(:,:), intent(in) :: tau, tau_sc
!  double precision, dimension(3,3), intent(in) :: at
!  integer, dimension(:), intent(in) :: itau
!  double complex, dimension(:,:,:,:), intent(out) :: phi
!
!  integer :: nat, natsc, nq
!  integer :: i, j, alpha, beta, ka, l, m, n, sup1, sup2, sup3
!  integer :: choose
!  double precision, dimension(3) :: latvec, vec, vect
!  double complex :: im, one, complex_number
!  double precision :: twopi, arg
!
!  double precision, dimension(:,:), allocatable :: list_latvec
!  double complex, dimension(:,:,:,:,:), allocatable :: phi_list
!  double precision, dimension(:,:,:,:,:,:,:), allocatable :: frc
!  double precision, dimension(:,:,:,:,:), allocatable :: phi_auxx 
!
!  one    = (1.0d0,0.0d0)
!  im     = (0.0d0,1.0d0)
!  twopi  = 6.283185307179586d0 
!
!  natsc = size(tau_sc(1,:))
!  nat   = size(tau(1,:))
!  nq    = natsc / nat
!
!  allocate(list_latvec(nq,3))
!  allocate(phi_list(nq,3,3,nat,nat))
!
!  ALLOCATE(frc(supercell_size(1),supercell_size(1),supercell_size(3),3,3,nat,nat))
!  allocate(phi_auxx(nat,nat,supercell_size(1),supercell_size(2),supercell_size(3)))
!
!  do alpha = 1, 3
!    do beta = 1, 3
!      do i = 1, natsc
!        do j = 1, nat
!          vect(:) = tau_sc(:,i)-tau(:,itau(i))
!          call asign_supercell_index_new(vect,at,supercell_size,l,m,n)
!          phi_auxx(itau(i),j,l,m,n) = phi_sc(alpha,beta,i,j)
!        end do
!      end do
!      do i = 1, nat 
!        do j = 1, nat
!          WRITE (*,'(4i4)') alpha,beta,i,j
!          do sup3 = 1, supercell_size(3) 
!            do sup2 = 1, supercell_size(2) 
!              do sup1 = 1, supercell_size(1)
!                frc(sup1,sup2,sup3,alpha,beta,i,j) = phi_auxx(i,j,sup1,sup2,sup3) 
!!                WRITE (*,'(3i4,2x,1pe18.11)')   &
!!                  sup1,sup2,sup3,  frc(sup1,sup2,sup3,alpha,beta,i,j) * 2.0d0
!              end do
!            end do
!          end do 
!        end do  
!      end do 
!    end do
!  end do
!
!  stop
!
!!  ka = 0
!!
!!  do bt = 1, natsc
!!    if (itau(bt) .ne. 1) cycle
!!    ka = ka + 1
!!    list_latvec(ka,:) = tau_sc(:,bt) - tau(:,1)
!!    print *, list_latvec(ka,:)
!!  end do
!!
!!  do i = 1, nat
!!    do j = 1, nat
!!      do alpha = 1, 3
!!        do beta = 1, 3
!!          do ka = 1, nq
!!            complex_number = (0.0d0,0.0d0) 
!!            do at = 1, nq 
!!              arg = twopi * dot_product(qtot(:,at),list_latvec(ka,:))
!!              complex_number = complex_number + &
!!                               exp( im * arg) * phitot(at,alpha,beta,i,j) / dble(nq)
!!            end do 
!!            phi_list(ka,alpha,beta,i,j) = complex_number
!!          end do
!!          phi(alpha,beta,i,j) = (0.0,0.0d0)
!!          do ka = 1, nq
!!            arg = twopi * dot_product(q(:),list_latvec(ka,:))
!!            phi(alpha,beta,i,j) = phi(alpha,beta,i,j) + &
!!                                  exp( - im * arg) * phi_list(ka,alpha,beta,i,j)  
!!          end do 
!!        end do
!!      end do
!!    end do
!!  end do
!
!!  do bt = 1, natsc
!!    if (itau(bt) .ne. 1) cycle
!!    ka = ka + 1
!!    list_latvec(ka,:) = tau_sc(:,bt) - tau(:,1)
!!    print *, list_latvec(ka,:)
!!    do i = 1, nat
!!      do j = 1, nat
!!        do alpha = 1, 3
!!          do beta = 1, 3
!!            do at = 1, natsc
!!              vec(:) = tau_sc(:,at)-(tau(:,j)+list_latvec(ka,:))
!!              if (dot_product(vec,vec) .lt. 1.0d-7) then
!!                choose = at
!!                exit
!!              end if
!!            end do
!!            print *, choose   
!!            phi_list(ka,alpha,beta,i,j) =  phi_sc(alpha,beta,i,choose)  
!!          end do
!!        end do
!!      end do
!!    end do
!!  end do
!!
!!  print *, ' ka ', ka
!!  print *, ' q =', q(:)
!!
!
!!  do i = 1, nat
!!    do j = 1, nat
!!      do alpha = 1, 3
!!        do beta = 1, 3
!!          complex_number = (0.0d0,0.0d0)
!!          do at = 1, natsc
!!            if (itau(at) .ne. i) cycle
!!            do bt = 1, natsc
!!              if (itau(bt) .ne. j) cycle
!!              latvec(:) = tau_sc(:,at) - tau(:,i) - tau_sc(:,bt) + tau(:,j)
!!              arg = twopi * dot_product(q,latvec)
!!              complex_number = complex_number + &
!!                               exp( - im * arg) * phi_sc(alpha,beta,at,bt) / dble(nq)  
!!
!!            end do  
!!          end do 
!!          phi(alpha,beta,i,j) = complex_number
!!        end do
!!      end do
!!    end do
!!  end do
!
!end subroutine dyn_from_fc

! This subroutine links the modes that are relevant in q
! space with those that are relevant in real space

subroutine wq_to_wr(wq,wr,mu_start,mu_end,mu_relevant)

  implicit none

  double precision, dimension(:,:), intent(in) :: wq
  double precision, dimension(:), intent(in) :: wr
  integer, intent(in) :: mu_start, mu_end
  logical, dimension(:), intent(out) :: mu_relevant

  integer :: nq, nmodes_r, nmodes_q, mur, muq, iq
  double precision :: tol

  tol = 0.001

  nq = size(wq(:,1))
  nmodes_r = size(wr)
  nmodes_q = nmodes_r / nq
  
  mu_relevant = .false.

  do mur = 1, nmodes_r
    do iq = 1, nq
      do muq = 1, nmodes_q
        if (abs((wr(mur)-wq(iq,muq))*13.6058d0*8065.5d0) .lt. tol &
            .and. muq .ge. mu_start .and. muq .le. mu_end) then
          mu_relevant(mur) = .true. 
        end if
      end do
    end do
  end do  

end subroutine wq_to_wr

! This subroutine Fourier transforms the generators to the
! real space in the supercell. The generators in real space
! are written in a 3*nat_sc x 3*nat_sc matrix format

subroutine ghr_supercell( ghr, nred, q, nqs, ityp, itau, tau, tau_sc, ghr_sc) 

  implicit none

  double complex, dimension(:,:,:,:,:,:), intent(in) :: ghr
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: nred, nqs
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(:,:,:,:), intent(out) :: ghr_sc

  integer :: nq, nat, natsc, nqirr
  integer :: i, j, alpha, beta, iq, qirr, sigma, lim1, lim2
  double precision, dimension(3) :: latvec
  double complex :: im, one, phase, z         
  double precision :: twopi

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0 

  nqirr = size(nred)
  natsc = size(ityp)
  nat   = size(tau(1,:))
  nq    = natsc / nat

  do qirr = 1, nqirr
    if (qirr .eq. 1) then
      lim1 = 1
      lim2 = nqs(qirr)
    else
      lim1 = sum(nqs(1:qirr-1)) + 1
      lim2 = sum(nqs(1:qirr))
    end if
    do sigma = 1, nred(qirr)
      do i = 1, natsc
        do j = 1, natsc
          latvec(:) = tau_sc(:,i) - tau(:,itau(i)) - tau_sc(:,j) + tau(:,itau(j))
          do alpha = 1, 3
            do beta = 1, 3
              z = (0.0d0,0.0d0)
              do iq = lim1, lim2
                phase = exp( im * twopi * dot_product(q(:,iq),latvec)) / dble(nq) 
                z = z + phase * ghr(iq,sigma,alpha,beta,itau(i),itau(j)) 
              end do 
              if (abs(aimag(z)) .gt. 1.0d-5) then
                print *, z             
                print *, ''
                print *, ' ERROR: There are force constants in the supercell that   '
                print *, '        are complex. This is not possible.              '
                print *, '        Stopping...                                     '
                print *, ' 3 '             
                stop
              end if
              ghr_sc(qirr,sigma,3*(i-1)+alpha,3*(j-1)+beta) = real(z)
            end do
          end do 
        end do
      end do
    end do
  end do

end subroutine ghr_supercell

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the symmetrized generators in all
! the q's with all the coefficients.

subroutine get_gamma_mu_nu ( ghr_sc, er, nred, ityp, amass, gamma_mu_nu) 

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: ghr_sc
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:,:,:), intent(in) :: er
  integer, dimension(:), intent(in) :: nred
  double precision, dimension(:), intent(in) :: amass
  double precision, dimension(:,:,:,:), intent(out) :: gamma_mu_nu 

  integer :: nqirr, natsc, nmodes
  integer :: i, j, alpha, beta, qirr, sigma, mu
  double precision, dimension(3) :: latvec
  double precision, dimension(:,:), allocatable :: e_mat, e_aux, prod, g, gr

  nqirr  = size(nred(:))
  natsc  = size(ityp)
  nmodes = natsc * 3
 
  allocate(e_mat(nmodes,nmodes))
  allocate(e_aux(nmodes,nmodes))
  allocate(prod(nmodes,nmodes))
  allocate(g(nmodes,nmodes))
  allocate(gr(nmodes,nmodes))

  do i = 1, natsc
    do alpha = 1, 3
      do mu = 1, nmodes
        e_mat(3*(i-1)+alpha,mu) = er(i,mu,alpha) / dsqrt(amass(ityp(i))) 
      end do
    end do
  end do

  e_aux = transpose(e_mat)

  do qirr = 1, nqirr    
    do sigma = 1, nred(qirr)
      g(:,:) = ghr_sc(qirr,sigma,:,:)
      call dgemm('N','N',nmodes,nmodes,nmodes,1.0d0,g, &
                 nmodes,e_mat,nmodes,0.0d0,prod,nmodes)
      call dgemm('N','N',nmodes,nmodes,nmodes,1.0d0,e_aux, &
                 nmodes,prod,nmodes,0.0d0,gr,nmodes)
      gamma_mu_nu(qirr,sigma,:,:) = gr(:,:)
       !prod = matmul(ghr_sc(qirr,sigma,:,:),e_aux)
       !gamma_mu_nu(qirr,sigma,:,:) = matmul(transpose(e_mat),prod)
    end do
  end do

  deallocate(e_mat, e_aux, prod, g, gr)

end subroutine get_gamma_mu_nu

! This subroutine calculates the force constant matrices in the
! real space in the supercell from the symmetrized generators in all
! the q's with all the coefficients.

subroutine get_gamma_mu_nu_new ( ghr_sc, er, nred, ityp, amass, qirr, & 
                                 sigma, gamma_mu_nu) 

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: ghr_sc
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:,:,:), intent(in) :: er
  integer, dimension(:), intent(in) :: nred
  double precision, dimension(:), intent(in) :: amass
  integer, intent(in) :: qirr, sigma
  double precision, dimension(:,:), intent(out) :: gamma_mu_nu 

  integer :: nqirr, natsc, nmodes
  integer :: i, j, alpha, beta, mu
  double precision, dimension(3) :: latvec
  double precision, dimension(:,:), allocatable :: e_mat, e_aux, prod, g, gr

  nqirr  = size(nred(:))
  natsc  = size(ityp)
  nmodes = natsc * 3
 
  allocate(e_mat(nmodes,nmodes))
  allocate(e_aux(nmodes,nmodes))
  allocate(prod(nmodes,nmodes))
  allocate(g(nmodes,nmodes))
  allocate(gr(nmodes,nmodes))

  do i = 1, natsc
    do alpha = 1, 3
      do mu = 1, nmodes
        e_mat(3*(i-1)+alpha,mu) = er(i,mu,alpha) / dsqrt(amass(ityp(i))) 
      end do
    end do
  end do

  e_aux = transpose(e_mat)

!  do qirr = 1, nqirr    
!    do sigma = 1, nred(qirr)
      g(:,:) = ghr_sc(qirr,sigma,:,:)
      call dgemm('N','N',nmodes,nmodes,nmodes,1.0d0,g, &
                 nmodes,e_mat,nmodes,0.0d0,prod,nmodes)
      call dgemm('N','N',nmodes,nmodes,nmodes,1.0d0,e_aux, &
                 nmodes,prod,nmodes,0.0d0,gr,nmodes)
      gamma_mu_nu(:,:) = gr(:,:)
       !prod = matmul(ghr_sc(qirr,sigma,:,:),e_aux)
       !gamma_mu_nu(qirr,sigma,:,:) = matmul(transpose(e_mat),prod)
!    end do
!  end do

  deallocate(e_mat, e_aux, prod, g, gr)

end subroutine get_gamma_mu_nu_new 

! Subroutine that creates the dynamical matrix from the force
! constants in the supercell and diagonalizes it. It gives as
! output the polarization vectors and the phonon frequencies.

subroutine diag_fc(fc,mass,type_atoms,w,e)  

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: fc
  double precision, dimension(:), intent(in) :: mass
  integer, dimension(:), intent(in) :: type_atoms
  double precision, dimension(:), intent(out) :: w
  double precision, dimension(:,:,:), intent(out) :: e

  integer :: nat, ntype, nmodes
  integer :: mu, nu, i, j, alpha, beta, info
  double precision, dimension(:,:), allocatable :: dyn
  double precision, dimension(:), allocatable :: work, eig

  nat    = size(type_atoms) 
  ntype  = size(mass) 
  nmodes = 3*nat

  allocate(dyn(nmodes,nmodes))
  allocate(work(3*nmodes-1))
  allocate(eig(nmodes))

  ! Create dynamical matrix
  do i = 1, nat
    do j = 1, nat
      do alpha = 1, 3
        do beta = 1, 3
          mu = 3*(i-1) + alpha
          nu = 3*(j-1) + beta  
          dyn(mu,nu) = fc(i,j,alpha,beta) / &
                       sqrt(mass(type_atoms(i))*mass(type_atoms(j)))
        end do
      end do
    end do 
  end do
 
  ! Diagonalize dynamical matrix
  call dsyev('V','U',nmodes,dyn,nmodes,eig,work,size(work),info)

  ! Assign polarization vectors
  do i = 1, nat
    do alpha = 1, 3
      mu = 3*(i-1) + alpha
      do nu = 1, nmodes
       ! e(i,nu,alpha) = dyn(nu,mu) 
        e(i,nu,alpha) = dyn(mu,nu) 
      end do
    end do
  end do

  ! Assign frequencies
  do i = 1, nmodes
    if (eig(i) < 0) then
      w(i) = - sqrt(-eig(i))
    else
      w(i) = sqrt(eig(i))
    end if
  end do

  deallocate(dyn)
  deallocate(work)
  deallocate(eig)
  
end subroutine diag_fc 

! Subroutine that reads the force constants in the supercell.
! The force constants are the dynamical matrix at Gamma point.    
! They are supposed to be given in a format of quantum-espressso 
! dynamical matrix. 

subroutine read_fc_supercell(file_fc,fc)

  character (len=*), intent(in) :: file_fc
  double precision, dimension(:,:,:,:), intent(out) :: fc

  integer :: natom
  integer :: i, j, k, l, atom1, atom2
  double precision :: kaka

  natom = size(fc(:,1,1,1))

  open (unit = 1, file = file_fc)

  do i = 1, natom
    do j = 1, natom
      read (unit = 1, fmt=*) atom1, atom2
      do k = 1, 3
        read (unit = 1, fmt=*) fc(atom1,atom2,k,1), kaka, &
                               fc(atom1,atom2,k,2), kaka, &
                               fc(atom1,atom2,k,3), kaka
      end do
    end do
  end do

  close (unit=1)

end subroutine read_fc_supercell 

! This subroutine creates the dynamical matrix from the 
! frequencies and the polarization vectors and diagonalizes it.
! This is useful to see whether the polarization vectors are
! correct.

subroutine test_dyn_mat(e_diag,w)

  double precision, dimension(:,:,:), intent(in) :: e_diag
  double precision, dimension(:), intent(in) :: w

  integer :: nmodes, natoms, i, j, mu, info
  double precision, dimension(:,:), allocatable :: dyn_mat
  double precision, dimension(:), allocatable :: eig    
  double complex, dimension(:,:), allocatable :: pol_e
  double precision, dimension(:), allocatable :: work

  natoms = size(e_diag(:,1,1))
  nmodes = 3*natoms

  allocate(dyn_mat(nmodes,nmodes))
  allocate(pol_e(nmodes,nmodes))
  allocate(eig(nmodes))
  allocate(work(3*nmodes-1))

  ! Read polarization vectors

  do i = 1, nmodes
    do j = 1, natoms
      pol_e(i,(3*(j-1)+1):(3*j)) = CMPLX(e_diag(j,i,:),kind=8)
    end do
  end do

  ! Create dynamical matrix

  do i = 1, nmodes
    do j = 1, nmodes
      dyn_mat(i,j) = 0.0d0
      do mu = 1, nmodes
        dyn_mat(i,j) = dyn_mat(i,j) + w(mu) * pol_e(mu,i) * pol_e(mu,j)
      end do
    end do
  end do

  ! Diagonalize dynamical matrix

  call dsyev('V','U',nmodes,dyn_mat,nmodes,eig,work,3*nmodes-1,info)

  do i = 1, nmodes
    write (unit=1,fmt=*) '  w = ', eig(i) * 219467.0943D0, ' cm^-1' 
    do j = 1, natoms
      !write (unit=1,fmt=*) '     at ', j, ':', dyn_mat(i,3*(j-1)+1:3*j) 
      write (unit=1,fmt=*) '     at ', j, ':', dyn_mat(3*(j-1)+1:3*j,i) 
    end do
  end do
 
  deallocate(dyn_mat)
  deallocate(pol_e)
  deallocate(eig)
  deallocate(work)

end subroutine test_dyn_mat

! This subroutine performs a check to see whether the polarization
! vectors are diagonal. However it does not modify the polarization
! vectors, only normalizes them

subroutine e_check(e)

  double precision, dimension(:,:,:), intent(in) :: e
 
  integer :: n_atom, n_type, n_mode
  integer :: i, j, k
  double precision, dimension(:,:), allocatable :: pol_e

  n_atom = size(e(:,1,1))
  n_mode = 3 * n_atom

  allocate(pol_e(n_mode,n_mode))

  ! Read polarization vectors

  do i = 1, n_mode
    do j = 1, n_atom
      pol_e(i,(3*(j-1)+1):(3*j)) = e(j,i,:)
    end do
  end do

  ! Check whether the eigenvectors satisfy the orthogonalization conditions with
  ! respect to the modes.

  write (unit=1,fmt=*) ''
  write (unit=1,fmt=*) '  CHECK 1: Are the eigenvectors orthogonal with respect to'
  write (unit=1,fmt=*) '           the mode index?                                '
  write (unit=1,fmt=*) ''
  do i = 1, 3*n_atom
    do j = 1, 3*n_atom
      write (unit=1,fmt='(a,i4,a,i4,a,f18.15)') '    e', i,' x e', j,' =  ', &
                                               dot_product(pol_e(i,:),pol_e(j,:)) 
    end do
  end do

  write (unit=1,fmt=*) ''
  write (unit=1,fmt=*) '  CHECK 2: Are the eigenvectors orthogonal with respect to'
  write (unit=1,fmt=*) '           the cartesian index?                           '
  write (unit=1,fmt=*) ''
  do i = 1, 3*n_atom
    do j = 1, 3*n_atom
      write (unit=1,fmt='(a,i4,a,i4,a,f18.15)') '    e', i,' x e', j,' =  ', &
                                               dot_product(pol_e(:,i),pol_e(:,j))
    end do
  end do

  deallocate(pol_e)

end subroutine e_check

! This subroutine performs a check to see whether the polarization
! vectors are diagonal. However it does not modify the polarization
! vectors, only normalizes them

subroutine e_check_complex(e)

  double complex, dimension(:,:,:), intent(in) :: e
 
  integer :: n_atom, n_type, n_mode
  integer :: i, j, k
  double complex, dimension(:,:), allocatable :: pol_e

  n_atom = size(e(1,:,1))
  n_mode = 3 * n_atom

  allocate(pol_e(n_mode,n_mode))

  ! Read polarization vectors

  do i = 1, n_mode
    do j = 1, n_atom
      pol_e(i,(3*(j-1)+1):(3*j)) = e(i,j,:)
    end do
  end do

  ! Check whether the eigenvectors satisfy the orthogonalization conditions with
  ! respect to the modes.

  write (unit=4,fmt=*) ''
  write (unit=4,fmt=*) '  CHECK 1: Are the eigenvectors orthogonal with respect to'
  write (unit=4,fmt=*) '           the mode index?                                '
  write (unit=4,fmt=*) ''
  do i = 1, 3*n_atom
    do j = 1, 3*n_atom
      write (unit=4,fmt='(a,i4,a,i4,a,2f16.12)') '    e', i,' x e', j,' =  ', &
                                               dot_product(pol_e(i,:),pol_e(j,:)) 
    end do
  end do

  write (unit=4,fmt=*) ''
  write (unit=4,fmt=*) '  CHECK 2: Are the eigenvectors orthogonal with respect to'
  write (unit=4,fmt=*) '           the cartesian index?                           '
  write (unit=4,fmt=*) ''
  do i = 1, 3*n_atom
    do j = 1, 3*n_atom
      write (unit=4,fmt='(a,i4,a,i4,a,2f16.12)') '    e', i,' x e', j,' =  ', &
                                               dot_product(pol_e(:,i),pol_e(:,j))
    end do
  end do

  deallocate(pol_e)

end subroutine e_check_complex 

! This function gets the phi phase that makes a complex 
! number real, so that a x e^(i phi) is real  

function phase(a) result(phase_a)            

  double complex, intent(in) :: a
  double precision :: halfpi
  double precision :: phase_a

  halfpi = 1.5707963267948966d0

  if (real(a) .eq. 0.0d0 .and. aimag(a) .gt. 0.0d0) then
    phase_a = - halfpi
  else if (real(a) .eq. 0.0d0 .and. aimag(a) .lt. 0.0d0) then
    phase_a = halfpi
  else if (real(a) .eq. 0.0d0 .and. aimag(a) .eq. 0.0d0) then
    phase_a = 0.0d0 
  else  
    phase_a = atan( - aimag(a) / real(a))   
  end if

end function phase

! This subroutine writes down the force constants in qe format

subroutine write_fc(flfrc,e_diag_over_mass,w,x_atoms_supercell,x_atoms_prim,&
                    primlatt_vec,supercell_size)

  character (len=20), intent(in) :: flfrc
  double complex, dimension(:,:,:), intent(in) :: e_diag_over_mass
  double precision, dimension(:), intent(in) :: w
  double precision, dimension(:,:), intent(in) :: x_atoms_supercell
  double precision, dimension(:,:), intent(in) :: x_atoms_prim
  double precision, dimension(3,3), intent(in) :: primlatt_vec
  integer, dimension(3), intent(in) :: supercell_size

  integer :: alpha, beta, i, j, natoms_supercell, natoms_prim, nmodes
  double precision, dimension(:,:,:,:), allocatable :: phi 
  integer :: s, l, m, n, this_atom, sup1, sup2, sup3
  double precision :: fc

  natoms_supercell = size(x_atoms_supercell(:,1))
  natoms_prim = size(x_atoms_prim(:,1))
  nmodes = size(w)
   
  allocate(phi(natoms_prim,supercell_size(1),supercell_size(2),supercell_size(3)))

  open (unit=13, file=flfrc)
  
  do alpha = 1, 3
    do beta = 1, 3
      do i = 1, natoms_prim
        phi = 0.0d0
        do j = 1, natoms_supercell
          if (x_atoms_supercell(j,1) .eq. x_atoms_prim(i,1) .and. &
              x_atoms_supercell(j,2) .eq. x_atoms_prim(i,2) .and. &
              x_atoms_supercell(j,3) .eq. x_atoms_prim(i,3) ) then
            this_atom = j
            exit
          end if 
        end do 
        do j = 1, natoms_supercell
          call asign_supercell_index(x_atoms_supercell(j,:),x_atoms_prim,&
                                     primlatt_vec,supercell_size,s,l,m,n)
          fc = sum(real(e_diag_over_mass(this_atom,:,alpha)) * &
                   real(e_diag_over_mass(j,:,beta)) * w(:)**2.0d0)
          phi(s,l,m,n) = fc
        end do 
        do j = 1, natoms_prim
          WRITE (13,'(4i4)') alpha,beta,i,j
          do sup3 = 1, supercell_size(1) 
            do sup2 = 1, supercell_size(2) 
              do sup1 = 1, supercell_size(3)
                WRITE (13,'(3i4,2x,1pe18.11)')   &
                  sup1,sup2,sup3,phi(j,sup1,sup2,sup3) * 2.0d0
              end do
            end do
          end do 
        end do  
      end do 
    end do
  end do

  close (unit=13)
  
  deallocate(phi)

end subroutine write_fc

! This subroutine writes down the force constants in qe format

subroutine write_fc_new(flfrc,phi_sc,tau_sc,tau,at,itau,supercell_size)

  character (len=20), intent(in) :: flfrc
  double precision, dimension(:,:,:,:), intent(in) :: phi_sc
  double precision, dimension(:,:), intent(in) :: tau_sc, tau
  double precision, dimension(3,3), intent(in) :: at
  integer, dimension(:), intent(in) :: itau
  integer, dimension(3), intent(in) :: supercell_size

  integer :: alpha, beta, i, j, nat_sc, nat
  double precision, dimension(:,:,:,:,:), allocatable :: phi 
  integer :: s, l, m, n, this_atom, sup1, sup2, sup3
  double precision :: fc
  double precision, dimension(3) :: vect

  nat_sc = size(tau_sc(1,:))
  nat    = size(tau(1,:))
   
  allocate(phi(nat,nat,supercell_size(1),supercell_size(2),supercell_size(3)))

  open (unit=13, file=flfrc)
  
  do alpha = 1, 3
    do beta = 1, 3
!----
      do i = 1, nat_sc
        do j = 1, nat
          vect(:) = tau_sc(:,i)-tau(:,itau(i))
          call asign_supercell_index_new(vect,at,l,m,n)
          fc = phi_sc(alpha,beta,i,j)
          phi(itau(i),j,l,m,n) = fc
        end do
      end do
      do i = 1, nat 
!----
!      do i = 1, nat
!        do j = 1, nat_sc
!          vect(:) = tau_sc(:,j)-tau(:,itau(j))
!          call asign_supercell_index_new(vect,at,supercell_size,l,m,n)
!          fc = phi_sc(alpha,beta,i,j)
!          phi(itau(j),l,m,n) = fc
!        end do 
        do j = 1, nat
          WRITE (13,'(4i4)') alpha,beta,i,j
          do sup3 = 1, supercell_size(3) 
            do sup2 = 1, supercell_size(2) 
              do sup1 = 1, supercell_size(1)
                WRITE (13,'(3i4,2x,1pe18.11)')   &
                  sup1,sup2,sup3,phi(i,j,sup1,sup2,sup3) * 2.0d0
!                  sup1,sup2,sup3,phi(j,i,sup1,sup2,sup3) * 2.0d0
              end do
            end do
          end do 
        end do  
      end do 
    end do
  end do

  close (unit=13)
  
  deallocate(phi)

end subroutine write_fc_new

! This subroutine writes the dynamical matrix on a file in quantum 
! espresso format using the generators and the coefficients

subroutine write_dyn_coeff (ghrtot, q, ntyp, nred, nqs, dyn_coeff, ityp, amass, &
                            population, fildyn_prefix, ibrav, celldm, tau, &
                            type_name, at, lrigid, epsil, zeu)
  use stochastic, only : population_idx
  implicit none

  double complex, dimension(:,:,:,:,:,:), intent(in) :: ghrtot
  double precision, dimension(:,:), intent(in) :: q
  double precision, dimension(:,:), intent(in) :: at
  integer, intent(in) :: ntyp
  integer, dimension(:), intent(in) :: nred, nqs
  double precision, dimension(:,:), intent(in) :: dyn_coeff
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:), intent(in) :: amass
  integer, intent(in) :: population
  character (len=50), intent(in) :: fildyn_prefix
  integer, intent(in) :: ibrav
  double precision, dimension(6), intent(in) :: celldm
  double precision, dimension(:,:), intent(in) :: tau 
  character (len=3), dimension(:), intent(in) :: type_name
  logical, intent(in) :: lrigid
  double precision, dimension(3,3), intent(in) :: epsil
  double precision, dimension(:,:,:), intent(in) :: zeu 

  integer :: nq, nqirr, nat, natsc, na, nt
  integer :: i, j, alpha, beta, iq, lim1, lim2, qtot, qirr, sigma 
  double complex, dimension(:,:,:,:), allocatable :: phi
  double complex, dimension(:,:), allocatable :: phi2
  character (len=200) :: fildyn
  double precision, dimension(:), allocatable :: w2 
  logical :: gamma_point

  character(len=6), EXTERNAL :: int_to_char
  CHARACTER(len=256),EXTERNAL :: trimcheck

  nq    = size(ghrtot(:,1,1,1,1,1))
  nqirr = size(nred(:))
  nat   = size(ghrtot(1,1,1,1,1,:))

  allocate(phi(3,3,nat,nat))
  allocate(phi2(3*nat,3*nat))
  allocate(w2(3*nat))

  do qirr = 1, nqirr
    ! Check if the point is equal to gamma
    gamma_point = .false.
    if (q(1,qirr) .eq. 0.0d0 .and. q(2,qirr) .eq. 0.0d0 .and. q(3,qirr) .eq. 0.0d0) then
      gamma_point = .true.
    end if    
    fildyn = trim(fildyn_prefix)//trim(int_to_char(qirr))//population_idx(population)
    open (unit=27,file=fildyn)
    ! Write header
    WRITE (unit=27,fmt='("Dynamical matrix file")')
    WRITE (unit=27,fmt='(a)') 'File generated with the SSCHA code by Ion Errea'
    WRITE (unit=27,fmt='(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
    IF(ibrav==0)THEN
      WRITE (27,'("Basis vectors")')
      WRITE (27,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
    ENDIF
    DO nt = 1, ntyp
      WRITE (unit=27,fmt=*) nt, ' ''', type_name(nt) , ' '' ', amass(nt)/2.0d0
    ENDDO
    DO na = 1, nat 
      WRITE (unit=27,fmt='(2i5,3f18.10)') na, ityp (na) , (tau (j, na) , j = 1, 3)
    ENDDO
    ! End header
    if (qirr .eq. 1) then
      lim1 = 1
      lim2 = nqs(qirr)
    else
      lim1 = sum(nqs(1:qirr-1)) + 1
      lim2 = sum(nqs(1:qirr))
    end if
    do qtot = lim1, lim2
      write (unit=27,fmt=9000) (q(i,qtot), i = 1, 3) 
      phi = (0.0d0,0.0d0)
      do i = 1, nred(qirr)
        phi = phi + dyn_coeff(qirr,i) * ghrtot(qtot,i,:,:,:,:)
      end do
      phi = 2.0d0 * phi
      if (qtot .eq. lim1) then
        do i = 1, nat
          do j = 1, nat
            do alpha = 1, 3
              do beta = 1, 3
                phi2(3*(i-1)+alpha,3*(j-1)+beta) = phi(alpha,beta,i,j)
              end do
            end do
          end do
        end do
      end if
      do i = 1, nat
        do j = 1, nat
          write (unit=27, fmt='(2i5)') i, j
          do alpha = 1, 3
            write (unit=27, fmt='(3(2f12.8,2x))') (phi(alpha,beta,i,j), beta=1,3)
          end do
        end do
      end do 
    end do
    ! Print effective charges in case there are     
    if (gamma_point .and. lrigid) then
      call write_epsilon_and_zeu (zeu, epsil, nat, 27)
    end if
    ! Diagonalize and print frequencies and polarization vectors
    call dyndia (q(:,lim2), 3*nat, nat, ntyp, ityp, amass/1822.8884482774586D0, &
                 27, phi2, w2)
    close (unit=27)
  end do

  deallocate(phi, phi2, w2)

 9000 format(/,5x,'Dynamical  Matrix in cartesian axes', &
         &       //,5x,'q = ( ',3f14.9,' ) ',/)

end subroutine write_dyn_coeff 

! This subroutine writes the dynamical matrix on a file in quantum 
! espresso format

subroutine write_dyn (phitot, q, ntyp, nqs, ityp, amass, &
                            fildyn_prefix, ibrav, celldm, tau, &
                            type_name, at, lrigid, epsil, zeu)

  implicit none

  double complex, dimension(:,:,:,:,:), intent(in) :: phitot
  double precision, dimension(:,:), intent(in) :: q
  double precision, dimension(:,:), intent(in) :: at
  integer, intent(in) :: ntyp
  integer, dimension(:), intent(in) :: nqs
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:), intent(in) :: amass
  character (len=50), intent(in) :: fildyn_prefix
  integer, intent(in) :: ibrav
  double precision, dimension(6), intent(in) :: celldm
  double precision, dimension(:,:), intent(in) :: tau 
  character (len=3), dimension(:), intent(in) :: type_name
  logical, intent(in) :: lrigid
  double precision, dimension(3,3), intent(in) :: epsil
  double precision, dimension(:,:,:), intent(in) :: zeu 

  integer :: nq, nqirr, nat, natsc, na, nt
  integer :: i, j, alpha, beta, iq, lim1, lim2, qtot, qirr, sigma 
  double complex, dimension(:,:,:,:), allocatable :: phi
  double complex, dimension(:,:), allocatable :: phi2
  character (len=200) :: fildyn
  double precision, dimension(:), allocatable :: w2 
  logical :: gamma_point

  character(len=6), EXTERNAL :: int_to_char
  CHARACTER(len=256),EXTERNAL :: trimcheck

  nq    = size(phitot(:,1,1,1,1))
  nqirr = size(nqs)
  nat   = size(phitot(1,1,1,1,:))

  allocate(phi(3,3,nat,nat))
  allocate(phi2(3*nat,3*nat))
  allocate(w2(3*nat))

  do qirr = 1, nqirr
    ! Check if the point is equal to gamma
    gamma_point = .false.
    if (q(1,qirr) .eq. 0.0d0 .and. q(2,qirr) .eq. 0.0d0 .and. q(3,qirr) .eq. 0.0d0) then
      gamma_point = .true.
    end if    
    fildyn = trim(fildyn_prefix) // int_to_char(qirr)
    open (unit=27,file=fildyn)
    ! Write header
    WRITE (unit=27,fmt='("Dynamical matrix file")')
    WRITE (unit=27,fmt='(a)') 'File generated with the SSCHA code by Ion Errea'
    WRITE (unit=27,fmt='(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
    IF(ibrav==0)THEN
      WRITE (27,'("Basis vectors")')
      WRITE (27,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
    ENDIF
    DO nt = 1, ntyp
      WRITE (unit=27,fmt=*) nt, ' ''', type_name(nt) , ' '' ', amass(nt)/2.0d0
    ENDDO
    DO na = 1, nat 
      WRITE (unit=27,fmt='(2i5,3f18.10)') na, ityp (na) , (tau (j, na) , j = 1, 3)
    ENDDO
    ! End header
    if (qirr .eq. 1) then
      lim1 = 1
      lim2 = nqs(qirr)
    else
      lim1 = sum(nqs(1:qirr-1)) + 1
      lim2 = sum(nqs(1:qirr))
    end if
    do qtot = lim1, lim2
      write (unit=27,fmt=9000) (q(i,qtot), i = 1, 3) 
      phi = 2.0d0 * phitot(qtot,:,:,:,:)
      if (qtot .eq. lim1) then
        do i = 1, nat
          do j = 1, nat
            do alpha = 1, 3
              do beta = 1, 3
                phi2(3*(i-1)+alpha,3*(j-1)+beta) = phi(alpha,beta,i,j)
              end do
            end do
          end do
        end do
      end if
      do i = 1, nat
        do j = 1, nat
          write (unit=27, fmt='(2i5)') i, j
          do alpha = 1, 3
            write (unit=27, fmt='(3(2f12.8,2x))') (phi(alpha,beta,i,j), beta=1,3)
          end do
        end do
      end do 
    end do
    ! Print effective charges in case there are     
    if (gamma_point .and. lrigid) then
      call write_epsilon_and_zeu (zeu, epsil, nat, 27)
    end if
    ! Diagonalize and print frequencies and polarization vectors
    call dyndia (q(:,lim2), 3*nat, nat, ntyp, ityp, amass/1822.8884482774586D0, &
                 27, phi2, w2)
    close (unit=27)
  end do

  deallocate(phi, phi2, w2)

 9000 format(/,5x,'Dynamical  Matrix in cartesian axes', &
         &       //,5x,'q = ( ',3f14.9,' ) ',/)

end subroutine write_dyn 

! This subroutine writes the header needed by import3py 
! 

subroutine write_mat2R (ntyp, ityp, amass, &
                            ibrav, celldm, tau, &
                            type_name, at, lrigid, epsil, zeu)

  implicit none

  double precision, dimension(:,:), intent(in) :: at
  integer, intent(in) :: ntyp
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:), intent(in) :: amass
  integer, intent(in) :: ibrav
  double precision, dimension(6), intent(in) :: celldm
  double precision, dimension(:,:), intent(in) :: tau 
  character (len=3), dimension(:), intent(in) :: type_name
  logical, intent(in) :: lrigid
  double precision, dimension(3,3), intent(in) :: epsil
  double precision, dimension(:,:,:), intent(in) :: zeu 

  integer :: nq, nqirr, nat, natsc, na, nt
  integer :: i, j, alpha, beta, iq, lim1, lim2, qtot, qirr, sigma 
  double complex, dimension(:,:,:,:), allocatable :: phi
  double complex, dimension(:,:), allocatable :: phi2
  character (len=200) :: fildyn
  double precision, dimension(:), allocatable :: w2 
  logical :: gamma_point

  character(len=6), EXTERNAL :: int_to_char
  CHARACTER(len=256),EXTERNAL :: trimcheck

  nat   = size(tau(1,:))

  fildyn = 'mat2R'

  open (unit=27,file=fildyn)
  ! Write header
  WRITE (unit=27,fmt='(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  WRITE (27,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  DO nt = 1, ntyp
    WRITE (unit=27,fmt=*) nt, ' ''', type_name(nt) , ' '' ', amass(nt)/2.0d0
  ENDDO
  DO na = 1, nat 
    WRITE (unit=27,fmt='(2i5,3f18.10)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  ENDDO
  write (unit=27,fmt=*) lrigid
  ! Print effective charges in case there are     
  if (lrigid) then
    call write_epsilon_and_zeu (zeu, epsil, nat, 27)
  end if
  close (unit=27)

end subroutine write_mat2R

! This subroutine takes as an input the dynamical matrix at a q
! point of the irreducible mesh and all the symmetry variables 
! needed. As an output the frequencies and the polarization vectors
! are given for the star of this q

subroutine phiqtophistar (nat, phi, ityp, amass, xq, at, bg, nsym, s, invs, &
                          nqs, sxq, imq, irt, rtau, isq, nqstot, phistar)
!                          nqs, sxq, imq, irt, rtau, isq, w, e)
  
  implicit none

  double complex, dimension(:,:,:,:), intent(in) :: phi
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(:), intent(in) :: amass
  double precision, dimension(3), intent(in) :: xq
  double precision, dimension(3,3), intent(in) :: at, bg
  integer, intent(in) :: nat, nsym, nqs, imq
  integer, dimension(3,3,48), intent(in) :: s
  integer, dimension(48), intent(in) :: invs, isq
  double precision, dimension(3,48), intent(in) :: sxq
  integer, dimension(:,:), intent(in) :: irt
  double precision, dimension(:,:,:), intent(in) :: rtau
!  double precision, dimension(nqs,3*nat), intent(out) :: w
!  double complex, dimension(nqs,3*nat,nat,3), intent(out) :: e 
  integer, intent(in) :: nqstot
  double complex, dimension(nqstot,3,3,nat,nat), intent(out) :: phistar

  double complex, dimension(:,:), allocatable :: d2 !, e_aux
  double complex, dimension(:,:,:,:,:), allocatable :: dynqstar
  integer :: i, j, na, nb, icar, jcar, nu 
!  double precision :: rydcm1
 
!  rydcm1 = 13.6058d0*8065.5d0

  ! repack phi to 3*nat,3*nat so that it can be repacked and then rerepacked again in q2qstar_ph

  ALLOCATE(d2(3*nat, 3*nat))
!  allocate(e_aux(3*nat,3*nat))
  allocate(dynqstar(nqstot,3,3,nat,nat))

  DO i = 1, 3 * nat
    na = (i - 1) / 3 + 1
    icar = i - 3 * (na - 1)
    DO j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcar = j - 3 * (nb - 1)
        d2 (i, j) = phi(icar, jcar, na, nb)
    ENDDO
  ENDDO
  !
print *, ' phiqtophistar nqs ', nqs
print *, ' phiqtophistar nqstot ', nqstot
  CALL q2qstar_out (d2, at, bg, nat, nsym, s, invs, irt, rtau, &
                    nqs, sxq, isq, imq, 1, nqstot, dynqstar)

  phistar = dynqstar

  deallocate(d2, dynqstar)

end subroutine phiqtophistar

! This subroutine diagonalizes the dynamical matrix at a given point
! and gives the polarization vectors and frequencies as output

subroutine dyndiag (nat, dyn, ityp, amass, w, e)

  implicit none

  integer, intent(in) :: nat
  double complex, dimension(3,3,nat,nat), intent(in) :: dyn
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(10), intent(in) :: amass
  double precision, dimension(3*nat), intent(out) :: w
  double complex, dimension(3*nat,nat,3), intent(out) :: e 

  integer :: i, j, icar, jcar, nu
  double precision :: rydcm1
  double complex, dimension(:,:), allocatable :: d2 , e_aux
 
  allocate(d2(3*nat,3*nat))
  allocate(e_aux(3*nat,3*nat))

  rydcm1 = 13.6058d0*8065.5d0

  do i = 1, nat
    do j = 1, nat
      do icar = 1, 3
        do jcar = 1, 3
          d2(3*(i-1)+icar,3*(j-1)+jcar) = dyn(icar,jcar,i,j) / &
                  dsqrt( amass(ityp(i)) * amass(ityp(j)))
        end do
      end do
    end do
  end do
  call cdiagh (3*nat,d2,3*nat,w,e_aux)
  do i = 1, 3*nat
    if (w(i) .lt. 0.0d0) then
      w(i) = - dsqrt( - w(i))
    else
      w(i) = dsqrt( w(i))
    end if
  end do
  do nu = 1, 3 * nat
    do i = 1, nat
      do icar = 1, 3
        e(nu,i,icar) = e_aux(3*(i-1)+icar,nu)
      end do   
    end do
  end do

  deallocate(d2,e_aux)

end subroutine dyndiag

! This subroutine diagonalizes the dynamical matrix in real space
! as it takes as an input a symmetric real matrix. Note that how
! the polarization vectors are defined is modified from the previous
! routine.              

subroutine dyndiag_real (nat, dyn, ityp, amass, w, e)

  implicit none

  integer, intent(in) :: nat
  double precision, dimension(3,3,nat,nat), intent(in) :: dyn
  integer, dimension(:), intent(in) :: ityp
  double precision, dimension(10), intent(in) :: amass
  double precision, dimension(3*nat), intent(out) :: w
  double precision, dimension(nat,3*nat,3), intent(out) :: e 

  integer :: i, j, icar, jcar, nu
  double precision :: rydcm1
  double precision, dimension(:,:), allocatable :: d2 , e_aux
 
  allocate(d2(3*nat,3*nat)) 
  allocate(e_aux(3*nat,3*nat))
 
  rydcm1 = 13.6058d0*8065.5d0

  do i = 1, nat
    do j = 1, nat
      do icar = 1, 3
        do jcar = 1, 3
          d2(3*(i-1)+icar,3*(j-1)+jcar) = dyn(icar,jcar,i,j) / &
                  dsqrt( amass(ityp(i)) * amass(ityp(j)))
        end do
      end do
    end do
  end do
  call cdiagh3 (3*nat,d2,3*nat,w,e_aux)
  do i = 1, 3*nat
    if (w(i) .lt. 0.0d0) then
      w(i) = - dsqrt( - w(i))
    else
      w(i) = dsqrt( w(i))
    end if
  end do
  do nu = 1, 3 * nat
    do i = 1, nat
      do icar = 1, 3
        e(i,nu,icar) = e_aux(3*(i-1)+icar,nu)
      end do   
    end do
  end do

  deallocate(d2,e_aux)

end subroutine dyndiag_real

! This subroutine assignes the translation modes

subroutine assign_trans (er,ityp_sc,amass,transmode)

  implicit none

  double precision, dimension(:,:,:), intent(in) :: er
  integer, dimension(:), intent(in) :: ityp_sc
  double precision, dimension(:), intent(in) :: amass
  logical, dimension(:), intent(out) :: transmode

  integer :: nat, nmodes
  integer :: i, alpha, mu, tot, check
  double precision :: tolerance, val
  logical :: notrans
  double precision, dimension(3) :: vect1, vect2

!  tolerance = 0.0001d0
  tolerance = 0.1d0

  nat = size(er(:,1,1))
  nmodes = size(transmode)

  transmode = .false.

  check = 0

  do mu = 1, nmodes
    notrans = .false.
    tot = 0
    do i = 1, nat
      if (notrans) exit
      vect1(:) = er(1,mu,:)/dsqrt(amass(ityp_sc(1))) - er(i,mu,:)/dsqrt(amass(ityp_sc(i)))
      vect2(:) = er(1,mu,:)/dsqrt(amass(ityp_sc(1))) 
      val = dot_product(vect1,vect1) / dot_product(vect2,vect2)
!      print *, ' v1 v1',  dot_product(vect1,vect1)
!      print *, ' v2 v2',  dot_product(vect2,vect2)
!      print *, ' val  ', val
      if ( val .lt. tolerance ) then
        tot = tot + 1
      else
        notrans = .true.
        exit   
      end if    
!      print *, ' tot ', tot
!      end do
!      do alpha = 1, 3
!        val = abs( (er(1,mu,alpha)/dsqrt(amass(ityp_sc(1))) - er(i,mu,alpha)/dsqrt(amass(ityp_sc(i)))) &
!                   /  er(1,mu,alpha)/dsqrt(amass(ityp_sc(1))) )
!        if ( val .lt. tolerance ) then
!          tot = tot + 1
!        else
!          notrans = .true.
!          exit   
!        end if    
!      end do
    end do
    if (tot .eq. nat) then
      transmode(mu) = .true.
      check = check + 1  
    end if
    if (check .eq. 3) exit
  end do

end subroutine assign_trans 


!
! GET HARMONIC ENERGY
! ===================
!
! This method returns the harmonic energy of the SCHA hamiltonian. It can be used to compute the SCHA pressure.
!
! This is done with the formula:
! .. math ::
!       < \mathcal H > = \sum_m \frac{\hbar \omega_m}{2 \tanh(\frac{\beta\hbar\omega_m}{2})}
!
!
! Parameters
! ----------
!
! - w : type (array of double precision)
!      The frequencies of the harmonic hoscillators.
! - transmode : type (array of logical)
!      A logical array of the same shape as w, it contains whether the corresponding mode is or not a translation.
! - T : type double precision
!      Temperature of the system. 
!
! Results
! ------
! - energy : double precision
!      The average energy of the system.
!
function get_harmonic_energy(w, transmode, T) result (energy)
  double precision, dimension(:), intent(in) :: w
  logical, dimension(:), intent(in) :: transmode
  double precision, intent(in) :: T
  double precision :: energy

  integer :: i, n
  double precision :: beta

  ! Convert boltzmann constant
  beta =  315774.65221921849D0 / T

  energy = 0.0d0
  n = size(w)
  do i = 1, n
     ! Skip the translational mode
     if (transmode(i)) continue

     if (T .eq. 0) then
        energy = energy + .5d0 * w(i)
     else
        energy = energy + w(i) / (2.0d0 * dtanh(beta * w(i) / 2.0d0))
     end if
  end do
end function get_harmonic_energy


!
! GET HARMONIC PRESSURE
! =====================
!
! According to the equipartition theorem the harmonic Pressure can be defined as
! ..math::
!     P = - \frac{<\matchal H>}{3\Omega}
!
! Where :math:`\Omega` is the volume of the unit cell. See Docs for details.
! The result is in untis of Ha / bohr^3
!
! Parameters
! ----------
!  - opt : type(options)
!     This type contains all the information about the structure. Used to compute the unit cell volume
!  - w : array of double precision
!     The harmonic frequencies
!  - transmode: array of logical
!     An array containing the info about which mode is a translation. It is of the same dimension as w.
!  - T : double precision
!     Temperature of the system.
!
! Results
! -------
!  - pressure : double precision
!     The harmonic pressure of the system.
!     (is the work made to expand the unit cell as the forces where armonic divided by the volume expansion)
!
function get_harmonic_pressure(opt, w, transmode, T) result (pressure)
  type(options) :: opt
  double precision, dimension(:), intent(in) :: w
  logical, dimension(:), intent(in) :: transmode
  double precision, intent(in) :: T
  double precision :: energy

  double precision :: Omega

  ! Compute the unit cell volume [bohr^3]
  Omega = opt%at(1,1)*(opt%at(2,2)*opt%at(3,3) - opt%at(2,3)*opt%at(3,2)) - &
       opt%at(2,1) *(opt%at(1,2)*opt%at(3,3) - opt%at(3,2)*opt%at(1,3)) + &
       opt%at(3,1) * (opt%at(1,2)*opt%at(2,3) - opt%at(2,2)*opt%at(1,3))
  
  ! Rescale from alat^3 to bohr^3
  Omega = Omega * opt%celldm(1) * opt%celldm(1) * opt%celldm(1)
  
  pressure = - get_harmonic_energy(w, transmode, T) / (3 * Omega)
end function get_harmonic_pressure

!
! GET HARMONIC STRESS TENSOR
! =====================
!
! This subroutine generalizes the previous one, but it computes the whole harmonic stress tensor.
!
! Parameters
! ----------
!  - opt : type(options)
!     This type contains all the information about the structure. Used to compute the unit cell volume
!  - w : array of double precision
!     The harmonic frequencies
!  - er_sc : double precision, dimension(nmodes, n_atoms_sc, 3)
!     The polarization vectors.
!  - transmode: array of logical
!     An array containing the info about which mode is a translation. It is of the same dimension as w.
!  - T : double precision
!     Temperature of the system.
!
! Results
! -------
!  - harmonic_stress : double precision, dimension(3,3)
!     The harmonic stress tensor of the system
!
subroutine get_harmonic_stress(opt, w, er_sc, transmode, T, harmonic_stress)
  type(options) :: opt
  double precision, dimension(:), intent(in) :: w
  double precision, dimension(:,:,:), intent(in) :: er_sc
  logical, dimension(:), intent(in) :: transmode
  double precision, intent(in) :: T
  double precision, dimension(3,3), intent(out) :: harmonic_stress
 
  double precision :: volume, beta
  integer i, j, mu

  ! Compute the volume [bohr^3]
  volume = opt%at_sc(1,1)*(opt%at_sc(2,2)*opt%at_sc(3,3) - opt%at_sc(2,3)*opt%at_sc(3,2)) - &
       opt%at_sc(2,1) *(opt%at_sc(1,2)*opt%at_sc(3,3) - opt%at_sc(3,2)*opt%at_sc(1,3)) + &
       opt%at_sc(3,1) * (opt%at_sc(1,2)*opt%at_sc(2,3) - opt%at_sc(2,2)*opt%at_sc(1,3))
  
  ! Rescale from alat^3 to bohr^3
  volume = volume * opt%celldm(1) * opt%celldm(1) * opt%celldm(1)
  
  ! Convert boltzmann constant
  beta =  315774.65221921849D0 / T
  
  ! Compute the stress tensor
  do i = 1, 3
     do j = i, 3
        harmonic_stress(i, j) = 0

        ! Sum over all the modes
        do mu = 1, 3 * opt%natoms_sc
           harmonic_stress(i, j) = harmonic_stress(i,j) - w(mu) / (2 * dtanh(beta * w(mu) / 2)) * &
                sum( er_sc(:, mu, i) * er_sc(:, mu, j))
        end do

        ! Get the pressure dividing by the volume
        harmonic_stress(i, j) = harmonic_stress(i, j) / volume
        
        ! Get the symmetric element
        if (i /= j) then
           harmonic_stress(j, i) = harmonic_stress(i, j)
        end if
     end do
  end do
end subroutine get_harmonic_stress

! This subroutine calculates a Fourier transform of a vector
! from real space to reciprocal space for a particular q point. There is a logical flag that
! allows to choose if a positive exponent is used or not.

subroutine fourier_r2q_vector ( q, latvec, tau_sc_latvec, plus, vectr, vectq )

  implicit none

  double precision, dimension(3), intent(in) :: q
  double precision, dimension(:,:), intent(in) :: latvec                                 
  integer, dimension(:,:), intent(in) :: tau_sc_latvec  
  logical, intent(in) :: plus 
  double precision, dimension(:,:), intent(in) :: vectr
  double complex, dimension(:,:), intent(out) :: vectq

  integer :: nat, nat_sc, nr
  double complex :: im, one, complex_number
  double precision :: twopi

  integer :: s, t

  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0

  ! Get number of atoms in unit cell and supercell, and the
  ! number of translation vectors

  nat    = size(vectq(:,1))
  nat_sc = size(vectr(:,1))

  nr = nat_sc / nat

  ! Perform Fourier transform

  vectq = (0.0d0,0.0d0)

  do s = 1, nat
    do t = 1, nr
      if (plus) then
        complex_number = exp (im * twopi * dot_product(q,latvec(t,:)))
      else if (.not. plus) then
        complex_number = exp (-im * twopi * dot_product(q,latvec(t,:)))
      else
        print *, ''
        print *, '  Problem in Fourier transform. Stopping...'
        print *, ''
        stop
      end if
      vectq(s,:) = vectq(s,:) + complex_number * vectr(tau_sc_latvec(s,t),:) 
    end do
  end do 

  vectq = vectq / sqrt(dble(nr))

end subroutine fourier_r2q_vector 

! This subroutine gets the phonon frequencies and polarization
! vectors in q space from the real spcae force constants matrix

subroutine weq_from_fc ( phitot_sc, q, tau, tau_sc, &
                         ityp, itau, amass, wq, eq )

  implicit none

  double precision, dimension(:,:,:,:), intent(in) :: phitot_sc
  double precision, dimension(:,:), intent(in) :: q, tau, tau_sc
  integer, dimension(:), intent(in) :: ityp, itau
  double precision, dimension(:), intent(in) :: amass
  double precision, dimension(:,:), intent(out) :: wq
  double complex, dimension(:,:,:,:), intent(out) :: eq
  
  double complex, dimension(:,:,:,:), allocatable :: dyn
  integer :: nq, nat, natsc
  integer :: i, j, k, alpha, beta, qtot, R
  integer :: ka
  double precision, dimension(:,:), allocatable :: latvec
  double precision, dimension(3) :: vecaux
  double complex :: im, one, complex_number
  double precision :: twopi, prec
  double complex, dimension(:,:,:), allocatable :: eaux
  double precision, dimension(:), allocatable :: waux


  one    = (1.0d0,0.0d0)
  im     = (0.0d0,1.0d0)
  twopi  = 6.283185307179586d0

  prec = 1.0d-6

  nq    = size(wq(:,1))
  nat   = size(ityp)
  natsc = nq * nat

  allocate(latvec(nq,3))
  allocate(dyn(3,3,nat,nat))

  allocate(eaux(3*nat,nat,3))
  allocate(waux(3*nat))

  ! Prepare list of lattice vectors

  ka = 0

  do i = 1, natsc
    if (itau(i) .ne. 1) cycle
    ka = ka + 1
    latvec(ka,:) = tau_sc(:,i) - tau(:,1)
  end do

  do qtot = 1, nq
    dyn = (0.0d0,0.0d0)
    do i = 1, nat
      do j = 1, nat
        do alpha = 1, 3
          do beta = 1, 3
            complex_number = (0.0d0,0.0d0)
            do R = 1, nq
              ! Check what the atom in the supercell is
              do k = 1, natsc
                vecaux = tau(:,j) + latvec(R,:) - tau_sc(:,k)
                if ( sqrt(dot_product(vecaux,vecaux)) .lt. prec ) then
                  complex_number = complex_number +  &
                                  exp( - im * twopi * dot_product(q(:,qtot),latvec(R,:))) * &
                                  phitot_sc(alpha,beta,i,k)
                end if
              end do
            end do
            dyn(alpha,beta,i,j) = complex_number
          end do
        end do
      end do
    end do
    call dyndiag (nat, dyn, ityp, amass, waux, eaux)
    wq(qtot,:)     = waux(:)
    do i = 1, nat
      do j = 1, 3*nat
        eq(qtot,i,j,:) = eaux(j,i,:)
      end do
    end do
  end do

  ! Deallocate stuff

  deallocate(latvec,dyn,eaux,waux)

end subroutine weq_from_fc

end module polarization
