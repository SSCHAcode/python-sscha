
! This subroutine calculates the L mat needed to get the average of the 
! third order derivatives. It is formed by four polarization vectors 
! times the mass^1/2 divided by the normal length.

subroutine get_odd_straight ( a, wr, er, transmode, amass, ityp_sc, T, v3, phi_sc_odd, &
    n_mode, nat_sc, ntyp)

  implicit none

  double precision, dimension(n_mode), intent(in) :: a, wr
  double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
  logical, dimension(n_mode), intent(in) :: transmode
  double precision, dimension(ntyp), intent(in) :: amass
  integer, dimension(nat_sc), intent(in) :: ityp_sc
  double precision, intent(in) :: T
  double precision, dimension(n_mode,n_mode,n_mode), intent(in) :: v3
  double precision, dimension(n_mode, n_mode), intent(out) :: phi_sc_odd

  integer :: nat_sc, n_mode, nl, ns, ntyp
  double precision, dimension(:,:), allocatable :: l, g, phi_aux, v1, v2, v32, maux
  double precision :: lsum 
  double precision, dimension(:), allocatable :: laux1, lres1, veclong
  double precision, dimension(:), allocatable :: laux2, lres2
 
  integer :: mu, nu, alpha
  integer :: ka, ja
  integer :: i, j, x, y, z, w

  real :: t1, t2
  logical, parameter :: debug = .true.

  ! Get integers

  !nat_sc = size(er(:,1,1))
  !n_mode = 3*nat_sc

  ns = n_mode
  nl = n_mode*n_mode

  ! Allocate stuff
  if (debug) then
    print *, "=== DEBUG ODD STRAIGHT ==="
    print *, "N_MODE:", n_mode
    print *, "NTYP:", ntyp 
    print *, "NAT_SC:", nat_sc
    call flush()
  end if

  allocate(l(n_mode,n_mode))
  allocate(g(n_mode,n_mode))
  !allocate(phi_aux(n_mode,n_mode))
  allocate(v1(n_mode,n_mode*n_mode))
  allocate(v2(n_mode,n_mode*n_mode))
  allocate(v32(n_mode,n_mode*n_mode))
  allocate(laux1(n_mode))
  allocate(lres1(n_mode))
  allocate(laux2(n_mode))
  allocate(lres2(n_mode))
  allocate(maux(n_mode,n_mode))
  allocate(veclong(nl))

  ! Allocate stuff
  if (debug) then
    print *, "=== DEBUG ODD STRAIGHT ==="
    print *, "AFTER ALLOCATION"
    call flush()
  end if


  ! Define the polarization vectors as a 3n x 3n matrix (n = nat_sc).
  ! The square root of the mass is also included in the new matrix
  ! and also the normal length 

  call get_emat ( er, a, amass, ityp_sc, .true., transmode, l, n_mode, nat_sc, ntyp)

  ! Calculate  the matrix g that will enter in the final equation

  call get_g (a, wr, transmode, T, g, n_mode)

  ! Allocate stuff
  if (debug) then
    print *, "=== DEBUG ODD STRAIGHT ==="
    print *, "AFTER G"
    call flush()
  end if
 
  ! Write third order force constants as rank 2

  ka = 0
   
  do x = 1, n_mode
    do y = 1, n_mode
      ka = ka + 1
      v32(:,ka) = v3(:,x,y)
    end do
  end do

  ! Calculate the auxiliary matrices

  ja = 0

!  do mu = 1, n_mode
!    do nu = 1, n_mode
!      veclong = 0.0d0
!      ka = 0
!      do x = 1, n_mode    
!        do y = 1, n_mode
!          ka = ka + 1
!          veclong(ka) = l(mu,x)*l(nu,y) 
!        end do
!      end do    
!      call dgemv ('N',ns,nl,1.0d0,v32,ns,veclong,1,0.0d0,lres1,1) 
!      ja = ja + 1
!      v2(:,ja) = lres1(:) * g(mu,nu)
!      v1(:,ja) = lres1(:) * 0.5d0
!    end do
!  end do

  do i = 1, n_mode
    maux(:,:) = v3(i,:,:)
    ja = 0 
    do mu = 1, n_mode 
      laux1 = l(mu,:)
      do nu = 1, n_mode         
        laux2 = l(nu,:)
        call dgemv ('N',ns,ns,1.0d0,maux,ns,laux1,1,0.0d0,lres1,1)
        ja = ja + 1
!        call dgemv ('T',ns,1,1.0d0,lres1,ns,laux2,1,0.0d0,lsum,1)
        v1(i,ja) = dot_product(lres1,laux2) 
        v2(i,ja) = v1(i,ja) * g(mu,nu)
        v1(i,ja) = v1(i,ja) * 0.5d0
!        v2(i,ja) = lsum * g(mu,nu)
!        v1(i,ja) = lsum * 0.5d0
      end do
    end do
  end do

!  do i = 1, n_mode
!    ka = 0
!    do mu = 1, n_mode
!      do nu = 1, n_mode
!        ka = ka +1
!        v1(i,ka) = 0.0d0
!        v2(i,ka) = 0.0d0
!        do x = 1, n_mode
!          do y = 1, n_mode
!            v1(i,ka) = v1(i,ka) + v3(i,x,y)*l(nu,x)*l(mu,y)*0.5d0
!            v2(i,ka) = v2(i,ka) + v3(i,x,y)*l(nu,x)*l(mu,y)*g(mu,nu)
!          end do
!        end do
!      end do
!    end do
!  end do
!
!  ! Make matrix product to get the odd correction
!
!  ns = n_mode
!  nl = n_mode*n_mode

  call dgemm('N','T',ns,ns,nl,1.0d0,v1,ns,v2,ns,0.0d0,phi_sc_odd,ns)

  ! Write the odd correction with four indexes

  !call twotofour_real (phi_aux,phi_sc_odd)

  ! Deallocate stuff

  deallocate(l,g,v1,v2,v32)

end subroutine get_odd_straight 
