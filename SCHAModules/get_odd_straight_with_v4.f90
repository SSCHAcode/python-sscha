
! This subroutine calculates the L mat needed to get the average of the 
! third order derivatives. It is formed by four polarization vectors 
! times the mass^1/2 divided by the normal length.

subroutine get_odd_straight_with_v4 ( a, wr, er, transmode, amass, ityp_sc, T, v3, v4, phi_sc_odd, &
  n_mode, nat_sc, ntyp)

  implicit none

  double precision, dimension(n_mode), intent(in) :: a, wr
  double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
  logical, dimension(n_mode), intent(in) :: transmode
  double precision, dimension(ntyp), intent(in) :: amass
  integer, dimension(nat_sc), intent(in) :: ityp_sc
  double precision, intent(in) :: T
  double precision, dimension(n_mode,n_mode,n_mode), intent(in) :: v3
  double precision, dimension(n_mode,n_mode,n_mode, n_mode), intent(in) :: v4
  double precision, dimension(n_mode, n_mode), intent(out) :: phi_sc_odd


  integer :: nat_sc, n_mode, nl, ns, ntyp
  double precision, dimension(:,:), allocatable :: l, g, phi_aux, v1, v2, v32, iden
  double precision :: lsum 
  double precision, dimension(:), allocatable :: laux1, lres1, veclong
  double precision, dimension(:), allocatable :: laux2, lres2
 
  double precision, dimension(:,:), allocatable :: lamat, v42, maux
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv
  integer :: info

  double precision, dimension(:), allocatable :: vv
  double precision, dimension(:), allocatable :: ww
  double precision, dimension(:,:), allocatable :: zz
  double precision, dimension(:,:), allocatable :: cf

  integer :: mu, nu, alpha
  integer :: ka, ja
  integer :: i, j, x, y, z, w

  real :: t1, t2

  logical, parameter :: debug = .true.

  ! Get integers

  if (debug) then
    print *, "=== DEBUG ODD STRAIGHT ==="
    print *, "N_MODE:", n_mode
    print *, "NTYP:", ntyp 
    print *, "NAT_SC:", nat_sc
    call flush()
  end if

  !nat_sc = size(er(:,1,1))
  !n_mode = 3*nat_sc

  ns = n_mode
  nl = n_mode*n_mode

  ! Allocate stuff

  allocate(lamat(nl,nl))
  allocate(v42(nl,nl))
  allocate(maux(nl,nl))
  allocate(ipiv(nl))
  allocate(work(nl))
  allocate(v32(n_mode,n_mode*n_mode))
  allocate(iden(nl,nl))

  allocate(cf(nl,ns))

  allocate(vv(nl*(nl+1)/2))
  allocate(ww(nl))
  allocate(zz(nl,nl))

  ! Get lambda matrix 

  call get_cmat ( a, wr, er, transmode, amass, ityp_sc, T, .true., lamat ) 

  ! Write third and fourth order force constants as rank 2

  ka = 0
   
  do x = 1, n_mode
    do y = 1, n_mode
      ka = ka + 1
      v32(:,ka) = v3(:,x,y)
      ja = 0
      do w = 1, n_mode
        do z = 1, n_mode
          ja = ja + 1
          v42(ja,ka) = v4(w,z,x,y)    
        end do
      end do
    end do
  end do

  ! Prepare identity matrix

  iden = 0.0d0

  do x = 1, nl
    iden(x,x) = 1.0d0
  end do

  ! Calculate ** iden - v4 lamat ** matrix

  maux = iden

  call dgemm('N','N',nl,nl,nl,-1.0d0,v42,nl,lamat,nl,1.0d0,maux,nl) 

  ! Invert ** iden - lamat v4 **

  call dgetrf ( nl, nl, maux, nl, ipiv, info )                                                                                                                                                      
  call dgetri ( nl, maux, nl, ipiv, work, nl, info )   

  ! Take product between lamat and the inverted matrix

  call dgemm('N','N',nl,nl,nl,1.0d0,lamat,nl,maux,nl,0.0d0,v42,nl)

  ! Calculate final matrix products and assign the correction matrix

  ! Calculate cf = ( 1 - lamat*v4)^-1 lamat *  v3
  call dgemm('N','T',nl,ns,nl,1.0d0,v42,nl,&
             v32,ns,0.0d0,cf,nl)

  ! Now get:
  ! v3 * ( 1 - lamat*v4)^-1 lamat *  v3
  call dgemm('N','N',ns,ns,nl,1.0d0,v32,ns,&
             cf,nl,0.0d0,phi_sc_odd,ns)


  !call get_odd_from_cmat_fu2 (v42, v32, phi_sc_odd)

  ! Deallocate stuff

  deallocate(lamat,v32,v42,maux,ipiv,work, cf)

end subroutine get_odd_straight_with_v4 