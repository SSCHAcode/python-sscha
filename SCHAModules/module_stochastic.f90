module stochastic

public :: random_gaussian
public :: check_gaussian 
public :: average_1D
public :: average_3D
public :: standard_deviation_1D
public :: get_gaussian_weight
public :: average_error_weight

contains

! This subroutine gives a random gaussian distributed
! number. The distribution is normalized so that
!
! (1 / 2 pi) * Exp [ -y^2 / 2 ]
!
! It has been more or less copied from the book 
! 'Numerical Recipes in C'

! function population_idx(pop_out, pop_in_)
!    implicit none 
!    integer,intent(in) :: pop_out
!    integer,optional :: pop_in_
!    character(len=32) :: population_idx
!    integer :: pop_in
!    character(len=6), EXTERNAL :: int_to_char

!    IF(present(pop_in_))THEN
!      pop_in = pop_in_
!    ELSE
!      pop_in = -1
!    ENDIF

!    IF(pop_in>0)THEN
!     population_idx="_population"//int_to_char(pop_in)
!    ELSE IF(pop_out>0)THEN
!     population_idx="_population"//int_to_char(pop_out)
!    ELSE
!     population_idx=''
!    ENDIF  
!    RETURN
! end function


subroutine random_gaussian(rand)
  double precision, intent(out) :: rand

  double precision :: x1, x2, v1, v2
  double precision :: rsq, fac 

  do

    call random_number(x1)
    call random_number(x2)

    v1 = 2.0d0 * x1 -1.0d0
    v2 = 2.0d0 * x2 -1.0d0

    rsq = v1**2 + v2**2

    if (rsq .le. 1.0d0 .and. rsq .gt. 0.0d0) exit

  end do

  fac = sqrt(-2.0d0 * log(rsq) / rsq)

  rand = v1 * fac

end subroutine random_gaussian

! subroutine random_sobol(method, n, m, rand)
!   USE constants,        ONLY : sqrt2
!   USE SOBOL,            ONLY : i8_sobol
!   implicit none
!   integer,intent(in)  :: n, m ! n configurations, m dimensions
!   real(8),intent(out) :: rand(n,m)
!   character(len=12)   :: method
!   !
!   integer ( kind = 8 ) :: seed, dim_num
!   real(8),allocatable    :: u1(:), u2(:)
!   integer :: i,j
!   !
!   LOGICAL,SAVE :: first = .true.
!   REAL(8),EXTERNAL :: dierfc2
!   !
!   ! just to be sure with integer types:
!   dim_num = m
!   allocate(u1(dim_num),u2(dim_num))
!   !
!   IF(first)THEN
!     if(method(7:7)=='c')then
!       read(method(8:len(method)),*,iostat=i) j
!       if(i/=0) j=0
!       seed = 2**(CEILING( j+ LOG(DBLE(dim_num))/LOG(2.d0)  ))
!     else if(method(7:7)=='p')then
!       read(method(8:len(method)),*,iostat=i) j
!       if(i/=0) j=0
!       seed = 2**j
!     else
!       read(method(7:len(method)),*,iostat=i) j
!       if(i/=0) j=0
!       seed = j
!     endif
!     print*, "Sobol method: "//method, seed
! !     seed = 0
!     first = .false.
!   ELSE
!     print*, "I did not consider the case where Sobol is initialized more than once"
!     stop 1
!   ENDIF
!   !
!   do i = 1,n
!     !
!     call i8_sobol ( dim_num, seed, u1(1:dim_num) )
!     IF(ANY(u1>=1.d0).or.ANY(u1<=-1.d0) ) THEN
!       print*, "WARNING Sobol: a point is on the border"
!     ENDIF
! !     WHERE(u1>=1.d0)  u1 = 1.d0-2**(-dim_num)
! !     WHERE(u1<=-1.d0) u1 = -1.d0+2**(-dim_num)
    
!     u1 = 2*u1-1.d0
! #if defined (__INTEL)
!     ! Next call is inverse Erf, only available in mkl.
!     CALL vderfinv(dim_num, u1, u2)
! #else 
!     ! Or use homebrew inverse erf
!     DO j = 1, dim_num
!       u2(j) = dierfc2(u1(j))
!     ENDDO
! #endif
!     !
!     ! Bad index ordering, but we need it like this:
!     rand(i,:) = -sqrt2 * u2
!   enddo
!   !
!   deallocate(u1,u2)
!   !
! end subroutine random_sobol
! !
! ! This subroutine checks whether the number is a given 
! ! array are formed according to a gaussian distribution.
! ! As an output it gives the file 'test_gaussian' 
! ! where the data is compared to a gaussian distribution.

subroutine check_gaussian(random)

  double precision, dimension(:), intent(in) :: random

  integer :: i1, kg, i2 
  double precision :: pi

  pi = 3.141592653589793238d0

  open (unit=1, file='test_gaussian')

  do i1 = 1, 100
    kg = 0
    do i2 = 1, size(random)
      if (random(i2) .ge. (-5 + (i1-1)*0.1) .and. & 
          random(i2) .lt. (-5 + i1*0.1) ) kg = kg + 1
    end do
    write (unit=1,fmt=*) &
           -5.05d0 + dble(i1)*0.1d0, dble(kg) / dble(size(random)/10), &
            exp(-0.5d0*(-5.05d0 + dble(i1)*0.1d0)**2) / &
            (sqrt(2.0d0 * pi))
  end do

  close (unit=1)

end subroutine check_gaussian

! This subroutine calculates the gaussian weight used for the  
! correlated sampling 

subroutine get_gaussian_weight(q,q_start,a,a_start,rho, n_config, n_modes)

  double precision, dimension(n_config, n_modes), intent(in) :: q_start, q
  double precision, dimension(n_modes), intent(in) :: a_start, a
  double precision, dimension(n_config), intent(out) :: rho
  integer, intent(in) :: n_config, n_modes

  integer :: i, j, N, nmodes

  N = n_config
  nmodes = n_modes

  do i = 1, N
    rho(i) = 1.0d0
    ! Skip the contribution of the acoustic modes
    do j = 1, nmodes
      rho(i) = rho(i) * (a_start(j) / a(j)) * dexp(-(q(i,j)**2)/(2.0d0 * a(j)**2)) / &
             dexp(-(q_start(i,j)**2)/(2.0d0 * a_start(j)**2))
    end do
  end do

end subroutine get_gaussian_weight

! This subroutine calculates the gaussian weight used for the  
! correlated sampling in case there are effective charges

subroutine get_gaussian_weight_eff_charge(q,q_start,a,a_start,mu_relevant,rho)

  double precision, dimension(:,:), intent(in) :: q_start, q
  double precision, dimension(:,:), intent(in) :: a_start, a
  logical, dimension(:), intent(in) :: mu_relevant 
  double precision, dimension(:), intent(out) :: rho

  integer :: i, j, N, nmodes

  N = size(q(:,1))
  nmodes = size(q(1,:))

  do i = 1, N
    rho(i) = 1.0d0
    ! Skip the contribution of the acoustic modes
    do j = 1, nmodes
      if (.not. mu_relevant(j)) cycle
      rho(i) = rho(i) * (a_start(i,j) / a(i,j)) * dexp(-(q(i,j)**2)/(2.0d0 * a(i,j)**2)) / &
             dexp(-(q_start(i,j)**2)/(2.0d0 * a_start(i,j)**2))
    end do
  end do

end subroutine get_gaussian_weight_eff_charge


! This subroutine calculates the average and the error of
! a correlated average of the following type:
!
!     \sum_i f(i) p(i) / \sum_i p(i)
!
! where f(i) is the function to average and p(i) is the 
! weight. The error is calculated using the formula for 
! the propagation of errors.
!

subroutine average_error_weight(f,rho,log_err,av_f,error_f)

  double precision, dimension(:), intent(in) :: f
  double precision, dimension(:), intent(in) :: rho
  character (len=10), intent(in) :: log_err
  double precision, intent(out) :: av_f
  double precision, intent(out) :: error_f

  integer :: nc
  double precision :: av_f1, av_rho, s_f, s_rho, s_f_rho

  double precision,allocatable :: rhof(:)
  logical :: correct_evenodd

  ! Get the number of configurations

  nc = size(f)
  allocate(rhof(nc))

  ! Calculate the auxiliar parameters

  ! profile 5 : avoid multiplying rho times f over and over again
  rhof = rho(:)*f(:)

  av_f1   = sum(rhof(:)) / dble(nc)
  av_rho  = sum(rho(:))  / dble(nc)
  
  s_f     = sum((rhof(:)-av_f1)**2) / dble(nc-1)
  s_rho   = sum((rho(:)-av_rho)**2) / dble(nc-1)
  s_f_rho = sum((rhof(:)-av_f1)*(rho(:)-av_rho)) / dble(nc-1)      

  ! Calculate the average and the error

  if (log_err .eq. 'err_yesrho') then
    av_f    = av_f1 / av_rho
    error_f = (1.0d0/dsqrt(dble(nc))) * abs(av_f1/av_rho)     *  &
              dsqrt( s_f/(av_f1**2) + s_rho/(av_rho**2) - &
                        2*s_f_rho/(av_rho * av_f1) )
  else if (log_err .eq. 'err_nonrho') then
    av_f    = av_f1
    error_f = (1.0d0/dsqrt(dble(nc))) * dsqrt( s_f )
 else
    stop "ERROR; log_err ill defined."
  end if 

  deallocate(rhof)

end subroutine average_error_weight

! This subroutine gives a random number distributed
! on a sphere of radious 1. 
! It has been more or less copied from 
! http://mathworld.wolfram.com/SpherePointPicking.html 

subroutine random_sphere(v_rand)

  double precision, dimension(3), intent(out) :: v_rand

  double precision :: v1, v2, u, theta
  double precision :: twopi

  twopi  = 6.2831853071795862d0 

  call random_number(v1)
  call random_number(v2)

  u     = 2.0 * v1 - 1.0d0
  theta = v2 * twopi
 
  v_rand(1) = dsqrt(1.0 - u**2) * cos(theta)
  v_rand(2) = dsqrt(1.0 - u**2) * sin(theta)
  v_rand(3) = u

end subroutine random_sphere

end module stochastic
