

! This subroutine calculates the stochastic average of the fourth order
! force constants based on fu2 and lmat. The output force constants
! are given with three indexes that correspond to atom-cartesian coordinates.
!

subroutine get_v4 ( a, er, transmode, amass, ityp_sc, f, u, rho, log_err, v4, &
    nat_sc, n_mode, n_random, ntyp)
    use omp_lib
    use stochastic

    implicit none
  
    double precision, dimension(n_mode), intent(in) :: a
    double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
    logical, dimension(n_mode), intent(in) :: transmode
    double precision, dimension(ntyp), intent(in) :: amass
    integer, dimension(nat_sc), intent(in) :: ityp_sc
    double precision, dimension(n_random,nat_sc,3), intent(in) :: f, u
    double precision, dimension(n_random), intent(in) :: rho
    character (len=10), intent(in) :: log_err
    double precision, dimension(n_mode,n_mode,n_mode, n_mode), intent(out) :: v4
  
    integer :: nat_sc, n_mode, n_mode2, n_random, ntyp

    ! Here the variables for this method
  
    double precision, dimension(:,:), allocatable :: e
    double precision, dimension(:), allocatable :: fun
    double precision, dimension(:,:), allocatable :: ur, u2, f2 
    double precision :: av, av_err 
    
    integer :: i, j, mu, nu, alpha, beta
    integer :: ka, ja
    integer :: x, y, z, w, thread_num

    logical, parameter :: debug = .true.
  
    real :: t1, t2
  
    ! Get integers
  
    !nat_sc = size(er(:,1,1))
    !n_mode = 3*nat_sc
    n_mode2 = n_mode*n_mode
    !n_random = size(u(:,1,1))


    if (debug) then
        print *, "=== V3 DEBUG ==="
        print *, "N_MODE:", n_mode
        print *, "N_RANDOM:", n_random
        print *, "NAT_SC:", nat_sc 
  
        print *, ""
        print *, "A displacemets:", a(:)
        print *, "Translations:", transmode(:)
        print *, "Weights:", rho(:)
        print *, "AMASS:", amass(:)
        call flush()
    end if
    
  
    ! Allocate stuff
  
    allocate(e(n_mode,n_mode))
    allocate(ur(n_random,n_mode))
    allocate(u2(n_random,n_mode))
    allocate(f2(n_random,n_mode))
    allocate(fun(n_random))
  
    ! Calculate e polarization vectors with lengths and masses
  
    call get_emat (er, a, amass, ityp_sc, .false., transmode, e) 
  
    ! Get displacements in a rank two tensor
  
    ka = 0
  
    do i = 1, nat_sc
      do alpha = 1, 3
        ka = ka + 1
        u2(:,ka) = u(:,i,alpha)
        f2(:,ka) = f(:,i,alpha)
      end do
    end do  
  
    ! Rotate displacementes
  
    do x = 1, n_mode
      ur(:,x) = 0.0d0
      do mu = 1, n_mode
        do y = 1, n_mode
          ur(:,x) = ur(:,x) + u2(:,y)*e(mu,x)*e(mu,y)
        end do
      end do
    end do
  
    ! Calculate the third order coefficients
  
    !thread_num = omp_get_max_threads ( )
  
    !$omp parallel SHARED ( v4, ur, f2, rho, n_mode, log_err ) PRIVATE ( x, y, z, w, fun, av, av_err )
    !$omp do COLLAPSE(4)
    do x = 1, n_mode
      do y = 1, n_mode
        do z = 1, n_mode 
          do w = 1, n_mode 
            fun(:) = - ur(:,x) * ur(:,y) * ur(:,z) * f2(:,w)
            call average_error_weight(fun,rho,log_err,av,av_err)
            v4(x,y,z,w) = av 
          end do
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  
    deallocate(e,ur,u2,f2,fun)
  
  end subroutine get_v4
  