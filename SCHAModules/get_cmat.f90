
! This subroutine gets the the C matrix needed for the odd term correction

subroutine get_cmat ( a, wr, er, transmode, amass, ityp_sc, T, v3_log, cmat, &
        n_mode, nat_sc, ntyp )

    implicit none
  
    double precision, dimension(n_mode), intent(in) :: a, wr
    double precision, dimension(nat_sc,n_mode,3), intent(in) :: er
    logical, dimension(n_mode), intent(in) :: transmode
    double precision, dimension(ntyp), intent(in) :: amass
    integer, dimension(nat_sc), intent(in) :: ityp_sc
    double precision, intent(in) :: T
    logical :: v3_log
    double precision, dimension(n_mode*n_mode, n_mode*n_mode), intent(out) :: cmat
  
    integer :: nat_sc, n_mode, n_mode2, ntyp
    double precision, dimension(:,:), allocatable :: e
    double precision, dimension(:,:), allocatable :: mat_w, mat_e, mat_et
  
    integer :: i, mu, nu, alpha
    integer :: ka, ja
    integer :: x, y, z, w
  
    real :: t1, t2

    logical, parameter :: debug = .true.
  
    ! Get integers
    if (debug) then
        print *, "=== DEBUG GET_CMAT ==="
        print *, "N_MODE:", n_mode
        print *, "NAT_SC:", nat_sc 
        print *, "NTYP:", ntyp
        call flush()
    endif
  
    !nat_sc = size(er(:,1,1))
    !n_mode = 3*nat_sc
    n_mode2 = n_mode*n_mode
  
    ! Allocate stuff

    print *, "CMAT: ALLOCATION"
    call flush()
  
  
    allocate(e(n_mode,n_mode))
    allocate (mat_w(n_mode,n_mode))
    allocate (mat_e(n_mode*n_mode,n_mode*n_mode))
    allocate (mat_et(n_mode*n_mode,n_mode*n_mode))

    print *, "CMAT: BEFORE EMAT"
    call flush()
  
    ! Define the polarization vectors as a 3n x 3n matrix (n = nat_sc).
    ! The square root of the mass is also included in the new matrix.
  
    call get_emat ( er, a, amass, ityp_sc, v3_log, transmode, e, n_mode, nat_sc, ntyp)

    print *, "CMAT: after get emat"
    call flush()
  
    ! Calculate the matrix based on the frequencies and lengths
  
    call get_g (a, wr, transmode, T, mat_w)

    print *, "CMAT: after get g"
    call flush()
  
    ! We get the e and et matrices
  
    call cpu_time(t1)
  
    ka = 0 ! counter for atom/cartesian index
  
    do x = 1, n_mode
      do y = 1, n_mode
        ka = ka + 1
        ja = 0 ! counter for mode index
        do mu = 1, n_mode 
          do nu = 1, n_mode 
            ja = ja + 1
            mat_e(ka,ja)  = e(nu,x) * e(mu,y)
            mat_et(ja,ka) = mat_e(ka,ja) * mat_w(mu,nu) * 0.5d0
          end do
        end do
      end do
    end do
  
    call cpu_time(t2)
   
    print *, ' Elapsed time calculating mat_e and mat_et = ', t2-t1
  
    ! Perfom matrix product to get c matrix
  
    call cpu_time(t1)
  
    call dgemm('N','N',n_mode2,n_mode2,n_mode2,1.0d0,mat_e,n_mode2,&
               mat_et,n_mode2,0.0d0,cmat,n_mode2)
  
    call cpu_time(t2)
   
    print *, ' Elapsed time calculating cmat = ', t2-t1
    
  !  ! Print c matrix
  !
  !  ka = 0
  !
  !  do x = 1, n_mode
  !    do y = 1, n_mode
  !      ka = ka + 1
  !      ja = 0
  !      do z = 1, n_mode
  !        do w = 1, n_mode
  !          ja = ja + 1
  !          print '(4i4,f24.10)', x, y, z, w,  cmat(ka,ja)  
  !        end do
  !      end do
  !    end do
  !  end do
  
    ! Deallocate objects
  
    deallocate(e)
    deallocate (mat_w)
    deallocate (mat_e)
    deallocate (mat_et)
  
  end subroutine get_cmat