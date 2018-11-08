module anharmonic
  use stochastic
  use thermodynamic
  
contains

  ! This subroutines computes df/da the first part of the
  ! Gradient for the SCHA code.
  !
  ! w <= Frequencies [Ha]
  ! w_harmonic <= Used only if stat-method = "stat_harmon"
  ! T <= Temperature [K]
  ! e <= Pol vectors (i,j) i is the coordinate, j the mode
  ! forces <= [Ha / bohr] (a,b,c) a is the random conf, b the atom, c the coord.
  ! q <= the q vector [bohr] (a,b) a is the config, b is the mode
  ! type_atoms <= same atoms have the same integer
  ! mass <= Mass (in Ha)
  ! stat_method <= string use "stat_normal"
  ! df_da => Result [Ha / bohr]
  subroutine get_df_da_nonav (w,w_harmonic,T,e,forces,q,mass,&
       stat_method,df_da, l, n, m)
 
    implicit none
    
    double precision, dimension(n), intent(in) :: w
    double precision, dimension(n), intent(in) :: w_harmonic
    double precision, intent(in) :: T
    double precision, dimension(3*m,n), intent(in) :: e
    double precision, dimension(l,m,3), intent(in) :: forces
    double precision, dimension(l,n), intent(in) :: q
    !double precision, dimension(n), intent(in) :: a
    double precision, dimension(m), intent(in) :: mass
    character (len=11), intent(in) :: stat_method
    double precision, dimension(l,n), intent(out) :: df_da

    integer, intent(in) ::  l, n, m
    
    integer :: i1, i2, i3, j, natoms, nrandom, nmodes, i3n
    double precision, dimension(:,:), allocatable :: f
    double precision, dimension(:,:), allocatable :: mat_a, mat_b, mat_c
    double precision, dimension(:,:), allocatable :: e_diag
    double precision, dimension(n) :: a, da_dw
    double precision, dimension(3*m) :: u_sqrtm 

    
    natoms  = m
    nrandom = l
    nmodes  = n
    
    allocate(mat_a(nmodes, 3 * natoms))
    allocate(mat_b(3 * natoms, nmodes))
    allocate(mat_c(nmodes,nmodes))
    allocate(f(nmodes,nrandom))
    allocate(e_diag(3*natoms, nmodes))

    ! Get a and da_dw
    call w_to_a(w, T, a, n)
    call w_to_da(w, T, da_dw, n)
    
    
!    ! Print all the info about the input parameters [DEBUG]
!    print *, ""
!    print *, "======== dF / dA INPUT VARIABLES ========="
!    print *, "w = ", w
!    print *, "w_harm = ", w_harmonic
!    print *, "a = ", a
!    print *, "T = ", T
!    print *, "Type atoms:", type_atoms
!    print *, "Mass = ", mass
!    print *, "Stat = ", stat_method
!    print *, "Pol vectors:"
!    do i1 = 1, nmodes
!        print "(I4, A3, 1000D16.8)", i1, ")", e(:, i1)
!    end do
!    print *, "Forces:"
!    do i1 = 1, nrandom
!        do i2 = 1, natoms
!            print "(A10, I4, A10, I4, A3, 3D16.8)", "Conf", i1, "Atom", i2, "=", forces(i1, i2, :)
!        end do
!    end do
!    print *, "q vectors:"
!    do i1 = 1, nrandom
!        print "(I4, A3, 1000D16.8)", i1, ")", q(i1, :)
!    end do
!    
    
    ! Create the auxiliar polarization vectors dividing by the
    ! the polarization vectors


    i3n = 0
    do i1 = 1, natoms
       do i3 = 1, 3
          i3n = i3n +1
          e_diag(i3n,:) = e(i3n,:) / dsqrt(mass(i1))
       end do
    end do

!    ! Print the quantity that should not change
!    print *, "Pol dot q:"
!    do i1 = 1, nrandom
!        u_sqrtm = 0.0d0
!        do i2 = 1, nmodes
!            u_sqrtm = u_sqrtm + e_diag(:, i2)* q(i1, i2)
!        end do
!        print "(A10, I4, A3, 1000D16.8)", "Conf", i1, "=", u_sqrtm
!    end do
!    print *, "====== dF / dA END INPUT VARIABLES ======="
!    print *, ""
    ! Create mat_a matrix

    do i1 = 1, nmodes
       ! if (.not. mu_relevant(i1)) then
       !    mat_a(i1,:) = 0.0d0 
       !    cycle
       ! end if
       i3n = 0
       do i2 = 1, natoms
          do i3 = 1, 3
             i3n = i3n + 1
             mat_a(i1,i3n) = e_diag(i3n, i1) / a(i1)
          end do
       end do
    end do

    ! Loop on the nmodes to build the function that we will
    ! use for the statistics.

    do j = 1, nrandom

       ! Build the mat_b matrix

       do i1 = 1, nmodes
          i3n = 0
          do i2 = 1, natoms
             do i3 = 1, 3
                i3n = i3n + 1
                mat_b(i3n,i1) = -forces(j,i2,i3)*q(j,i1)
             end do
          end do
       end do

       ! Call lapack subroutine to calculate matrix product

       call dgemm('N','N',nmodes,nmodes,3*natoms,1.0d0,mat_a,nmodes,&
            mat_b,3*natoms,0.0d0,mat_c,nmodes)
            
       ! Print the new matrix
!       
!       print *, ""
!       print *, "MAT A  Conf", j
!       do i1 = 1, 3*natoms
!           print *, mat_a(:, i1)
!       end do 
!       print *, "MAT B  Conf", j
!       do i1 = 1, nmodes
!           print *, mat_b(:, i1)
!       end do 
!       print *, "MAT C  Conf", j
!       do i1 = 1, nmodes
!           print *, mat_c(:, i1)
!       end do 

       ! Create the gradient according to the different possible statistics

       do i1 = 1, nmodes
          f(i1,j) = mat_c(i1,i1)
       end do

       ! Include in the gradients the analytical part

       do i1 = 1, nmodes
          if (stat_method .eq. 'stat_normal') then
             df_da(j,i1) = f(i1,j) + dW_f0_u0(w(i1),T) / da_dw(i1)
          else if (stat_method .eq. 'stat_harmon') then
             df_da(j,i1) = f(i1,j) + dW_f0_u0(w(i1),T) / da_dw(i1) &
                  + w_harmonic(i1)**2.0d0 * a(i1)
          else if (stat_method .eq. 'stat_schapp') then
             df_da(j,i1) = f(i1,j) !+ dW_f0_u0(w(i1),T) / w_to_da(w(i1),1.0d0,T) &
             !+ w(i1)**2.0d0 * a(i1)
          else
             stop "ERROR, STAT_METHOD NOT VALID"
          end if
       end do
    end do

    deallocate(mat_a)
    deallocate(mat_b)
    deallocate(mat_c)
    deallocate(f)
    deallocate(e_diag)

  end subroutine get_df_da_nonav
  
  
  ! The following subroutine computes the last derivatives
  ! That are 
  
  subroutine get_da_dcr_and_de_dcr (wr, er, T, mass, x_i, y_i, da_dcr, de_dcr, n, m ) 

    implicit none

    double precision, dimension(n), intent(in) :: wr
    double precision, dimension(3*m,n), intent(in) :: er
    double precision, intent(in) :: T
    double precision, dimension(m), intent(in) :: mass
    integer, intent(in) :: x_i, y_i
    double precision, dimension(n), intent(out) :: da_dcr
    double precision, dimension(3*m, n), intent(out) :: de_dcr
    
    ! Optional integer
    integer, intent(in) :: n, m

    integer :: nmodes, natsc
    integer :: i, mu, nu, alpha, atm_x, atm_y
    double precision, dimension(:,:), allocatable :: g, gamma_mu_nu
    double precision, dimension(n) :: da_dw

    natsc  = m
    nmodes = n
    
    atm_x = (x_i - 1) / 3 + 1
    atm_y = (y_i -1 ) / 3 + 1

    allocate(g(nmodes,nmodes))
    allocate(gamma_mu_nu(nmodes, nmodes))
    
    ! Prepare gamma_mu_nu
    do mu = 1, nmodes
        do nu = 1, nmodes
            gamma_mu_nu(mu, nu) = er(x_i, mu) * er(y_i, nu) /&
                dsqrt(mass(atm_x) * mass(atm_y))
        end do
    end do

    ! Take the derivative of w_to_da
    call w_to_da(wr, T, da_dw, nmodes)

    do mu = 1, nmodes
        da_dcr(mu) = 0.5d0 * da_dw(mu) * gamma_mu_nu(mu,mu) / wr(mu)
    end do

    ! Prepare matrices for product

    do nu = 1, nmodes
       do mu = 1, nmodes
          if (abs((wr(mu)-wr(nu))/wr(nu)) .lt. 0.00001d0) then
             g(mu,nu) = 0.0d0
          else
             g(mu,nu) = gamma_mu_nu(mu,nu) / (wr(nu)**2 - wr(mu)**2) 
          end if
       end do
    end do

    ! Make matrix product

    call dgemm('N','N',nmodes,nmodes,nmodes,1.0d0,er, &
         nmodes,g,nmodes,0.0d0,de_dcr,nmodes)
    deallocate(g)

  end subroutine get_da_dcr_and_de_dcr
  
  
  ! This subroutine gets the derivative of the Free energy with respect
  ! to the minimum number of coefficients that decompose the dynamical 
  ! matrices. In this case the statistics are performed here. Thus,
  ! the error of the gradoient is also calculated in this subroutine.
  ! In this new version the value of dF/de is calculated in this 
  ! subroutine. This is done to solve memory problems. 
  
  subroutine get_df_dcoeff_av_new (df, da, forces, q, mass, dedc, rho, &
       minim_polvec, log_err, dfdc, delta_dfdc, r, s, t)

    implicit none

    integer, intent(in) :: r,s,t
    double precision, dimension(r,s), intent(in) :: df, q
    double precision, dimension(s), intent(in) :: da
    double precision, dimension(t), intent(in) ::  mass
    double precision, dimension(r,t,3), intent(in) :: forces
    double precision, dimension(3*t,s), intent(in) :: dedc
    double precision, dimension(r), intent(in) :: rho
    logical, intent(in) :: minim_polvec
    !  character (len=11), intent(in) :: stat_method
    character (len=10), intent(in) :: log_err
    double precision, intent(out) :: dfdc, delta_dfdc
    

    double precision, dimension(:), allocatable :: f
    integer :: nat, nmodes, nrandom
    double precision :: aux
    integer :: i, alpha, mu, n, i3n
    double precision :: sqrtm1mass(t)


    sqrtm1mass = 1.d0/DSQRT(mass)
    nmodes = s
    nat = t
    nrandom = r
    
    allocate(f(nrandom))

    ! Compute df/da * da/dcr
    f = 0.d0
    do mu = 1,nmodes
       do n = 1, nrandom
          f(n) = f(n)+df(n,mu)*da(mu)
       end do
    enddo


    ! Compute df/de * de/dcr
    if (minim_polvec) then
      !
      do mu = 1, nmodes
        i3n = 1
        do i = 1, nat
            do alpha = 1, 3
               aux =  sqrtm1mass(i) * dedc(i3n,mu)  
               i3n = i3n + 1
               !do n = 1, nrandom
               FORALL(n=1:nrandom) 
                  f(n) = f(n) - aux * forces(n,i,alpha) * q(n,mu)
               END FORALL
            end do 
        end do
      end do
      !
    endif

    ! Calculate the statistics
    call average_error_weight(f,rho,log_err,dfdc,delta_dfdc)

    deallocate(f)

  end subroutine get_df_dcoeff_av_new

end module anharmonic
