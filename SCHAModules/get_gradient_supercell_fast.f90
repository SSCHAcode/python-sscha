! By LORENZO MONACELLI
!
!
! Get the gradient of the system in the supercell using the Raffaello Bianco's equation
! EQ. A18b of PRB paper (2017) Second order structural....
!
! The equation is modified to exploit the noise reduction as:
!
!
! dF/dPhi = Lambda  Upsilon < u (f - f_scha) > 
!
! The Lambda matrix is the Hessian of the minimization Phi variable (Monacelli et al 2018 PRB)
! Upsilon matrix is the inverse of the covariance matrix Upsilon^{-1}_ab = <u_a u_b>
!
! An optional argumet can be passed, if true the Lambda matrix multiplication is not
! performed. In this way the gradient is automatically preconditioned.
  
subroutine get_gradient_supercell_fast( n_random, natsc, n_modes, ntyp_sc, rho, u_disp, eforces, &
     wr_sc, epols_sc, trans, T, mass, ityp_sc, log_err, grad, grad_err, preconditioned)


  use stochastic
  
  implicit none
  
  integer, intent(in) :: n_random
  ! Number of random configurations
  integer, intent(in) :: natsc
  ! Number of atoms in the supercell
  integer, intent(in) :: n_modes
  ! Number of modes in the supercell (usually 3 * nat_sc)
  integer, intent(in) :: ntyp_sc
  ! Number of different type of atoms in the structure

  double precision, dimension(n_random), intent(in) :: rho
  ! The importance sampling weights of the configurations

  double precision, dimension(n_random, natsc, 3), intent(in) :: u_disp, eforces
  !
  ! The displacement vectors and forces - sscha forces for each configuration.
  ! Displacements are required in bohr, while forces in Ha/bohr
  ! Most importantly, forces and wr unit of measurement must coincide:
  ! If wr is in Ha and u_disp in bohr then forces must be Ha/bohr
  ! if wr is in Ry and u_disp in angstrom then forces must be in Ry/angstrom
  !
  
  double precision, dimension(n_modes), intent(in) :: wr_sc
  double precision, dimension(3 * natsc, n_modes), intent(in) :: epols_sc
  !
  ! Frequencies (in Ha or look above) and polarization vectors for each mode
  ! Polarization vectors are contracted so that the x coord of the i-th atom
  ! of the j-th vector is accessed by epols_sc(x + 3*(i-1), j)
  ! this allows an easy wrapping to the subroutine with python code.
  !
  
  logical, dimension(n_modes), intent(in) :: trans
  ! This is a vector that is True if the mode is a translation, False otherwise
  ! It is needed to generate the Upsilon matrix (as it has a diverging eigenvalue for w->0)

  double precision, intent(in) :: T
  ! Temperature of the system (also this is used to compute the upsilon matrix)

  double precision, dimension(ntyp_sc),  intent(in) :: mass
  ! Mass of the atoms
  ! The unit of measure must be in masses of the electron. (the mass saved in the QE dynamical matrix file
  ! must be multiplied by 2 to match this unit).

  integer, dimension(natsc), intent(in) :: ityp_sc
  ! Type of the atoms in the supercell
  
  character(len=10), intent(in) :: log_err
  ! A flag used to chose the kind of average to be performed 

  double precision, dimension(3*natsc, 3*natsc), intent(out) :: grad, grad_err
  !
  ! The output gradient and its error in real space.
  ! Note that it needs to be symmetrized.
  !
  logical, optional, intent(in) :: preconditioned
  ! If true (default) the gradient is computed without the Lambda matrix mutliplication,
  ! That is as if it was already preconditioned

  
  ! ---------------------------------- END OF INPUT DEFINITION ------------------------------------
  integer i, j, alpha, beta, ical, jcal
  double precision, dimension(3*natsc, 3*natsc) :: uf_mat, err_uf_mat, ups_mat, tmp
  double precision, dimension(3*natsc) :: v_aux1, v_aux2
  double precision, dimension(3*natsc, n_random) :: v_disp, f_aux
  double precision t1, t2
  logical precond
  logical, parameter :: print_input = .false.

  ! Setup the default value of the preconditioned variable
  if (present(preconditioned)) then
     precond = preconditioned
  else
     precond = .true.
  end if

  print *, "AJFNKJ"
  call flush()
  ! DEBUG INPUT
  if (print_input) then
     print *, "GET_GRADIENT_SUPERCELL INPUT VALUES:"
     print *, "NRAND:", n_random
     print *, "NAT_SC:", natsc
     print *, "N_MODES:", n_modes
     print *, "NTYP_SC:", ntyp_sc
     print *, "RHO:", rho(:)
     print *, "W:", wr_sc(:)
     print *, "TRANS:", trans(:)
     print *, "T:", T
     print *, "MASS:", mass
     print *, "ITYP_SC:", ityp_sc
     print *, "LOG_ERR:", log_err

     print *, ""
     print *, "POL VECTORS:"
     do i = 1, 3*natsc
        print "(A10, I8, A10, 1000E13.4)", "MODE", i, "VECTOR", epols_sc(:, i) 
    end do
    
    print *, ""
    print *, "INPUT U and F"
    do i = 1, n_random
        ! Create the auxiliary vector
        do j = 1, natsc
            do alpha = 1, 3
                v_aux1(alpha + 3* (j-1)) = eforces(i, j, alpha)
                v_aux2(alpha + 3*(j-1)) = u_disp(i, j, alpha)
            end do
        end do
        
        print "(A10, I8, A10, 1000E13.4)", "CONF", i, "FORCE", v_aux1(:)
        print "(A10, I8, A10, 1000E13.4)", "CONF", i, "DISP", v_aux2(:)
    end do
  end if
  
  ! If the gradient is preconditioned then use the 
  ! fast code.
  if (.not. precond) then
      print *, "COMPUTING PRECOND!"
      call flush()
      ! Prepare u and f resizing the arrays
      do i = 1, n_random
          do j= 1, natsc
              do alpha = 1, 3
                  v_disp(alpha + 3*(j-1), i) = u_disp(i, j, alpha)
                  f_aux(alpha + 3*(j-1), i) = eforces(i, j, alpha)
              enddo
          enddo
      enddo
      print *, "CALLING THE FUNCTION"
      call flush()
      call get_unpreconditioned_gradient(v_disp, f_aux, epols_sc, wr_sc, T, mass, ityp_sc, &
              rho, log_err, grad, grad_err, 3*natsc, n_random, n_modes, ntyp_sc)
      print *, "EXITING THE FUNCTION"
      call flush()
      return
  endif
  
  ! Compute the <uf> matrix in the supercell
  uf_mat = 0.0d0
  err_uf_mat = 0.0d0
  
  ! Compute the upsilon matrix in the supercell
  call cpu_time(t1)
  call get_upsilon_matrix(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, mass, ityp_sc, T, ups_mat)
  call cpu_time(t2)
  
  print *, " --> [FORT] Time to compute upsilon:", t2 -t1
  call flush()
  
  call cpu_time(t1)
  do i = 1, n_random
    ! Create the auxiliary vector
    do j = 1, natsc
        do alpha = 1, 3
            !v_aux1(alpha + 3* (j-1)) = eforces(i, j, alpha)
            v_aux2(alpha + 3*(j-1)) = u_disp(i, j, alpha)
        end do
    end do
    
    ! Get the new v_disp 
    ! Take into account the -1 in the gradient equation
    !print *, "BEFORED", i
    call flush()
    call dgemv("N", 3*natsc, 3*natsc, 1.0d0, ups_mat, 3*natsc, v_aux2, 1, 0.0d0, v_disp(:, i), 1)
    !print *, "AFTERD", i
    call flush()
  enddo
  call cpu_time(t2)
  
  print *, " --> [FORT] Time to prepare the <(-Yu)f> computation:", t2 -t1
  call flush()
  
  ! First of all compute Yu for all the vectors

  ! print *, "DISP FORCES:"
  ! do i = 1, size(u_disp(:, 1,1))
  !    print *, "CONFIG", i, "RHO:", rho(i)
  !    do j = 1, natsc
  !       print *, "U:", u_disp(i,j,:), "F:", eforces(i,j,:)
  !    end do
  ! end do
  
  
  call cpu_time(t1)
  do alpha = 1, 3
     do beta = 1, 3
        do i = 1, natsc
           do j = 1, natsc
              ical = alpha + (i-1)*3
              jcal = beta + (j-1)*3
              call average_error_weight(v_disp(ical, :) * eforces(:,  j, beta), &
                   rho, log_err, grad(ical, jcal), grad_err(ical, jcal))
              !print *, "Terra di mezzo:", ical, jcal, "UF:", uf_mat(ical, jcal), &
              !     sum(u_disp(:, i, alpha) * eforces(:, j, beta)) / sum(rho)
           end do
        end do
     end do
  end do

  ! Impose the hermitianity
  ! do ical = 1, 3*natsc
  !    do jcal = ical, 3*natsc
  !       uf_mat(ical, jcal) = 0.5d0 * (uf_mat(ical, jcal) + uf_mat(jcal, ical))
  !       err_uf_mat(ical, jcal) = 0.5d0 *dsqrt(err_uf_mat(ical, jcal)**2 + err_uf_mat(ical, jcal)**2)
  !       if (ical /= jcal) then
  !          uf_mat(jcal, ical) = uf_mat(ical, jcal)
  !          err_uf_mat(jcal, ical) = err_uf_mat(ical, jcal)
  !       end if
  !    end do
  ! end do
  call cpu_time(t2)
  
  print *, " Time to compute <(Yu)f> in real space: ", t2 - t1
  call flush()

  !print *, "Time to compute the upsilon matrix:", t2 - t1

  ! Perform the matrix multiplication
  ! Grad = - Upsilon . <u(f - fscha)>
  ! Note the '-' sign is applied only in the gradient not in the error

  if (print_input) then
     print *, "======== UF MAT ========="
     do ical = 1, 3*natsc
        print "(10000E16.5)", uf_mat(:, ical)
     end do
     
     print *, ""
     print *, "========== UPS MAT ========"
     do ical = 1, 3*natsc
        print "(10000E16.5)", ups_mat(:, ical)
     end do
     print *, ""
  end if
     

  ! Symmetrize the gradient
  ! In fact the product of symmetric matrices is not symmetric!!!!
  do ical = 1, 3*natsc-1
     do jcal = ical + 1, 3*natsc
        grad(ical, jcal) = 0.5d0 * (grad(ical, jcal) + grad(jcal, ical))
        grad(jcal,ical) = grad(ical, jcal)
     end do
  end do

  ! ! Print the gradient
  ! print *, "======= GRADIENT NOW ======="
  ! do ical = 1, 3*natsc
  !    print "(1000E16.5)", grad(:, ical)
  ! end do
  ! call flush()
  
  ! Perform the inverse preconditioning if required:
  if (.not. precond) then
     call cpu_time(t1)
     print *, "Computing the inverse preconditioning..."
     call multiply_lambda_tensor(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, &
          mass, ityp_sc, T, grad, tmp, .false.)
     ! The lambda matrix is negative defined so multiply it by -1
     grad = -tmp
     !print *, "Setting the error..."
     !call flush()
     call multiply_lambda_tensor(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, &
          mass, ityp_sc, T, grad_err, tmp, .false.)
     grad_err = tmp
     call cpu_time(t2)
     print *, " --> [FORT] Time to compute the inverse preconditioning:", t2 -t1
  end if
  !print *, "Exiting..."
  call flush()
  
end subroutine get_gradient_supercell_fast



!
! Get the standard gradient
! Without doing a computationally expensive 4-rank matrix product calculation
!
! Parameters
! ----------
!  u : displacement for each config (3*natoms vector, n_random)
!  f : forces for each config(f - fscha) (3*natoms vector, n_random)
!  e_pols : polarization vectors : (3*natoms, nmodes)
!  w : frequencies (Ha) : (nmodes)
!  T : temperature (K)
!  gradient : matrix(3*natoms x 3*natoms)
subroutine get_unpreconditioned_gradient(u, f, e_pols, w, T, mass, &
        ityp, rho, log_err, gradient, err_gradient, n, nrand, nmodes, ntyp)  
    use stochastic
  
    implicit none
    integer, intent(in) :: n, nmodes, nrand, ntyp

    double precision, dimension(n, nrand), intent(in) :: u, f
    double precision, dimension(n, nmodes), intent(in) :: e_pols
    double precision, dimension(nmodes), intent(in) :: w
    double precision, intent(in) :: T
    double precision, dimension(ntyp), intent(in) :: mass
    integer, dimension(n/3), intent(in) :: ityp
    double precision, dimension(nrand), intent(in) :: rho
    character(len=10), intent(in) :: log_err
    double precision, dimension(n,n), intent(out) :: gradient, err_gradient
    
    
    ! Now the standard variables
    integer i, mu, nu, alpha, beta
    double precision, dimension(nmodes, nmodes) :: fmunu
    double precision, dimension(nmodes, n) :: fe_aux
    double precision, dimension(nmodes) :: ynu
    double precision, dimension(nmodes, nrand) :: u_new, f_new
    double precision, dimension(nrand) :: av_aux
    
    print *, "IN THE FUNCTION"
    call flush()
    
    ! Get the fmunu matrix and the ynu vector once for all
    do mu = 1, nmodes
        call get_ymu(w(mu), T, ynu(mu))
        do nu = mu, nmodes
            call get_fmunu(w(mu), w(nu), T, fmunu(mu, nu))
            fmunu(mu,nu) = fmunu(mu,nu)
        enddo
    enddo
    
    print *, "GET THE FE"
    call flush()
    
    ! Get feaux
    do alpha = 1, n
        do nu = 1, nmodes
            fe_aux(nu, alpha) = 0.0d0
            
            do mu = 1, nmodes
                fe_aux(nu, alpha) = fe_aux(nu, alpha) + fmunu(nu, mu) * e_pols(alpha, mu)
            enddo
        enddo
    enddo
    print *, "GET THE F and U"
    call flush()
    
    ! Get the displacements and the forces in the polarization vector basis
    do i = 1, nrand
        do nu = 1, nmodes
            u_new(nu, i) = 0.0d0
            f_new(nu, i) = 0.0d0
            do alpha = 1, n
                u_new(nu, i) = u_new(nu, i) + &
                    e_pols(alpha, nu)*u(alpha, i)*dsqrt(mass(ityp(1+(alpha -1)/3)))
                f_new(nu, i) = u_new(nu, i) + &
                    e_pols(alpha, nu)*u(alpha, i)/dsqrt(mass(ityp(1+(alpha -1)/3)))
            enddo
        enddo
    enddo
    
    print *, "START TO GET THE GRADIENT"
    call flush()
    
    ! Get the gradient
    do alpha = 1, n
        do beta = alpha, n
            av_aux(:) = 0.0d0
            do i = 1, nrand
                av_aux(i) = 0.5*sum(fe_aux(:, alpha)*ynu(:)*u_new(:, i)*f_new(:, i)*e_pols(beta,:))
                av_aux(i) = av_aux(i) + &
                    0.5*sum(fe_aux(:, beta)*ynu(:)*u_new(:, i)*f_new(:, i)*e_pols(alpha,:))
            enddo
            
            ! There is a -1 in front of the gradient equation (TO BE ADDED)
            call average_error_weight(av_aux, rho, log_err,&
                gradient(alpha,beta), err_gradient(alpha,beta))
            gradient(beta,alpha) = gradient(alpha,beta)
            err_gradient(beta,alpha) = err_gradient(alpha,beta)
        enddo
    enddo          
       
    print *, "GOT THE GRADIENT"
    call flush()
end subroutine get_unpreconditioned_gradient


subroutine get_ymu(w, T, y_mu) 
  implicit none
  double precision, intent(out) :: y_mu
  double precision, intent(in) :: w, T

  double precision :: n_mu
  double precision, parameter :: epsilon = 1d-6
  double precision, parameter :: K_TO_HA =  6.336857346553283d-06 / 2
  
  ! Get the occupation numbers
  n_mu = 0.0d0
  if (T .gt. epsilon) then
     n_mu = 1 / (dexp(w / (K_to_Ha * T)) - 1)
  end if
  
  if (w .lt. epsilon) then
      y_mu = 0.0d0
  else
      y_mu = 2 * w / (1 + 2 * n_mu)
  endif
end subroutine get_ymu
