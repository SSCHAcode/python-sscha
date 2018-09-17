! By LORENZO MONACELLI
!
!
! Get the gradient of the system in the supercell using the Raffaello Bianco's equation
! EQ. A18b of PRB paper (2017) Second order structural....
!
! The equation is modified to exploit the noise reduction as:
!
!
! dF/dPhi = - Lambda  Upsilon < u (f - f_scha) > 
!
! The Lambda matrix is the Hessian of the minimization Phi variable (Monacelli et al 2018 PRB)
! Upsilon matrix is the inverse of the covariance matrix Upsilon^{-1}_ab = <u_a u_b>
!
! An optional argumet can be passed, if true the Lambda matrix multiplication is not
! performed. In this way the gradient is automatically preconditioned.
  
subroutine get_gradient_supercell( n_random, natsc, n_modes, ntyp_sc, rho, u_disp, eforces, &
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

  double precision, dimension(natsc), intent(in) :: ityp_sc
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
  double precision, dimension(3*natsc, 3*natsc) :: uf_mat, err_uf_mat, ups_mat
  double precision t1, t2
  logical precond

  ! Setup the default value of the preconditioned variable
  if (present(preconditioned)) then
     precond = preconditioned
  else
     precond = .true.
  end if
  
  ! Compute the <uf> matrix in the supercell
  uf_mat = 0.0d0
  err_uf_mat = 0.0d0
  call cpu_time(t1)
  do alpha = 1, 3
     do beta = 1, 3
        do i = 1, natsc
           do j = 1, natsc
              ical = alpha + (i-1)*3
              jcal = beta + (j-1)*3
              call average_error_weight(u_disp(:, i, alpha) * eforces(:,  j, beta), &
                  rho, log_err, uf_mat(ical, jcal), err_uf_mat(ical, jcal))
           end do
        end do
     end do
  end do

  ! Impose the hermitianity
  do ical = 1, 3*natsc
     do jcal = ical, 3*natsc
        uf_mat(ical, jcal) = 0.5d0 * (uf_mat(ical, jcal) + uf_mat(jcal, ical))
        err_uf_mat(ical, jcal) = 0.5d0 *dsqrt(err_uf_mat(ical, jcal)**2 + err_uf_mat(ical, jcal)**2)
        if (ical /= jcal) then
           uf_mat(jcal, ical) = uf_mat(ical, jcal)
           err_uf_mat(jcal, ical) = err_uf_mat(ical, jcal)
        end if
     end do
  end do
  call cpu_time(t2)
  
  print *, " Time to compute <uf> in real space: ", t2 - t1

  ! Compute the upsilon matrix in the supercell
  call get_upsilon_matrix(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, mass, ityp_sc, T, ups_mat)


  ! Perform the matrix multiplication
  ! Grad = - Upsilon . <u(f - fscha)>
  ! Note the '-' sign is applied only in the gradient not in the error
  call dgemm("N", "N", 3*natsc, 3*natsc, -1.0d0, ups_mat, 3*natsc,  uf_mat, 3*natsc, 0.0d0, grad, 3*natsc)
  call dgemm("N", "N", 3*natsc, 3*natsc, 1.0d0, ups_mat, 3*natsc,  err_uf_mat, 3*natsc, 0.0d0, grad_err, 3*natsc)

  !TODO : It could be that now the sum rule must be applied here.
  ! Perform the inverse preconditioning if required:
  if (.not. precond) then
     call multiply_lambda_tensor(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, &
          mass, ityp_sc, T, grad, grad, .false.)
     call multiply_lambda_tensor(n_modes, natsc, ntyp_sc, wr_sc, epols_sc, trans, &
          mass, ityp_sc, T, grad_err, grad_err, .false.)
  end if
  
end subroutine get_gradient_supercell
