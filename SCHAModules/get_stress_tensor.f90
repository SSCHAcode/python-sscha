!
! Stress tensors are requested and returned in Ha / [x]^3
! where [x]^3 is the unit of measure of the volume variable passed in input
!
!
! Parameters
! ----------
!      volume : double precision
!          The volume of the unit cell
!
subroutine get_stress_tensor(volume, forces_e, u_vector, abinit_stress_tensors, wr, er, T, rho, log_err, &
        stress_tensor, error_on_stress_tensor, nrandom, natoms_sc, nmodes)
  use stochastic
  !use polarization
  
  implicit none
  
  
  double precision, intent(in) :: volume, T
  !type(dynamical_matrix), intent(in) :: dyn
  double precision, dimension(3,3,nrandom), intent(in) :: abinit_stress_tensors
  double precision, dimension(3,3), intent(out) :: stress_tensor
  double precision, dimension(3,3), intent(out) :: error_on_stress_tensor
  double precision, dimension(nmodes), intent(in) :: wr
  double precision, dimension(nrandom), intent(in) :: rho
  double precision, dimension(natoms_sc, nmodes, 3), intent(in) :: er
  double precision, dimension(nrandom,natoms_sc, 3), intent(in) :: forces_e, u_vector
  character (len=10), intent(in) :: log_err
  integer, intent(in) :: nrandom, natoms_sc, nmodes
  

  ! Additional variables
  double precision, dimension(3,3) :: anharmonic_stress, anharmonic_error, harmonic_stress, abinitio_stress, abinitio_error
  double precision, dimension(:,:,:), allocatable :: f_dot_r!, debug_f_dot_r, f_scha, fscha_dot_r
  !double precision, dimension(:,:), allocatable :: f_atom
  !double precision :: volume, scale_factor, iso_pressure, iso_p_err, tmp
  integer :: i, j, k
  logical :: pon

 ! pon = .false.
  
  ! Manage the optional argument print_on_screen
  ! if (present(print_on_screen)) then
  !    pon = print_on_screen
  ! else
  !    pon = .false.
  ! end if

  !scale_factor = opt%celldm(1)

  ! Initialize the f dot r 
  allocate(f_dot_r(3,3,nrandom))
  !allocate(fscha_dot_r(3,3,nrandom))
  !allocate(f_scha(nrandom, natoms_sc, 3))
  !if (opt%debug_stress) allocate(debug_f_dot_r(3,3,nrandom))

  ! If the debugging flag is active, open a file named 'debug_stress_evaluation.dat'
  ! if (opt%debug_stress) then
  !    open (unit = 23, file="debug_stress_evaluation.dat", status="replace")
  !    write (23, *) "# config, 3P abinitio, 3P with ab-forces, 3P harmonic [Ha/bohr^3]"
  ! end if

  ! Compute the average of the force on each atom to suppress the noise
  ! if (opt%stress_suppress_noise) then
  !    allocate(f_atom(opt%natoms_sc, 3))

  !    do i = 1, opt%natoms_sc
  !       do j = 1,3
  !          call average_error_weight(opt%forces(:, i, j), opt%rho, opt%log_err, &
  !               f_atom(i, j), tmp)
  !       end do
  !    end do
  ! end if

  ! Extract the harmonic sscha forces
  ! if (opt%lrigid .and. .not. opt%single_q) then
  !    do k = 1, opt%nrandom
  !       call get_harmonic_force_from_fc(dyn%phi_sc + dyn%phi_sc_nonanal, opt%u_vector(k,:,:), &
  !            f_scha(k,:,:))
  !    end do
  ! else
  !    do k= 1, opt%nrandom
  !       call get_harmonic_force_from_fc(dyn%phi_sc, opt%u_vector(k,:,:), &
  !            f_scha(k,:,:))
  !    end do
  ! end if


  ! Get the scalar product between forces and atomic positions
  do k = 1, nrandom
     do i = 1, 3
        do j = 1, 3
           f_dot_r(i, j, k) = sum(forces_e(k, :, i) * u_vector(k, :, j))
           ! fscha_dot_r(i, j, k) = sum(f_scha(k, :, i) * u_vector(k, :, j))
           ! if (.not. stress_suppress_noise) then 
           !    f_dot_r(i, j, k) = sum(forces_e(k, :, i) * u_vector(k, :, j))
           ! else
           !    f_dot_r(i, j, k) = sum(( forces(k, :, i)) * u_vector(k, :, j))                
           !    ! f_dot_r(i, j, k) = sum(( opt%forces(k, :, i) - f_atom(:, i)) * opt%u_vector(k, :, j))
           ! end if

           ! if (opt%debug_stress) then
           !    debug_f_dot_r(i,j,k) = sum(forces(k, :, i) * u_vector(k, :, j))
           
           ! end if 
        end do
     end do
 
!     
!     print *, ""
!     print *, "CONF", k
!     do i = 1, natoms_sc
!         print "(A8, I8, A10, 3D16.8, A10, 3D16.8)", "AT:", i, "FORCE:", forces_e(k, i,  :), "DISP:", u_vector(k, i, :)
!     end do
!     
!     print *, "F DOT R:"    
!     print "(3D16.8)", f_dot_r(1, :, k)
!     print "(3D16.8)", f_dot_r(2, :, k)
!     print "(3D16.8)", f_dot_r(3, :, k)

     ! If the debugging flag is active, print the trace of the two stress tensors
     ! if (opt%debug_stress) write (23, "(I10, D20.12, D20.12, D20.12)") k, &
     !      abinit_stress_tensors(1, 1, k) + abinit_stress_tensors(2,2,k) + abinit_stress_tensors(3,3,k), &
     !      - (debug_f_dot_r(1,1,k) + debug_f_dot_r(2,2,k) + debug_f_dot_r(3,3,k)) / volume, &
     !      - (fscha_dot_r(1,1,k) +  fscha_dot_r(2,2,k) + fscha_dot_r(3,3,k)) / volume
  end do

  ! Close the debugging file
  !if (opt%debug_stress) close(23)

  ! Average the Ab-Initio stress
  do i = 1, 3
     do j = 1, 3
        ! if (.not. stress_suppress_noise) then
        !    call average_error_weight(abinit_stress_tensors(i, j, :), opt%rho, opt%log_err, &
        !         abinitio_stress(i, j), abinitio_error(i, j))
        !    call average_error_weight(-f_dot_r(i, j, :) / volume, opt%rho, opt%log_err, &
        !         anharmonic_stress(i, j), anharmonic_error(i, j))
        ! else
           call average_error_weight(abinit_stress_tensors(i, j, :) - f_dot_r(i, j, : ) / volume, rho, &
                log_err, stress_tensor(i, j), error_on_stress_tensor(i, j))


           ! Delete the stress tensor
           ! if (opt%is_there_a_stress_matrix_offset) then
           !    stress_tensor(i,j) = stress_tensor(i,j) - opt%stress_offset_matrix(i,j)
        !    end if
        !  end if
     end do
  end do




  ! TODO ---- Implement symmetry operations from quantum espresso to the STRESS TENSOR -------
  ! For now implemented only the fact that the tensor must be symmetric 

  ! If the old formula is used also the harmonic tensor must be computed.
  ! if (.not. stress_suppress_noise) then
  !    do i = 1, 3
  !       do j = i +1, 3
  !          abinitio_stress(i, j) = 0.5d0 * (abinitio_stress(i,j) + abinitio_stress(j, i))
  !          abinitio_stress(j, i) = abinitio_stress(i, j)

  !          abinitio_error(i, j) = 0.5d0 * dsqrt( abinitio_error(i, j)**2 + abinitio_error(j, i)**2)
  !          abinitio_error(j, i) = abinitio_error(i, j)

  !          anharmonic_stress(i, j) = 0.5d0 * (anharmonic_stress(i, j) + anharmonic_stress(j, i))
  !          anharmonic_stress(j, i) = anharmonic_stress(i, j)

  !          anharmonic_error(i, j) = 0.5d0 * dsqrt( anharmonic_error(i, j)**2 + anharmonic_error(j, i)**2)
  !          anharmonic_error(j, i) = anharmonic_error(i, j)
  !       end do
  !    end do

  !    ! ! Add the convergence offset on the abinitio stress tensor diagonal elements
  !    ! do i = 1, 3
  !    !    abinitio_stress(i,i) = abinitio_stress(i,i) + opt%stress_offset
  !    ! end do

  !    ! ! Delete the stress tensor
  !    ! if (opt%is_there_a_stress_matrix_offset) then
  !    !    abinitio_stress(:,:) = abinitio_stress(:,:) - opt%stress_offset_matrix(:,:)
  !    ! end if

  !    ! Get the harmonic pressure
  !    call get_harmonic_stress(volume, wr_sc, er_sc, T, harmonic_stress)

  !    ! Compute the total stress tensor (and its error)
  !    stress_tensor = abinitio_stress - harmonic_stress + anharmonic_stress
  !    error_on_stress_tensor = dsqrt(abinitio_error**2 + anharmonic_error**2)
  ! end if


!
!  if (pon) then
!     ! Compute also the isotropic pressure
!     iso_pressure = 0
!     iso_p_err = 0
!     do i = 1, 3
!        iso_pressure = iso_pressure + stress_tensor(i, i)
!        iso_p_err = iso_p_err + error_on_stress_tensor(i, i)**2       
!     end do
!     iso_pressure = iso_pressure / 3.0d0
!     iso_p_err = dsqrt(iso_p_err) / 3.0d0
!
!     ! Print details on stdout
!     print *, ""
!     print *, "* * * * * * * * * * * * * *"
!     print *, "*                         *"
!     print *, "*          STRESS         *"
!     print *, "*                         *"    
!     print *, "* * * * * * * * * * * * * *"
!     print *, ""
!     print *, "Input unit of measure: Ry/bohr^3 (Pay attention that the input stress matches this unit)"
!
!     ! If all the components of the stress are evaluated, print them in output
!     if (.not. opt%stress_suppress_noise) then
!        print *, "Ab-initio stress:"
!        do i = 1, 3
!           print "(3D18.10)", (abinitio_stress(i, j) * 2.0d0, j = 1, 3)
!        end do
!        print *, "+-"
!        do i = 1, 3
!           print "(3D18.10)", (abinitio_error(i, j) * 2.0d0, j = 1, 3)
!        end do
!        print *, ""
!        print *, "Anharmonic stress:"
!        do i = 1, 3
!           print "(3D18.10)", (anharmonic_stress(i, j) * 2.0d0, j = 1, 3)
!        end do
!        print *, "+-"
!        do i = 1, 3
!           print "(3D18.10)", (anharmonic_error(i, j) * 2.0d0, j = 1, 3)
!        end do
!        print *, ""
!        print *, "Harmonic stress:"
!        do i = 1, 3
!           print "(3D18.10)", (harmonic_stress(i, j) * 2.0d0, j = 1, 3)
!        end do
!        print *, ""
!     end if
!
!     ! Print the total stress
!     print *, "-- TOTAL STRESS [Ry/bohr^3] --"
!     do i = 1, 3
!        print "(3D18.10)", (stress_tensor(i, j) * 2.0d0, j = 1, 3)
!     end do
!     print *, "+-"
!     do i = 1, 3
!        print "(3D18.10)", (error_on_stress_tensor(i, j) * 2.0d0, j = 1, 3)
!     end do
!     print *, ""
!     print *, "-- TOTAL STRESS [GPa] --"
!     do i = 1, 3
!        print "(3D18.10)", (stress_tensor(i, j) * 29421.6144000, j = 1, 3)
!     end do
!     print *, "+-"
!     do i = 1, 3
!        print "(3D18.10)", (error_on_stress_tensor(i, j) * 29421.6144000, j = 1, 3)
!     end do
!     print *, ""
!     print *, "Isotropic pressure (1/3 of the trace) "
!     print "(A10, D18.10, A4,D18.10, A15)", " P = ",  iso_pressure * 2.0d0, " +- ", iso_p_err*2.d0, "[Ry/bohr^3]"
!     print "(A10, D18.10, A4,D18.10, A15)", " P = ",  iso_pressure*29421.6144000, " +- ", &
!          iso_p_err*29421.6144000, "[GPa]"
!
!     print *, ""
!     print *, "* * * * * END STRESS * * * * *"
!     print *, ""
!  end if

  deallocate(f_dot_r)
end subroutine get_stress_tensor
