!
! By LORENZO MONACELLI
!
!
! Multiply the given matrix for the Lambda tensor, and overwrite the result


subroutine get_fmunu(w_mu, w_nu, T, f_munu) 
  implicit none
  double precision, intent(out) :: f_munu
  double precision, intent(in) :: w_mu, w_nu, T

  double precision :: n_mu, n_nu, dn_dw
  double precision, parameter :: epsilon = 1d-5
  double precision, parameter :: K_TO_HA =  6.336857346553283d-06 / 2
  
  ! Get the occupation numbers
  n_mu = 0.0d0
  n_nu = 0.0d0
  dn_dw = 0.0d0
  if (T /= 0) then
     n_mu = 1 / (dexp(w_mu / (K_to_Ha * T)) - 1)
     n_nu = 1 / (dexp(w_nu / (K_to_Ha * T)) - 1)
     dn_dw = -(K_TO_HA / T) / (2 * dcosh(w_mu * K_TO_HA / T) - 1)
   end if
     
  
  ! If w_mu != w_nu
  if ((abs(w_mu - w_nu) / dsqrt(w_mu* w_nu)) .gt. epsilon) then
     f_munu = (n_mu + n_nu +1) / (4*(w_mu + w_nu)) - (n_mu - n_nu) / (4*(w_mu - w_nu))
  else
     ! Degenerate case
     f_munu = (2 * n_mu + 1) / (8*w_mu) - dn_dw/4
  end if
  f_munu = - f_munu / (w_mu * w_nu)
end subroutine get_fmunu
        
  

subroutine multiply_lambda_tensor(nmodes, nat, ntyp, wr, pols, trans, &
     mass, ityp, T, input_matrix, output_matrix, inverse)
  implicit none
  

  integer, intent(in) :: nmodes, nat, ntyp
  !
  ! The number of modes, atoms and different types in the current structure
  !

  double precision, dimension(nmodes), intent(in) :: wr
  !
  ! Frequencies (in Ha) of the current dynamical matrix
  !

  double precision, dimension(3*nat, nmodes), intent(in) :: pols
  !
  ! Polarization vectors of the current dynamical matrix
  !

  logical, dimension(nmodes), intent(in) :: trans
  !
  ! True if the given mode is a translation, false otherwise
  !

  double precision, dimension(ntyp), intent(in) :: mass
  !
  ! Mass of each atomic type (in bohr^2 / Ha)
  !

  integer, dimension(nat), intent(in) :: ityp
  !
  ! Type of each atom
  !

  double precision, intent(in) :: T
  !
  ! Temperature of the system (in K)
  !

  double precision, dimension(3*nat, 3*nat), intent(in) :: input_matrix
  !
  ! The matrix that must be multiplied by Lambda tensor (and also the result)
  !
  
  double precision, dimension(3*nat, 3*nat), intent(out) :: output_matrix
  !
  ! The result of the multiplication
  !
  logical, intent(in) :: inverse
  !
  ! If true the inverse of Lambda tensor is applied
  
  ! -------------------------------- END OF INPUT DEFINITION ---------------------------------------

  integer :: nu, mu, i, j
  double precision, dimension(3*nat) :: v_aux1
  double precision, dimension(3*nat, nmodes) :: epols_aux ! Polarization vectors renormalized by masses
  double precision :: matrix_munu, fm

  ! Get the polarization vectors renormalized by the masses
  do i = 1, 3*nat
     if (inverse) then
        epols_aux(i, :) = pols(i, :) * dsqrt(mass(ityp( 1+ (i -1)/3)))
     else
        epols_aux(i, :) = pols(i, :) / dsqrt(mass(ityp(1 + (i-1) / 3)))
     end if
  end do
     
  
  
  output_matrix = 0.0d0
  do nu = 1, nmodes
     if (trans(nu)) cycle

     ! Start by computing v_aux1 = matrix |e_nu>
     !print *, "MODE", nu
     !call flush()
     call dgemv("N",  3*nat, 3*nat, 1.0d0, input_matrix, 3*nat, epols_aux(:, nu), 1, 0.0d0, v_aux1, 1)
     do mu = 1, nmodes
        if (trans(mu)) cycle
        ! Compute the matrix element <e_mu | matrix | e_nu>
        matrix_munu =  dot_product(v_aux1, epols_aux(:, mu))
        
        call get_fmunu(wr(mu), wr(nu), T, fm)

        !print *, "WMU:", wr(mu), "WNU:", wr(nu), "FMUNU:", fm, "<emu| M |enu>:", matrix_munu
        !print *, "VAUX:", v_aux1(:)
        call flush()
        ! Multiply the f_munu
        if (inverse) then
           matrix_munu = matrix_munu / fm
        else
           matrix_munu = matrix_munu * fm
        end if

        
        ! Write the output matrix
        do i = 1, 3*nat
           do j = 1, 3*nat
              output_matrix(i, j) = output_matrix(i, j) + matrix_munu * epols_aux(i, mu) * epols_aux(j, nu)
           end do
        end do
     end do
  end do
end subroutine multiply_lambda_tensor
