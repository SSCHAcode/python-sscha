! by LORENZO MONACELLI
! 
!
! This subroutine computes the upsilon matrix
! It is the inverse of the 
!
subroutine get_upsilon_matrix(nmodes, nat, ntyp, wr, epols, trans, mass, ityp, T, ups)
  implicit none
  integer, intent(in) :: nmodes, nat, ntyp
  !
  ! The number of modes
  ! The number of atoms
  ! The number of different types
  !

  double precision, dimension(nmodes), intent(in) :: wr
  !
  ! The frequencies (in Ha) used to compute upsilon matrix
  !

  double precision, dimension(3*nat, nmodes), intent(in) :: epols
  !
  ! Polarization vectors
  !

  logical, dimension(nmodes), intent(in) :: trans
  ! True if the mode is a translation, false otherwise

  double precision, dimension(ntyp), intent(in) :: mass
  ! Mass of each type of atoms in the structure
  
  integer, dimension(nat), intent(in) :: ityp
  ! Type of each atom in the structure
  
  double precision, intent(in) :: T
  ! Temperature in K

  double precision, dimension(3*nat, 3*nat), intent(out) :: ups
  !
  ! The output upsilon matrix in bohr^-2
  !

  ! --------------------------------- END OF INPUT DEFINITION -----------------------------
  double precision, parameter :: K_to_Ha =  6.336857346553283d-06 / 2
  integer :: i, j, k
  double precision :: eigenvalue, nb
  
  ups = 0.0d0
  do i = 1, nmodes
     ! Avoid translational modes
     if (trans(i)) cycle
     
     eigenvalue = 2 * wr(i)

     if ( T /= 0) then
        nb = 1 / (dexp(wr(i) / (T * K_to_Ha)) - 1)
        eigenvalue = eigenvalue / (1 + 2*nb)
     end if

     do j = 1, 3*nat
        do k= 1, 3*nat
           ups(j,k) = ups (j, k) + eigenvalue * epols(j,i) * epols(k,i) * &
               dsqrt(mass(ityp(1+ (j-1) / 3)) * mass(ityp(1 + (k-1)/3)))
        end do
     end do
  end do
end subroutine get_upsilon_matrix
