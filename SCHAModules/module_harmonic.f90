module harmonic  

public :: get_harmonic_force
public :: get_harmonic_force_from_fc
public :: get_harmonic_energy_from_fc

contains

! Subroutine to calculate the force on each atom from the
! force constants. The force constants can be obtained calculating
! the dynamical matrix at G point in a supercell. This subroutine
! reads a file that is given as an input which is written in
! dynamical matrix format. As an input a given displacement pattern
! for each atom needs to be given.

subroutine get_harmonic_force(file_fc,u,forces)

  implicit none

  character (len=*), intent(in) :: file_fc
  double precision, dimension(:,:), intent(in) :: u
  double precision, dimension(:,:), intent(out) :: forces 

  double precision, dimension(:,:,:,:), allocatable :: fc
  integer :: natom
  integer :: i, j, k, l, atom1, atom2 
  double precision :: kaka

  natom = size(u(:,1))

  allocate(fc(natom,natom,3,3))
 
  open (unit = 1, file = file_fc)

  do i = 1, natom
    do j = 1, natom
      read (unit = 1, fmt=*) atom1, atom2
      do k = 1, 3
        read (unit = 1, fmt=*) fc(atom1,atom2,k,1), kaka, &
                               fc(atom1,atom2,k,2), kaka, &
                               fc(atom1,atom2,k,3), kaka
      end do
    end do
  end do

  close (unit=1)

  do i = 1, natom
    do j = 1, 3
      forces(i,j) = 0.0d0
      do k = 1, natom
        do l = 1, 3     
          forces(i,j) = forces(i,j) - fc(i,k,j,l) * u(k,l)
        end do
      end do
    end do
  end do

  deallocate(fc)
        
end subroutine get_harmonic_force 

! Subroutine to calculate the force on each atom from the
! force constants. The force constants are given as input.        

subroutine get_harmonic_force_from_fc(fc,u,natom,forces)

  implicit none

  double precision, dimension(3*natom,3*natom), intent(in) :: fc
  double precision, dimension(3*natom), intent(in) :: u
  double precision, dimension(natom*3), intent(out) :: forces 

  integer,intent(in) :: natom
  integer :: i, j, k, l,n 

  double precision, parameter :: conv = 3.571064313502028121 !a_to_bohr^2

!  natom = size(u(:,1))
  n =natom *3
  do i=1,n
    forces(i)=0
    do j=1, n
      forces(i)=forces(i)-fc(i,j)*u(j)
    end do
  end do
  
  forces = forces * conv
!  do i = 1, natom
!    do j = 1, 3
!      forces(i,j) = 0.0d0
!      do k = 1, natom
!       do l = 1, 3     
!          forces(i,j) = forces(i,j) - fc(i,j,k,l) * u(k,l)
!        end do
!      end do
!    end do
!  end do    

end subroutine get_harmonic_force_from_fc

! Subroutine to calculate the harmonic energy from the
! force constants. The force constants are given as input.        

subroutine get_harmonic_energy_from_fc(fc,u,v)

  implicit none

  double precision, dimension(:,:), intent(in) :: fc
  double precision, dimension(:), intent(in) :: u
  double precision, intent(out) :: v
  
  integer :: natom
  integer :: i, j, alpha, beta 

  double precision, parameter :: conv = 3.571064313502028121 !a_to_bohr^2

  natom = size(u(:))



  v = 0.0d0
  do i=1, natom
    do j=1,natom
      v = v+ 0.50d0 *u(i)*fc(i,j)*u(j)
    end do
  end do

  v = v*conv
  
!  do i = 1, natom
!    do j = 1, natom
!      do alpha = 1, 3
!        do beta = 1, 3
!          v = v + 0.5d0 * u(i,alpha) * fc(i,alpha,j,beta) * u(j,beta)
!        end do
!      end do
!    end do
!  end do
        
end subroutine get_harmonic_energy_from_fc


subroutine update(fc,u,natom,n_random,forces,v)

  use omp_lib 

  implicit none

  double precision, dimension(3*natom,3*natom), intent(in) :: fc
  double precision, dimension(n_random,3*natom), intent(in) :: u
  double precision, dimension(n_random,natom*3), intent(out) :: forces 
  double precision, dimension(n_random), intent(out) :: v
  
  integer, intent(in) :: natom,n_random

  integer :: i
  
  !$OMP PARALLEL DO COLLAPSE(1) PRIVATE(i) 
  do i=1,n_random
   call get_harmonic_force_from_fc(fc,u(i,:),natom,forces(i,:))
   call get_harmonic_energy_from_fc(fc,u(i,:),v(i))
  end do
  !$OMP END PARALLEL DO

end subroutine update



end module harmonic
