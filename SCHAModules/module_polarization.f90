module polarization

  implicit none
  
  public :: get_harmonic_energy
  public :: get_harmonic_pressure



contains


  !
  ! GET HARMONIC ENERGY
  ! ===================
  !
  ! This method returns the harmonic energy of the SCHA hamiltonian. It can be used to compute the SCHA pressure.
  !
  ! This is done with the formula:
  ! .. math ::
  !       < \mathcal H > = \sum_m \frac{\hbar \omega_m}{2 \tanh(\frac{\beta\hbar\omega_m}{2})}
  !
  !
  ! Parameters
  ! ----------
  !
  ! - w : type (array of double precision)
  !      The frequencies of the harmonic hoscillators.
  ! - transmode : type (array of logical)
  !      A logical array of the same shape as w, it contains whether the corresponding mode is or not a translation.
  ! - T : type double precision
  !      Temperature of the system. 
  !
  ! Results
  ! ------
  ! - energy : double precision
  !      The average energy of the system.
  !
  function get_harmonic_energy(w, transmode, T) result (energy)
    double precision, dimension(:), intent(in) :: w
    logical, dimension(:), intent(in) :: transmode
    double precision, intent(in) :: T
    double precision :: energy

    integer :: i, n
    double precision :: beta

    ! Convert boltzmann constant
    beta =  315774.65221921849D0 / T

    energy = 0.0d0
    n = size(w)
    do i = 1, n
       ! Skip the translational mode
       if (transmode(i)) continue

       if (T .eq. 0) then
          energy = energy + .5d0 * w(i)
       else
          energy = energy + w(i) / (2.0d0 * dtanh(beta * w(i) / 2.0d0))
       end if
    end do
  end function get_harmonic_energy

  !
  ! GET HARMONIC PRESSURE
  ! =====================
  !
  ! According to the equipartition theorem the harmonic Pressure can be defined as
  ! ..math::
  !     P = - \frac{<\matchal H>}{3\Omega}
  !
  ! Where :math:`\Omega` is the volume of the unit cell. See Docs for details.
  ! The result is in untis of Ha / bohr^3
  !
  ! Parameters
  ! ----------
  !  - opt : type(options)
  !     This type contains all the information about the structure. Used to compute the unit cell volume
  !  - w : array of double precision
  !     The harmonic frequencies
  !  - transmode: array of logical
  !     An array containing the info about which mode is a translation. It is of the same dimension as w.
  !  - T : double precision
  !     Temperature of the system.
  !
  ! Results
  ! -------
  !  - pressure : double precision
  !     The harmonic pressure of the system.
  !     (is the work made to expand the unit cell as the forces where armonic divided by the volume expansion)
  !
  function get_harmonic_pressure(volume, w, transmode, T) result (pressure)
    double precision, intent(in) :: volume
    double precision, dimension(:), intent(in) :: w
    logical, dimension(:), intent(in) :: transmode
    double precision, intent(in) :: T
    double precision :: energy

    pressure = - get_harmonic_energy(w, transmode, T) / (3 * Omega)
  end function get_harmonic_pressure

  !
  ! GET HARMONIC STRESS TENSOR
  ! =====================
  !
  ! This subroutine generalizes the previous one, but it computes the whole harmonic stress tensor.
  !
  ! Parameters
  ! ----------
  !  - opt : type(options)
  !     This type contains all the information about the structure. Used to compute the unit cell volume
  !  - w : array of double precision
  !     The harmonic frequencies
  !  - er_sc : double precision, dimension(nmodes, n_atoms_sc, 3)
  !     The polarization vectors.
  !  - transmode: array of logical
  !     An array containing the info about which mode is a translation. It is of the same dimension as w.
  !  - T : double precision
  !     Temperature of the system.
  !
  ! Results
  ! -------
  !  - harmonic_stress : double precision, dimension(3,3)
  !     The harmonic stress tensor of the system
  !
  subroutine get_harmonic_stress(volume, w, er_sc, T, harmonic_stress)
    double precision, intent(in) :: volume
    double precision, dimension(:), intent(in) :: w
    double precision, dimension(:,:,:), intent(in) :: er_sc
    double precision, intent(in) :: T
    double precision, dimension(3,3), intent(out) :: harmonic_stress

    double precision :: volume, beta
    integer i, j, mu, nmodes

    nmodes = size(w)

    ! ! Compute the volume [bohr^3]
    ! volume = opt%at_sc(1,1)*(opt%at_sc(2,2)*opt%at_sc(3,3) - opt%at_sc(2,3)*opt%at_sc(3,2)) - &
    !      opt%at_sc(2,1) *(opt%at_sc(1,2)*opt%at_sc(3,3) - opt%at_sc(3,2)*opt%at_sc(1,3)) + &
    !      opt%at_sc(3,1) * (opt%at_sc(1,2)*opt%at_sc(2,3) - opt%at_sc(2,2)*opt%at_sc(1,3))

    ! ! Rescale from alat^3 to bohr^3
    ! volume = volume * opt%celldm(1) * opt%celldm(1) * opt%celldm(1)

    ! Convert boltzmann constant
    beta =  315774.65221921849D0 / T

    ! Compute the stress tensor
    do i = 1, 3
       do j = i, 3
          harmonic_stress(i, j) = 0

          ! Sum over all the modes
          do mu = 1, nmodes
             harmonic_stress(i, j) = harmonic_stress(i,j) - w(mu) / (2 * dtanh(beta * w(mu) / 2)) * &
                  sum( er_sc(:, mu, i) * er_sc(:, mu, j))
          end do

          ! Get the pressure dividing by the volume
          harmonic_stress(i, j) = harmonic_stress(i, j) / volume

          ! Get the symmetric element
          if (i /= j) then
             harmonic_stress(j, i) = harmonic_stress(i, j)
          end if
       end do
    end do
  end subroutine get_harmonic_stress

end module polarization
