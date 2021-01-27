module thermodynamic

  PUBLIC :: w_to_a
  PUBLIC :: w_to_da
  
  ! The following subroutine translates the frequency in
  ! the a distance
  ! Note w must be in Ha units and T in K, then a is in Bohr.
contains
  subroutine w_to_a(w,T, a, n)
    double precision, intent(in) :: T
    double precision, dimension(n), intent(in) :: w
    double precision, dimension(n), intent(out) :: a
    integer, intent(in) :: n
    
    a = 0.0d0
    if (T .eq. 0.0D0) then
       a(:) = dsqrt(1.0d0 / (2.0d0 * w(:)))
    else
       a(:) = dsqrt((1.0d0 / dtanH(0.5d0 * w * &
            315774.65221921849D0 / T )) / &
            (2.0d0 * w(:)) )
    end if
  end subroutine w_to_a


  ! Computes da/dw starting from the w frequencies
  ! w unit <= Ha
  ! T unit <= K
  ! a unit => Bohr
  subroutine w_to_da(w,T, da, n)
    double precision, intent(in) :: T
    double precision, dimension(n), intent(in) :: w
    double precision, dimension(n), intent(out) :: da


    integer, intent(in) :: n
    double precision, dimension(n) :: a, b
    double precision :: beta


    if (T .eq. 0.0D0) then
       da = - dsqrt(1.0d0 / (8.0d0 *  w**3.0d0))
    else
       !    result_w_to_da = - dsqrt((1.0d0 / dtanH(0.5d0 * w * &
       !                      315774.65221921849D0 / T )) / &
       !                     (8.0d0 * m * w**3.0d0) )     * &
       !                     (1.0d0 +  (w * 315774.65221921849D0 / T) * &
       !                     (1.0d0 / dsinH(0.5d0 * w * &
       !                      315774.65221921849D0 / T )))
       beta =  315774.65221921849D0 / T
       a = w * beta + dsinH(w * beta)
       b = dsqrt(1.0d0 / (32.0d0 *  (w**3.0d0) * &
            (dsinH(0.5d0 * w * beta)**3.0) *  dcosH(0.5d0 * w * beta)))  
       da = - a * b
    end if
  end subroutine w_to_da

  
  ! This function calculates the value of the derivative of the
  ! total free energy of the harmonic oscillator minus the potential
  ! of the harmonic oscillator with respect to the arbitrary frequency:
  !
  ! d [ F_0 - 1/2 m W^2 <u^2>_0 ] / dW   
  ! 
  ! The frequency needs to
  ! be given in Ha, the temperature in K and the mass in Ha atomic
  ! units.
  
  function dW_f0_u0(w,T) result(result_dW_f0_u0)

    double precision, intent(in) :: w, T
    double precision :: result_dW_f0_u0

    if (T .eq. 0.0d0) then
       result_dW_f0_u0 = 0.25d0
    else 
       result_dW_f0_u0 = 0.25d0 * (2.0d0 * nb(w,T) + 1.0d0 +    & 
            2.0d0 * (315774.65221921849D0 / T) * w * &
            dexp(w * 315774.65221921849D0 / T) *      &
            nb(w,T)**2.0d0)
    end if
  end function dW_f0_u0

  
  ! This function calculates the Bose-Einstein distribution function for
  ! temperature T given in K and frequency given in Ha
  function nb(w,T) result(result_nb)
    double precision, intent(in) :: w, T
    double precision :: result_nb
    !  double precision :: T          !MODIFIED
    !  double complex   :: w          !MODIFIED
    !  double complex   :: result_nb  !MODIFIED

    if (T .eq. 0.0D0) then
       result_nb = 0.0D0
    else
       result_nb = 1.0D0 / (dexp(w * 315774.65221921849D0 / T )  - 1.0D0)
    end if

  end function nb

end module thermodynamic
