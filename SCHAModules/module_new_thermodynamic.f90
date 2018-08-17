module thermodynamic

  PUBLIC :: w_to_a
  
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
end module thermodynamic
