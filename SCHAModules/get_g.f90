
subroutine get_g (a, wr, transmode, T, g, n_mode)

    use thermodynamic
    implicit none
  

    double precision, dimension(n_mode), intent(in) :: a, wr
    logical, dimension(n_mode), intent(in) :: transmode
    double precision, intent(in) :: T
    double precision, dimension(n_mode,n_mode), intent(out) :: g                 
  
    integer :: mu, nu
    integer :: n_mode
    double precision, allocatable, dimension(:) :: da

    logical, parameter :: debug = .true.


    if (debug) then
      print *, "=== DEBUG GET_G ==="
      print *, "N_MODE:", n_mode
      call flush() 
    end if

    allocate(da(n_mode))

    call w_to_da(wr, T, da, n_mode)
  
    ! Calculate  the matrix g that will enter in the final equation
  
    do mu = 1, n_mode
      do nu = 1, n_mode
        if (transmode(mu) .or. transmode(nu)) then
          g(mu,nu) = 0.0d0
        else if (mu .eq. nu) then
          g(mu,nu) = da(mu) / wr(mu) / a(mu)**3.0d0
        else if ( mu .ne. nu .and. abs((wr(mu)-wr(nu))/wr(nu)) .lt. 0.00001d0) then
          !g(mu,nu) = 0.0d0
          g(mu,nu) = da(mu) / wr(mu) / a(mu)**3.0d0
        else
          g(mu,nu) = (a(mu)**2.0d0 - a(nu)**2.0d0) / (wr(mu)**2.0d0 - wr(nu)**2.0d0)  / &
                     ( a(nu)**2.0d0 * a(mu)**2.0d0 )
        end if
      end do
    end do
  
  end subroutine get_g
  