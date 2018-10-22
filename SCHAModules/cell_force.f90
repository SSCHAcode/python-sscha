
subroutine cell_force( fcell, ainv, stress, omega, press, wmassIN )
    DOUBLE PRECISION, intent(out) :: fcell(3,3)
    DOUBLE PRECISION, intent(in) :: stress(3,3), ainv(3,3)
    DOUBLE PRECISION, intent(in) :: omega, press
    DOUBLE PRECISION, intent(in), optional :: wmassIN
    integer        :: i, j
    DOUBLE PRECISION :: wmass

    IF (.not. present(wmassIN)) THEN
        wmass = 1.0
    ELSE
        wmass = wmassIN
    END IF

    print *, "omega:", omega
    print *, "press:", press
    print *, "wmass:", wmass
    call flush()

    do j=1,3
        do i=1,3
        fcell(i,j) = ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) + ainv(j,3)*stress(i,3)
        end do
    end do
    do j=1,3
        do i=1,3
        fcell(i,j) = fcell(i,j) - ainv(j,i) * press
        end do
    end do
    fcell = omega * fcell / wmass
    return
end subroutine cell_force
