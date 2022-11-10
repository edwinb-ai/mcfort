module movement
    use types
    use parameters
    use energies

    implicit none

    public mcmove, adjust
contains
    ! This subroutine displace the system to a new configuration
    subroutine mcmove(x, y, z, ener, nattemp, nacc, del)
    real(dp), intent(in) :: del
    real(dp), intent(inout) :: ener
    integer, intent(inout) :: nattemp, nacc
    real(dp), intent(inout) :: x(:), y(:), z(:)

    ! Local variables
    integer :: no
    real(dp) :: xo, yo, zo, enero, enern, dener
    real(dp) :: rng

    nattemp = nattemp + 1

    call random_number(rng)
    no = int(rng*np)+1
    call denergy(x, y, z, no, enero)

    xo = x(no)
    yo = y(no)
    zo = z(no)

    call random_number(rng)
    x(no) = x(no)+(rng-0.5_dp)*del
    call random_number(rng)
    y(no) = y(no)+(rng-0.5_dp)*del
    call random_number(rng)
    z(no) = z(no)+(rng-0.5_dp)*del

    ! periodic boundary conditions
    x(no) = x(no)-boxl*nint(x(no)/boxl)
    y(no) = y(no)-boxl*nint(y(no)/boxl)
    z(no) = z(no)-boxl*nint(z(no)/boxl)

    call denergy(x, y, z, no, enern)

    dener = enern - enero
    call random_number(rng)
    if (rng < exp(-dener / ktemp)) then
        ener = ener + dener
        nacc = nacc + 1
    else
        x(no) = xo
        y(no) = yo
        z(no) = zo
    end if
    end subroutine mcmove

    subroutine adjust(nattemp, nacc, del, tol)
        ! This subroutine adjusts the displacement of particles
        integer, intent(in) :: nattemp, nacc
        real(dp), intent(in) :: tol
        real(dp), intent(inout) :: del
        ! Local variables
        real(dp) :: ratio, half_box

        ratio = real(nacc, dp)/real(nattemp, dp)
        if (ratio > tol) then
            del = del*1.05_dp
        else
            del = del*0.95_dp
        end if

        half_box = boxl * 0.5_dp
        if (del > half_box) then
            del = half_box
        end if
    end subroutine adjust
end module movement