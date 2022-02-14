module utils
    use types
    use parameters

    implicit none

    public iniconfig, parse_input

contains
    subroutine iniconfig(x, y, z, d)
    ! defining three vector of mp dimension, it indicate that only are out variables
        real(dp), intent(inout) :: x(:), y(:), z(:)
        real(dp), intent(in) :: d

        ! Local variables
        integer :: i ! it is neccesary for restart i
        real(dp) :: half_box

        half_box = boxl / 2.0_dp

        x(1) = -half_box + (d / 2.0_dp)
        y(1) = -half_box + (d / 2.0_dp)
        z(1) = -half_box + (d / 2.0_dp)

        do i = 1, np-1
            x(i+1) = x(i)+d
            y(i+1) = y(i)
            z(i+1) = z(i)

            if (x(i+1) > half_box) then
                x(i+1) = -half_box + (d / 2.0_dp)
                y(i+1) = y(i+1)+d
                z(i+1) = z(i)

                if (y(i+1) > half_box) then
                    x(i+1) = -half_box + (d / 2.0_dp)
                    y(i+1) = -half_box + (d / 2.0_dp)
                    z(i+1) = z(i+1)+d
                end if
            end if
        end do
    end subroutine iniconfig ! note that the out is x, y, z vector

    subroutine parse_input(filein, limg)
        character(len = *), intent(in) :: filein
        integer, intent(inout) :: limg

        ! Variables locales
        integer :: u

        open(newunit=u, file=filein, status='old')
        ! Leer todas las variables
        read(u, *) rho
        read(u, *) ktemp
        read(u, *) np
        read(u, *) nvq
        read(u, *) mr
        read(u, *) limg
        read(u, *) stfac
        close(u)
        
    end subroutine parse_input
end module utils