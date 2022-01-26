module energies
    use ieee_arithmetic, only: ieee_positive_inf, ieee_value
    use types
    use parameters
    use potentials
    use omp_lib
    
    implicit none
    
    public energy, denergy

contains
    ! This configuration calculates the energy of a given configuration
    subroutine energy(x, y, z, ener)
        ! Input variables
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: ener

        ! Local variables
        integer :: i, j
        real(dp) :: rij, xij, yij, zij, uij
        
        ! Variable initialization
        ener = 0.0_dp

        !$omp parallel default(shared) private(i,j,uij,xij,yij,zij,rij)
        !$omp do reduction(+:ener)
        do i = 1, np
            do j = 1, np
                if (i == j) cycle
                uij = 0.0_dp

                xij = x(j) - x(i)
                yij = y(j) - y(i)
                zij = z(j) - z(i)

                ! Minimum image convention
                xij = xij-boxl*nint(xij/boxl)
                yij = yij-boxl*nint(yij/boxl)
                zij = zij-boxl*nint(zij/boxl)

                rij = norm2([xij, yij, zij])

                if (rij < rc) then
                    call pseudohs(rij, uij)
                    ener = ener + uij
                end if
            end do
        end do
        !$omp end do
        !$omp end parallel

        ener = ener / 2.0_dp
    end subroutine energy

    ! This subroutine calculates the difference in energy when a particle is displaced
    subroutine denergy(x, y, z, no, dener)
        ! Input variables
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: dener
        integer, intent(in) :: no

        ! Local variables
        integer :: i
        real(dp) :: rij, xij, yij, zij, uij

        ! Variable initialization
        dener = 0.0_dp

        !$omp parallel do default(shared) private(i,xij,yij,zij,rij,uij) &
        !$omp reduction(+:dener)
        do i = 1, np
            if ( i == no ) cycle

            xij = x(no)-x(i)
            yij = y(no)-y(i)
            zij = z(no)-z(i)
            
            ! Minimum image convention
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            
            rij = norm2([xij, yij, zij])

            if (rij < rc) then
                call pseudohs(rij, uij)
                dener = dener + uij
            end if
        end do
        !$omp end parallel do
    end subroutine denergy
end module energies