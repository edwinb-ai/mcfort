module observables
    use types, only: dp
    use parameters
    implicit none
    save
    public rdf, sq
contains
    ! This subroutine calculates the pair potential between particles i & j
    subroutine rdf(x, y, z, dr, g, pbc)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(inout) :: g(:)
    real(dp), intent(in) :: dr
    integer, intent(in) :: pbc

    ! Local variables
    integer :: i, j, nbin
    real(dp) :: rij, xij, yij, zij

    do i = 1, np - 1
        do j = i + 1, np
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            if (pbc == 1) then
                xij = xij - boxl*nint(xij/boxl)
                yij = yij - boxl*nint(yij/boxl)
                zij = zij - boxl*nint(zij/boxl)
            end if
            rij = norm2([xij, yij, zij])

            nbin = int(rij/dr) + 1
            if (nbin <= mr) then
                g(nbin) = g(nbin)+2.0_dp
            end if
        end do
    end do
    end subroutine rdf ! out g

    subroutine sq(x, y, z, s, pbc)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: pbc
    real(dp), intent(inout) :: s(:)
    
    ! Local variables
    real(dp) :: auxc(nvq), auxs(nvq)
    integer :: i, k, j
    real(dp) :: xaux, yaux, zaux, rij, arg, sum, parti, auxsq

    do i = 2, mr
        auxc = 0.0_dp
        auxs = 0.0_dp

        parti = 0.0_dp
        do k = 1, np
            if (pbc == 1) then
                xaux = x(k)-boxl*nint(x(k)/boxl)
                yaux = y(k)-boxl*nint(y(k)/boxl)
                zaux = z(k)-boxl*nint(z(k)/boxl)
            else
                xaux = x(k)
                yaux = y(k)
                zaux = z(k)
            end if
            rij = norm2([xaux, yaux, zaux])

            if (rij < rc) then
                parti = parti + 1.0_dp

                do j = 1, nvq
                arg = qx(i,j)*xaux + qy(i,j)*yaux + qz(i,j)*zaux
                auxc(j) = auxc(j) + cos(arg)
                auxs(j) = auxs(j) + sin(arg)
                end do
            end if
        end do

        sum = 0.0_dp
        do j = 1, nvq
            sum = sum + auxc(j)**2 + auxs(j)**2
        end do

        auxsq = sum / (nvq * parti)
        s(i) = s(i) + auxsq
    end do
    end subroutine sq
end module observables