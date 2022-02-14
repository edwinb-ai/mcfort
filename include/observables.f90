module observables
    use types, only: dp
    use parameters
    use omp_lib
    implicit none
    save
    public rdf, sq, init_sq
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
    real(dp) :: gbins(mr)

    !$omp parallel do default(shared) private(i,j,xij,yij,zij,rij,nbin,gbins)
    do i = 1, np
        do j = 1, np
            if (i == j) cycle
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            if (pbc == 1) then
                xij = xij - boxl*nint(xij/boxl)
                yij = yij - boxl*nint(yij/boxl)
                zij = zij - boxl*nint(zij/boxl)
            end if
            rij = norm2([xij, yij, zij])

            nbin = nint(rij/dr) + 1
            if (nbin <= mr) then
                gbins(nbin) = gbins(nbin) + 2.0_dp
            end if
        end do
    end do
    !$omp end parallel do

    !$omp critical
    do i = 1, mr
        g(i) = g(i) + gbins(i)
    end do
    !$omp end critical
    end subroutine rdf ! out g

    subroutine init_sq(q, qx, qy, qz, filein)
    character(len=*), intent(in) :: filein
    real(dp), intent(inout) :: q(:)
    real(dp), intent(inout) :: qx(:,:), qy(:,:), qz(:,:)

    open(newunit=u, file=filein, status='unknown')

    ! Gives values to q vector
    do i=1, mr
        q(i) = (i-1)*dq
    end do

    ! Create random values in the unit sphere
    do i = 1, mr
        do j = 1, nvq
            call random_number(rng)
            dt2 = pi*rng
            call random_number(rng)
            dphi = 2.0_dp * pi * rng

            qx(i,j) = q(i)*cos(dphi)*sin(dt2)
            qy(i,j) = q(i)*sin(dphi)*sin(dt2)
            qz(i,j) = q(i)*cos(dt2)

            ncq = ncq + 1
            qs = norm2([qx(i,j), qy(i,j), qz(i,j)])

            write(unit=u, fmt='(2f15.7)') real(ncq, dp), qs
        end do
    end do
    close(u)
    end subroutine init_sq

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