module observables
    use types, only: dp
    use parameters
    use omp_lib
    implicit none
    save
    public rdf, sq, init_sq, normalize_sq, normalize_gr
contains
    ! This subroutine calculates the pair potential between particles i & j
    subroutine rdf(x, y, z, dr, g)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(inout) :: g(:)
    real(dp), intent(in) :: dr

    ! Local variables
    integer :: i, j, nbin
    real(dp) :: rij, xij, yij, zij
    real(dp) :: gbins(mr)

    !$omp parallel default(shared) private(i,j,xij,yij,zij,rij,nbin,gbins)
    gbins = 0.0_dp

    !$omp do
    do i = 1, np
        do j = 1, np
            if (i == j) cycle
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            xij = xij - boxl*nint(xij/boxl)
            yij = yij - boxl*nint(yij/boxl)
            zij = zij - boxl*nint(zij/boxl)
            rij = norm2([xij, yij, zij])

            nbin = nint(rij/dr) + 1
            if (nbin <= mr) then
                gbins(nbin) = gbins(nbin) + 2.0_dp
            end if
        end do
    end do
    !$omp end do
    
    !$omp critical
    do i = 1, mr
        g(i) = g(i) + (gbins(i) / 2.0_dp)
    end do
    !$omp end critical

    !$omp end parallel
    end subroutine rdf ! out g

    subroutine normalize_gr(g, r, dr, naveg, filein)
    !! Arguments
    character(len=*), intent(in) :: filein
    integer, intent(in) :: naveg
    real(dp), intent(in) :: dr
    real(dp), intent(inout) :: g(:), r(:)
    !! Local variables
    integer :: u, i
    real(dp) :: dv

    open(newunit=u, file=filein, status='unknown')
    do i = 2, mr
        r(i) = (i-1)*dr
        dv = (4.0_dp * pi * r(i)*r(i) * dr) * rho
        g(i) = g(i) / (np*naveg*dv)
        write(unit=u, fmt='(2f15.8)') r(i), g(i)
    end do
    close(u)
    end subroutine normalize_gr
    
    subroutine init_sq(q, qx, qy, qz, dq, filein)
    !! Arguments
    character(len=*), intent(in) :: filein
    real(dp), intent(in) :: dq
    real(dp), intent(inout) :: q(:)
    real(dp), intent(inout) :: qx(:,:), qy(:,:), qz(:,:)

    !! Local variables
    real(dp) :: dphi, qs, dt2, rng
    integer :: i, j, ncq, u
    ncq = 0

    open(newunit=u, file=filein, status='unknown')

    ! Gives values to q vector
    do i=1, mr
        q(i) = (i-1)*dq
    end do

    ! Create random values in the unit sphere
    do i = 1, mr
        do j = 1, nvq
            call random_number(rng)
            dt2 = pi * rng
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

    subroutine sq(x, y, z, qx, qy, qz, s)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(in) :: qx(:,:), qy(:,:), qz(:,:)
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
            xaux = x(k)-boxl*nint(x(k)/boxl)
            yaux = y(k)-boxl*nint(y(k)/boxl)
            zaux = z(k)-boxl*nint(z(k)/boxl)
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

    subroutine normalize_sq(s, q, naveg, filein)
    !! Arguments
    character(len=*), intent(in) :: filein
    integer, intent(in) :: naveg
    real(dp), intent(inout) :: s(:), q(:)
    !! Local variables
    integer :: u, i

    !! This is the structure factor from the definition
    open(newunit=u, file=filein, status='unknown')
    do i = 3, mr
        s(i) = s(i) / naveg
        if (q(i) < 40.0_dp) then
            write(unit=u, fmt='(2f15.7)') q(i), s(i)
        end if
    end do
    close(u)
    end subroutine normalize_sq

end module observables