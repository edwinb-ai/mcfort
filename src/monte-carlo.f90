program main
    use types, only: dp
    use parameters
    use utils
    use energies
    use movement
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    ! Local variables, note that somes variables are initialized
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: r(:), g(:), q(:), s(:)
    real(dp) :: del = 0.1_dp, ener, dv, dt2, dphi, qs
    real(dp) :: rng, d, dr, dq
    integer :: nattemp = 0
    integer :: nacc = 1, nacco, nav, i, j, ncq = 0
    integer :: ng = 0, naveg = 0
    integer, parameter :: limT = 2e7
    integer :: limG, u
    ! Condiciones periódicas a la frontera
    integer :: pbc = 1

    ! Inicializar el RNG
    call random_seed()

    ! Read an input file that contains all the necessary information
    call parse_input('input.in', limG)
    ! Update the simulation parameters with this information
    boxl = (np / rho)**(1.0_dp/3.0_dp)
    ! rc = boxl * 0.5_dp
    rc= 3.0_dp
    d = (1.0_dp/rho)**(1.0_dp/3.0_dp)
    dr = rc / mr
    dq = pi / rc

    print*, 'rc = ', rc
    print*, 'dr = ', dr
    print*, 'dq = ', dq, 'boxl =', boxl
    print*, 'Mean interparticle distance: ', d
    print*, rho

    ! Allocate memory for arrays
    allocate(x(np), y(np), z(np))
    allocate(r(mr), g(mr), s(mr), q(mr))

    ! Gives values to q vector
    do i=1, mr
        q(i) = (i-1)*dq
    end do

    ! Randomly initialize the q-vectors
    allocate( qx(mr, nvq), qy(mr, nvq), qz(mr, nvq) )
    open(newunit=u, file='qvectors.dat', status='unknown')
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

    ! initial configuration and initial energy
    call iniconfig(x, y, z, d)
    call energy(x, y, z, ener)

    print*, 'E/N for the initial configuration:', ener/np

    ! MC cycle to thermalize the system
    open(newunit=u, file='energy.dat', status='unknown')
    do i = 1, limT
        call mcmove(x, y, z, ener, nattemp, nacc, del)
        call adjust(nattemp, nacc, del, 0.35_dp)
        
        if (mod(i, 1000) == 0) then
            write(unit=u, fmt='(2f15.7)') real(i, dp), ener/np
        end if
        
        if (mod(i, 1000000) == 0) then
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N'
            write(unit=output_unit, fmt='(3f14.10)') real(i, dp), del, ener/np
        end if
    end do

    print*, 'The system has thermalized'
    close(u)
    ! write the final configuration and the energy
    open(newunit=u, file='finalconf.dat', status='unknown')
    do i = 1, np
        write(unit=u, fmt='(3f15.7)') x(i), y(i), z(i)
    end do
    close(u)

    !MC cycle to calculate the g(r)
    nacco = nacc
    g = 0.0_dp

    do i = 1, limG
        call average(x, y, z, g, s, ener, nattemp, nacc, ng, naveg, del, dr, pbc)
        call adjust(nattemp, nacc, del, 0.35_dp)

        if (mod(i, 2500) == 0) then
            write(unit=output_unit, fmt='(a)') 'calculating g(r) and S(q)'
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N'
            write(unit=output_unit, fmt='(3f14.10)') real(i, dp), del, ener/np
        end if
    end do

    ! True value for averaging
    nav = nacc - nacco

    ! This is the radial distribution function
    open(newunit=u, file='gr.dat', status='unknown')
    do i = 2, mr
        r(i) = (i-1)*dr
        dv = (4.0_dp * pi * r(i)**2 * dr) * rho
        g(i) = g(i) / (np*naveg*dv)
        write(unit=u, fmt='(2f15.8)') r(i), g(i)
    end do
    close(u)

    ! This is the structure factor from the definition
    open(newunit=u, file='sq.dat', status='unknown')
    do i = 3, mr
        s(i) = s(i) / naveg
        if (q(i) < 40.0_dp) then
            write(unit=u, fmt='(2f15.7)') q(i), s(i)
        end if
    end do
    close(u)

    deallocate(x, y, z, r, g, s, q)

end program main