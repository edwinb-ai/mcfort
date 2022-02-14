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
    real(dp), allocatable :: qx(:, :), qy(:, :), qz(:, :)
    real(dp) :: del, ener, dv, dt2, dphi, qs
    real(dp) :: rng, d, dr, dq
    integer :: nacc, i, j, ncq, nattemp, avefreq
    integer :: limT, limG, u, ng, naveg

    ! PRNG initialization
    call random_seed()

    !! Scalar variable initialization
    ! Read an input file that contains all the necessary information
    call parse_input('input.in', limG)
    ! Update the simulation parameters with this information
    boxl = (np / rho)**(1.0_dp/3.0_dp)
    rc = boxl * 0.5_dp
    d = (1.0_dp/rho)**(1.0_dp/3.0_dp)
    dr = rc / mr
    dq = pi / rc
    limT = 1e8 ! Thermalization steps
    ng = 0
    naveg = 0
    ncq = 0
    nacc = 1
    nattemp = 0
    del = 0.1_dp ! Particle displacement
    avefreq = 10000 ! Average frequency

    !! Array initialization
    ! Allocate memory for arrays
    allocate(x(np), y(np), z(np))
    allocate(r(mr), g(mr), s(mr), q(mr))
    allocate(qx(mr, nvq), qy(mr, nvq), qz(mr, nvq))
    ! Give initial values
    g = 0.0_dp
    q = 0.0_dp

    print*, 'rc = ', rc
    print*, 'dr = ', dr
    print*, 'dq = ', dq, 'boxl =', boxl
    print*, 'Mean interparticle distance: ', d
    print*, rho

    ! Gives values to q vector
    if (stfac == .true.) then
        init_sq(q, qx, qy, qz, 'qvectors.dat')
    end if

    ! initial configuration and initial energy
    call iniconfig(x, y, z, d)
    call energy(x, y, z, ener)

    print*, 'E/N for the initial configuration:', ener/np

    ! MC cycle to thermalize the system
    open(newunit=u, file='energy.dat', status='unknown')
    do i = 1, limT
        call mcmove(x, y, z, ener, nattemp, nacc, del)
        call adjust(nattemp, nacc, del, 0.35_dp)
        
        if (mod(i, 100000) == 0) then
            write(unit=u, fmt='(i15.10,f15.10)') i, ener/real(np,dp)
        end if
        
        if (mod(i, 1000000) == 0) then
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N'
            print*, i, del, ener/np
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
    do i = 1, limG
        call mcmove(x, y, z, ener, nattemp, nacc, del)
        call adjust(nattemp, nacc, del, 0.35_dp)

        if (mod(i, avefreq) == 0) then
            naveg = naveg + 1
            write(unit=output_unit, fmt='(a)') 'calculating g(r) and S(q)'
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N'
            print*, i, del, ener/np
            ! Accumulation step for the RDF
            call rdf(x, y, z, dr, g, pbc)
            ! Accumulation step for the structure factor
            ! call sq(x, y, z, s, pbc)
        end if
    end do

    ! This is the radial distribution function
    open(newunit=u, file='gr.dat', status='unknown')
    do i = 2, mr
        r(i) = (i-1)*dr
        dv = (4.0_dp * pi * r(i)*r(i) * dr) * rho
        g(i) = g(i) / (np*naveg*dv)
        write(unit=u, fmt='(2f15.8)') r(i), g(i)
    end do
    close(u)

    ! This is the structure factor from the definition
    ! open(newunit=u, file='sq.dat', status='unknown')
    ! do i = 3, mr
    !     s(i) = s(i) / naveg
    !     if (q(i) < 40.0_dp) then
    !         write(unit=u, fmt='(2f15.7)') q(i), s(i)
    !     end if
    ! end do
    ! close(u)

    deallocate(x, y, z, r, g, s, q)

end program main