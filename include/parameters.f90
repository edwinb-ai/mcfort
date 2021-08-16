module parameters
    use types, only: dp
    implicit none
    save
! CONSTANT PARAMETERS
    ! pot arguments
    real(dp), parameter :: dlr = 50.0_dp, dT = 1.4737_dp
    real(dp), parameter :: dla = 49.0_dp
    real(dp), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
    real(dp), parameter :: bpot = (dlr/dla)**(1.0_dp/(dlr-dla))

    ! constant values
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: diam = 1.0_dp

    ! User dependant parameters
    real(dp) :: phi, rho, boxl, rc, ktemp
    integer :: np, nvq, mr
    real(dp), allocatable :: qx(:, :), qy(:, :), qz(:, :)
end module parameters