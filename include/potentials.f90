module potentials

use ieee_arithmetic, only: ieee_positive_inf, ieee_value
use types
use parameters

implicit none

! Intrinsic pseudohs potential arguments that need to be computed
! only once
real(dp), parameter, private :: dlr = 50.0_dp, dT = 1.4737_dp
real(dp), parameter, private :: dla = 49.0_dp
real(dp), parameter, private :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
real(dp), parameter, private :: bpot = (dlr/dla)**(1.0_dp/(dlr-dla))

! Intrinsic parameters for the composition of square well, hard sphere, square
! shoulder potential
real(dp), parameter :: eps1 = -1.0_dp
real(dp), parameter :: eps2 = 0.25_dp
real(dp), parameter :: eps3 = 0.0_dp
real(dp), parameter :: k1 = 1000.0_dp
real(dp), parameter :: delta1 = pi / k1
real(dp), parameter :: k2 = 1000.0_dp
real(dp), parameter :: delta2 =pi / k2
real(dp), parameter :: kappa = 1000000.0_dp
real(dp), parameter :: lambda1 = 1.15_dp
real(dp), parameter :: lambda2 = 1.25_dp
real(dp), parameter :: lambda3 = 0.0_dp
real(dp), parameter :: Asw = lambda1 - 0.5 * delta1
real(dp), parameter :: Ass = lambda2 - 0.5 * delta2
real(dp), parameter :: ao = 0.4675_dp
real(dp), parameter :: bo = 0.1074_dp
real(dp), parameter :: doo = 0.2049_dp
real(dp), parameter :: cparam = -1.4733_dp
real(dp), parameter :: alpha = ao + bo / (doo + abs(eps3)**cparam)
real(dp), parameter :: delta3 = sqrt(-log(alpha) / kappa)
real(dp), parameter :: A3sw = lambda3 - delta3

! Parameters for the true square well potential
real(dp), parameter :: lambdasw = 1.5_dp

! Export all the potentials
! Make them protected so that their parameters can't be modified accidentally
protected pseudohs, hardsphere, lennardjones, hertzian, &
softsphere, yukawa_attr, gaussian, smooth_sw, square_well

contains
    subroutine smooth_sw(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < bpot) then
            call pseudohs(rij, uij)
            uij = uij + eps1
        elseif ((bpot <= rij) .and. (rij <= Asw)) then
            uij = eps1
        elseif ((Asw < rij) .and. (rij <= (Asw + delta1))) then
            uij = (eps1 + eps2) / 2.0_dp
            uij = uij - ((abs(eps1) + abs(eps2)) / 2.0_dp) * (cos(k1 * (rij - Asw)))
        elseif (((Asw + delta1) < rij) .and. (rij <= Ass)) then
            uij = eps2
        elseif ((Ass < rij) .and. (rij <= (Ass + delta2))) then
            uij = (eps2 + eps3) / 2.0_dp
            uij = uij + ((abs(eps2) + abs(eps3)) / 2.0_dp) * (cos(k2 * (rij - Ass)))
        elseif (((Ass + delta2) < rij) .and. (rij <= A3sw)) then
            uij = eps3
        elseif (rij > A3sw) then
            uij = eps3 * exp(-kappa * (rij - A3sw)**2)
        else
            uij = 0.0_dp
        end if
    end subroutine smooth_sw
    
    subroutine pseudohs(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < bpot) then
            uij = (a2/dT)*((1.0_dp/rij)**dlr-(1.0_dp/rij)**dla)
            uij = uij + 1.0_dp/dT
        else
            uij = 0.0_dp
        end if

    end subroutine pseudohs

    subroutine hardsphere(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        real(dp) :: rinf

        if (rij < 1.0_dp) then
            uij = ieee_value(rinf, ieee_positive_inf)
        else
            uij = 0.0_dp
        end if

    end subroutine hardsphere

    subroutine square_well(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        real(dp) :: rinf

        if (rij < 1.0_dp) then
            uij = ieee_value(rinf, ieee_positive_inf)
        elseif ((1.0_dp < rij) .and. (rij < lambdasw)) then
            uij = -1.0_dp
        else
            uij = 0.0_dp
        end if
    end subroutine square_well

    subroutine lennardjones(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        real(dp) :: temp

        temp = 1.0_dp / (rij**6)

        uij = 4.0_dp * (temp**2 - temp)
        
    end subroutine lennardjones

    subroutine softsphere(rij, uij, n)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        integer, intent(in) :: n

        if (rij <= 1.0) then
            uij = (1.0_dp / rij)**n
        else
            uij = 0.0_dp
        end if

    end subroutine softsphere

    subroutine yukawa_attr(rij, uij, k)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij, k

        real(dp) :: rinf

        if (rij < 1.0) then
            uij = ieee_value(rinf, ieee_positive_inf)
        else
            uij = -(1.0_dp / rij) * exp(-k * (rij - 1.0_dp))
        end if

    end subroutine yukawa_attr

    subroutine gaussian(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        uij = exp(-0.5_dp * rij**2)
    end subroutine gaussian

    subroutine hertzian(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < 1.0_dp) then
            uij = (1.0_dp - rij)**(5.0_dp / 2.0_dp)
        else
            uij = 0.0_dp
        end if

    end subroutine hertzian

end module potentials