module Parameters
    implicit none
    real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510
    real, parameter :: YerhoE2a = 1.52588e-4 ! Converts electron fraction times density (g/cc) times neutrino energy (GeV) to eV^2
    real, parameter :: eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4
end module Parameters

! Probability_Matter_LBL calculates all nine oscillation probabilities including
! the matter effect in an optimized, fast, and efficient way. The precision can
! be controlled with N_Newton. For many applications N_Newton=0 may be enough,
! but many years of DUNE or HK-LBL may require N_Newton=1. This code may be
! suitable for atmospheric neutrinos. The code is standalone with the Parameters
! module for several constants.
!
! Inputs:
!   mixing angles (usual parameterization)
!   phase (usual parameterization) make Dmsq31 positive/negative for the NO/IO
!   Delta msq's (eV^2)
!   L (km)
!   E (GeV) positive for neutrinos, negative for antineutrinos
!   rho (g/cc)
!   Ye: electron fraction, typically around 0.5
!   N_Newton: number of Newton's method iterations to do. should be zero, one, two (or higher)
! Outputs:
!   probs_returned is all nine oscillation probabilities: e.g. probs_returned(2,1) is mu->e
subroutine Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)
    use Parameters
    implicit none
    real, intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye
    integer, intent(in) :: N_Newton
    real, intent(out) :: probs_returned(3, 3)

    real :: c13sq, sind, cosd, Jrr, Jmatter, Dmsqee, Amatter
    real :: Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq
    real :: A, B, C
    real :: See, Tee, Smm, Tmm
    real :: xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32
    real :: Xp2, Xp3, PiDlambdaInv
    real :: Lover4E, D21, D32
    real :: sinD21, sinD31, sinD32
    real :: sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin
    real :: Pme_CPC, Pme_CPV, Pmm, Pee
    integer :: i

    ! --------------------------------------------------------------------- !
    ! First calculate useful simple functions of the oscillation parameters !
    ! --------------------------------------------------------------------- !
    c13sq = 1 - s13sq

    ! Ueisq's
    Ue2sq = c13sq * s12sq
    Ue3sq = s13sq

    ! Umisq's, Utisq's and Jvac     
    Um3sq = c13sq * s23sq
    ! Um2sq and Ut2sq are used here as temporary variables, will be properly defined later     
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)
      
    Jrr = sqrt(Um2sq * Ut2sq)
    sind = sin(delta)
    cosd = cos(delta)

    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jmatter = 8 * Jrr * c13sq * sind
    Amatter = Ye * rho * E * YerhoE2a
    Dmsqee = Dmsq31 - s12sq * Dmsq21

    ! calculate A, B, C, See, Tee, and part of Tmm
    A = Dmsq21 + Dmsq31 ! temporary variable
    See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq
    Tmm = Dmsq21 * Dmsq31 ! using Tmm as a temporary variable      
    Tee = Tmm * (1 - Ue3sq - Ue2sq)
    C = Amatter * Tee
    A = A + Amatter

    ! ---------------------------------- !
    ! Get lambda3 from lambda+ of MP/DMP !
    ! ---------------------------------- !
    xmat = Amatter / Dmsqee
    lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + sqrt((1 - xmat) ** 2 + 4 * s13sq * xmat))

    ! ---------------------------------------------------------------------------- !
    ! Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) !
    ! ---------------------------------------------------------------------------- !
    B = Tmm + Amatter * See ! B is only needed for N_Newton >= 1
    do i = 1, N_Newton
        lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B) ! this strange form prefers additions to multiplications
    enddo ! i, 1, N_Newton

    ! ------------------- !
    ! Get  Delta lambda's !
    ! ------------------- !
    Dlambda21 = sqrt((A - lambda3) ** 2 - 4 * C / lambda3)
    lambda2 = 0.5 * (A - lambda3 + Dlambda21)
    Dlambda32 = lambda3 - lambda2
    Dlambda31 = Dlambda32 + Dlambda21

    ! ----------------------- !
    ! Use Rosetta for Veisq's !
    ! ----------------------- !
    ! denominators      
    PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21)
    Xp3 = PiDlambdaInv * Dlambda21
    Xp2 = -PiDlambdaInv * Dlambda31

    ! numerators
    Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3
    Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2

    Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq
    Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A)

    Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3
    Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2

    ! ------------- !
    ! Use NHS for J !
    ! ------------- !
    Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv

    ! ----------------------- !
    ! Get all elements of Usq !
    ! ----------------------- !
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq

    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq

    ! ----------------------- !
    ! Get the kinematic terms !
    ! ----------------------- !
    Lover4E = eVsqkm_to_GeV_over4 * L / E

    D21 = Dlambda21 * Lover4E
    D32 = Dlambda32 * Lover4E
      
    sinD21 = sin(D21)
    sinD31 = sin(D32 + D21)
    sinD32 = sin(D32)

    triple_sin = sinD21 * sinD31 * sinD32

    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32

    ! ------------------------------------------------------------------- !
    ! Calculate the three necessary probabilities, separating CPC and CPV !
    ! ------------------------------------------------------------------- !
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 &
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 &
            + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    Pme_CPV = -Jmatter * triple_sin

    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 &
                 + Um3sq * Um1sq * sinsqD31_2 &
                 + Um3sq * Um2sq * sinsqD32_2)

    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 &
                 + Ue3sq * Ue1sq * sinsqD31_2 &
                 + Ue3sq * Ue2sq * sinsqD32_2)

    ! ---------------------------- !
    ! Assign all the probabilities !
    ! ---------------------------- !
    probs_returned(1, 1) = Pee                              ! Pee
    probs_returned(1, 2) = Pme_CPC - Pme_CPV                ! Pem
    probs_returned(1, 3) = 1 - Pee - probs_returned(1, 2)   ! Pet

    probs_returned(2, 1) = Pme_CPC + Pme_CPV                ! Pme
    probs_returned(2, 2) = Pmm                              ! Pmm
    probs_returned(2, 3) = 1 - probs_returned(2, 1) - Pmm   ! Pmt

    probs_returned(3, 1) = 1 - Pee - probs_returned(2, 1)
    probs_returned(3, 2) = 1 - probs_returned(1, 2) - Pmm
    probs_returned(3, 3) = 1 - probs_returned(1, 3) - probs_returned(2, 3)
end subroutine Probability_Matter_LBL

subroutine Probability_Vacuum_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, probs_returned)
    use Parameters
    implicit none
    real, intent(in) :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E
    real, intent(out) :: probs_returned(3, 3)

    real :: c13sq, sind, cosd, Jrr, Jvac
    real :: Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq
    real :: Lover4E, D21, D31
    real :: sinD21, sinD31, sinD32
    real :: sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin
    real :: Pme_CPC, Pme_CPV, Pmm, Pee

    ! --------------------------------------------------------------------- !
    ! First calculate useful simple functions of the oscillation parameters !
    ! --------------------------------------------------------------------- !
    c13sq = 1 - s13sq

    ! Ueisq's
    Ue3sq = s13sq
    Ue2sq = c13sq * s12sq

    ! Umisq's, Utisq's and Jvac     
    Um3sq = c13sq * s23sq
    ! Um2sq and Ut2sq are used here as temporary variables, will be properly defined later     
    Ut2sq = s13sq * s12sq * s23sq
    Um2sq = (1 - s12sq) * (1 - s23sq)
      
    Jrr = sqrt(Um2sq * Ut2sq)
    sind = sin(delta)
    cosd = cos(delta)
    Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd
    Jvac = 8 * Jrr * c13sq * sind
    
    ! ----------------------- !
    ! Get all elements of Usq !
    ! ----------------------- !
    Ue1sq = 1 - Ue3sq - Ue2sq
    Um1sq = 1 - Um3sq - Um2sq

    Ut3sq = 1 - Um3sq - Ue3sq
    Ut2sq = 1 - Um2sq - Ue2sq
    Ut1sq = 1 - Um1sq - Ue1sq

    ! ----------------------- !
    ! Get the kinematic terms !
    ! ----------------------- !
    Lover4E = eVsqkm_to_GeV_over4 * L / E

    D21 = dmsq21 * Lover4E
    D31 = dmsq31 * Lover4E
      
    sinD21 = sin(D21)
    sinD31 = sin(D31)
    sinD32 = sin(D31-D21)

    triple_sin = sinD21 * sinD31 * sinD32

    sinsqD21_2 = 2 * sinD21 * sinD21
    sinsqD31_2 = 2 * sinD31 * sinD31
    sinsqD32_2 = 2 * sinD32 * sinD32

    ! ------------------------------------------------------------------- !
    ! Calculate the three necessary probabilities, separating CPC and CPV !
    ! ------------------------------------------------------------------- !
    Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2 &
            + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2 &
            + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2
    
    Pme_CPV = -Jvac * triple_sin

    Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2 &
                 + Um3sq * Um1sq * sinsqD31_2 &
                 + Um3sq * Um2sq * sinsqD32_2)

    Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2 &
                 + Ue3sq * Ue1sq * sinsqD31_2 &
                 + Ue3sq * Ue2sq * sinsqD32_2)

    ! ---------------------------- !
    ! Assign all the probabilities !
    ! ---------------------------- !
    probs_returned(1, 1) = Pee                              ! Pee
    probs_returned(1, 2) = Pme_CPC - Pme_CPV                ! Pem
    probs_returned(1, 3) = 1 - Pee - probs_returned(1, 2)   ! Pet

    probs_returned(2, 1) = Pme_CPC + Pme_CPV                ! Pme
    probs_returned(2, 2) = Pmm                              ! Pmm
    probs_returned(2, 3) = 1 - probs_returned(2, 1) - Pmm   ! Pmt

    probs_returned(3, 1) = 1 - Pee - probs_returned(2, 1)
    probs_returned(3, 2) = 1 - probs_returned(1, 2) - Pmm
    probs_returned(3, 3) = 1 - probs_returned(1, 3) - probs_returned(2, 3)
end subroutine Probability_Vacuum_LBL

program NuFast
    use Parameters
    implicit none
    external :: Probability_Matter_LBL, Probability_Vacuum_LBL, Speed
    real :: s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye
    real :: probs_returned(3, 3)
    integer :: N_Newton, alpha, beta

    ! ------------------------------- !
    ! Set the experimental parameters !
    ! ------------------------------- !
    L = 1300 ! km
    E = 2.5 ! GeV
    rho = 3 ! g/cc
    Ye = 0.5

    ! --------------------------------------------------------------------- !
    ! Set the number of Newton-Raphson iterations which sets the precision. !
    ! 0 is close to the single precision limit and is better than DUNE/HK   !
    ! in the high statistics regime. Increasing N_Newton to 1,2,... rapidly !
    ! improves the precision at a modest computational cost                 !
    ! --------------------------------------------------------------------- !
    N_Newton = 0

    ! ------------------------------------- !
    ! Set the vacuum oscillation parameters !
    ! ------------------------------------- !
    s12sq = 0.31
    s13sq = 0.02
    s23sq = 0.55
    delta = 0.7 * PI
    Dmsq21 = 7.5e-5 ! eV^2
    Dmsq31 = 2.5e-3 ! eV^2

    ! ------------------------------------------ !
    ! Calculate all 9 oscillations probabilities !
    ! ------------------------------------------ !
    call Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, probs_returned)

    ! --------------------------- !
    ! Print out the probabilities !
    ! --------------------------- !
    write (*,"(3(A,F10.2))") "L = ", L, " E = ", E, " rho = ", rho
    write (*,*) "Probabilities:"
    write (*,*) "alpha beta P(nu_alpha -> nu_beta)"
    do alpha = 1, 3
        do beta = 1, 3
            write (*,"(2I3,F10.6)") alpha, beta, probs_returned(alpha, beta)
        end do ! beta, 1, 3
    end do ! alpha, 1, 3

end program NuFast
