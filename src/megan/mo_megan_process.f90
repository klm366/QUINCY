!>
!> Contains structures and methods for the topmodel for calculating innundation fractions
!>
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-24
!>

MODULE mo_megan_process
#ifndef __NO_JSBACH__
USE mo_kind, only: lp

PUBLIC calc_gamma_pfd, calc_gamma_temp_iso, calc_gamma_temp_other, calc_gamma_age, calc_gamma_sm, calc_gamma_co2
PUBLIC update_megan_runmean

CONTAINS
!------------------------------------------------------------------------------------
  !-----------------------------------calculate adjustment factors: gamma -------------------------
  ! (from original MEGAN 2.10 source code)
  !     Scientific algorithm
  !
  !           Emission = [EF][GAMMA][RHO]
  !         where [EF]    = emission factor (ug/m2h)
  !               [GAMMA] = emission activity factor (non-dimension)
  !               [RHO]   = production and loss within plant canopies
  !                         (non-dimensino)
  !               Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
  !
  !           GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
  !         where [GAMMA_CE]  = canopy correction factor
  !               [GAMMA_age] = leaf age correction factor
  !               [GAMMA_SM]  = soil moisture correction factor
  !               Assumption: [GAMMA_SM]  = 1 (11/27/06)
  !
  !           GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
  !         where [GAMMA_LAI] = leaf area index factor
  !               [GAMMA_P]   = PPFD emission activity factor
  !               [GAMMA_T]   = temperature response factor
  !
  !           Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
  !      Derivation:
  !           Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
  !           Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
  !           Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
  !         where LDF = light dependent function (non-dimension)
  !
  !   For ISOPRENE
  !               Assumption: LDF = 1 for isoprene            (11/27/06)
  !
  !      Final Equation
  !           Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
  !
  !   For NON-ISOPRENE
  !      Final Equation
  !           Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
  !                      { (1-LDF) + [LDF][GAMMA_P] }
  !
  !------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
  PURE SUBROUTINE calc_gamma_lai(lai, gamma)
    ! ----- leaf area index adjustment
    REAL(lp), INTENT(in)  :: lai
    REAL(lp), INTENT(out) :: gamma

    ! ### assumes that LAI is zero where not valid ! ###
    gamma = (0.49_lp * lai) / SQRT(1.0_lp + 0.2_lp * (lai**2._lp))

  END SUBROUTINE calc_gamma_lai

!------------------------------------------------------------------------------------
  ! PURE
  SUBROUTINE calc_gamma_pfd(zenith, pfd, tpfd, avg_pfd, gamma)
    USE mo_jsb_math_constants, ONLY: rad2deg
    ! ----- photosynthetically active radiation adjustment
    !             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
    !
    !             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
    !                                   0<a<180, sin(a) > 0.0
    !           where PHI    = above canopy PPFD transmission (non-dimension)
    !                 Pdaily = daily average above canopy PPFD (umol/m2s)
    !                 a      = solar angle (degree)
    !                 PPFD   = photosynthetic flux density
    !
    !                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
    !                       BBB = (1+0.0005(Pdaily-400))
    !                       GAMMA_P = sin(a)*AAA
    !
    !                       Pac
    !             PHI = -----------
    !                   sin(a)*Ptoa
    !           where Pac  = above canopy PFD (umol/m2s)
    !                 Ptoa = PFD at the top of atmosphere (umol/m2s)
    !
    !             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
    !
    !             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
    !           where DOY = day of year

    REAL(lp), INTENT(in)  :: zenith
    REAL(lp), INTENT(in)  :: pfd, tpfd, avg_pfd
    REAL(lp), INTENT(out) :: gamma

    REAL(lp)              :: sinbeta, PHI
    REAL(lp)              :: AAA

    REAL(lp), PARAMETER :: conv_pfd = 4.766_lp          ! factor to convert avg_pfd in W/m2 to umol/m2/s

    ! ### assumes zenith angle in degrees
    IF (zenith < 1.0e-10_lp) THEN   !zenith = SIN(solar angle)
      gamma = 0.0_lp
    ELSE

      !  ### PFD and PFD at the top of the atmosphere are adjusted for solar angle
      IF(tpfd <= 0.0_lp) THEN
        PHI = 0.0_lp
      ELSEIF(tpfd > 0.0_lp) THEN
        PHI = pfd/tpfd
      ENDIF

      BBB = 1._lp + 0.0005_lp * (avg_pfd*conv_pfd - 400._lp)

      AAA = ( 2.46_lp * BBB * PHI ) - (0.9_lp * PHI**2._lp)

      gamma = zenith * AAA
    ENDIF

    IF ((ASIN(zenith) * rad2deg < 1.0_lp) .AND. (gamma > 0.1_lp)) THEN
      gamma = 0.0_lp
    ENDIF

  END SUBROUTINE calc_gamma_pfd

!------------------------------------------------------------------------------------
  PURE SUBROUTINE calc_gamma_temp_iso(temp, temp_relax, gamma)
    ! ----- temperature adjustment for isoprene
    !                          Eopt*CT2*exp(CT1*x)
    !             GAMMA_T =  ------------------------
    !                        [CT2-CT1*(1-exp(CT2*x))]
    !           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
    !                 Eopt   = 1.75*exp(0.08(Tdaily-297)
    !                 CT1    = 80
    !                 CT2    = 200
    !                 Thr    = hourly average air temperature (K)
    !                 Tdaily = daily average air temperature (K)
    !                 Topt   = 313 + 0.6(Tdaily-297)
    !
    !                 Note: AAA = Eopt*CT2*exp(CT1*x)
    !                       BBB = [CT2-CT1*(1-exp(CT2*x))]
    !                       GAMMA_T = AAA/BBB

    USE mo_physical_constants,          ONLY: argas

    REAL(lp), INTENT(in)  :: temp, temp_relax ! current and average air temperature
                                              ! at canopy level (original: "hourly" instead of current)
    REAL(lp), INTENT(out) :: gamma

    REAL(lp), PARAMETER   :: CT1 = 80._lp, CT2 = 200._lp
    REAL(lp), PARAMETER   :: invR = 1.E3_lp/argas          ! (1/R) expressed in (mol K / kJ)
    REAL(lp)              :: Eopt, Topt, X
    REAL(lp)              :: AAA, BBB

    Eopt = 1.75_lp * EXP(0.08_lp * (temp_relax - 297.0_lp))
    Topt = 313.0_lp +    (0.6_lp * (temp_relax - 297.0_lp))
    X = ((1.0_lp / Topt) - (1.0_lp / temp)) * invR

    AAA = Eopt * CT2 * EXP(CT1 * X)
    BBB = (CT2 - CT1 * (1.0_lp - EXP(CT2 * X)))

    gamma = AAA / BBB

  END SUBROUTINE calc_gamma_temp_iso


!------------------------------------------------------------------------------------
  PURE SUBROUTINE calc_gamma_temp_other(temp, beta, gamma)
    ! ----- temperature adjustment for non-isoprene species
    !             GAMMA_T =  exp[BETA*(T-Ts)]
    !           where BETA   = temperature dependent parameter
    !                 Ts     = standard temperature (normally 303K, 30C)

    REAL(lp), INTENT(in)  :: temp         ! current canopy level air temperature
    REAL(lp), INTENT(in)  :: beta         ! species dependent parameter
    REAL(lp), INTENT(out) :: gamma

    REAL(lp), PARAMETER   :: Ts = 303._lp

    gamma = EXP(beta * (temp - Ts))

  END SUBROUTINE calc_gamma_temp_other

!------------------------------------------------------------------------------------
  !----- adjustment for leaf age
  !### need to verify !
  ! (unchanged from v2.04)
  PURE SUBROUTINE calc_gamma_age(kclass, lai, laip, tempd, gamma)
    !
    !             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
    !           where Fnew = new foliage fraction
    !                 Fgro = growing foliage fraction
    !                 Fmat = mature foliage fraction
    !                 Fold = old foliage fraction
    !                 Anew = relative emission activity for new foliage
    !                 Agro = relative emission activity for growing foliage
    !                 Amat = relative emission activity for mature foliage
    !                 Aold = relative emission activity for old foliage
    !
    !
    !             For foliage fraction
    !             Case 1) LAIc = LAIp
    !             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
    !
    !             Case 2) LAIp > LAIc
    !             Fnew = 0.0  , Fgro = 0.0
    !             Fmat = 1-Fold
    !             Fold = (LAIp-LAIc)/LAIp
    !
    !             Case 3) LAIp < LAIc
    !             Fnew = 1-(LAIp/LAIc)                       t <= ti
    !                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
    !
    !             Fmat = LAIp/LAIc                           t <= tm
    !                  = (LAIp/LAIc) +
    !                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
    !
    !             Fgro = 1 - Fnew - Fmat
    !             Fold = 0.0
    !
    !           where
    !             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
    !                = 2.9                                  Tt >  303
    !             tm = 2.3*ti
    !
    !             t  = length of the time step (days)
    !             ti = number of days between budbreak and the induction of
    !                  emission
    !             tm = number of days between budbreak and the initiation of
    !                  peak emissions rates
    !             Tt = average temperature (K) near top of the canopy during
    !                  current time period (daily ave temp for this case)
    !             LAIc = current month's LAI
    !             LAIp = previous month's LAI
    !
    !             Relative emission activity
    !             Class           Anew       Agro       Amat       Aold
    !             Constant         1.0        1.0        1.0        1.0
    !             Monoterpenes     2.0        1.8        0.95       1.0
    !             Sesquiterpenes   0.4        0.6        1.075      1.0
    !             Methanol         3.0        2.6        0.85       1.0
    !             Isoprene         0.05       0.6        1.125      1.0
    !
    !             Species:
    !             Constant:  'ACTO','ACTA','FORM','CH4','NO','CO'
    !             Monoterpenes: 'MYRC','SABI','LIMO','3CAR','OCIM','BPIN','APIN','OMTP'
    !             Sesquiterpenes: 'FARN','BCAR','OSQT'
    !             Methanol: 'MEOH'
    !             Isoprene: 'C5H8'
    !-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kclass      ! compound class
    REAL(lp), INTENT(in)  :: lai, laip   ! LAI of current and previous month
    REAL(lp), INTENT(in)  :: tempd       ! daily average above canopy temperature
    REAL(lp), INTENT(out) :: gamma


!    INTEGER, PARAMETER    :: NCLASS = 6          ! constant,isoprene,monoterpenes,sesquiterpenes,MBO,methanol
!                                                 ! updated from Guenther et al. 2012
!    INTEGER, PARAMETER    :: NAGE = 4            ! new, grown, mature, old

!     ! Fortran 90/95 reference: For 2D matrices, values are listed column wise.
!     REAL(lp), PARAMETER   :: afac(NCLASS, NAGE) = &
!                             RESHAPE(SOURCE = (/  1.0_lp,  0.05_lp,  2.0_lp,   0.4_lp,   0.05_lp,  3.5_lp,     &
!                                                  1.0_lp,  0.6_lp,   1.8_lp,   0.6_lp,   0.6_lp,   3.0_lp,     &
!                                                  1.0_lp,  1.25_lp,  1.0_lp,   1.0_lp,   1.0_lp,   1.0_lp,     &
!                                                  1.0_lp,  1.0_lp,   1.05_lp,  0.95_lp,  0.9_lp,   1.2_lp /),  &
! !                                                  1.0_lp,  1.0_lp,   1.0_lp,   1.0_lp,   1.0_lp,   1.0_lp,     &
! !                                                  1.0_lp,  0.9_lp,   1.05_lp,  0.95_lp,  0.9_lp,   1.2_lp /),  &
!                                     SHAPE = (/NCLASS, NAGE/))
    REAL(lp), PARAMETER   :: dt = 30._lp  ! original time step in MEGAN [days] (set to TSTLEN = 30.)
                                          ! ### this is unclear! See formulae for fnew and fmat below

    REAL(lp)              :: Fnew, Fgro, Fmat, Fold     ! fraction of new, grown, mature and aged leafs
    REAL(lp)              :: Anew, Agro, Amat, Aold

    REAL(lp)              :: dti         ! number of days between budbreak
                                         ! and the induction of emission
    REAL(lp)              :: dtm         ! number of days between budbreak
                                         ! and the initiation of peak
                                         ! emissions rates
    REAL(lp)              :: fac_lai     ! Relation between previous and current lai

    !  This way of assigning parameters is less elegant, but secure...
    IF (kclass == 2) THEN
      Anew = 0.05_lp
      Agro = 0.6_lp
      Amat = 1.125_lp
      Aold = 1.0_lp
    ELSE
      Anew = 1.0_lp
      Agro = 1.0_lp
      Amat = 1.0_lp
      Aold = 1.0_lp
    END IF

    ! TODO: This code may break if either lai or laip is == 0
    ! Good enough fix? Needs to be checked!
    IF (lai > 1.0e-10_lp) THEN
      fac_lai = laip / lai
    ELSE
      fac_lai = 1.0_lp
    ENDIF

    !...  Calculate foliage fraction
    IF (laip == lai) THEN      ! mature state
      Fnew = 0.0_lp
      Fgro = 0.1_lp
      Fmat = 0.8_lp
      Fold = 0.1_lp

    ELSEIF (laip > lai) THEN  ! senescing state
      Fnew = 0.0_lp
      Fgro = 0.0_lp
      Fold = (laip - lai) / laip
      Fmat = 1.0_lp - Fold

    ELSE ! laip < lai
      !  Calculate ti and tm
      IF (tempd <= 303.0_lp) THEN
        dti = 5.0_lp + 0.7_lp*(300._lp - tempd)
      ELSEIF (tempd > 303.0_lp) THEN
        dti = 2.9_lp
      ENDIF
      dtm = 2.3_lp*dti

      !      Calculate Fnew and Fmat, then Fgro and Fold
      ! ### questionable algorithm! This seems to implicitly assume that time step for lai and "t" is the same (= 1 month)!
      ! ### also: if t <= tm (by definition tm>ti) then there will never be grown leaves, because Fnew+Fmat = 1
      !      Fnew
      IF (dt <= dti) THEN
        Fnew = 1.0_lp - fac_lai
      ELSEIF (dt > dti) THEN
        Fnew = (dti / dt) * (1.0_lp - fac_lai)
      ENDIF

      ! Fmat
      IF (dt <= dtm) THEN
        Fmat = fac_lai
      ELSEIF (dt > dtm) THEN
        Fmat = fac_lai + ((dt - dtm) / dt) * (1._lp - fac_lai)
      ENDIF

      Fgro = 1.0_lp - Fnew - Fmat
      Fold = 0.0_lp
    ENDIF

    !...  Calculate GAMMA_A
    gamma = Fnew * Anew + Fgro * Agro + Fmat * Amat + Fold * Aold
!    original way of assigning parameters -- buggy due to bug in matrix above
!    pgamma(ji) = Fnew*afac(kclass,1) + Fgro*afac(kclass,2) + Fmat*afac(kclass,3) + Fold*afac(kclass,4)

  END SUBROUTINE calc_gamma_age

!------------------------------------------------------------------------------------
  PURE SUBROUTINE calc_gamma_sm(soilm, gamma)
  ! --- adjustment for soil moisture
  !  gamma_sm = 1 in standard version of MEGAN v2.10
  !> From Guenther et al., 2012, MEGANv2.1
  !> soil moisture dependence algorithm for isoprene emission:
  !>
  !> GAM_SM = 1              IF SM > SM1
  !> GAM_SM = (SM-WP)/DSM1   IF WP < SM < SM1
  !> GAM_SM = 0              IF SM < WP
  !>
  !> where SM = soil moisture (m3/m3)
  !>       WP = wilting point (m3/m3)
  !>       DSM1 = 0.04 empirical parameter (from Pegoraro et al., 2004; modified in MEGANv2.1)
  !>       SM1 = WP + DSM1
  !----------------------------------------------------------------------------------------------
    REAL(lp), INTENT(in)  :: soilm
    REAL(lp), INTENT(out) :: gamma

    REAL(lp), PARAMETER :: delta1 = 0.06_lp

    REAL(lp) :: wiltp
    REAL(lp) :: soilm1

    wiltp = 0.35_lp
    soilm1 = wiltp + delta1


    IF (soilm > soilm1) THEN
      gamma = 1._lp
    ELSE
      gamma = (soilm - wiltp) / delta1
    ENDIF
    IF (soilm < wiltp) THEN
      gamma = 0._lp
    ENDIF

  END SUBROUTINE calc_gamma_sm

!------------------------------------------------------------------------------------
  PURE SUBROUTINE calc_gamma_co2(co2,gamma)
    USE mo_methane_constants, ONLY: CO2_consts
    USE mo_methane_util,      ONLY: concentration_mmr_2_vmr

    REAL(lp), INTENT(in)  :: co2            ! CO2 mmr
    REAL(lp), INTENT(out) :: gamma

    REAL(lp), PARAMETER :: Ismax =   1.344_lp
    REAL(lp), PARAMETER :: h     =   1.4614_lp
    REAL(lp), PARAMETER :: Cstar = 585.0_lp

    REAL(lp) :: Ci
    REAL(lp) :: CO2_ppm

    CO2_ppm = concentration_mmr_2_vmr(CO2_consts, co2) * 1.e6 ! conversion of mass mixing ratio to volume
                                                              ! mixing ratio and then ppmv
    Ci = 0.7_lp * CO2_ppm
    gamma = Ismax - ((Ismax * Ci**h) / (Cstar**h + Ci**h))

  END SUBROUTINE calc_gamma_co2

  ! Update needed running means of temp, lai and pfd
  PURE SUBROUTINE update_megan_runmean(linit, dt,t_short,t_long,              &
    &                                            temp,temp_relax_short,                 &
    &                                            lai,  lai_relax_short, lai_relax_long, &
    &                                            pfd,  pfd_relax_short, pfd_relax_long)
    USE mo_top_util, ONLY: time_relaxation
    LOGICAL,  INTENT(in)    :: linit            ! Initialize means?
    REAL(lp), INTENT(in)    :: dt               ! Time step for which instantaneous data are assumed valid [s]
    REAL(lp), INTENT(in)    :: t_short          ! E-fold time for short term averages [s]
    REAL(lp), INTENT(in)    :: t_long           ! E-fold time for long  term averages [s]
    REAL(lp), INTENT(in)    :: temp             ! Instantaneous temperature [K]
    REAL(lp), INTENT(in)    :: lai              ! Instantaneous LAI
    REAL(lp), INTENT(in)    :: pfd              ! Instantaneous pfd
    REAL(lp), INTENT(inout) :: temp_relax_short ! Short term mean of temp
    REAL(lp), INTENT(inout) :: lai_relax_short  ! Short term mean of LAI
    REAL(lp), INTENT(inout) :: pfd_relax_short  ! Short term mean of pfd
    REAL(lp), INTENT(inout) :: lai_relax_long   ! Long  term mean of LAI
    REAL(lp), INTENT(inout) :: pfd_relax_long   ! Long  term mean of pfd

    REAL(lp) :: nstep_short, nstep_long ! #time steps for short and long averages respectively (ideally integer numbers)

    nstep_short = t_short / dt
    nstep_long  = t_long  / dt

    CALL time_relaxation(temp,temp_relax_short,nstep=nstep_short, linit=linit)
    CALL time_relaxation(lai , lai_relax_short,nstep=nstep_short, linit=linit)
    CALL time_relaxation(pfd , pfd_relax_short,nstep=nstep_short, linit=linit)
    CALL time_relaxation(lai , lai_relax_long ,nstep=nstep_long,  linit=linit)
    CALL time_relaxation(pfd , pfd_relax_long ,nstep=nstep_long,  linit=linit)


  END SUBROUTINE update_megan_runmean

#endif

END MODULE mo_megan_process
