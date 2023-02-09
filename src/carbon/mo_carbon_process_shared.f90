!>
!! @brief Contains carbon process routines shared between the standard carbon usage and the methane version
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag -> allow for changing max_LAI
!! Julia Nabel (2021-05-28): start of some refactoring (more advisable)
!! Stiig Wilkenskjeld, MPI-M (2023-02-08): Separated from main process module for clear separation of standard and methane versions

!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_carbon_process_shared
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_carbon_constants, ONLY: days_per_year, sec_per_year,                        &
                               & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol,    &
                               & i_lctlib_nonsoluble, i_lctlib_humus,                &
                               & i_pool_acid_ag, i_pool_water_ag, i_pool_ethanol_ag, &
                               & i_pool_nonsoluble_ag,                               &
                               & i_pool_acid_bg, i_pool_water_bg, i_pool_ethanol_bg, &
                               & i_pool_nonsoluble_bg, i_pool_humus,                 &
                               & i_flux_resp_total, i_flux_2_humus_ag,               &
                               & i_flux_2_humus_bg, i_flux_from_humus_diag,          &
                               & i_flux_resp_litter_ag, i_flux_resp_litter_bg,       &
                               & i_rel_flux_AWEN, i_flux_from_humus, i_flux_litter_2_humus
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: yasso, add_litter_to_yasso_pool, get_per_tile

CONTAINS

  PURE SUBROUTINE yasso (Yasso_in_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, fract_aboveground, &
                         NPP_2_rootExudates, redFact_Nlimit, is_anoxic, is_peat,                                  &
                         A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H)
    !! DSG: 11.03.2013
    !! modified Yasso soil C model using an Euler scheme and readable coding style
    !! AWEN (Acid, Water, Ethanol, Nonsoluble) pools are the yasso pools without the humus pools.
    !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !!
    !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !! litter.
    !!
    !! The routine needs as input 1. Yasso pools
    !!                            2. Weather: precipitation and Temperature
    !!                            3. daily litter flux and NPP_2_rootExudates
    !!                            4. PFT specific parameters
    !!
    !! Even when this subroutine is called every day, the dimensions of the variables refere to years!
    !! output pools: AWEN + humus
    !! The routine must be called twice, first for non-woody litter, then for woody litter
    !!-------------------------------------------------------------------------------------------------------------------------
    !! Yasso - Soil carbon model (Jari Liski)
    !!
    !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is
    !! distributed as part of the MPI Earth System Model under the MPI-M license conditions as stated in the header of this module.
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Raeisaenen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    USE mo_jsb_physical_constants, ONLY: tmelt
    IMPLICIT NONE

    REAL(wp), DIMENSION(9),  INTENT(IN)  :: Yasso_in_pools     ! Yasso pools in [mol(C) m-2]
    REAL(wp),                INTENT(IN)  :: WoodLitterSize     ! Size of woody litter;
                                                               !   is zero when the subroutine is called for non-woody litter
    REAL(wp), DIMENSION(2),  INTENT(IN)  :: Weather            ! climatic drivers: air temperature and precipitation, 15-d means
                                                               ! Weather(1) in [K] ; Weather(2) in [mm s-1]
    REAL(wp),                INTENT(IN)  :: litter             ! fresh litter inputs  (above and belowground) [mol(C) m-2 day-1]
    REAL(wp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV          ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(wp),                INTENT(IN)  :: fract_aboveground  ! parameter to seperate litter influx into above- and belowground
                                                               !     part                                                   [ ]
    REAL(wp),                INTENT(IN)  :: NPP_2_rootExudates ! root exudates [daily]  [mol(C)/m2/day]

    REAL(wp), DIMENSION(18), INTENT(OUT) :: Yasso_out          ! updated Y pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                               ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more
                                                               ! information

    REAL(wp), OPTIONAL,      INTENT(IN)  ::  redFact_Nlimit
    LOGICAL,  OPTIONAL,      INTENT(IN)  ::  is_anoxic
    LOGICAL,  OPTIONAL,      INTENT(IN)  ::  is_peat

    ! Possibilities to overwrite (some of) the default exchange coefficients (AWEN2H for backward compatibility)
    REAL(wp), OPTIONAL,      INTENT(IN)  :: A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H

    ! Yasso parameters (Tuomi '11)
    !   reference decompositation rates for T = 0C and unlimited water availability  in [a-1]...
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !  ... acid-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !  ... water-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !  ... ethanol-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !  ... nonsoluble-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !  ... humus
    !   parameters for temperature dependence
    REAL(wp), PARAMETER   :: temp_p1                    =  0.095_wp   ! first  order temperature dependence [C-1]
    REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  ! second order temperature dependence [C-1]
    !   parameter for precipitation dependence
    REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    ! first  order precipitation dependence [a m-1]
    !   parameters for size dependence
    REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    ! first  order size dependence [cm-1]
    REAL(wp), PARAMETER   :: size_p2                    =  0.86_wp    ! second order size dependence [cm-1]
    REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   ! size dependence power [ ]

    ! locals

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(wp)   :: A_2_W
    REAL(wp)   :: A_2_E
    REAL(wp)   :: A_2_N
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(wp)   :: W_2_A
    REAL(wp)   :: W_2_E
    REAL(wp)   :: W_2_N
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(wp)   :: E_2_A
    REAL(wp)   :: E_2_W
    REAL(wp)   :: E_2_N
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(wp)   :: N_2_A
    REAL(wp)   :: N_2_W
    REAL(wp)   :: N_2_E
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(wp)   :: AWEN_2_H
    ! Alternative pathway for N pool which enters the humus pool [ ]
    REAL(wp)   :: N_2_H
    REAL(wp)   :: AWE_2_H

    ! Yasso pools   [mol(C) m-2]
    !   Aboveground
    REAL(wp)     :: c_acid_ag
    REAL(wp)     :: c_water_ag
    REAL(wp)     :: c_ethanol_ag
    REAL(wp)     :: c_nonsoluble_ag
    !   Belowground
    REAL(wp)     :: c_acid_bg
    REAL(wp)     :: c_water_bg
    REAL(wp)     :: c_ethanol_bg
    REAL(wp)     :: c_nonsoluble_bg
    REAL(wp)     :: c_humus

    ! Carbon fluxes [mol(C) m-2 a-1]
    !   Aboveground
    REAL(wp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool
    REAL(wp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool
    REAL(wp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool
    REAL(wp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool
                                                                                                         
    REAL(wp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool
    REAL(wp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool
    REAL(wp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool
    REAL(wp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool
                                                                                                          
    !   Belowground                                                                                         
    REAL(wp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool
    REAL(wp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool
    REAL(wp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool
    REAL(wp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool
    REAL(wp)     :: Cflx_from_humus                ! Loss flux of humus pool
                                                                                                         
    REAL(wp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool
    REAL(wp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool
    REAL(wp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool
    REAL(wp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool

    !   Humus
    REAL(wp)     :: Cflx_2_humus                   ! Gain flux of humus pool
    REAL(wp)     :: Cflx_2_humusAG                 ! Gain flux of humus pool; from ag only
    REAL(wp)     :: Cflx_2_humusBG                 ! Gain flux of humus pool; from bg only
    REAL(wp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool

    ! Respiration is output of Yasso
    REAL(wp)     :: Cflx_respiration               ! Total flux to the atmosphere
    REAL(wp)     :: Cflx_resp_from_litter_ag       ! Flux to the atmosphere; from AWEN ag only
    REAL(wp)     :: Cflx_resp_from_litter_bg       ! Flux to the atmosphere; from AWEN bg only

    ! Decomposition rates [a-1]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus
    !   For methane relevant processes
    REAL(wp)     :: d_anoxic_ag             ! decomposition rate modifier for ag pools, if conditions are anoxic
    REAL(wp)     :: d_anoxic_bg             ! decomposition rate modifier for bg pools, if conditions are anoxic
    REAL(wp)     :: d_peat_acrotelm         ! decomposition rate modifier for peatland C pools acrotelm
    REAL(wp)     :: d_peat_catotelm         ! decomposition rate modifier for peatland C pools catotelm

    REAL(wp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(wp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(wp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition

    ! Decomposition rate for N litter green (needed for N cycle)
    REAL(wp)     :: d_litter_green          ! the rate with which the sum of AWEN pools decays
    REAL(wp)     :: Pseudo_litter_green     ! Sum of all AWEN pools

    ! climatic drivers
    REAL(wp)     :: precip                  ! 15 runnning mean of precipitation      [m a-1]
    REAL(wp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(wp)     :: redFactor  ! dummy for redFact_Nlimit

    ! anoxic & peat decomposition modifier
    LOGICAL      :: anoxic_flag, peat_flag

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates.
    REAL(wp)     :: years_per_day

    years_per_day = 1.0_wp/days_per_year

    ! Treatment of optional parameters
    !   Set exchange coefficients
    IF (PRESENT(AWEN2H)) THEN
      AWE_2_H = AWEN2H
      N_2_H   = AWEN2H
    ENDIF
    AWE_2_H  = 0.0045_wp; IF (PRESENT(AWE2H))  AWE_2_H  = AWE2H
    N_2_H    = 0.0045_wp; IF (PRESENT(N2H))    N_2_H    = N2H
    A_2_W    = 0.99_wp;   IF (PRESENT(A2W))    A_2_W    = A2W
    A_2_E    = 0.00_wp;   IF (PRESENT(A2E))    A_2_E    = A2E
    A_2_N    = 0.00_wp;   IF (PRESENT(A2N))    A_2_N    = A2N
    W_2_A    = 0.48_wp;   IF (PRESENT(W2A))    W_2_A    = W2A
    W_2_E    = 0.00_wp;   IF (PRESENT(W2E))    W_2_E    = W2E
    W_2_N    = 0.015_wp;  IF (PRESENT(W2N))    W_2_N    = W2N
    E_2_A    = 0.01_wp;   IF (PRESENT(E2A))    E_2_A    = E2A
    E_2_W    = 0.00_wp;   IF (PRESENT(E2W))    E_2_W    = E2W
    E_2_N    = 0.95_wp;   IF (PRESENT(E2N))    E_2_N    = E2N
    N_2_A    = 0.83_wp;   IF (PRESENT(N2A))    N_2_A    = N2A
    N_2_W    = 0.01_wp;   IF (PRESENT(N2W))    N_2_W    = N2W
    N_2_E    = 0.02_wp;   IF (PRESENT(N2E))    N_2_E    = N2E

    anoxic_flag    = .FALSE.; IF (PRESENT(is_anoxic)) anoxic_flag    = is_anoxic
    peat_flag      = .FALSE.; IF (PRESENT(is_peat))   peat_flag      = is_peat
    redFactor      =  1.0_wp; 
    IF (PRESENT(redFact_Nlimit)) THEN
      ! nitrogen cycle active & we know already the Nlimitation factor
      redFactor = redFact_Nlimit
    ENDIF

    ! set up anoxic decomposition rate modifiers
    d_anoxic_ag = 1._wp
    d_anoxic_bg = 1._wp
    IF (anoxic_flag) THEN
      d_anoxic_ag = 0.35_wp ! same as Rita Wania
      d_anoxic_bg = 0.35_wp ! same as Rita Wania
    ENDIF

    ! set up peat decomposition rate modifiers
    d_peat_acrotelm = 1._wp
    d_peat_catotelm = 1._wp
    IF (peat_flag) THEN
      d_peat_acrotelm = 1._wp                   ! no decomposition in litter pools modifier -- we reduce GPP instead
      ! decomposition timescale in LPJ is 7.e-5 -> modifier is applied to (ref_decomp_rate_humus * d_anoxic_bg)
      d_peat_catotelm = 0.125_wp                ! decomposition in catotelm is reduced  -> 0.0016*0.35*0.125 = 7e-5
    ENDIF

    ! 0. Preparations

    !   Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a
    temp   = Weather(1) - Tmelt                ! K -> C

    !   initialize yasso pools (AWEN + humus)
    c_acid_ag        = Yasso_in_pools(i_pool_acid_ag)
    c_water_ag       = Yasso_in_pools(i_pool_water_ag)
    c_ethanol_ag     = Yasso_in_pools(i_pool_ethanol_ag)
    c_nonsoluble_ag  = Yasso_in_pools(i_pool_nonsoluble_ag)

    c_acid_bg        = Yasso_in_pools(i_pool_acid_bg)
    c_water_bg       = Yasso_in_pools(i_pool_water_bg)
    c_ethanol_bg     = Yasso_in_pools(i_pool_ethanol_bg)
    c_nonsoluble_bg  = Yasso_in_pools(i_pool_nonsoluble_bg)

    c_humus          = Yasso_in_pools(i_pool_humus)

    !   store sum of the AWEN pools before decomposition
    Pseudo_litter_green = c_acid_ag          &
                        + c_water_ag         &
                        + c_ethanol_ag       &
                        + c_nonsoluble_ag    &
                        + c_acid_bg          &
                        + c_water_bg         &
                        + c_ethanol_bg       &
                        + c_nonsoluble_bg


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    !   Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0
    d_size   = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)
    !   Temperature dependence of decomposition
    d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_wp)
    !   Precipitation dependence of decomposition (not for peat)
    IF (peat_flag) THEN
      d_precip = 1.0_wp
    ELSE
      d_precip = 1.0_wp - EXP(precip_p1*precip)
    ENDIF

    !   decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation
    d_acid       =  redFactor * ref_decomp_rate_acid       * d_temp * d_precip * d_size * d_peat_acrotelm
    d_water      =  redFactor * ref_decomp_rate_water      * d_temp * d_precip * d_size * d_peat_acrotelm
    d_ethanol    =  redFactor * ref_decomp_rate_ethanol    * d_temp * d_precip * d_size * d_peat_acrotelm
    d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble * d_temp * d_precip * d_size * d_peat_acrotelm

    !   the decomposition of humus is not limited by nutrients; thus no redFactor
    d_humus      =              ref_decomp_rate_humus      * d_temp * d_precip * d_peat_catotelm ! no size effect on humus

    ! 2. Calculate fluxes

    !   loss fluxes (negative values):
    Cflx_from_acid_ag         = c_acid_ag       * d_acid       * d_anoxic_ag
    Cflx_from_water_ag        = c_water_ag      * d_water      * d_anoxic_ag
    Cflx_from_ethanol_ag      = c_ethanol_ag    * d_ethanol    * d_anoxic_ag
    Cflx_from_nonsoluble_ag   = c_nonsoluble_ag * d_nonsoluble * d_anoxic_ag

    Cflx_from_acid_bg         = c_acid_bg       * d_acid       * d_anoxic_bg
    Cflx_from_water_bg        = c_water_bg      * d_water      * d_anoxic_bg
    Cflx_from_ethanol_bg      = c_ethanol_bg    * d_ethanol    * d_anoxic_bg
    Cflx_from_nonsoluble_bg   = c_nonsoluble_bg * d_nonsoluble * d_anoxic_bg

    Cflx_from_humus           = c_humus         * d_humus      * d_anoxic_bg

    !   gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_2_acid_ag           = ABS(Cflx_from_water_ag      * W_2_A  &    ! returns positive fluxes
                                 + Cflx_from_ethanol_ag    * E_2_A  &
                                 + Cflx_from_nonsoluble_ag * N_2_A)

    Cflx_2_water_ag          = ABS(Cflx_from_acid_ag       * A_2_W  &
                                 + Cflx_from_ethanol_ag    * E_2_W  &
                                 + Cflx_from_nonsoluble_ag * N_2_W)

    Cflx_2_ethanol_ag        = ABS(Cflx_from_acid_ag       * A_2_E  &
                                 + Cflx_from_water_ag      * W_2_E  &
                                 + Cflx_from_nonsoluble_ag * N_2_E)

    Cflx_2_nonsoluble_ag     = ABS(Cflx_from_acid_ag       * A_2_N  &
                                 + Cflx_from_water_ag      * W_2_N  &
                                 + Cflx_from_ethanol_ag    * E_2_N)

    Cflx_2_acid_bg           = ABS(Cflx_from_water_bg      * W_2_A  &
                                 + Cflx_from_ethanol_bg    * E_2_A  &
                                 + Cflx_from_nonsoluble_bg * N_2_A)

    Cflx_2_water_bg          = ABS(Cflx_from_acid_bg       * A_2_W  &
                                 + Cflx_from_ethanol_bg    * E_2_W  &
                                 + Cflx_from_nonsoluble_bg * N_2_W)

    Cflx_2_ethanol_bg        = ABS(Cflx_from_acid_bg       * A_2_E  &
                                 + Cflx_from_water_bg      * W_2_E  &
                                 + Cflx_from_nonsoluble_bg * N_2_E)

    Cflx_2_nonsoluble_bg     = ABS(Cflx_from_acid_bg       * A_2_N  &
                                 + Cflx_from_water_bg      * W_2_N  &
                                 + Cflx_from_ethanol_bg    * E_2_N)

    Cflx_2_humusAG           = ABS(Cflx_from_acid_ag                &
                                 + Cflx_from_water_ag               &
                                 + Cflx_from_ethanol_ag             &
                                 ) * AWE_2_H                        & 
                                 +                                  &
                               ABS(Cflx_from_nonsoluble_ag)         &
                                   * N_2_H
                             
    Cflx_2_humusBG           = ABS(Cflx_from_acid_bg                &
                                 + Cflx_from_water_bg               &
                                 + Cflx_from_ethanol_bg             &
                                 ) * AWE_2_H                        &
                                 +                                  &
                               ABS(Cflx_from_nonsoluble_bg)         &
                                   * N_2_H
                             
    Cflx_2_humus             = Cflx_2_humusAG + Cflx_2_humusBG

    !   the remaining fractions of the loss fluxes enter the atmosphere as respiration
!    Cflx_respiration        =  (Cflx_from_acid_ag + Cflx_from_acid_bg)             &
!                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H)          &
!                             + (Cflx_from_water_ag + Cflx_from_water_bg)           &
!                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H)          &
!                             + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &
!                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H)          &
!                             + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &
!                             * (1.0_wp - N_2_A - N_2_W - N_2_E - N_2_H)            &
!                             + Cflx_from_humus

    !   litter ag & bg respiration (needed for N cycle)
    Cflx_resp_from_litter_ag =  (Cflx_from_acid_ag )                          &
                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H)    &
                              + (Cflx_from_water_ag )                         &
                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H)    &
                              + (Cflx_from_ethanol_ag )                       &
                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H)    &
                              + (Cflx_from_nonsoluble_ag )                    &
                              * (1.0_wp - N_2_A - N_2_W - N_2_E - N_2_H)

    Cflx_resp_from_litter_bg =  (Cflx_from_acid_bg )                          &
                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H)    &
                              + (Cflx_from_water_bg )                         &
                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H)    &
                              + (Cflx_from_ethanol_bg )                       &
                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H)    &
                              + (Cflx_from_nonsoluble_bg )                    &
                              * (1.0_wp - N_2_A - N_2_W - N_2_E - N_2_H)

    Cflx_respiration         =  Cflx_resp_from_litter_ag                      &
                              + Cflx_resp_from_litter_bg                      &
                              + Cflx_from_humus

    ! 3. update Yasso pools

    !   3a. redistribution among yasso pools
    c_acid_ag           = MAX(0.0_wp,                                         &
                             c_acid_ag                                        & ! old pool
                             + (Cflx_from_acid_ag                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_ag) * years_per_day)
    c_water_ag          = MAX(0.0_wp,                                         & ! and so on .....
                             c_water_ag                                       &
                             + (Cflx_from_water_ag                            &
                               + Cflx_2_water_ag) * years_per_day)
    c_ethanol_ag        = MAX(0.0_wp,                                         &
                             c_ethanol_ag                                     &
                             + (Cflx_from_ethanol_ag                          &
                               + Cflx_2_ethanol_ag) * years_per_day)
    c_nonsoluble_ag     = MAX(0.0_wp,                                         &
                             c_nonsoluble_ag                                  &
                             + (Cflx_from_nonsoluble_ag                       &
                               + Cflx_2_nonsoluble_ag) * years_per_day)

    c_acid_bg           = MAX(0.0_wp,                                         &
                             c_acid_bg                                        & ! old pool
                             + (Cflx_from_acid_bg                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_bg) * years_per_day)
    c_water_bg          = MAX(0.0_wp,                                         & ! ...
                             c_water_bg                                       &
                             + (Cflx_from_water_bg                            &
                               + Cflx_2_water_bg) * years_per_day             &
                             + NPP_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg        = MAX(0.0_wp,                                         &
                             c_ethanol_bg                                     &
                             + (Cflx_from_ethanol_bg                          &
                               + Cflx_2_ethanol_bg) * years_per_day)
    c_nonsoluble_bg     = MAX(0.0_wp,                                         &
                             c_nonsoluble_bg                                  &
                             + (Cflx_from_nonsoluble_bg                       &
                               + Cflx_2_nonsoluble_bg) * years_per_day)

    c_humus          = MAX(0.0_wp,                            &
                             c_humus                          &
                             + (Cflx_from_humus               &
                               + Cflx_2_humus) * years_per_day)

    !   3b. distribute litter
    !    ag
    CALL distribute_yasso_litter(c_acid_ag, c_water_ag, c_ethanol_ag, c_nonsoluble_ag, c_humus, &
      &                          litter, fract_aboveground, Lit_coefV)
    !    bg
    CALL distribute_yasso_litter(c_acid_bg, c_water_bg, c_ethanol_bg, c_nonsoluble_bg, c_humus, &
      &                          litter, 1.0_wp - fract_aboveground, Lit_coefV)

    ! Cflx_litter_2_humus is not only required for distribution, but possibly also as output for the eq script
    Cflx_litter_2_humus = Lit_coefV(i_lctlib_humus) * litter * days_per_year ! in [a-1] as all other Cflx-variables (litter [d-1])

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green > 0.0_wp) THEN
       d_litter_green = (  c_acid_ag          &
                         + c_water_ag         &
                         + c_ethanol_ag       &
                         + c_nonsoluble_ag    &
                         + c_acid_bg          &
                         + c_water_bg         &
                         + c_ethanol_bg       &
                         + c_nonsoluble_bg)   &
                         / Pseudo_litter_green - 1._wp
    ELSE
      d_litter_green = 0.0_wp
    ENDIF

    ! 4.1 Write pools into the output array
    Yasso_out(i_pool_acid_ag)       = c_acid_ag
    Yasso_out(i_pool_water_ag)      = c_water_ag
    Yasso_out(i_pool_ethanol_ag)    = c_ethanol_ag
    Yasso_out(i_pool_nonsoluble_ag) = c_nonsoluble_ag
    Yasso_out(i_pool_acid_bg)       = c_acid_bg
    Yasso_out(i_pool_water_bg)      = c_water_bg
    Yasso_out(i_pool_ethanol_bg)    = c_ethanol_bg
    Yasso_out(i_pool_nonsoluble_bg) = c_nonsoluble_bg
    Yasso_out(i_pool_humus)         = c_humus

    ! 4.2 Write respiration into the output array
    Yasso_out(i_flux_resp_total)    = Cflx_respiration * years_per_day ! convert back to daily

    IF (.NOT. PRESENT(redFact_Nlimit)) THEN
      ! the routine might have been called to diagnose the nutrient demand, so we give out
      ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
      Yasso_out(i_flux_2_humus_ag)      =  Cflx_2_humusAG * years_per_day ! daily ! the flux is positive
      Yasso_out(i_flux_2_humus_bg)      =  Cflx_2_humusBG * years_per_day ! daily ! the flux is positive
      ! 2. flux for mineralisation:  this is the flux from humus to atmosphere (Cflx_slow_2_atmos)
      Yasso_out(i_flux_from_humus_diag) = -Cflx_from_humus * years_per_day ! daily ! sign CHECKED!
      ! ... as well as all fluxes from AWEN pools to atmosphere
      Yasso_out(i_flux_resp_litter_ag)  = -Cflx_resp_from_litter_ag * years_per_day !daily ! the flux is negative
      Yasso_out(i_flux_resp_litter_bg)  = -Cflx_resp_from_litter_bg * years_per_day !daily ! the flux is negative
      ! 3. decomposition rate for Npool litter green
      Yasso_out(i_rel_flux_AWEN)        =  d_litter_green
    ENDIF

    ! 4.3 Write diagnostics to determine equilibrium of humus pools
    Yasso_out(i_flux_from_humus)     = -Cflx_from_humus * years_per_day
    Yasso_out(i_flux_litter_2_humus) = (Cflx_2_humus + Cflx_litter_2_humus) * years_per_day

  END SUBROUTINE yasso

  ! ====================================================================================================== !
  !
  !> Distributes incomming litter to yasso pools (i.e. proportional to litter coefficients and passed factors)
  !
  PURE SUBROUTINE distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
    &                          litter, fraction, coefficients)
    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT)             :: c_acid        !< yasso acid pool
    REAL(wp), INTENT(INOUT)             :: c_water       !< yasso water pool
    REAL(wp), INTENT(INOUT)             :: c_ethanol     !< yasso ethanol pool
    REAL(wp), INTENT(INOUT)             :: c_nonsoluble  !< yasso nonsoluble pool
    REAL(wp), INTENT(INOUT)             :: c_humus       !< yasso humus pool
    REAL(wp), INTENT(IN)                :: litter        !< litter input (above and belowground) [mol(C) m-2 day-1]
    REAL(wp), INTENT(IN)                :: fraction      !< fraction of litter put into yasso pools []
    REAL(wp), DIMENSION(5), INTENT(IN)  :: coefficients  !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !

    CALL add_litter_to_yasso_pool(c_acid,       litter, fraction, coefficients(i_lctlib_acid))
    CALL add_litter_to_yasso_pool(c_water,      litter, fraction, coefficients(i_lctlib_water))
    CALL add_litter_to_yasso_pool(c_ethanol,    litter, fraction, coefficients(i_lctlib_ethanol))
    CALL add_litter_to_yasso_pool(c_nonsoluble, litter, fraction, coefficients(i_lctlib_nonsoluble))
    CALL add_litter_to_yasso_pool(c_humus,      litter, fraction, coefficients(i_lctlib_humus))

  END SUBROUTINE distribute_yasso_litter

  ! ====================================================================================================== !
  !
  !> Adds incomming litter share (i.e. proportional to litter coefficient and factor)
  !
  ELEMENTAL SUBROUTINE add_litter_to_yasso_pool(this_pool, litter, fraction, coefficient)
    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT) :: this_pool     !< target yasso pool
    REAL(wp), INTENT(IN)    :: litter        !< litter input  (above and belowground) [mol(C) m-2 day-1]
    REAL(wp), INTENT(IN)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(wp), INTENT(IN)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_wp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  ELEMENTAL SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,                 & ! out
    & Cpool_per_canopy,               & ! in
    & veg_fract_correction,           & ! in
    & fract_fpc_max)                    ! in

    !$claw model-data
    REAL(wp), INTENT(IN) ::  fract_fpc_max, veg_fract_correction
    REAL(wp), INTENT(IN) ::  Cpool_per_canopy
    REAL(wp), INTENT(OUT)::  Cpool_per_tile
    !$claw end model-data

    !$claw sca

    Cpool_per_tile = Cpool_per_canopy * veg_fract_correction * fract_fpc_max

  END SUBROUTINE get_per_tile

#endif
END MODULE mo_carbon_process_shared
