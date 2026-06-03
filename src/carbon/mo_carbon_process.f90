!> Contains the routines for the carbon processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_carbon_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control,      ONLY: acc_stream
  USE mo_kind,             ONLY: wp
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
                               & days_per_year, sec_per_year, sec_per_day,        &
                               & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, &
                               & i_lctlib_nonsoluble, i_lctlib_humus,             &
                               & molarMassCO2_kg
  USE mo_carbon_shared,    ONLY: yasso, distribute_yasso_litter

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_soil_c, relocate_carbon_fire, relocate_carbon_damage

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools.
#ifndef _OPENACC
  PURE ELEMENTAL &
#endif
  SUBROUTINE calc_soil_c( &
    ! General inputs
    & dt,                         &
    ! Input coefficients (fixed)
    & LeafLit_coef_acid,          &
    & LeafLit_coef_water,         &
    & LeafLit_coef_ethanol,       &
    & LeafLit_coef_nonsoluble,    &
    & LeafLit_coef_humus,         &
    & WoodLit_coef_acid,          &
    & WoodLit_coef_water,         &
    & WoodLit_coef_ethanol,       &
    & WoodLit_coef_nonsoluble,    &
    & WoodLit_coef_humus,         &
    & WoodLitterSize,             &
    ! Variable inputs
    & temp,                       &
    & precip,                     &
    & Cflx_wood_2_litter,         &
    & Cflx_green_2_litter,        &
    & Cflx_2_root_exudates,       &
    ! Yasso pools. Used on input and updated on output
    & c_acid_ag1,                 &
    & c_water_ag1,                &
    & c_ethanol_ag1,              &
    & c_nonsoluble_ag1,           &
    & c_acid_bg1,                 &
    & c_water_bg1,                &
    & c_ethanol_bg1,              &
    & c_nonsoluble_bg1,           &
    & c_humus_1,                  &
    & c_acid_ag2,                 &
    & c_water_ag2,                &
    & c_ethanol_ag2,              &
    & c_nonsoluble_ag2,           &
    & c_acid_bg2,                 &
    & c_water_bg2,                &
    & c_ethanol_bg2,              &
    & c_nonsoluble_bg2,           &
    & c_humus_2,                  &
    ! Output fluxes
    & c_decomp_humus_1_sum,       &
    & c_decomp_humus_2_sum,       &
    & c_into_humus_1_sum,         &
    & c_into_humus_2_sum,         &
    & Cflx_soil_respiration       &
    & )

    !$ACC ROUTINE SEQ
    ! Arguments

    REAL(wp), INTENT(in)    :: dt                         ! Time step for carbon processes
    !
    ! Fixed input parameters for Yasso
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: LeafLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_humus         ! fraction going to the humus pool
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: WoodLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_humus         ! fraction going to the humus pool
    REAL(wp), INTENT(in)    :: WoodLitterSize             ! size of coarse debris

    ! Carbon inputs [mol m-2 s-1]:
    REAL(wp), INTENT(inout) :: Cflx_wood_2_litter
    REAL(wp), INTENT(inout) :: Cflx_green_2_litter
    REAL(wp), INTENT(inout) :: Cflx_2_root_exudates

    ! Input meteorology for Yasso
    REAL(wp), INTENT(in)    :: temp              ! multi-day mean temperature [K]
    REAL(wp), INTENT(in)    :: precip            ! multi-mean precipitation [mm s-1]

    ! Yasso pools [mol(C) m-2(canopy)], input and updated on output
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag1        ! Acid soluble litter pool
    REAL(wp), INTENT(inout) :: c_water_ag1       ! Water soluble litter pool
    REAL(wp), INTENT(inout) :: c_ethanol_ag1     ! Ethanol soluble litter pool
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1  ! Non-soluble litter pool
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg1        ! Acid soluble litter pool
    REAL(wp), INTENT(inout) :: c_water_bg1       ! Water soluble litter pool
    REAL(wp), INTENT(inout) :: c_ethanol_bg1     ! Ethanol soluble litter pool
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1  ! Non-soluble litter pool
    REAL(wp), INTENT(inout) :: c_humus_1         ! Humus litter pool
    ! Size class 2; woody litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag2        ! Acid soluble litter pool
    REAL(wp), INTENT(inout) :: c_water_ag2       ! Water soluble litter pool
    REAL(wp), INTENT(inout) :: c_ethanol_ag2     ! Ethanol soluble litter pool
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2  ! Non-soluble litter pool
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg2        ! Acid soluble litter pool
    REAL(wp), INTENT(inout) :: c_water_bg2       ! Water soluble litter pool
    REAL(wp), INTENT(inout) :: c_ethanol_bg2     ! Ethanol soluble litter pool
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2  ! Non-soluble litter pool
    REAL(wp), INTENT(inout) :: c_humus_2         ! Humus litter pool

    ! Decomposistion fluxes
    REAL(wp), INTENT(inout) :: c_decomp_humus_1_sum   ! annual cflux sum of humus decomposed (leaf)
    REAL(wp), INTENT(inout) :: c_decomp_humus_2_sum   ! annual cflux sum of humus decomposed (wood)
    REAL(wp), INTENT(inout) :: c_into_humus_1_sum     ! annual cflux sum into humus (leaf)
    REAL(wp), INTENT(inout) :: c_into_humus_2_sum     ! annual cflux sum into humus (wood)

    REAL(wp), INTENT(out)   :: Cflx_soil_respiration ! Soil (=heterotrophic) respiration rate  [mol(C) m-2 s-1]
                                                     ! Note: this is a loss, therefore negative!
    !REAL(wp), INTENT(out)   :: soilResp_rate_pot     ! soil respiration without N limitation

    ! local variables for Yasso
    REAL(wp), DIMENSION(18) :: Yasso_io_pools
    REAL(wp), DIMENSION(2)  :: Weather
    REAL(wp), DIMENSION(18) :: Yasso_out
    REAL(wp), DIMENSION(5)  :: LeafLit_coefV
    REAL(wp), DIMENSION(5)  :: WoodLit_coefV
    REAL(wp)                :: leafLitter      ! amount of leafLitter
    REAL(wp)                :: sizeZero

    !REAL(wp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot
    !REAL(wp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    !REAL(wp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    !REAL(wp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    !REAL(wp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    !REAL(wp) :: N2O_emissions_litter_pot

    !REAL(wp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    REAL(wp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    !REAL(wp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    !REAL(wp) :: redFact_Nmobile
    REAL(wp) :: redFact_Nlimitation
    !REAL(wp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    !REAL(wp) :: Nflx_crop_harvest_2_SMINN

    !LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated
    !LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used
    !REAL(wp), PARAMETER :: N2O_ef_grazing = 0._wp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_wp

    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    Cflx_soil_respiration     = 0.0_wp
    !soilResp_rate_pot        = 0.0_wp

    ! Preparations
    !-------------------------------------------

    ! Yasso pools:
    Yasso_io_pools(1)  = c_acid_ag1
    Yasso_io_pools(2)  = c_water_ag1
    Yasso_io_pools(3)  = c_ethanol_ag1
    Yasso_io_pools(4)  = c_nonsoluble_ag1
    Yasso_io_pools(5)  = c_acid_bg1
    Yasso_io_pools(6)  = c_water_bg1
    Yasso_io_pools(7)  = c_ethanol_bg1
    Yasso_io_pools(8)  = c_nonsoluble_bg1
    Yasso_io_pools(9)  = c_humus_1
    Yasso_io_pools(10) = c_acid_ag2
    Yasso_io_pools(11) = c_water_ag2
    Yasso_io_pools(12) = c_ethanol_ag2
    Yasso_io_pools(13) = c_nonsoluble_ag2
    Yasso_io_pools(14) = c_acid_bg2
    Yasso_io_pools(15) = c_water_bg2
    Yasso_io_pools(16) = c_ethanol_bg2
    Yasso_io_pools(17) = c_nonsoluble_bg2
    Yasso_io_pools(18) = c_humus_2

    Weather(1) = temp
    Weather(2) = precip

    LeafLit_coefV(i_lctlib_acid)       = LeafLit_coef_acid
    LeafLit_coefV(i_lctlib_water)      = LeafLit_coef_water
    LeafLit_coefV(i_lctlib_ethanol)    = LeafLit_coef_ethanol
    LeafLit_coefV(i_lctlib_nonsoluble) = LeafLit_coef_nonsoluble
    LeafLit_coefV(i_lctlib_humus)      = LeafLit_coef_humus
    WoodLit_coefV(i_lctlib_acid)       = WoodLit_coef_acid
    WoodLit_coefV(i_lctlib_water)      = WoodLit_coef_water
    WoodLit_coefV(i_lctlib_ethanol)    = WoodLit_coef_ethanol
    WoodLit_coefV(i_lctlib_nonsoluble) = WoodLit_coef_nonsoluble
    WoodLit_coefV(i_lctlib_humus)      = WoodLit_coef_humus

    Yasso_out(1:18)     = 0.0_wp           ! array to store the Yasso pools/fluxes
    sizeZero            = 0.0_wp           ! Green litter by definition has size 0
    redFact_Nlimitation = 1.0_wp           ! Presently no nutrient limitation is used

    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    Cflx_soil_respiration = 0._wp


    !! call yasso for leaf litter
    CALL yasso(dt,Yasso_io_pools(1:9), Weather, Cflx_green_2_litter * dt, LeafLit_coefV,                &
      &        sizeZero, Yasso_out, fract_green_aboveGround, Cflx_2_root_exudates * dt, redFact_Nlimitation)

    c_acid_ag1        = Yasso_out(1)
    c_water_ag1       = Yasso_out(2)
    c_ethanol_ag1     = Yasso_out(3)
    c_nonsoluble_ag1  = Yasso_out(4)
    c_acid_bg1        = Yasso_out(5)
    c_water_bg1       = Yasso_out(6)
    c_ethanol_bg1     = Yasso_out(7)
    c_nonsoluble_bg1  = Yasso_out(8)
    c_humus_1         = Yasso_out(9)

    c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out(17)
    c_into_humus_1_sum   = c_into_humus_1_sum + Yasso_out(18)

    ! rescale respiration (from leaf litter decomposition) from "per day" to "per second"
    Cflx_soil_respiration = Yasso_out(10) / dt

    !! call yasso for woody litter
    CALL yasso(dt,Yasso_io_pools(10:18), Weather, Cflx_wood_2_litter * dt, WoodLit_coefV,       &
      & WoodLitterSize, Yasso_out, fract_wood_aboveGround, 0.0_wp, redFact_Nlimitation) ! no exudates

    c_acid_ag2        = Yasso_out(1)
    c_water_ag2       = Yasso_out(2)
    c_ethanol_ag2     = Yasso_out(3)
    c_nonsoluble_ag2  = Yasso_out(4)
    c_acid_bg2        = Yasso_out(5)
    c_water_bg2       = Yasso_out(6)
    c_ethanol_bg2     = Yasso_out(7)
    c_nonsoluble_bg2  = Yasso_out(8)
    c_humus_2         = Yasso_out(9)

    !! Yasso total respiration
    ! add respiration from wood litter decomposition and rescale respiration from "per day" to "per second"
    Cflx_soil_respiration = CFlx_soil_respiration + Yasso_out(10) / dt

    c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out(17)
    c_into_humus_2_sum   = c_into_humus_2_sum + Yasso_out(18)

  END SUBROUTINE calc_soil_c

  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (c_litter_green_ag,
  ! c_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (c_green_st, c_reserve_st, c_woods_st) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire( &
    & burned_fract,                      & ! in from disturbed_frac
    & fire_fract_wood_2_atmos,           & ! in
    & LeafLit_coef,                      & ! in
    & WoodLit_coef,                      & ! in
    & fuel_moisture,                     & ! in
    & c_green,                           & ! inout
    & c_woods,                           & ! inout
    & c_reserve,                         & ! inout
    & cflux_dist_green_2_soil,           & ! inout
    & cflux_dist_woods_2_soil,           & ! inout
    & c_acid_ag1,                        & ! inout
    & c_water_ag1,                       & ! inout
    & c_ethanol_ag1,                     & ! inout
    & c_nonsoluble_ag1,                  & ! inout
    & c_acid_ag2,                        & ! inout
    & c_water_ag2,                       & ! inout
    & c_ethanol_ag2,                     & ! inout
    & c_nonsoluble_ag2,                  & ! inout
    & c_acid_bg1,                        & ! inout
    & c_water_bg1,                       & ! inout
    & c_ethanol_bg1,                     & ! inout
    & c_nonsoluble_bg1,                  & ! inout
    & c_acid_bg2,                        & ! inout
    & c_water_bg2,                       & ! inout
    & c_ethanol_bg2,                     & ! inout
    & c_nonsoluble_bg2,                  & ! inout
    & c_humus_1,                         & ! inout
    & c_humus_2,                         & ! inout
    & co2flux_fire_all_2_atm             & ! out
    & )

    ! Arguments
    REAL(wp), INTENT(in)    :: burned_fract(:)              ! Fraction of the vegetated area of each tile burned till the last
                                                            ! call of this routine
    REAL(wp), INTENT(in)    :: fire_fract_wood_2_atmos      ! Fraction of above ground wood immediately emitted to atmosphere
                                                            ! by fire
    REAL(wp), INTENT(in)    :: fuel_moisture(:)             !
    REAL(wp), INTENT(in)    :: LeafLit_coef(:)              ! Factor to spread non woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(in)    :: WoodLit_coef(:)              ! Factor to spread woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(inout) :: c_green(:)                   ! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods(:)                   ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve(:)                 ! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    ! CO2 immediately released by fire per tile [kg(CO2)m-2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_green_2_soil(:)   ! Amount of carbon relocated by wind and fire damage
                                                            ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   ! Amount of carbon relocated by wind damage and fire
                                                            ! .. to the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction, bfe ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)
      bfe = burned_fract(ic) * MIN(0.945,1._wp-fuel_moisture(ic))

      ! Amount of carbon released from the green and reserve pools entering the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._wp - fract_green_aboveGround) * bfe / sec_per_day
      ! Amount of carbon released from the wood pool entering the woody litter pools
      cflux_dist_woods_2_soil(ic)  = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * bfe / sec_per_day

      ! Determine the amount of carbon released from the wood pool, the living tissue pools, and the
      ! litter pools going to the atmosphere
      co2flux_fire_all_2_atm(ic) = ( &
          ! above ground fraction of living biomass
        & non_woody_litter * fract_green_aboveGround                          &
          ! add the above ground YASSO pools
        & + c_water_ag1(ic) + c_acid_ag1(ic) + c_ethanol_ag1(ic)              &
        & + c_nonsoluble_ag1(ic)                                              &
        & + c_water_ag2(ic) + c_acid_ag2(ic) + c_ethanol_ag2(ic)              &
        & + c_nonsoluble_ag2(ic)                                              &
          ! add the fire affected fraction of the above ground wood pools
        & + c_woods(ic) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
          ! scale with the burned fraction
        & ) * burned_fract(ic) * bfe * molarMassCO2_kg / sec_per_day

      ! Reduce the above ground yasso pools accordingly
      fraction = 1._wp - bfe
      c_acid_ag1      (ic) = fraction * c_acid_ag1       (ic)
      c_water_ag1     (ic) = fraction * c_water_ag1      (ic)
      c_ethanol_ag1   (ic) = fraction * c_ethanol_ag1    (ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1 (ic)

      c_acid_ag2      (ic) = fraction * c_acid_ag2       (ic)
      c_water_ag2     (ic) = fraction * c_water_ag2      (ic)
      c_ethanol_ag2   (ic) = fraction * c_ethanol_ag2    (ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2 (ic)

      ! Add below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._wp-fract_green_aboveGround) * bfe
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remaings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._wp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * bfe
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - bfe
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (c_green_st, c_reserve_st and
  ! c_woods_st) is partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_damage(  &
    & damaged_fract,                  & ! in
    & LeafLit_coef,                   & ! in
    & WoodLit_coef,                   & ! in
    & c_green,                        & ! inout
    & c_woods,                        & ! inout
    & c_reserve,                      & ! inout
    & cflux_dist_green_2_soil,        & ! inout
    & cflux_dist_woods_2_soil,        & ! inout
    & c_acid_ag1,                     & ! inout
    & c_water_ag1,                    & ! inout
    & c_ethanol_ag1,                  & ! inout
    & c_nonsoluble_ag1,               & ! inout
    & c_acid_ag2,                     & ! inout
    & c_water_ag2,                    & ! inout
    & c_ethanol_ag2,                  & ! inout
    & c_nonsoluble_ag2,               & ! inout
    & c_acid_bg1,                     & ! inout
    & c_water_bg1,                    & ! inout
    & c_ethanol_bg1,                  & ! inout
    & c_nonsoluble_bg1,               & ! inout
    & c_acid_bg2,                     & ! inout
    & c_water_bg2,                    & ! inout
    & c_ethanol_bg2,                  & ! inout
    & c_nonsoluble_bg2,               & ! inout
    & c_humus_1,                      & ! inout
    & c_humus_2                       & ! inout
    & )

    ! Arguments
    real(wp),intent(in)    :: damaged_fract(:)         ! Fraction of the vegetated area of each tile damaged since the last
                                                       ! call of this routine
    real(wp),intent(in)    :: LeafLit_coef(:)          ! Factor to spread non woody litter into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)          ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: c_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(wp),intent(inout) :: cflux_dist_green_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
    real(wp),intent(inout) :: cflux_dist_woods_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
                                                                  ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_1(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_2(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc

      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Amount of carbon transfered from damaged wood pool to the wood litter pools
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + c_woods(ic) * damaged_fract(ic) / sec_per_day

      ! Amount of carbon transfered from the damaged non-woody vegetation pools to the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + non_woody_litter * damaged_fract(ic) / sec_per_day

      ! Add the carbon of damaged non-woody vegetation pools to the yasso soil pools
      !  Above ground
      fraction = fract_green_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._wp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Add the carbon of damaged woody vegetation pool to the yasso soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_damage

#endif
END MODULE mo_carbon_process
