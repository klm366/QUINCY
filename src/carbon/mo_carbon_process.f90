!> Contains the routines for the carbon processes
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
MODULE mo_carbon_process
#ifndef __NO_JSBACH__

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
  PURE ELEMENTAL SUBROUTINE calc_soil_c( &
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
    REAL(wp), INTENT(in)    :: Cflx_wood_2_litter
    REAL(wp), INTENT(in)    :: Cflx_green_2_litter
    REAL(wp), INTENT(in)    :: Cflx_2_root_exudates

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

    ! locals
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
    REAL(wp), DIMENSION(SIZE(burned_fract))  :: non_woody_litter, fraction ! Helper
    INTEGER                                  :: ind, nc

    nc = SIZE(burned_fract)
    non_woody_litter(:) = c_green(:) + c_reserve(:)

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    cflux_dist_green_2_soil(:) = cflux_dist_green_2_soil(:) + &
      & (c_green(:) + c_reserve(:)) * (1._wp - fract_green_aboveGround) * burned_fract(:) / sec_per_day
    cflux_dist_woods_2_soil(:)  = cflux_dist_woods_2_soil(:) + &
      & c_woods(:) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(:) / sec_per_day

    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere
    co2flux_fire_all_2_atm(:) = ( &
        ! take living parts of plants above ground:
      & non_woody_litter(:) * fract_green_aboveGround                      &
        ! now add the YASSO pools above ground:
      & + c_water_ag1(:) + c_acid_ag1(:) + c_ethanol_ag1(:)                &
      & + c_nonsoluble_ag1(:)                                              &
      & + c_water_ag2(:) + c_acid_ag2(:) + c_ethanol_ag2(:)                &
      & + c_nonsoluble_ag2(:)                                              &
        ! now add the fire affected fraction of wood pools above ground:
      & + c_woods(:) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
      & ) * burned_fract(:) * molarMassCO2_kg / sec_per_day

    ! Lower down density of aboveground litter pools
    fraction(:) = 1._wp - burned_fract(:)
    c_acid_ag1      (:) = fraction(:) * c_acid_ag1       (:)
    c_water_ag1     (:) = fraction(:) * c_water_ag1      (:)
    c_ethanol_ag1   (:) = fraction(:) * c_ethanol_ag1    (:)
    c_nonsoluble_ag1(:) = fraction(:) * c_nonsoluble_ag1 (:)

    c_acid_ag2      (:) = fraction(:) * c_acid_ag2       (:)
    c_water_ag2     (:) = fraction(:) * c_water_ag2      (:)
    c_ethanol_ag2   (:) = fraction(:) * c_ethanol_ag2    (:)
    c_nonsoluble_ag2(:) = fraction(:) * c_nonsoluble_ag2 (:)

    ! Add green and reserve vegetation carbon to belowground green litter pools
    fraction(:) = (1._wp-fract_green_aboveGround) * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg1(ind), c_water_bg1(ind), c_ethanol_bg1(ind), c_nonsoluble_bg1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO

    ! Add the remaings of burned wood:
    ! AG
    fraction(:) = (1._wp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag2(ind), c_water_ag2(ind), c_ethanol_ag2(ind), c_nonsoluble_ag2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp-fract_wood_aboveGround) * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg2(ind), c_water_bg2(ind), c_ethanol_bg2(ind), c_nonsoluble_bg2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO

    ! Lower down the carbon density of living plant pools
    fraction(:) = 1._wp - burned_fract(:)
    c_green(:)   = c_green(:)   * fraction(:)
    c_woods(:)   = c_woods(:)   * fraction(:)
    c_reserve(:) = c_reserve(:) * fraction(:)

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
    REAL(wp), DIMENSION(SIZE(damaged_fract)) :: non_woody_litter, fraction  ! Helper
    INTEGER                                  :: ind, nc

    ! Initializations
    nc = SIZE(damaged_fract)
    non_woody_litter(:) = c_green(:) + c_reserve(:)

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool
    cflux_dist_woods_2_soil(:) = cflux_dist_woods_2_soil(:) + c_woods(:) * damaged_fract(:) / sec_per_day
    cflux_dist_green_2_soil(:) = cflux_dist_green_2_soil(:) + non_woody_litter(:) * damaged_fract(:) / sec_per_day

    ! Add damaged green and reserve
    ! AG
    fraction(:) = fract_green_aboveGround * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag1(ind), c_water_ag1(ind), c_ethanol_ag1(ind), c_nonsoluble_ag1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp - fract_green_aboveGround) * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg1(ind), c_water_bg1(ind), c_ethanol_bg1(ind), c_nonsoluble_bg1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO

    ! Add damaged wood
    ! AG
    fraction(:) = fract_wood_aboveGround  * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag2(ind), c_water_ag2(ind), c_ethanol_ag2(ind), c_nonsoluble_ag2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp - fract_wood_aboveGround)  * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg2(ind), c_water_bg2(ind), c_ethanol_bg2(ind), c_nonsoluble_bg2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO

    ! Lower down the carbon density of living plant pools
    fraction(:) = 1._wp - damaged_fract(:)
    c_green(:)   = c_green(:)   * fraction(:)
    c_woods(:)   = c_woods(:)   * fraction(:)
    c_reserve(:) = c_reserve(:) * fraction(:)

  END SUBROUTINE relocate_carbon_damage

  PURE SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, fract_aboveground, &
                         NPP_2_rootExudates, redFact_Nlimit)
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
    !! distributed as part of the ICON Earth System Model under the license conditions as stated in the header of this module.
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Raisanen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    USE mo_jsb_physical_constants, ONLY: tmelt
    IMPLICIT NONE

    REAL(wp), DIMENSION(9),  INTENT(IN)  :: Yasso_io_pools     ! Yasso pools IN [mol(c)/m2]
    REAL(wp),                INTENT(IN)  :: WoodLitterSize     ! Size of woody litter;
                                                               !   is zero when the subroutine is called for non-woody litter
    REAL(wp), DIMENSION(2),  INTENT(IN)  :: Weather            ! climatic drivers: air temperature and precipitation, 15-d means
                                                               ! Weather(1) in [K] ; Weather(2) in [mm/s]
    REAL(wp),                INTENT(IN)  :: litter             ! fresh litter inputs  (above and belowground) [mol(C)/m2/day]
    REAL(wp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV          ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(wp),                INTENT(IN)  :: fract_aboveground  ! parameter to seperate litter influx into above- and belowground
                                                               !     part                                                   [ ]
    REAL(wp),                INTENT(IN)  :: NPP_2_rootExudates ! root exudates [daily]  [mol(C)/m2/day]

    REAL(wp), DIMENSION(18), INTENT(OUT) :: Yasso_out          ! updated Y pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                               ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more
                                                               ! information

    REAL(wp), OPTIONAL, INTENT(IN)       ::  redFact_Nlimit

    !local
    REAL(wp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool

    ! Yasso pools   [mol(c)/m2]
    REAL(wp)     :: c_acid_ag    ! Aboveground
    REAL(wp)     :: c_water_ag
    REAL(wp)     :: c_ethanol_ag
    REAL(wp)     :: c_nonsoluble_ag
    REAL(wp)     :: c_acid_bg    ! Belowground
    REAL(wp)     :: c_water_bg
    REAL(wp)     :: c_ethanol_bg
    REAL(wp)     :: c_nonsoluble_bg

    REAL(wp)     :: c_humus

    ! Aboveground
    REAL(wp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool         [mol(c)/m2/a]

    REAL(wp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool         [mol(c)/m2/a]

    ! Belowground
    REAL(wp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_humus                ! Loss flux of humus pool               [mol(c)/m2/a]

    REAL(wp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humus                   ! Gain flux of humus pool               [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusAG                 ! Gain flux of humus pool; from ag only [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusBG                 ! Gain flux of humus pool; from bg only [mol(c)/m2/a]

    ! Respiration is output of Yasso therefore a daily flux
    REAL(wp)     :: soilResp_rateYasso             ! Flux to the atmosphere                [mol(c)/m2/d]
    REAL(wp)     :: soilResp_rateLitterAG          ! Flux to the atmosphere; from AWEN ag only  [mol(c)/m2/d]
    REAL(wp)     :: soilResp_rateLitterBG          ! Flux to the atmosphere; from AWEN bg only  [mol(c)/m2/d]

    ! Decomposition rates [1/a]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus

    ! Dcomposition rate for N litter green (needed for N cycle)
    REAL(wp)     :: d_litter_green              ! the rate with which the sum of AWEN pools decays
    REAL(wp)     :: Pseudo_litter_green         ! Sum of all AWEN pools

    ! Yasso parameters (Tuomi '11)
    ! reference decompositation rates for T = 0C and unlimited water availability  in [1/yr]...
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !  ... acid-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !  ... water-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !  ... ethanol-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !  ... nonsoluble-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !  ... humus
    ! parameters for temperature dependence
    REAL(wp), PARAMETER   :: temp_p1                    = 0.095_wp    ! first  order temperature dependence [1/C]
    REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  ! second order temperature dependence [1/C]
    ! parameter for precipitation dependence
    REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    ! first  order precipitation dependence [a/m]
    ! parameters for size dependence
    REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    ! first  order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p2                    = 0.86_wp     ! second order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   ! size dependence power [ ]

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(wp), PARAMETER   :: A_2_W                      =  0.99_wp
    REAL(wp), PARAMETER   :: A_2_E                      =  0.00_wp
    REAL(wp), PARAMETER   :: A_2_N                      =  0.00_wp
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(wp), PARAMETER   :: W_2_A                      = 0.48_wp
    REAL(wp), PARAMETER   :: W_2_E                      = 0.00_wp
    REAL(wp), PARAMETER   :: W_2_N                      = 0.015_wp
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(wp), PARAMETER   :: E_2_A                      = 0.01_wp
    REAL(wp), PARAMETER   :: E_2_W                      = 0.00_wp
    REAL(wp), PARAMETER   :: E_2_N                      = 0.95_wp
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(wp), PARAMETER   :: N_2_A                      = 0.83_wp
    REAL(wp), PARAMETER   :: N_2_W                      = 0.01_wp
    REAL(wp), PARAMETER   :: N_2_E                      = 0.02_wp
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(wp), PARAMETER   :: AWEN_2_H                   = 0.0045_wp

    ! climatic drivers
    REAL(wp)     :: precip                  ! 15 runnning mean of precipitation      [m/a]
    REAL(wp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(wp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(wp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(wp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition

    REAL(wp)     :: redFactor  ! dummy for redFact_Nlimit

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates.
    REAL(wp)     :: dt
    dt = 1.0_wp/days_per_year

    IF(PRESENT(redFact_Nlimit)) THEN ! nitrogen cycle active & we know already the Nlimitation factor
        redFactor  = redFact_Nlimit
    ELSE
        redFactor  = 1.0_wp
    ENDIF

   ! 0. Preparations
    ! initialize yasso pools (AWEN + humus)
    c_acid_ag        = Yasso_io_pools(1)
    c_water_ag       = Yasso_io_pools(2)
    c_ethanol_ag     = Yasso_io_pools(3)
    c_nonsoluble_ag  = Yasso_io_pools(4)

    c_acid_bg        = Yasso_io_pools(5)
    c_water_bg       = Yasso_io_pools(6)
    c_ethanol_bg     = Yasso_io_pools(7)
    c_nonsoluble_bg  = Yasso_io_pools(8)

    c_humus          = Yasso_io_pools(9)

    ! get the sum of the AWEN pools before decomposition
    Pseudo_litter_green = c_acid_ag          &
                        + c_water_ag         &
                        + c_ethanol_ag       &
                        + c_nonsoluble_ag    &
                        + c_acid_bg          &
                        + c_water_bg         &
                        + c_ethanol_bg       &
                        + c_nonsoluble_bg

    ! Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a
    temp   = Weather(1) - Tmelt                ! K -> C

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    ! Temperature dependence of decomposition
    d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_wp)
    ! Precipitation dependence of decomposition
    d_precip = 1.0_wp - EXP(precip_p1*precip)
    ! Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0
    d_size   = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)

     ! decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation
     d_acid       =  redFactor * ref_decomp_rate_acid       * d_temp * d_precip * d_size
     d_water      =  redFactor * ref_decomp_rate_water      * d_temp * d_precip * d_size
     d_ethanol    =  redFactor * ref_decomp_rate_ethanol    * d_temp * d_precip * d_size
     d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble * d_temp * d_precip * d_size

     ! the decomposition of humus is not limited by nutrients; thus no redFactor
     d_humus      =              ref_decomp_rate_humus      * d_temp * d_precip        ! no size effect on humus

    ! 2. Calculate fluxes

    ! loss fluxes (negative values):
    Cflx_from_acid_ag         = c_acid_ag       * d_acid
    Cflx_from_water_ag        = c_water_ag      * d_water
    Cflx_from_ethanol_ag      = c_ethanol_ag    * d_ethanol
    Cflx_from_nonsoluble_ag   = c_nonsoluble_ag * d_nonsoluble

    Cflx_from_acid_bg         = c_acid_bg       * d_acid
    Cflx_from_water_bg        = c_water_bg      * d_water
    Cflx_from_ethanol_bg      = c_ethanol_bg    * d_ethanol
    Cflx_from_nonsoluble_bg   = c_nonsoluble_bg * d_nonsoluble

    Cflx_from_humus           = c_humus         * d_humus

    ! gain fluxes (positive values):
    ! fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
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

    Cflx_2_humusAG        = ABS(Cflx_from_acid_ag                &
                              + Cflx_from_water_ag               &
                              + Cflx_from_ethanol_ag             &
                              + Cflx_from_nonsoluble_ag          &
                              ) * AWEN_2_H

    Cflx_2_humusBG        = ABS(Cflx_from_acid_bg                &
                              + Cflx_from_water_bg               &
                              + Cflx_from_ethanol_bg             &
                              + Cflx_from_nonsoluble_bg          &
                              ) * AWEN_2_H

    Cflx_2_humus          = Cflx_2_humusAG + Cflx_2_humusBG

    ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
    soilResp_rateYasso    =  (Cflx_from_acid_ag + Cflx_from_acid_bg)               &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_ag + Cflx_from_water_bg)           &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         &
                             + Cflx_from_humus

    ! litter ag & bg respiration (needed for N cycle)
    soilResp_rateLitterAG   =  (Cflx_from_acid_ag )                                &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_ag )                               &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_ag )                             &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_ag )                          &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

    soilResp_rateLitterBG   =  (Cflx_from_acid_bg )                                &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_bg )                               &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_bg )                             &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_bg )                          &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

    ! 3. update Yasso pools

    ! 3a. redistribution among yasso pools
    c_acid_ag           = MAX(0.0_wp,                                         &
                             c_acid_ag                                        & ! old pool
                             + (Cflx_from_acid_ag                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_ag) *dt)
    c_water_ag          = MAX(0.0_wp,                                         & ! and so on .....
                             c_water_ag                                       &
                             + (Cflx_from_water_ag                            &
                               + Cflx_2_water_ag) *dt)
    c_ethanol_ag        = MAX(0.0_wp,                                         &
                             c_ethanol_ag                                     &
                             + (Cflx_from_ethanol_ag                          &
                               + Cflx_2_ethanol_ag) *dt)
    c_nonsoluble_ag     = MAX(0.0_wp,                                         &
                             c_nonsoluble_ag                                  &
                             + (Cflx_from_nonsoluble_ag                       &
                               + Cflx_2_nonsoluble_ag) *dt)

    c_acid_bg           = MAX(0.0_wp,                                         &
                             c_acid_bg                                        & ! old pool
                             + (Cflx_from_acid_bg                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_bg) *dt)
    c_water_bg          = MAX(0.0_wp,                                         & ! ...
                             c_water_bg                                       &
                             + (Cflx_from_water_bg                            &
                               + Cflx_2_water_bg) *dt                         &
                             + NPP_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg        = MAX(0.0_wp,                                         &
                             c_ethanol_bg                                     &
                             + (Cflx_from_ethanol_bg                          &
                               + Cflx_2_ethanol_bg) *dt)
    c_nonsoluble_bg     = MAX(0.0_wp,                                         &
                             c_nonsoluble_bg                                  &
                             + (Cflx_from_nonsoluble_bg                       &
                               + Cflx_2_nonsoluble_bg) *dt)

    c_humus          = MAX(0.0_wp,                            &
                             c_humus                          &
                             + (Cflx_from_humus               &
                               + Cflx_2_humus) *dt)

    ! 3b. distribute litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag, c_water_ag, c_ethanol_ag, c_nonsoluble_ag, c_humus, &
      &                          litter, fract_aboveground, Lit_coefV)
    ! bg
    CALL distribute_yasso_litter(c_acid_bg, c_water_bg, c_ethanol_bg, c_nonsoluble_bg, c_humus, &
      &                          litter, 1.0_wp - fract_aboveground, Lit_coefV)

    ! Cflx_litter_2_humus is not only required for distribution, but possibly also as output for the eq script
    Cflx_litter_2_humus = Lit_coefV(i_lctlib_humus) * litter

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green .GT. 0.0_wp) THEN
       d_litter_green = (c_acid_ag            &
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
    Yasso_out(1)           = c_acid_ag
    Yasso_out(2)           = c_water_ag
    Yasso_out(3)           = c_ethanol_ag
    Yasso_out(4)           = c_nonsoluble_ag
    Yasso_out(5)           = c_acid_bg
    Yasso_out(6)           = c_water_bg
    Yasso_out(7)           = c_ethanol_bg
    Yasso_out(8)           = c_nonsoluble_bg
    Yasso_out(9)           = c_humus
    ! 4.2 Write respiration into the output array
    Yasso_out(10)           = soilResp_rateYasso * dt ! convert back to daily

    IF (.NOT.PRESENT(redFact_Nlimit)) THEN
      ! the routine might have been called to diagnose the nutrient demand, so we give out
      ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
      Yasso_out(11)           = Cflx_2_humusAG * dt ! daily ! the flux is positive
      Yasso_out(12)           = Cflx_2_humusBG * dt ! daily ! the flux is positive
      ! 2. flux for mineralisation:  this is the flux from humus to atmosphere (Cflx_slow_2_atmos)
      Yasso_out(13)           =  - Cflx_from_humus * dt ! daily ! sign CHECKED!
      ! ... as well as all fluxes from AWEN pools to atmosphere
      Yasso_out(14)           = -soilResp_rateLitterAG * dt !daily ! the flux is negative
      Yasso_out(15)           = -soilResp_rateLitterBG * dt !daily ! the flux is negative
      ! 3. decomposition rate for Npool litter green
      Yasso_out(16)           = d_litter_green
    ENDIF

    ! 4.3 Write diagnostics to determine equilibrium of humus pools (according to jsb3 svn rev. )
    Yasso_out(17)           = -1._wp * Cflx_from_humus * dt
    Yasso_out(18)           = (Cflx_2_humus * dt) + Cflx_litter_2_humus

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
    REAL(wp), INTENT(IN)                :: litter        !< litter input  (above and belowground) [mol(C)/m2/day]
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
    REAL(wp), INTENT(IN)    :: litter        !< litter input  (above and belowground) [mol(C)/m2/day]
    REAL(wp), INTENT(IN)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(wp), INTENT(IN)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_wp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,                 & ! out
    & Cpool_per_canopy,               & ! in
    & veg_fract_correction,           & ! in
    & fract_fpc_max)                    ! in

    REAL(wp), INTENT(IN), DIMENSION(:) :: &
      & Cpool_per_canopy, &
      & fract_fpc_max,    &
      & veg_fract_correction
    REAL(wp), INTENT(OUT), DIMENSION(:) :: &
      & Cpool_per_tile

    INTEGER :: ic

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, SIZE(fract_fpc_max)
      Cpool_per_tile(ic) = Cpool_per_canopy(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE get_per_tile

