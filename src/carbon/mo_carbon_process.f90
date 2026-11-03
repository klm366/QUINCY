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
  USE mo_kind,             ONLY: lp
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround,    &
    & Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS, ABOVE_GROUND, BELOW_GROUND, &
    & BELOW_GROUND_SL, molarMassCO2_kg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_plant_c, calc_soil_c, relocate_carbon_fire, relocate_carbon_damage, &
    & get_per_tile, add_litter_to_yasso_pool, sum_up_carbon_layers

  INTERFACE calc_soil_c
    MODULE PROCEDURE calc_soil_c_single_layer
    MODULE PROCEDURE calc_soil_c_multi_layer
  END INTERFACE calc_soil_c

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  !
  !> Allocate NPP to plant carbon pools and calculate respiration and litter fluxes.
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_plant_c(        &
    & delta_time_C,               &
    & LAI,                        &
    & LAI_previous,               &
    & max_LAI,                    &
    & cflx_npp_pot_yday_mean,     &
    & fract_npp_2_woodPool,       &
    & fract_npp_2_reservePool,    &
    & fract_npp_2_exudates,       &
    & fract_green_2_herbivory,    &
    & tau_c_woods,                &
    & LAI_shed_constant,          &
    & c_woods_max,                &
    & specific_leaf_area_C,       &
    & reserveC2leafC,             &
    & pheno_type,                 &
    & c_green,                    &
    & c_woods,                    &
    & c_reserve,                  &
    & c_crop_harvest,             &
    & cflx_npp_correction,        &
    & cflx_npp_excess,            &
    & cflx_npp_actual,            &
    & cflx_to_root_exudates,      &
    & cflx_vegetation_to_litter,  &
    & cflx_green_to_litter,       &
    & cflx_wood_to_litter,        &
    & cflx_herbivory,             &
    & cflx_herbivory_to_litter,   &
    & cflx_herbivory_to_atmosphere, &
    & cflx_crop_harvest_to_atmosphere &
    & )

    !$ACC ROUTINE SEQ

    !! Carbon fluxes use the cflx_ prefix and are rates [mol(C) m-2 s-1].
    !! Carbon pools use the c_ prefix and are amounts [mol(C) m-2].
    !! Intermediate cflx_*_step variables are transfers integrated over delta_time_C.
    !! All area-based quantities are per m2 of canopy area.

    ! Input arguments
    REAL(lp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: LAI                     !< Yesterday's mean LAI
    REAL(lp), INTENT(in)    :: LAI_previous            !< Mean LAI from the day before yesterday
    REAL(lp), INTENT(in)    :: max_LAI                 !< Maximum value of LAI
    REAL(lp), INTENT(in)    :: cflx_npp_pot_yday_mean  !< Yesterday's mean potential NPP rate [mol(C) m-2 s-1]
    REAL(lp), INTENT(in)    :: fract_npp_2_woodPool    !< Fraction of NPP allocated to the wood pool
    REAL(lp), INTENT(in)    :: fract_npp_2_reservePool !< Fraction of NPP allocated to the reserve pool
    REAL(lp), INTENT(in)    :: fract_npp_2_exudates    !< Fraction of NPP allocated to root exudates
    REAL(lp), INTENT(in)    :: fract_green_2_herbivory !< Fraction of the green pool removed by herbivory
    REAL(lp), INTENT(in)    :: tau_c_woods             !< Wood-pool turnover time, converted from years to days [d]
    REAL(lp), INTENT(in)    :: LAI_shed_constant       !< Leaf shedding at a constant rate for evergreens etc. [d-1]
    REAL(lp), INTENT(in)    :: c_woods_max             !< Maximum wood-pool carbon content (from lctlib) [mol(C) m-2]
    REAL(lp), INTENT(in)    :: specific_leaf_area_C    !< Specific leaf area (from lctLib) [m2(leaf) mol(C)-1]
    REAL(lp), INTENT(in)    :: reserveC2leafC          !< Ratio of maximum reserve-pool carbon to maximum leaf carbon
    INTEGER,  INTENT(in)    :: pheno_type              !< Phenology type (checks for crop == 5)
    ! Plant carbon pools
    REAL(lp), INTENT(inout) :: c_green                 !< Green carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_woods                 !< Wood carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_reserve               !< Reserve carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_crop_harvest          !< Crop-harvest carbon pool [mol(C) m-2]
    ! Output arguments
    REAL(lp), INTENT(out)   :: cflx_npp_correction             !< NPP correction rate [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_npp_excess                 !< NPP rate not stored in plant pools [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_npp_actual                 !< Actual NPP rate after excess carbon is removed [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_to_root_exudates           !< Carbon flux into root exudates [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_vegetation_to_litter       !< Vegetation-to-litter carbon flux [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_green_to_litter            !< Green-litter carbon flux [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_wood_to_litter             !< Wood-litter carbon flux [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_herbivory                  !< Carbon removed from the green pool by herbivory [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_herbivory_to_litter        !< Herbivory carbon flux to green litter [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_herbivory_to_atmosphere    !< Signed herbivory carbon flux to the atmosphere [mol(C) m-2 s-1]
    REAL(lp), INTENT(out)   :: cflx_crop_harvest_to_atmosphere !< Signed crop-harvest decay flux to atmosphere [mol(C) m-2 s-1]

    ! Parameters
    REAL(lp), PARAMETER :: days_per_year        = 365.0_lp !< Number of days used to convert annual turnover times to days
    REAL(lp), PARAMETER :: greenC2leafC         = 4.0_lp   !< ratio of carbon green pool (leaves, fine roots, starches,
                                                           !<   sugars) to leaf carbon
    REAL(lp), PARAMETER :: fract_c_crop_harvest = 0.5_lp   !< Fraction of eligible litter carbon routed to crop harvest
    REAL(lp), PARAMETER :: fract_c_faeces2_LG   = 0.3_lp   !< Fraction of faecal carbon returned to the green-litter pool
    REAL(lp), PARAMETER :: tau_c_reserve        = 1.0_lp * days_per_year !< Annual reserve-pool turnover time [d]
    REAL(lp), PARAMETER :: tau_c_crop_harvest   = 1.0_lp * days_per_year !< Annual crop-harvest-pool decay time [d]

    LOGICAL  :: is_crop                 ! logical crop mask

    REAL(lp) :: c_green_max, c_reserve_optimal
    REAL(lp) :: cflx_npp_to_green_step, cflx_npp_to_wood_step, cflx_npp_to_reserve_step
    REAL(lp) :: cflx_npp_to_root_exudates_step, cflx_green_to_herbivores_step
    REAL(lp) :: cflx_green_to_litter_step, cflx_wood_to_litter_step, cflx_reserve_to_litter_step
    REAL(lp) :: c_excess
    REAL(lp) :: leaf_shedding

    REAL(lp) :: c_green_pot, c_woods_pot         !! Candidate pool values during carbon allocation
    REAL(lp) :: c_reserve_pot

    REAL(lp) :: cflx_to_crop_harvest_step
    REAL(lp) :: cflx_faeces_to_litter_step, cflx_faeces_to_atmosphere_step

    IF (pheno_type == 5) THEN
      is_crop = .TRUE.
    ELSE
      is_crop = .FALSE.
    END IF

    ! Initializations
    cflx_to_crop_harvest_step = 0.0_lp

    c_green_max       = greenC2leafC * LAI / specific_leaf_area_C
    c_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

    !! ============== CARBON ALLOCATION AND POOL TURNOVER ================================================================
    !!
    !! Apply leaf and wood shedding, reserve turnover, herbivory, and NPP allocation.
    !! ----------------------------------------------------------------------------
    !! Leaf shedding transfers carbon from the green pool to green litter.

    leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI)  !! Apply the minimum constant rate, or the larger LAI decrease.
    cflx_green_to_litter_step =  &                                   !! Transfer no more than the available green-pool carbon
      & MIN(greenC2leafC * leaf_shedding / specific_leaf_area_C, &   !!   according to leaf shedding and specific leaf area.
      & c_green)
    c_green = c_green - cflx_green_to_litter_step
    c_excess = MAX(0.0_lp,c_green - c_green_max)                     !! If the green pool exceeds its LAI-dependent maximum,
    cflx_green_to_litter_step = cflx_green_to_litter_step + c_excess !!   shed the excess to green litter.
    c_green = c_green - c_excess

    IF (is_crop) THEN
      !! Route part of above-ground green litter to the crop-harvest pool.
      cflx_to_crop_harvest_step = &
        & fract_green_aboveGround * fract_c_crop_harvest * cflx_green_to_litter_step
    END IF

    !! ---------------------------------------------------------------------------------
    !! Herbivory removes a fraction of the remaining green pool.

    cflx_green_to_herbivores_step = c_green * fract_green_2_herbivory

    c_green = c_green - cflx_green_to_herbivores_step
    cflx_faeces_to_litter_step = cflx_green_to_herbivores_step * fract_c_faeces2_LG
    cflx_faeces_to_atmosphere_step = cflx_green_to_herbivores_step * (1.0_lp - fract_c_faeces2_LG)
    !! ---------------------------------------------------------------------------------
    !! Wood-pool turnover transfers carbon to the wood-litter pools.
    cflx_wood_to_litter_step = c_woods / tau_c_woods
    c_woods = c_woods - cflx_wood_to_litter_step

    !! Reserve-pool turnover transfers carbon to green litter.

    cflx_reserve_to_litter_step = c_reserve / tau_c_reserve
    IF (is_crop) THEN
      cflx_to_crop_harvest_step = cflx_to_crop_harvest_step + &
        & fract_green_aboveGround * fract_C_crop_harvest * cflx_reserve_to_litter_step
    END IF
    c_reserve      = c_reserve      - cflx_reserve_to_litter_step
    c_crop_harvest = c_crop_harvest + cflx_to_crop_harvest_step

    !!
    !! ============== DETERMINE POTENTIAL PLANT-CARBON ALLOCATION ========================================================

    c_green_pot   = c_green
    c_woods_pot   = c_woods
    c_reserve_pot = c_reserve

    !! Determine carbon allocation.

    cflx_npp_correction = 0.0_lp
    IF (cflx_npp_pot_yday_mean <= 0.0_lp) THEN
      cflx_npp_excess = 0.0_lp

      cflx_npp_to_reserve_step = cflx_npp_pot_yday_mean * delta_time_C
      c_reserve = c_reserve + cflx_npp_to_reserve_step
      IF (c_reserve < 0.0_lp) THEN
        cflx_npp_to_reserve_step = cflx_npp_to_reserve_step - c_reserve
        cflx_npp_correction = -c_reserve / delta_time_C
        c_reserve = 0.0_lp
      END IF

      !! No growth allocation is applied for zero or negative NPP.
      cflx_npp_to_reserve_step = 0.0_lp
      cflx_npp_to_green_step   = 0.0_lp
      cflx_npp_to_wood_step    = 0.0_lp
      cflx_to_root_exudates = 0.0_lp
      c_reserve_pot        = c_reserve

    ELSE

      !! Allocate NPP according to the configured fractions.
      cflx_npp_to_wood_step     = delta_time_C * fract_npp_2_woodPool    * cflx_npp_pot_yday_mean
      cflx_npp_to_reserve_step  = delta_time_C * fract_npp_2_reservePool * cflx_npp_pot_yday_mean
      cflx_npp_to_root_exudates_step = delta_time_C * fract_npp_2_exudates * cflx_npp_pot_yday_mean
      cflx_npp_to_green_step    = delta_time_C * &
        & (1.0_lp-fract_npp_2_woodPool-fract_npp_2_reservePool-fract_npp_2_exudates) * cflx_npp_pot_yday_mean

      !! Cap wood growth at the configured maximum and redirect excess to green growth.

      c_woods_pot       = c_woods_pot + cflx_npp_to_wood_step
      c_excess          = MAX(0.0_lp,c_woods_pot-c_woods_max)
      IF (c_excess > 0._lp) THEN
        cflx_npp_to_green_step = cflx_npp_to_green_step + c_excess
        cflx_npp_to_wood_step  = cflx_npp_to_wood_step  - c_excess
      END IF

      !! Fill the reserve pool up to its target; redirect any excess to green growth.
      IF (c_reserve_pot < c_reserve_optimal) THEN
        c_reserve_pot     = c_reserve_pot     + cflx_npp_to_reserve_step
        c_excess          =                                         &
          & MAX(0.0_lp,c_reserve_pot - c_reserve_optimal)
        IF (c_excess > 0._lp) THEN
          c_reserve_pot       = c_reserve_pot       - c_excess
          cflx_npp_to_green_step   = cflx_npp_to_green_step   + c_excess
          cflx_npp_to_reserve_step = cflx_npp_to_reserve_step - c_excess
        END IF
      ELSE
        cflx_npp_to_green_step   = cflx_npp_to_green_step + cflx_npp_to_reserve_step
        cflx_npp_to_reserve_step = 0.0_lp
      END IF

      !! Cap green growth at the LAI-dependent maximum and use available reserve capacity.

      c_green_pot       = c_green_pot + cflx_npp_to_green_step
      c_excess          = MAX(0.0_lp,c_green_pot - c_green_max)
      IF (c_excess > 0._lp) THEN
        cflx_npp_to_green_step = cflx_npp_to_green_step - c_excess
        c_green_pot       = c_green_pot - c_excess
      END IF
      IF (c_reserve_pot < c_reserve_optimal) THEN
        IF (c_excess > 0._lp) THEN
          c_reserve_pot       = c_reserve_pot       + c_excess
          cflx_npp_to_reserve_step = cflx_npp_to_reserve_step + c_excess
        END IF
        c_excess = MAX(0.0_lp,c_reserve_pot - c_reserve_optimal)
        IF (c_excess > 0._lp) THEN
          c_reserve_pot       = c_reserve_pot       - c_excess
          cflx_npp_to_reserve_step = cflx_npp_to_reserve_step - c_excess
        END IF
      END IF

      cflx_npp_excess       = c_excess                       / delta_time_C
      cflx_to_root_exudates = cflx_npp_to_root_exudates_step / delta_time_C

    END IF

    !! ============== APPLY NPP ALLOCATION ====================================================================

    !! ============== UPDATE PLANT CARBON POOLS ===============================================================

    c_green   = c_green   + cflx_npp_to_green_step
    c_woods   = c_woods   + cflx_npp_to_wood_step
    c_reserve = c_reserve + cflx_npp_to_reserve_step

    !! ============== COMPUTE ACTUAL NPP RATE =================================================================

    IF (cflx_npp_pot_yday_mean > 0.0_lp) THEN
      cflx_npp_actual = (cflx_npp_to_root_exudates_step + &
        & (cflx_npp_to_green_step + cflx_npp_to_wood_step) + cflx_npp_to_reserve_step) / delta_time_C
    ELSE
      cflx_npp_actual = cflx_npp_pot_yday_mean + cflx_npp_correction
    END IF

    !! ============== CROP-HARVEST POOL DECAY =================================================================

    cflx_crop_harvest_to_atmosphere = c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest = c_crop_harvest - cflx_crop_harvest_to_atmosphere

    cflx_crop_harvest_to_atmosphere = -cflx_crop_harvest_to_atmosphere / delta_time_C

    !! ============== HERBIVORY AND LITTER FLUX DIAGNOSTICS ====================================================

    cflx_herbivory_to_atmosphere = -(cflx_faeces_to_atmosphere_step) / delta_time_C

    cflx_vegetation_to_litter = &
      & (cflx_green_to_litter_step + cflx_wood_to_litter_step + cflx_reserve_to_litter_step) / delta_time_C
    cflx_green_to_litter = cflx_vegetation_to_litter + &
      & (-cflx_wood_to_litter_step+cflx_faeces_to_litter_step-cflx_to_crop_harvest_step) / delta_time_C
    cflx_wood_to_litter = cflx_wood_to_litter_step / delta_time_C
    cflx_herbivory = cflx_green_to_herbivores_step / delta_time_C
    cflx_herbivory_to_litter = cflx_faeces_to_litter_step / delta_time_C

  END SUBROUTINE calc_plant_c

  !
  !> This routine handles the impact of fire on the carbon pools.
  !>
  !> We assume that the above ground litter carbon (yasso ag1 and ag2 pools), is completely released
  !> to the atmosphere while the carbon of the living plants (c_green, c_reserve, c_woods) is partly
  !> released to the atmosphere and partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire(       &
    & delta_time_C,                      & ! in
    & burned_fract,                      & ! in
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
    REAL(lp), INTENT(in)    :: delta_time_C                 !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: burned_fract(:)              !< Fraction of the vegetated area of each tile burned till the last
                                                            !< call of this routine
    REAL(lp), INTENT(in)    :: fire_fract_wood_2_atmos      !< Fraction of above ground wood immediately emitted to atmosphere
                                                            !< in case of fire
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)              !< Factor to spread non woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)              !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                   !< Amount of green carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                   !< Amount of wood carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                 !< Amount of reserve carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)   !< Carbon relocated by disturbances from non-woody biomass
                                                            !< to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   !< Carbon relocated by disturbances from woody biomass to
                                                            !< the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)             !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)             !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1(:)                 !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)             !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)             !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2(:)                 !< Humus pool from woody litter [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    !< CO2 released by fire [kg(CO2)m-2(canopy)s-1]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Carbon released due to disturbances from the green and reserve pools, that enters the green litter and
      ! soil pools. We here add the contribution of fires.
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._lp - fract_green_aboveGround) * burned_fract(ic) / delta_time_C

      ! Carbon released due to disturbances from the woody vegetation pool, that enters the woody litter and soil
      ! pools. We here add the contribution of fires.
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._lp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(ic) / delta_time_C

      ! Determine the carbon flux to the atmosphere due to fire: We assume that the above ground non-woody
      ! biomass and above ground litter are burnt completely, while parts of the above ground woody biomass
      ! and below ground matter resists the flames.
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
        & ) * burned_fract(ic) * molarMassCO2_kg / delta_time_C

      ! Reduce the above ground yasso pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_acid_ag1(ic)       = fraction * c_acid_ag1(ic)
      c_water_ag1(ic)      = fraction * c_water_ag1(ic)
      c_ethanol_ag1(ic)    = fraction * c_ethanol_ag1(ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1(ic)

      c_acid_ag2(ic)       = fraction * c_acid_ag2(ic)
      c_water_ag2(ic)      = fraction * c_water_ag2(ic)
      c_ethanol_ag2(ic)    = fraction * c_ethanol_ag2(ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2(ic)

      ! Add below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._lp - fract_green_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remaings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._lp - fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_fire

  !
  !> This routine handles the impact of vegetation damage on the carbon pools.
  !>
  !< In this routine we model the impact of vegetation damage (i.e. wind through) on the carbon pools:
  !< Depending on the disturbed area fraction carbon from the living vegetation pools (c_green,
  !< c_reserve and c_woods) is relocated  to the litter and soil pools.
  !
  SUBROUTINE relocate_carbon_damage(  &
    & delta_time_C,                   & ! in
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
    REAL(lp), INTENT(in)    :: delta_time_C                !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: damaged_fract(:)            !< Fraction of the vegetated area of each tile damaged since the last
                                                           !< call of this routine
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)             !< Factor to spread non woody litter into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)             !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                  !< Amount of carbon in green pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                  !< Amount carbon in wood pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                !< Amount of carbon in reserve pool [mol(C)/m^(canopy)2]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)  !< Disturbance-induced carbon flux of non-woody living biomass
                                                           !< to the green litter and soil pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)  !< Disturbance-induced carbon flux of woody living biomass to
                                                           !< the woody litter and soil pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag1(:)             !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg1(:)             !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_1(:)                 !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag2(:)             !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg2(:)             !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_2(:)                 !< Humus pool from woody litter [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Carbon transferred from the wood pool to the woody litter and soil pools due to disturbances;
      ! we here add the contribution of wind throw.
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + c_woods(ic) * damaged_fract(ic) / delta_time_C

      ! Carbon transferred from the non-woody vegetation pools to the green litter and soil pools due to
      ! disturbances; we here add the contribution of wind throw.
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + non_woody_litter * damaged_fract(ic) / delta_time_C

      ! Distribute the carbon of damaged non-woody vegetation to the yasso litter and soil pools
      !  Above ground
      fraction = fract_green_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), &
        &                          c_nonsoluble_ag1(ic), c_humus_1(ic),                &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._lp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), &
        &                          c_nonsoluble_bg1(ic), c_humus_1(ic),                &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the carbon of damaged woody vegetation to the yasso litter and soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_damage


  !
  !> Compute soil carbon dynamics for single layer below ground pools.
  !>
  !> Carbon fluxes from litter and root exudates are added to the respective carbon pools,
  !> which are then decomposed into humus according to the single-layer Yasso model.
  !
  SUBROUTINE calc_soil_c_single_layer(                                &
    & nc, nsoil, delta_time_C, LeafLit_coef, WoodLit_coef, MiscPar,   &
    & t_air, precip,                                                  &
    & Cflx_wood_2_litter, Cflx_green_2_litter, Cflx_2_root_exudates,  &
    & c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1,       &
    & c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2,       &
    & c_acid_bg1, c_water_bg1, c_ethanol_bg1, c_nonsoluble_bg1,       &
    & c_acid_bg2, c_water_bg2, c_ethanol_bg2, c_nonsoluble_bg2,       &
    & c_humus_1,  c_humus_2,                                          &
    & c_into_humus_1_sum, c_decomp_humus_1_sum,                       &
    & c_into_humus_2_sum, c_decomp_humus_2_sum, Cflx_soil_respiration &
    & )

    USE mo_carbon_constants, ONLY: WOODLITSIZE, ROOT_DIST_BETA

    INTEGER,  INTENT(in)    :: nc                       !< Number of grid cells []
    INTEGER,  INTENT(in)    :: nsoil                    !< Number of soil layers in vertical yasso scheme []
    REAL(lp), INTENT(in)    :: delta_time_C             !< Carbon time step [s]
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)          !< Leaf litter fraction going to specific pools
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)          !< Woody litter fraction going to specific pools
    REAL(lp), INTENT(in)    :: MiscPar(:)               !< miscellaneous parameters
    ! Input meteorology for Yasso
    REAL(lp), INTENT(in)    :: t_air(:)                 !< 15-day mean air temperature [K]
    REAL(lp), INTENT(in)    :: precip(:)                !< 15-day mean precipitation [kg m-2 s-1]
    ! Carbon inputs [mol m-2 s-1]:
    REAL(lp), INTENT(in)    :: Cflx_wood_2_litter(:)    !< Carbon flux from woody litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(in)    :: Cflx_green_2_litter(:)   !< Carbon flux from green litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(in)    :: Cflx_2_root_exudates(:)  !< Carbon flux from root exudates to soil [mol(C) m-2 s-1]
    ! Soil carbon pools [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)         !< Above-ground acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)        !< Above-ground water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)      !< Above-ground ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)   !< Above-ground nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)         !< Above-ground acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)        !< Above-ground water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)      !< Above-ground ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)   !< Above-ground nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)         !< Below-ground acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)        !< Below-ground water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)      !< Below-ground ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)   !< Below-ground nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)         !< Below-ground acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)        !< Below-ground water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)      !< Below-ground ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)   !< Below-ground nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1(:)          !< Below-ground humus carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2(:)          !< Below-ground humus carbon pools of green woody [mol(C) m-2]

    ! Decomposistion fluxes
    REAL(lp), INTENT(inout) :: c_into_humus_1_sum(:)    !< Annual sum of carbon decomposed into humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_decomp_humus_1_sum(:)  !< Annual sum of carbon of humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_into_humus_2_sum(:)    !< Annual sum of carbon decomposed into humus (wood) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_decomp_humus_2_sum(:)  !< Annual sum of carbon of humus (wood) [mol(C) m-2 yr-1]

    REAL(lp), INTENT(out)   :: Cflx_soil_respiration(:) !< Soil (=heterotrophic) respiration rate  [mol(C) m-2 s-1]
                                                        ! Note: this is a loss, therefore negative!

    ! local variables for Yasso
    INTEGER  :: ic                           !< Grid cell index []
    REAL(lp) :: LeafLitterSize               !< Leaf litter size (always 0 -- > size has no effect on decomposition)
    REAL(lp) :: WoodLitterSize               !< Wood litter size (PFT dependent)
    REAL(lp) :: Cflx_soil_respiration_pool   !< Soil respiration flux from carbon pool [mol(C) m-2 d-1]
    REAL(lp) :: Cflx_humus_gain              !< Humus pool carbon gain (litter decomposition) [mol(C) m-2 d-1]
    REAL(lp) :: Cflx_humus_loss              !< Humus pool carbon loss (humus decomposition) [mol(C) m-2 d-1]
    REAL(lp) :: d_temp                       !< surface or soil temperature impact on litter decomposition []
    REAL(lp) :: d_moist                      !< precipitation or soil moisture impact on litter decomposition []
    REAL(lp) :: d_size_leaf                  !< green (leaf) litter size impact on litter decomposition []
    REAL(lp) :: d_size_wood                  !< wood litter size impact on litter decomposition []
    REAL(lp) :: d_depth                      !< soil deptj impact on litter decomposition []
    REAL(lp) :: f_roots                      !< root distribution factor for at soil layer []


    ! Set globally uniform values
    LeafLitterSize = 0._lp
    WoodLitterSize = MiscPar(WOODLITSIZE)
    ! Litter size dependence of decomposition -- no effect if LitterSize = 0.0 (e.g. for leafs)
    d_size_leaf    = get_modifier_size(LeafLitterSize)
    d_size_wood    = get_modifier_size(WoodLitterSize)
    ! As the single layer yasso does not account for depth dependence, these parameters are set to one
    d_depth = 1._lp
    f_roots = 1._lp

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss) &
    !$ACC   PRIVATE(d_temp, d_moist)
    DO ic = 1, nc
      ! Initialize and set decomposition parameters
      !-------------------------------------------
      Cflx_soil_respiration(ic) = 0._lp
      d_temp  = get_modifier_temp(t_air(ic))
      d_moist = get_modifier_precip(precip(ic))

      ! ============== COMPUTE SOIL RESPIRATION ===================================================================
      !
      ! Compute yasso with single layer below ground pools
      ! Above ground green (leaf) litter
      CALL yasso(ABOVE_GROUND, delta_time_C, d_temp, d_moist, d_size_leaf, d_depth, f_roots,              &
        &        Cflx_green_2_litter(ic) * delta_time_C, LeafLit_coef(:), fract_green_aboveGround,        &
        &        c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
      c_into_humus_1_sum(ic)    = c_into_humus_1_sum(ic)    + Cflx_humus_gain
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C

      ! Above ground woody litter
      CALL yasso(ABOVE_GROUND, delta_time_C, d_temp, d_moist, d_size_wood, d_depth, f_roots,              &
        &        Cflx_wood_2_litter(ic) * delta_time_C, WoodLit_coef(:), fract_wood_aboveGround,          &
        &        c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
      c_into_humus_2_sum(ic)    = c_into_humus_2_sum(ic)    + Cflx_humus_gain
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C

      ! Below ground green (leaf) litter, root exudates are used only here
      CALL yasso(BELOW_GROUND, delta_time_C, d_temp, d_moist, d_size_leaf, d_depth, f_roots,              &
        &        Cflx_green_2_litter(ic) * delta_time_C, LeafLit_coef(:), fract_green_aboveGround,        &
        &        c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss,                            &
        &        Cflx_2_root_exudates(ic) * delta_time_C)
      c_into_humus_1_sum(ic)    = c_into_humus_1_sum(ic)    + Cflx_humus_gain
      c_decomp_humus_1_sum(ic)  = c_decomp_humus_1_sum(ic)  + Cflx_humus_loss
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C

      ! Below ground woody litter
      CALL yasso(BELOW_GROUND, delta_time_C, d_temp, d_moist, d_size_wood, d_depth, f_roots,              &
        &        Cflx_wood_2_litter(ic) * delta_time_C, WoodLit_coef(:), fract_wood_aboveGround,          &
        &        c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
      c_into_humus_2_sum(ic)    = c_into_humus_2_sum(ic)    + Cflx_humus_gain
      c_decomp_humus_2_sum(ic)  = c_decomp_humus_2_sum(ic)  + Cflx_humus_loss
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_soil_c_single_layer
  !
  !> Compute soil carbon dynamics.
  !>
  !> Carbon fluxes from litter and root exudates are added to the respective carbon pools,
  !> which are then decomposed into humus according to the Yasso model.
  !
  SUBROUTINE calc_soil_c_multi_layer( &
    & nc, nsoil, delta_time_C, LeafLit_coef, WoodLit_coef, MiscPar,         &
    & soil_depth, soil_depth_sl, max_depth, max_depth_sl,                   &
    & layer_depth_sl, layer_thick, thaw_depth, permafrost_period_sl,           &
    & t_air, precip, t_soil_sl, wtr_soil_rel_sl, n_thaw_cycle,              &
    & Cflx_wood_2_litter, Cflx_green_2_litter, Cflx_2_root_exudates,        &
    & c_acid_ag1,    c_water_ag1,    c_ethanol_ag1,    c_nonsoluble_ag1,    &
    & c_acid_ag2,    c_water_ag2,    c_ethanol_ag2,    c_nonsoluble_ag2,    &
    & c_acid_bg1,    c_water_bg1,    c_ethanol_bg1,    c_nonsoluble_bg1,    &
    & c_acid_bg2,    c_water_bg2,    c_ethanol_bg2,    c_nonsoluble_bg2,    &
    & c_acid_bg1_sl, c_water_bg1_sl, c_ethanol_bg1_sl, c_nonsoluble_bg1_sl, &
    & c_acid_bg2_sl, c_water_bg2_sl, c_ethanol_bg2_sl, c_nonsoluble_bg2_sl, &
    & c_humus_1,     c_humus_2,      c_humus_1_sl,     c_humus_2_sl,        &
    & c_into_humus_1_sum, c_decomp_humus_1_sum,                             &
    & c_into_humus_2_sum, c_decomp_humus_2_sum, Cflx_soil_respiration       &
    & )

    USE mo_carbon_constants, ONLY: WOODLITSIZE, ROOT_DIST_BETA

    INTEGER,  INTENT(in)    :: nc                    !< Number of grid cells []
    INTEGER,  INTENT(in)    :: nsoil                 !< Number of soil layers in vertical yasso scheme []
    REAL(lp), INTENT(in)    :: delta_time_C          !< Carbon time step [s]
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)       !< Leaf litter fraction going to specific pools
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)       !< Woody litter fraction going to specific pools
    REAL(lp), INTENT(in)    :: MiscPar(:)            !< miscellaneous parameters
    ! Input meteorology for Yasso
    REAL(lp), INTENT(in)    :: soil_depth(:)            !< Total soil depth until bedrock [m]
    REAL(lp), INTENT(in)    :: soil_depth_sl(:,:)       !< Soil depth within each layer [m]
    REAL(lp), INTENT(in)    :: max_depth(:)             !< maximum available depth [m], either soil or thaw
    REAL(lp), INTENT(in)    :: max_depth_sl(:,:)        !< maximum available depth [m], either soil or thaw, per layer
    REAL(lp), INTENT(in)    :: layer_depth_sl(:)        !< Depths of soil layer center [m]
    REAL(lp), INTENT(in)    :: layer_thick(:)           !< Thickness of soil layer [m]
    REAL(lp), INTENT(in)    :: thaw_depth(:)            !< thaw depth [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:,:)  !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(in)    :: t_air(:)                 !< 15-day mean air temperature [K]
    REAL(lp), INTENT(in)    :: precip(:)                !< 15-day mean precipitation [kg m-2 s-1]
    REAL(lp), INTENT(in)    :: t_soil_sl(:,:)           !< soil layer temperatures [K]
    REAL(lp), INTENT(in)    :: wtr_soil_rel_sl(:,:)     !< relative soil moisture []
    REAL(lp), INTENT(in)    :: n_thaw_cycle(:)          !< number of freeze-thaw cycles in previous year
    ! Carbon inputs [mol m-2 s-1]:
    REAL(lp), INTENT(in)    :: Cflx_wood_2_litter(:)    !< Carbon flux from woody litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(in)    :: Cflx_green_2_litter(:)   !< Carbon flux from green litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(in)    :: Cflx_2_root_exudates(:)  !< Carbon flux from root exudates to soil [mol(C) m-2 s-1]
    ! Soil carbon pools [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)          !< Above-ground acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)         !< Above-ground water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)       !< Above-ground ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)    !< Above-ground nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)          !< Above-ground acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)         !< Above-ground water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)       !< Above-ground ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)    !< Above-ground nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)          !< Below-ground acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)         !< Below-ground water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)       !< Below-ground ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)    !< Below-ground nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)          !< Below-ground acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)         !< Below-ground water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)       !< Below-ground ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)    !< Below-ground nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg1_sl(:,:)       !< Soil layer acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg1_sl(:,:)      !< Soil layer water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1_sl(:,:)    !< Soil layer ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1_sl(:,:) !< Soil layer nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg2_sl(:,:)       !< Soil layer acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg2_sl(:,:)      !< Soil layer water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2_sl(:,:)    !< Soil layer ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2_sl(:,:) !< Soil layer nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1(:)           !< Below-ground humus carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2(:)           !< Below-ground humus carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:,:)      !< Soil layer humus carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:,:)      !< Soil layer humus carbon pools of green woody [mol(C) m-2]
    ! Decomposistion fluxes
    REAL(lp), INTENT(inout) :: c_into_humus_1_sum(:)    !< Annual sum of carbon decomposed into humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_decomp_humus_1_sum(:)  !< Annual sum of carbon of humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_into_humus_2_sum(:)    !< Annual sum of carbon decomposed into humus (wood) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_decomp_humus_2_sum(:)  !< Annual sum of carbon of humus (wood) [mol(C) m-2 yr-1]
    ! Soil respiration flux
    REAL(lp), INTENT(out)   :: Cflx_soil_respiration(:) !< Soil (=heterotrophic) respiration rate  [mol(C) m-2 s-1]
                                                        ! Note: this is a loss, therefore negative!

    ! local variables for Yasso
    INTEGER  :: ic                         !< Grid cell index []
    INTEGER  :: is                         !< Soil layer index []
    REAL(lp) :: root_dist_sl(nc,nsoil)     !< Root distribution factor for each soil layer []
    REAL(lp) :: diff_coeff(nc,nsoil)       !< Vertical carbon flux diffusion coefficient [m2 s-1]
    REAL(lp) :: LeafLitterSize             !< Leaf litter size (always 0 -- > size has no effect on decomposition)
    REAL(lp) :: WoodLitterSize             !< Wood litter size (PFT dependent)
    REAL(lp) :: Cflx_soil_respiration_pool !< Soil respiration flux from carbon pool [mol(C) m-2 d-1]
    REAL(lp) :: Cflx_humus_gain            !< Humus pool carbon gain (litter decomposition) [mol(C) m-2 d-1]
    REAL(lp) :: Cflx_humus_loss            !< Humus pool carbon loss (humus decomposition) [mol(C) m-2 d-1]
    REAL(lp) :: d_temp                     !< surface or soil temperature impact on litter decomposition []
    REAL(lp) :: d_moist                    !< precipitation or soil moisture impact on litter decomposition []
    REAL(lp) :: d_size_leaf                !< green (leaf) litter size impact on litter decomposition []
    REAL(lp) :: d_size_wood                !< wood litter size impact on litter decomposition []
    REAL(lp) :: d_depth                    !< soil deptj impact on litter decomposition []
    REAL(lp) :: f_roots                    !< root distribution factor for at soil layer []

    !$ACC DATA CREATE(root_dist_sl, diff_coeff)

    ! Set globally uniform values
    LeafLitterSize = 0._lp
    WoodLitterSize = MiscPar(WOODLITSIZE)
    ! Litter size dependence of decomposition -- no effect if LitterSize = 0.0 (e.g. for leafs)
    d_size_leaf    = get_modifier_size(LeafLitterSize)
    d_size_wood    = get_modifier_size(WoodLitterSize)

    ! ============== COMPUTE SOIL RESPIRATION for above ground pools ============================================
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss) &
    !$ACC   PRIVATE(d_temp, d_moist, d_depth, f_roots)
    DO ic = 1, nc
      ! Initialize and set decomposition parameters
      !-------------------------------------------
      Cflx_soil_respiration(ic) = 0._lp
      d_temp  = get_modifier_temp(t_air(ic))
      d_moist = get_modifier_precip(precip(ic))
      ! Above ground pools do not account for depth dependence, these parameters are set to one
      d_depth = 1._lp
      f_roots = 1._lp

      ! Above ground green (leaf) litter
      CALL yasso(ABOVE_GROUND, delta_time_C, d_temp, d_moist, d_size_leaf, d_depth, f_roots,                   &
        &        Cflx_green_2_litter(ic) * delta_time_C, LeafLit_coef(:), fract_green_aboveGround,             &
        &        c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1_sl(ic,1), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
      c_into_humus_1_sum(ic)    = c_into_humus_1_sum(ic)    + Cflx_humus_gain
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C

      ! Above ground woody litter
      CALL yasso(ABOVE_GROUND, delta_time_C, d_temp, d_moist, d_size_wood, d_depth, f_roots,                   &
        &        Cflx_wood_2_litter(ic) * delta_time_C, WoodLit_coef(:), fract_wood_aboveGround,               &
        &        c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2_sl(ic,1), &
        &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
      c_into_humus_2_sum(ic)    = c_into_humus_2_sum(ic)    + Cflx_humus_gain
      Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C
    END DO
    !$ACC END PARALLEL LOOP

    ! ============== COMPUTE SOIL RESPIRATION for below ground pools ============================================
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      root_dist_sl(ic,:) = determine_root_dist(max_depth_sl(ic,:), layer_depth_sl(:), MiscPar(ROOT_DIST_BETA), &
        &                                   t_soil_sl(ic,:), permafrost_period_sl(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    !$ACC   PRIVATE(Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss) &
    !$ACC   PRIVATE(d_temp, d_moist, d_depth, f_roots)
    DO is = 1, nsoil
      DO ic = 1, nc
        ! Initialize and set decomposition parameters
        !-------------------------------------------
        d_temp  = get_modifier_temp(t_soil_sl(ic,is))
        d_moist = get_modifier_water(wtr_soil_rel_sl(ic,is))
        d_depth = get_modifier_depth(layer_depth_sl(is))
        f_roots = root_dist_sl(ic,is)
        ! Below ground green (leaf) litter, root exudates are used only here
        CALL yasso(BELOW_GROUND_SL, delta_time_C, d_temp, d_moist, d_size_leaf, d_depth, f_roots,           &
          &        Cflx_green_2_litter(ic) * delta_time_C, LeafLit_coef(:), fract_green_aboveGround,        &
          &        c_acid_bg1_sl(ic,is), c_water_bg1_sl(ic,is), c_ethanol_bg1_sl(ic,is),                    &
          &        c_nonsoluble_bg1_sl(ic,is), c_humus_1_sl(ic,is),                                         &
          &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss,                            &
          &        Cflx_2_root_exudates(ic) * delta_time_C)
        c_into_humus_1_sum(ic)    = c_into_humus_1_sum(ic)    + Cflx_humus_gain
        c_decomp_humus_1_sum(ic)  = c_decomp_humus_1_sum(ic)  + Cflx_humus_loss
        Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C

        ! Below ground woody litter
        CALL yasso(BELOW_GROUND_SL, delta_time_C, d_temp, d_moist, d_size_wood, d_depth, f_roots,           &
          &        Cflx_wood_2_litter(ic) * delta_time_C, WoodLit_coef(:), fract_wood_aboveGround,          &
          &        c_acid_bg2_sl(ic,is), c_water_bg2_sl(ic,is), c_ethanol_bg2_sl(ic,is),                    &
          &        c_nonsoluble_bg2_sl(ic,is), c_humus_2_sl(ic,is),                                         &
          &        Cflx_soil_respiration_pool, Cflx_humus_gain, Cflx_humus_loss)
        c_into_humus_2_sum(ic)    = c_into_humus_2_sum(ic)    + Cflx_humus_gain
        c_decomp_humus_2_sum(ic)  = c_decomp_humus_2_sum(ic)  + Cflx_humus_loss
        Cflx_soil_respiration(ic) = Cflx_soil_respiration(ic) + Cflx_soil_respiration_pool / delta_time_C
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! ============== VERTICAL MIXING OF CARBON POOLS ==========================================================
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      ! Move carbon downwards that exceed the maximum possible carbon content of a layer
      CALL move_ag_2_bg(nsoil, max_depth(ic), layer_thick(:), permafrost_period_sl(ic,:),               &
        & c_acid_ag1(ic),      c_water_ag1(ic),      c_ethanol_ag1(ic),      c_nonsoluble_ag1(ic),      &
        & c_acid_ag2(ic),      c_water_ag2(ic),      c_ethanol_ag2(ic),      c_nonsoluble_ag2(ic),      &
        & c_acid_bg1_sl(ic,:), c_water_bg1_sl(ic,:), c_ethanol_bg1_sl(ic,:), c_nonsoluble_bg1_sl(ic,:), &
        & c_acid_bg2_sl(ic,:), c_water_bg2_sl(ic,:), c_ethanol_bg2_sl(ic,:), c_nonsoluble_bg2_sl(ic,:), &
        & c_humus_1_sl(ic,:),  c_humus_2_sl(ic,:))
      ! Compute diffusivity for vertical carbon transport
      CALL set_diffusivity(nsoil, soil_depth(ic), soil_depth_sl(ic,:), thaw_depth(ic), layer_depth_sl(:), &
        & t_soil_sl(ic,:), wtr_soil_rel_sl(ic,:), n_thaw_cycle(ic), permafrost_period_sl(ic,:), diff_coeff(ic,:))
      ! Compute vertical transport (cryoturbation and bioturbation) for each below ground carbon pool
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_acid_bg1_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_water_bg1_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_ethanol_bg1_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_nonsoluble_bg1_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_humus_1_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_acid_bg2_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_water_bg2_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_ethanol_bg2_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_nonsoluble_bg2_sl(ic,:))
      CALL vertical_C_transport(max_depth_sl(ic,:), diff_coeff(ic,:), delta_time_C, layer_thick(:), &
        & permafrost_period_sl(ic,:), c_humus_2_sl(ic,:))
      ! Remove any C that might have ended up below the bedrock boundary
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_acid_bg1_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_water_bg1_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_ethanol_bg1_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_nonsoluble_bg1_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_humus_1_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_acid_bg2_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_water_bg2_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_ethanol_bg2_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_nonsoluble_bg2_sl(ic,:))
      CALL remove_C_below_bedrock(soil_depth(ic), layer_thick(:), c_humus_2_sl(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

    ! Provide all variables as sum, too, not just on levels
    CALL sum_up_carbon_layers(nc, nsoil, c_acid_bg1_sl(:,:)      , c_acid_bg1(:)       )
    CALL sum_up_carbon_layers(nc, nsoil, c_water_bg1_sl(:,:)     , c_water_bg1(:)      )
    CALL sum_up_carbon_layers(nc, nsoil, c_ethanol_bg1_sl(:,:)   , c_ethanol_bg1(:)    )
    CALL sum_up_carbon_layers(nc, nsoil, c_nonsoluble_bg1_sl(:,:), c_nonsoluble_bg1(:) )
    CALL sum_up_carbon_layers(nc, nsoil, c_acid_bg2_sl(:,:)      , c_acid_bg2(:)       )
    CALL sum_up_carbon_layers(nc, nsoil, c_water_bg2_sl(:,:)     , c_water_bg2(:)      )
    CALL sum_up_carbon_layers(nc, nsoil, c_ethanol_bg2_sl(:,:)   , c_ethanol_bg2(:)    )
    CALL sum_up_carbon_layers(nc, nsoil, c_nonsoluble_bg2_sl(:,:), c_nonsoluble_bg2(:) )
    CALL sum_up_carbon_layers(nc, nsoil, c_humus_1_sl(:,:)       , c_humus_1(:)        )
    CALL sum_up_carbon_layers(nc, nsoil, c_humus_2_sl(:,:)       , c_humus_2(:)        )

  END SUBROUTINE calc_soil_c_multi_layer


#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE yasso (yasso_case, delta_time_C, d_temp, d_moist, d_size, d_depth, f_roots, &
    &               c_2_litter, Lit_coef, fract_aboveground,                             &
    &               c_acid, c_water, c_ethanol, c_nonsoluble, c_humus,                   &
    &               Cflx_soil_respiration, Cflx_humus_gain, Cflx_humus_loss,             &
    &               c_2_rootExudates)
    !$ACC ROUTINE SEQ

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
    USE mo_carbon_constants, ONLY: sec_per_year, ref_decomp_rate_acid, ref_decomp_rate_water, &
      & ref_decomp_rate_ethanol, ref_decomp_rate_nonsoluble, ref_decomp_rate_humus,           &
      & A_2_W, A_2_E, A_2_N, W_2_A, W_2_E, W_2_N, E_2_A, E_2_W, E_2_N, N_2_A, N_2_W, N_2_E,   &
      & AWE_2_H, N_2_H

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: yasso_case                !< Yasso case: ABOVE_GROUND or BELOW_GROUND
    REAL(lp), INTENT(in)  :: delta_time_C              !< Carbon time step [s]
    REAL(lp), INTENT(in)  :: d_temp                    !< Surface or soil temperature impact on litter decomposition []
    REAL(lp), INTENT(in)  :: d_moist                   !< Precipitation or soil moisture impact on litter decomposition []
    REAL(lp), INTENT(in)  :: d_size                    !< Litter size impact on litter decomposition []
    REAL(lp), INTENT(in)  :: d_depth                   !< Soil depth impact on litter decomposition []
    REAL(lp), INTENT(in)  :: f_roots                   !< Root distribution factor for at soil layer []
    REAL(lp), INTENT(in)  :: c_2_litter                !< Fresh litter inputs  (above and belowground) [mol(C) m-2 d-1]
    REAL(lp), INTENT(in)  :: Lit_coef(:)               !< Fractions to seperate incoming litter to the yasso pool []
    REAL(lp), INTENT(in)  :: fract_aboveground         !< Parameter to seperate litter influx into above- and belowground parts []

    REAL(lp), INTENT(inout) :: c_acid                  !< Acid soluble carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water                 !< Water soluble carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol               !< Ethanol soluble carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble            !< Nonsoluble carbon pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus                 !< Humus carbon pool [mol(C) m-2]

    REAL(lp), INTENT(out)   :: Cflx_soil_respiration   !< Soil respiration flux [mol(C) m-2 d-1]
    REAL(lp), INTENT(out)   :: Cflx_humus_gain         !< Carbon flux into humus pool [mol(C) m-2 d-1]
    REAL(lp), INTENT(out)   :: Cflx_humus_loss         !< Carbon flux from humus pool [mol(C) m-2 d-1]

    REAL(lp), INTENT(in), OPTIONAL :: c_2_rootExudates !< root exudates flux (below ground leaf litter only) [mol(C) m-2 d-1]

    ! Carbon loss and gain fluxes [mol(C) m-2 yr-1]
    REAL(lp)   :: Cflx_yr_from_acid                 !< Loss flux of acid soluble pool
    REAL(lp)   :: Cflx_yr_from_water                !< Loss flux of water soluble pool
    REAL(lp)   :: Cflx_yr_from_ethanol              !< Loss flux of ethanol soluble pool
    REAL(lp)   :: Cflx_yr_from_nonsoluble           !< Loss flux of non soluble pool
    REAL(lp)   :: Cflx_yr_from_humus                !< Loss flux of humus pool

    REAL(lp)   :: Cflx_yr_2_acid                    !< Gain flux of acid soluble pool
    REAL(lp)   :: Cflx_yr_2_water                   !< Gain flux of water soluble pool
    REAL(lp)   :: Cflx_yr_2_ethanol                 !< Gain flux of ethanol soluble pool
    REAL(lp)   :: Cflx_yr_2_nonsoluble              !< Gain flux of non soluble pool
    REAL(lp)   :: Cflx_yr_2_humus                   !< Gain flux of humus pool

    REAL(lp)   :: Cflx_yr_litter_2_humus            !< Litter influx to humus soluble pool

    ! Respiration is output of Yasso
    REAL(lp)   :: Cflx_yr_respiration               !< Total flux to the atmosphere
    REAL(lp)   :: Cflx_yr_resp_from_litter          !< Flux to the atmosphere; from AWEN ag only

    REAL(lp)   :: fract_litter            !< litter influx fraction into above- or belowground pools [ ]

    ! Decomposition rates [yr-1]
    REAL(lp)   :: d_acid
    REAL(lp)   :: d_water
    REAL(lp)   :: d_ethanol
    REAL(lp)   :: d_nonsoluble
    REAL(lp)   :: d_humus

    ! Climate dependencies
    REAL(lp)   :: d_pool                  !< term which accounds for the pool type (above,below) for humus decomposition

    ! Decomposition rate for N litter green (needed for N cycle)
    REAL(lp)   :: d_litter_green          !< the rate with which the sum of AWEN pools decays
    REAL(lp)   :: Pseudo_litter_green     !< Sum of all AWEN pools

    ! time stepping: yasso is parameterized based on annual time
    ! steps, thus the fluxes have to be scaled up to annual rates.
    REAL(lp)   :: years_per_ts, ts_per_year

    ts_per_year  = REAL(sec_per_year, lp) / delta_time_C
    years_per_ts = 1.0_lp / ts_per_year

    ! Effects of above vs below ground pools
    IF (yasso_case == ABOVE_GROUND) THEN
      fract_litter = fract_aboveground
      d_pool       = 0._lp
    ELSE
      fract_litter = MAX(0._lp, 1.0_lp - fract_aboveground) * f_roots
      d_pool       = 1._lp
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Decomposition rates accounting for temperature, precipitation, litter size or pool type, and nutrient limitation
    d_acid       = ref_decomp_rate_acid       * d_temp * d_moist * d_size * d_depth
    d_water      = ref_decomp_rate_water      * d_temp * d_moist * d_size * d_depth
    d_ethanol    = ref_decomp_rate_ethanol    * d_temp * d_moist * d_size * d_depth
    d_nonsoluble = ref_decomp_rate_nonsoluble * d_temp * d_moist * d_size * d_depth
    d_humus      = ref_decomp_rate_humus      * d_temp * d_moist * d_pool * d_depth

    ! 2. Calculate carbon fluxes

    !   annual loss fluxes (negative values):
    Cflx_yr_from_acid       = c_acid       * d_acid
    Cflx_yr_from_water      = c_water      * d_water
    Cflx_yr_from_ethanol    = c_ethanol    * d_ethanol
    Cflx_yr_from_nonsoluble = c_nonsoluble * d_nonsoluble
    Cflx_yr_from_humus      = c_humus      * d_humus

    !   annual gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_yr_2_acid         = ABS(Cflx_yr_from_water        * W_2_A  &
      &                        + Cflx_yr_from_ethanol      * E_2_A  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_A)

    Cflx_yr_2_water        = ABS(Cflx_yr_from_acid         * A_2_W  &
      &                        + Cflx_yr_from_ethanol      * E_2_W  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_W)

    Cflx_yr_2_ethanol      = ABS(Cflx_yr_from_acid         * A_2_E  &
      &                        + Cflx_yr_from_water        * W_2_E  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_E)

    Cflx_yr_2_nonsoluble   = ABS(Cflx_yr_from_acid         * A_2_N  &
      &                        + Cflx_yr_from_water        * W_2_N  &
      &                        + Cflx_yr_from_ethanol      * E_2_N)

    Cflx_yr_2_humus        = ABS(Cflx_yr_from_acid                  &
      &                        + Cflx_yr_from_water                 &
      &                        + Cflx_yr_from_ethanol               &
      &                          ) * AWE_2_H                        &
      &                    +                                        &
      &                      ABS(Cflx_yr_from_nonsoluble)           &
      &                            * N_2_H

    !   litter respiration
    Cflx_yr_resp_from_litter =                                                 &
      &   Cflx_yr_from_acid       * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      & + Cflx_yr_from_water      * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      & + Cflx_yr_from_ethanol    * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      & + Cflx_yr_from_nonsoluble * (1.0_lp - N_2_A - N_2_W - N_2_E -   N_2_H)

    ! 3. update Yasso pools

    !   redistribution among yasso pools
    c_acid        = MAX(0.0_lp,                                      &
      &                 c_acid                                       & ! old pool
      &                 + (Cflx_yr_from_acid                         & ! incoming flux from AWEN pools
      &                 + Cflx_yr_2_acid) * years_per_ts)
    c_ethanol     = MAX(0.0_lp,                                      & ! and so on .....
      &                 c_ethanol                                    &
      &                 + (Cflx_yr_from_ethanol                      &
      &                 + Cflx_yr_2_ethanol) * years_per_ts)
    c_nonsoluble  = MAX(0.0_lp,                                      &
      &                 c_nonsoluble                                 &
      &                 + (Cflx_yr_from_nonsoluble                   &
      &                 + Cflx_yr_2_nonsoluble) * years_per_ts)
    IF (yasso_case == ABOVE_GROUND) THEN
      ! Above ground only consideres decomposition into humus
      c_humus     = MAX(0.0_lp,                                      &
        &               c_humus                                      &
        &               + Cflx_yr_2_humus * years_per_ts)
    ELSE
      ! Below ground also considers decomposition of humus
      c_humus     = MAX(0.0_lp,                                      &
        &               c_humus                                      &
        &               + (Cflx_yr_from_humus                        &
        &               + Cflx_yr_2_humus) * years_per_ts)
    END IF
    IF (PRESENT(c_2_rootExudates)) THEN
      ! Root exudates are only considered carbonhydrates and only for below ground leaf litter
      c_water   = MAX(0.0_lp,                                      &
        &             c_water                                      &
        &             + (Cflx_yr_from_water                        &
        &             + Cflx_yr_2_water) * years_per_ts            &
        &             + c_2_rootExudates * f_roots)
    ELSE
      c_water     = MAX(0.0_lp,                                      &
        &               c_water                                      &
        &               + (Cflx_yr_from_water                        &
        &               + Cflx_yr_2_water) * years_per_ts)
    END IF

    ! 4.0 Distribute litter to the pools
    CALL distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
      &                          c_2_litter, fract_litter, Lit_coef)

    Cflx_yr_litter_2_humus = Lit_coef(Y_HUMUS) * c_2_litter * ts_per_year

    ! 4.2 Compute the total soil respiration flux (from AWEN and humus pools)
    Cflx_yr_respiration   = Cflx_yr_resp_from_litter + Cflx_yr_from_humus
    Cflx_soil_respiration = Cflx_yr_respiration * years_per_ts ! convert back to daily

    ! 4.3 Write diagnostics to determine equilibrium of humus pools
    Cflx_humus_gain = (Cflx_yr_2_humus + Cflx_yr_litter_2_humus) * years_per_ts
    Cflx_humus_loss = -Cflx_yr_from_humus * years_per_ts

  END SUBROUTINE yasso

  ! ====================================================================================================== !
  !
  !> Distributes incomming litter to yasso pools (i.e. proportional to litter coefficients and passed factors)
  !
#ifndef _OPENACC
  PURE &
#endif
    SUBROUTINE distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
    &                                  litter, fraction, coefficients)

    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp), INTENT(inout)             :: c_acid        !< yasso acid pool
    REAL(lp), INTENT(inout)             :: c_water       !< yasso water pool
    REAL(lp), INTENT(inout)             :: c_ethanol     !< yasso ethanol pool
    REAL(lp), INTENT(inout)             :: c_nonsoluble  !< yasso nonsoluble pool
    REAL(lp), INTENT(inout)             :: c_humus       !< yasso humus pool
    REAL(lp), INTENT(in)                :: litter        !< litter input (above and belowground) [mol(C) m-2 day-1]
    REAL(lp), INTENT(in)                :: fraction      !< fraction of litter put into yasso pools []
    REAL(lp), DIMENSION(5), INTENT(in)  :: coefficients  !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !

    CALL add_litter_to_yasso_pool(c_acid,       litter, fraction, coefficients(Y_ACID))
    CALL add_litter_to_yasso_pool(c_water,      litter, fraction, coefficients(Y_WATER))
    CALL add_litter_to_yasso_pool(c_ethanol,    litter, fraction, coefficients(Y_ETHANOL))
    CALL add_litter_to_yasso_pool(c_nonsoluble, litter, fraction, coefficients(Y_NONSOLUBLE))
    CALL add_litter_to_yasso_pool(c_humus,      litter, fraction, coefficients(Y_HUMUS))

  END SUBROUTINE distribute_yasso_litter

  ! ====================================================================================================== !
  !
  !> Adds incomming litter share (i.e. proportional to litter coefficient and factor)
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE add_litter_to_yasso_pool(this_pool, litter, fraction, coefficient)
    !$ACC ROUTINE SEQ
    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp), INTENT(inout) :: this_pool     !< target yasso pool
    REAL(lp), INTENT(in)    :: litter        !< litter input  (above and belowground) [mol(C) m-2 day-1]
    REAL(lp), INTENT(in)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(lp), INTENT(in)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_lp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,       & ! out
    & Cpool_per_canopy,     & ! in
    & veg_fract_correction, & ! in
    & fract_fpc_max)          ! in

    REAL(lp), INTENT(in), DIMENSION(:) :: &
      & Cpool_per_canopy, &
      & fract_fpc_max,    &
      & veg_fract_correction
    REAL(lp), INTENT(out), DIMENSION(:) :: &
      & Cpool_per_tile

    INTEGER :: ic

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, SIZE(fract_fpc_max)
      Cpool_per_tile(ic) = Cpool_per_canopy(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE get_per_tile

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> We here determine the distribution of plant roots in the soil (based on Jackson 1996), assuming an
  !> exponential decay of root material with depth.
  !>
  !> Note: The root distribution assumed here is only used to determine the vertical profile of root litter components.
  !>       It has no direct affect on the plant hydrology.
  !
  FUNCTION determine_root_dist(max_thick, layer_depth, root_dist_beta, t_soil_sl, permafrost_period_sl) RESULT(root_dist)
    USE mo_carbon_constants,       ONLY: permafrost_threshold
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(lp), INTENT(in)  :: max_thick(:)            !< Active layer thickness [m]
    REAL(lp), INTENT(in)  :: layer_depth(:)          !< Model layer depth [m]
    REAL(lp), INTENT(in)  :: root_dist_beta          !< PFT dependent factor for root decay with depth
    REAL(lp), INTENT(in)  :: t_soil_sl(:)            !< Soil temperature [K]
    REAL(lp), INTENT(in)  :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]


    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':determine_root_dist'
    REAL(lp) :: root_dist(size(max_thick))     ! root distribution in each layer
    REAL(lp) :: root_dist_cum(size(max_thick)) ! cumulative root distribution for layers
    INTEGER  :: is, nsoil, nsoil_p

    nsoil = size(max_thick)
    nsoil_p = nsoil
    root_dist(:) = 0._lp
    root_dist_cum(:) = 0._lp

    ! Determine which soil layers belong to the active layer depth limited by either permafrost
    ! diagnostics or temperatures below zero.
    DO is = 1, nsoil
      IF ((permafrost_period_sl(is) > permafrost_threshold) .OR. (t_soil_sl(is) < tmelt)) THEN
        nsoil_p = MAX(is-1, 1)
        EXIT
      END IF
    END DO

    IF (root_dist_beta > 0._lp) THEN
      DO is = 1, nsoil_p
        root_dist_cum(is) = 1._lp - (root_dist_beta ** (layer_depth(is) * 100._lp))
          IF (is == 1) THEN
            root_dist(is) = root_dist_cum(is)
          ELSE
            root_dist(is) = root_dist_cum(is) - root_dist_cum(is-1)
          END IF
      END DO
    END IF

    ! No root distribution is computed for permafrost affected soil layers
    DO is = 1, nsoil
      IF (permafrost_period_sl(is) > permafrost_threshold) THEN
        root_dist(is) = 0._lp
      END IF
    END DO

    ! normalise root distribution -> SUM(root_dist) = 1
    IF (SUM(root_dist) > 0._lp) THEN
      root_dist(:) = root_dist(:) / SUM(root_dist(:))
    ELSE
      root_dist(:) = 0._lp
    END IF

  END FUNCTION determine_root_dist

  ! ---------------------------------------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  FUNCTION get_modifier_temp(temp_soil, use_Q10) RESULT(d_temp)
    USE mo_carbon_constants,       ONLY: Q10, temp_p1, temp_p2
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(lp), INTENT(in)          :: temp_soil
    LOGICAL, OPTIONAL, INTENT(in) :: use_Q10

    LOGICAL :: do_Q10
    REAL(lp):: d_temp

    IF (PRESENT(use_Q10)) THEN
      do_Q10 = use_Q10
    ELSE
      do_Q10 = .FALSE.
    END IF

    IF (.NOT. do_Q10) THEN       ! use YASSO formulation
      d_temp  = EXP(temp_p1*(temp_soil - tmelt) + &                 ! YASSO formulation for temperature dependence
        &           temp_p2*(temp_soil - tmelt)**2.0_lp)
    ELSE                        ! use alternative Q10 formulation
      d_temp  = Q10**((temp_soil - (tmelt + 10._lp))/10._lp)        ! Temperature modifier as CLM
    END IF
  END FUNCTION get_modifier_temp

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_water(water_soil) RESULT(d_water)
    USE mo_carbon_constants,  ONLY: moist_p_a, moist_p_b, moist_p_d

    REAL(lp), INTENT(in)  :: water_soil
    REAL(lp)              :: d_water

    d_water = ( ( (water_soil - moist_p_b) / (moist_p_a-moist_p_b) )**  &
      &         ( moist_p_d*( (moist_p_b-moist_p_a) / moist_p_a ) ) ) * &
      &       ( (water_soil / moist_p_a)**moist_p_d )

  END FUNCTION get_modifier_water

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_depth(layer_depth) RESULT(d_depth)
    USE mo_carbon_constants,  ONLY: ref_dep

    REAL(lp), INTENT(in)  :: layer_depth
    REAL(lp) :: d_depth

    d_depth = exp(-layer_depth/ref_dep)                       ! depth modifier as CLM

  END FUNCTION get_modifier_depth

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_size(WoodLitterSize) RESULT(d_size)
    USE mo_carbon_constants,       ONLY: size_p1, size_p2, size_p3

    REAL(lp), INTENT(in)  :: WoodLitterSize
    REAL(lp) :: d_size

    d_size = MIN(1.0_lp,(1.0_lp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_lp)**size_p3)

  END FUNCTION get_modifier_size

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_precip(precip) RESULT(d_precip)
    USE mo_carbon_constants,       ONLY: precip_p1, sec_per_year

    REAL(lp), INTENT(in)  :: precip
    REAL(lp) :: d_precip

    d_precip = 1.0_lp - EXP(precip_p1*(precip*sec_per_year/1000._lp))      ! YASSO formulation for precipitation dependence

  END FUNCTION get_modifier_precip

  ! ---------------------------------------------------------------------------------------------------------
  SUBROUTINE vertical_C_transport(act_depth_sl, diff_coeff, delta_time_C, dz, permafrost_period_sl, cpool_bg_sl)
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    REAL(lp), INTENT(in)    :: act_depth_sl(:)         !< actual layer thickness, considering soil and active layer depth [m]
    REAL(lp), INTENT(in)    :: diff_coeff(:)           !< diffusivity
    REAL(lp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: dz(:)                   !< layer thickness [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout) :: cpool_bg_sl(:)          !< below ground C pools on soil layers

    ! layer counters, # soil layers
    INTEGER :: is, nsoil, nsoil_loc
    REAL(lp), PARAMETER :: implicit_contrib = 0.5_lp !< weighting factor for Crank-Nicholson scheme

    ! Diffusion constants averaged over adjacent layers.
    REAL(lp) :: D_below, D_above

    ! Layer-dependent weighting factors in Crank-Nicholson scheme.
    REAL(lp) :: dz_factor, W_below, W_above
    REAL(lp) :: dz_here, dz_above, dz_below
    REAL(lp) :: C_here, C_above, C_below

    ! Diagonal, left and right subdiagonals for Crank-Nicholson matrix
    REAL(lp) :: left(size(act_depth_sl))
    REAL(lp) :: diag(size(act_depth_sl))
    REAL(lp) :: right(size(act_depth_sl))

    ! Requires carbon concentrations for diffusion algorithm
    REAL(lp) :: C_conc(size(act_depth_sl))
    REAL(lp) :: src_C_conc(size(act_depth_sl))

    ! Right hand side vector for Crank-Nicholson scheme equations
    ! This is overwritten with solution by dgtsv tridiagonal solver
    REAL(lp) :: rhs(size(act_depth_sl))
    INTEGER  :: info             ! dgtsv status variable

    REAL(lp) :: C_before, C_after

    nsoil = size(act_depth_sl)
    nsoil_loc = nsoil
    src_C_conc(:) = 0._lp

    C_before = SUM(cpool_bg_sl(:))

    IF (C_before > 0._lp) THEN
      C_conc(:) = 0._lp
      DO is = 1, nsoil
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_period_sl(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          C_conc(is) = cpool_bg_sl(is) / act_depth_sl(is)
        ELSE
          nsoil_loc = MAX(is - 1, 1)
          EXIT
        END IF
      END DO

      !------------------------------------------------------------
      ! build tridiagonal matrix and known right hand side

      ! Process the local layers
      DO is = 1, nsoil_loc

        ! Extract sensible values to use for concentration and
        ! thickness of the current layer, the layer above and below
        dz_here = dz(is)
        C_here = C_conc(is)
        IF (is == 1) THEN
          dz_above = 0._lp
          C_above = 0._lp
        ELSE
          dz_above = dz(is-1)
          C_above = C_conc(is-1)
        END IF
        IF (is == nsoil_loc) THEN
          dz_below = 0._lp
          C_below = 0._lp
        ELSE
          dz_below = dz(is+1)
          C_below = C_conc(is+1)
        END IF

        ! Determine diffusion coefficients to the layers above and below the current layer.
        IF (is == 1) THEN
          ! upwards diffusivity in uppermost layer is zero
          D_above = 0._lp
        ELSE
          ! for the layers below, upwards diffusivity is the weighted mean of the layer
          ! diffusivity in layer and layer above
          D_above = (dz(is-1)*diff_coeff(is-1) + &
            &        dz(is  )*diff_coeff(is))  / &
            &       (dz(is-1) + dz(is))
        END IF
        IF (is == nsoil_loc) THEN
          ! downwards diffusivity in the bottom layer is set to zero
          ! to impose no flux boundary condition
          D_below = 0._lp
        ELSE
          ! for the layers above, downwards diffusivity is the weighted mean of the
          ! diffusivity in the layer and the layer below
          D_below = (dz(is  )*diff_coeff(is  )  + &
            &        dz(is+1)*diff_coeff(is+1)) / &
            &       (dz(is) + dz(is+1))
        END IF

        ! Calculate layer-dependent weighting factors.
        dz_factor = 0.25_lp * (dz_below + 2._lp * dz_here + dz_above)
        W_below = D_below * (delta_time_C / dz_factor) * (2._lp / (dz_below + dz_here))
        W_above = D_above * (delta_time_C / dz_factor) * (2._lp / (dz_here + dz_above))

        ! boundary conditions for surface and bottom layers: no flux through top / bottom
        IF (is == 1) THEN
          W_above = 0._lp
        END IF
        IF (is == nsoil_loc) THEN
          W_below = 0._lp
        END IF

        ! Fill in matrix diagonal and off-diagonal elements.
        diag(is) = 1.0_lp + (implicit_contrib * (W_below + W_above))
        left(is) = -implicit_contrib * W_above
        right(is) = -implicit_contrib * W_below

        ! Calculate right hand side vector values.
        rhs(is) = (1.0_lp - (1.0_lp - implicit_contrib)*(W_below + W_above)) * C_here + &
          &        (1.0_lp - implicit_contrib)*W_below*C_below +                         &
          &        (1.0_lp - implicit_contrib)*W_above*C_above +                         &
          &        src_C_conc(is)
      END DO

      ! Semi-implicit case -- whenever implicit_contrib > 0.1,
      ! fully explicit case is solved already ( = rhs)
      IF (implicit_contrib >= 0.1_lp) THEN
        !------------------------------------------------------------
        ! solve tridiagonal system using dgtsv from NETLIB
        CALL dgtsv(nsoil_loc, 1, left(2:nsoil_loc), diag, &
          &        right(1:nsoil_loc-1), rhs, nsoil_loc, info)
  !       IF (PRESENT(lissue)) lissue = info /= 0
  !      IF (info /= 0) THEN
  !        WRITE(message_text,*) 'dgtsv info ', info
  !        CALL message(routine, message_text, all_print=.TRUE., level=em_warn)
  !      END IF
      END IF

      ! solution is ALWAYS written into the rhs vector, also in the fully explicit case
      C_conc = rhs

      DO is = 1, nsoil_loc
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_period_sl(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          cpool_bg_sl(is) = C_conc(is) * act_depth_sl(is)
        END IF
      END DO

      WHERE(cpool_bg_sl(:) < 0._lp)
        cpool_bg_sl(:) = 0._lp
      END WHERE

      C_after = SUM(cpool_bg_sl(:))

      ! re-normalise to conserve carbon
      IF (C_after > 0._lp) THEN
        cpool_bg_sl(:) = cpool_bg_sl(:)/C_after*C_before
      END IF
    END IF    ! C_before > 0

  END SUBROUTINE vertical_C_transport

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Computes the diffusivity for the vertical carbon transport through the soil column. Permafrost
  !> regions are dominated by cryoturbation and the diffusivity computation is limited to the
  !> active layer depth. In regions not affected by permafrost, the vertical transport is
  !> facilitated by bioturbation.
  !>
  SUBROUTINE set_diffusivity(nsoil, soil_depth, soil_depth_sl, thaw_depth, layer_depth, t_soil_sl, &
    &                        wtr_rel_sl, n_thaw_cycle, permafrost_period_sl, diff_coeff)

    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    INTEGER,  INTENT(in)    :: nsoil                 !< Number of soil layers in vertical yasso scheme []
    REAL(lp), INTENT(in)    :: soil_depth
    REAL(lp), INTENT(in)    :: soil_depth_sl(:)
    REAL(lp), INTENT(in)    :: thaw_depth
    REAL(lp), INTENT(in)    :: layer_depth(:)
    REAL(lp), INTENT(in)    :: t_soil_sl(:)
    REAL(lp), INTENT(in)    :: wtr_rel_sl(:)
    REAL(lp), INTENT(in)    :: n_thaw_cycle
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(out)   :: diff_coeff(:)           !< Vertical carbon mixing coefficient [m2 s-1]
    INTEGER :: is, isp1

    REAL(lp), PARAMETER :: cryo_max = 4.7e-12_lp !< Maximum diffusivity for cryoturbation [m^2 s-1]
    REAL(lp), PARAMETER :: bio_max  = 4.7e-13_lp !< Maximum diffusivity for bioturbation [m^2 s-1]

    diff_coeff(:) = 0._lp

    IF (thaw_depth < soil_depth) THEN    ! active layer depth < soil depth --> permafrost
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF (  (permafrost_period_sl(is)   < permafrost_threshold) .AND. &
          &   (permafrost_period_sl(isp1) < permafrost_threshold) .AND. &
          &   (t_soil_sl(is)          > tmelt)) THEN
          diff_coeff(is) = cryo_max *                               &
            &             MIN(1._lp, wtr_rel_sl(is)) *              &
            &              MIN(1._lp, (n_thaw_cycle / 40._lp))*     &
            &              (1._lp - layer_depth(is)/thaw_depth)
        ELSE
          diff_coeff(is) = 0._lp
        END IF
      END DO
    ELSE
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF ((soil_depth_sl(is) > minimum_layer_thickness) .AND. (t_soil_sl(is) > tmelt)) THEN
          diff_coeff(is) = MAX(0._lp, bio_max *(1._lp - (layer_depth(is)/3._lp)))
        ELSE
          diff_coeff(is) = 0._lp
        END IF
      END DO
    END IF

  END SUBROUTINE set_diffusivity

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> This subroutine moves carbon from the above ground litter pools to the below ground litter pools,
  !> if the above ground litter pools exceed a certain threshold. The threshold is defined as a fraction
  !> of the maximum carbon content in the above ground layer. The excess carbon is moved to the uppermost
  !> below ground layer. Similarily, if the below ground litter pools exceed a certain threshold, the
  !> excess carbon is moved to the next below ground layer.
  !>
    SUBROUTINE move_ag_2_bg(nsoil, max_depth, layer_thick, permafrost_period_sl, &
      &    c_acid_ag1,    c_water_ag1,    c_ethanol_ag1,    c_nonsoluble_ag1,    &
      &    c_acid_ag2,    c_water_ag2,    c_ethanol_ag2,    c_nonsoluble_ag2,    &
      &    c_acid_bg1_sl, c_water_bg1_sl, c_ethanol_bg1_sl, c_nonsoluble_bg1_sl, &
      &    c_acid_bg2_sl, c_water_bg2_sl, c_ethanol_bg2_sl, c_nonsoluble_bg2_sl, &
      &    c_humus_1_sl,  c_humus_2_sl)

    USE mo_carbon_constants, ONLY: C_2_biomass_ratio, rho_org_ag, rho_org_bg, molarMassC_kg, &
      &                            permafrost_threshold

    INTEGER,  INTENT(in)    :: nsoil                   !< Number of soil layers []
    REAL(lp), INTENT(in)    :: max_depth               !< Soil depth until bedrock or permafrost [m]
    REAL(lp), INTENT(in)    :: layer_thick(:)          !< Thickness of soil layers [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout) :: c_acid_ag1              !< Above-ground acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag1             !< Above-ground water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1           !< Above-ground ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1        !< Above-ground nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_ag2              !< Above-ground acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag2             !< Above-ground water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2           !< Above-ground ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2        !< Above-ground nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg1_sl(:)        !< Soil layer acid soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg1_sl(:)       !< Soil layer water soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1_sl(:)     !< Soil layer ethanol soluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1_sl(:)  !< Soil layer nonsoluble carbon pools of green litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg2_sl(:)        !< Soil layer acid soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg2_sl(:)       !< Soil layer water soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2_sl(:)     !< Soil layer ethanol soluble carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2_sl(:)  !< Soil layer nonsoluble carbon pools of green woody [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:)         !< Soil layer humus carbon pools of woody litter [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:)         !< Soil layer humus carbon pools of green woody [mol(C) m-2]

    REAL(lp),  PARAMETER :: thresh_frac_org    = 1.25_lp  !< Threshold organic fraction above which C is being moved downwards
    REAL(lp),  PARAMETER :: ag_layer_depth     = 0.13_lp  !< Assumed depth of above ground litter layer [m]
    REAL(lp) :: cpool_tot
    REAL(lp) :: depth_upper
    REAL(lp) :: rfrac
    REAL(lp) :: tmp

    REAL(lp) :: max_C_ag
    REAL(lp) :: max_C_bg(nsoil)
    INTEGER  :: is, isp1

    ! determine max C amount per layer
    ! C [kgC] = height [m] x A [m2] x rho [kg/m3] x C_2_bio [kgC/kg]
    ! finally converted to mol
    max_C_ag    = ag_layer_depth * rho_org_ag * C_2_biomass_ratio / molarMassC_kg
    max_C_bg(:) = layer_thick(:) * rho_org_bg * C_2_biomass_ratio / molarMassC_kg

    ! Check ag soil
    cpool_tot = c_acid_ag1 + c_acid_ag2       +      &
      &         c_water_ag1 + c_water_ag2     +      &
      &         c_ethanol_ag1 + c_ethanol_ag2 +      &
      &         c_nonsoluble_ag1 + c_nonsoluble_ag2

    IF (cpool_tot > max_C_ag) THEN
      ! Rescale above ground carbon pools such, that maximum amount is no longer exceeded
      ! and move excess carbon into the uppermost below ground layer.

      rfrac = max_C_ag / cpool_tot

      tmp                    = c_acid_ag1
      c_acid_ag1             = c_acid_ag1 * rfrac
      c_acid_bg1_sl(1)       = c_acid_bg1_sl(1) + (tmp - c_acid_ag1)
      tmp                    = c_acid_ag2
      c_acid_ag2             = c_acid_ag2 * rfrac
      c_acid_bg2_sl(1)       = c_acid_bg2_sl(1) + (tmp - c_acid_ag2)
      tmp                    = c_water_ag1
      c_water_ag1            = c_water_ag1 * rfrac
      c_water_bg1_sl(1)      = c_water_bg1_sl(1) + (tmp - c_water_ag1)
      tmp                    = c_water_ag2
      c_water_ag2            = c_water_ag2 * rfrac
      c_water_bg2_sl(1)      = c_water_bg2_sl(1) + (tmp - c_water_ag2)
      tmp                    = c_ethanol_ag1
      c_ethanol_ag1          = c_ethanol_ag1 * rfrac
      c_ethanol_bg1_sl(1)    = c_ethanol_bg1_sl(1) + (tmp - c_ethanol_ag1)
      tmp                    = c_ethanol_ag2
      c_ethanol_ag2          = c_ethanol_ag2 * rfrac
      c_ethanol_bg2_sl(1)    = c_ethanol_bg2_sl(1) + (tmp - c_ethanol_ag2)
      tmp                    = c_nonsoluble_ag1
      c_nonsoluble_ag1       = c_nonsoluble_ag1 * rfrac
      c_nonsoluble_bg1_sl(1) = c_nonsoluble_bg1_sl(1) + (tmp - c_nonsoluble_ag1)
      tmp                    = c_nonsoluble_ag2
      c_nonsoluble_ag2       = c_nonsoluble_ag2 * rfrac
      c_nonsoluble_bg2_sl(1) = c_nonsoluble_bg2_sl(1) + (tmp - c_nonsoluble_ag2)
    END IF

    ! check bg soil layers, move C downwards, if possible
    depth_upper = 0._lp
    DO is = 1,nsoil-1
      isp1 = MIN(is+1, nsoil)
      IF ((depth_upper < max_depth) .AND. (permafrost_period_sl(is) < permafrost_threshold) &
        & .AND. (permafrost_period_sl(isp1) < permafrost_threshold)) THEN
        depth_upper = depth_upper + layer_thick(is)
        cpool_tot = c_acid_bg1_sl(is)       + c_acid_bg2_sl(is)       +      &
          &         c_water_bg1_sl(is)      + c_water_bg2_sl(is)      +      &
          &         c_ethanol_bg1_sl(is)    + c_ethanol_bg2_sl(is)    +      &
          &         c_nonsoluble_bg1_sl(is) + c_nonsoluble_bg2_sl(is) +      &
          &         c_humus_1_sl(is)        + c_humus_2_sl(is)

        IF (cpool_tot > max_C_bg(is)) THEN
          ! Rescale this layers carbon pools such, that maximum amount is no longer exceeded
          ! and move excess carbon into the layer below.

          rfrac = max_C_bg(is) / cpool_tot

          tmp                       = c_acid_bg1_sl(is)
          c_acid_bg1_sl(is)         = c_acid_bg1_sl(is) * rfrac
          c_acid_bg1_sl(is+1)       = c_acid_bg1_sl(is+1) + (tmp - c_acid_bg1_sl(is))
          tmp                       = c_acid_bg2_sl(is)
          c_acid_bg2_sl(is)         = c_acid_bg2_sl(is) * rfrac
          c_acid_bg2_sl(is+1)       = c_acid_bg2_sl(is+1) + (tmp - c_acid_bg2_sl(is))
          tmp                       = c_water_bg1_sl(is)
          c_water_bg1_sl(is)        = c_water_bg1_sl(is) * rfrac
          c_water_bg1_sl(is+1)      = c_water_bg1_sl(is+1) + (tmp - c_water_bg1_sl(is))
          tmp                       = c_water_bg2_sl(is)
          c_water_bg2_sl(is)        = c_water_bg2_sl(is) * rfrac
          c_water_bg2_sl(is+1)      = c_water_bg2_sl(is+1) + (tmp - c_water_bg2_sl(is))
          tmp                       = c_ethanol_bg1_sl(is)
          c_ethanol_bg1_sl(is)      = c_ethanol_bg1_sl(is) * rfrac
          c_ethanol_bg1_sl(is+1)    = c_ethanol_bg1_sl(is+1) + (tmp - c_ethanol_bg1_sl(is))
          tmp                       = c_ethanol_bg2_sl(is)
          c_ethanol_bg2_sl(is)      = c_ethanol_bg2_sl(is) * rfrac
          c_ethanol_bg2_sl(is+1)    = c_ethanol_bg2_sl(is+1) + (tmp - c_ethanol_bg2_sl(is))
          tmp                       = c_nonsoluble_bg1_sl(is)
          c_nonsoluble_bg1_sl(is)   = c_nonsoluble_bg1_sl(is) * rfrac
          c_nonsoluble_bg1_sl(is+1) = c_nonsoluble_bg1_sl(is+1) + (tmp - c_nonsoluble_bg1_sl(is))
          tmp                       = c_nonsoluble_bg2_sl(is)
          c_nonsoluble_bg2_sl(is)   = c_nonsoluble_bg2_sl(is) * rfrac
          c_nonsoluble_bg2_sl(is+1) = c_nonsoluble_bg2_sl(is+1) + (tmp - c_nonsoluble_bg2_sl(is))
          tmp                       = c_humus_1_sl(is)
          c_humus_1_sl(is)          = c_humus_1_sl(is) * rfrac
          c_humus_1_sl(is+1)        = c_humus_1_sl(is+1) + (tmp - c_humus_1_sl(is))
          tmp                       = c_humus_2_sl(is)
          c_humus_2_sl(is)          = c_humus_2_sl(is) * rfrac
          c_humus_2_sl(is+1)        = c_humus_2_sl(is+1) + (tmp - c_humus_2_sl(is))
        END IF
      END IF
    END DO

  END SUBROUTINE move_ag_2_bg

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Remove any carbon that is below the bedrock depth and move it to the lowest layer above bedrock.
  !> This should not happen in reality, but can occur due to numerical issues.
  !> This routine is a safeguard against that.
  !>
  SUBROUTINE remove_C_below_bedrock(soil_depth, layer_thick, c_pool_bg)
    REAL(lp), INTENT(in)     :: soil_depth
    REAL(lp), INTENT(in)     :: layer_thick(:)
    REAL(lp), INTENT(inout)  :: c_pool_bg(:)

    INTEGER  :: is, nsoil, max_layer_above
    REAL(lp) :: soil_above, C_below_bedrock, upper_bound

    nsoil = size(layer_thick)

    soil_above = 0._lp
    C_below_bedrock = 0._lp
    max_layer_above = 0
    upper_bound = 0._lp

    ! Check for C below bedrock depth
    DO is = 1, nsoil
      IF (soil_depth > upper_bound) THEN
        max_layer_above = is
        soil_above = soil_above + layer_thick(is)
      ELSE
        IF (c_pool_bg(is) > 0._lp) THEN
          C_below_bedrock = C_below_bedrock + c_pool_bg(is)
        END IF
      END IF
      upper_bound = upper_bound + layer_thick(is)
    END DO

    ! Move C to lowest layer above bedrock, the lowest layer with a soil fraction
    IF ((C_below_bedrock > 0._lp) .AND. (soil_above > 0._lp) .AND. (max_layer_above > 0)) THEN
      c_pool_bg(max_layer_above) = c_pool_bg(max_layer_above) + C_below_bedrock
      DO is = MIN(max_layer_above+1,nsoil), nsoil
        c_pool_bg(is) = 0._lp
      END DO
    END IF

  END SUBROUTINE remove_C_below_bedrock

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Sum up layered carbon storages for single below ground pools
  !>
  !> This routine accumulates the vertically distinct carbon fields of the multi-layer
  !> yasso scheme over all soil layers to fill the arrays used for the single
  !> layer yasso.
  !>
  SUBROUTINE sum_up_carbon_layers(nc, nsoil, cpool_sl, cpool_bg)

    INTEGER,  INTENT(in)     :: nc            !< Number of grid cells [ ]
    INTEGER,  INTENT(in)     :: nsoil         !< Number of soil layers [ ]
    REAL(lp), INTENT(in)     :: cpool_sl(:,:) !< Soil carbon on layers [mol(C) m-2]
    REAL(lp), INTENT(out)    :: cpool_bg(:)   !< Below ground soil carbon pool [mol(C) m-2]

    INTEGER  :: ic                            !< Grid cell index
    INTEGER  :: is                            !< Soil layer index

    ! First set below ground carbon pool to zero
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      cpool_bg(ic) = 0._lp
    END DO
    !$ACC END PARALLEL LOOP

    ! Add carbon from each layer to below ground pool
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO ic = 1, nc
      DO is = 1, nsoil
        cpool_bg(ic) = cpool_bg(ic) + cpool_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE sum_up_carbon_layers
#endif
END MODULE mo_carbon_process
