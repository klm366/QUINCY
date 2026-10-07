!> Contains carbon process routines for vertically stratified version of YASSO
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
MODULE mo_carbon_process_vertical
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_vertical_soil_c, depth_to_layers, relocate_carbon_vertical_fire, relocate_carbon_vertical_damage, &
    &       determine_root_dist

  CHARACTER(*), PARAMETER :: modname = 'mo_carbon_process_vertical'

CONTAINS

  SUBROUTINE calc_vertical_soil_c(            &
       ! General inputs
    &  delta_time_C, nsoil,                   &
    &  layer_depth, layer_thick,              &
       ! Input coefficients (fixed)
    &  LeafLit_coefV, WoodLit_coefV,          &
    &  MiscPar,                               &
       ! Climate inputs
    &  temp,                                  &
    &  precip,                                &
    &  t_soil_sl,                             &
    &  wtr_soil_rel_sl,                       &
    &  n_thaw_cycle,                          &
       ! Topmodel inputs
    &  frac_inun,                             &
       ! Soil inputs
    &  soil_depth,                            &
    &  soil_depth_sl,                         &
    &  max_depth,                             &
    &  max_depth_sl,                          &
    &  thaw_depth,                            &
    &  root_dist_sl,                          &
    &  permafrost_layer,                      &
       ! Litter inputs
    &  woodLitter,                            &
    &  leafLitter,                            &
    &  rootExudates,                          &
       ! Yasso pools. Used on input and updated on output
    &  c_acid_ag1,                            &
    &  c_water_ag1,                           &
    &  c_ethanol_ag1,                         &
    &  c_nonsoluble_ag1,                      &
    &  c_acid_ag2,                            &
    &  c_water_ag2,                           &
    &  c_ethanol_ag2,                         &
    &  c_nonsoluble_ag2,                      &
    &  c_acid_bg1_sl,                         &
    &  c_water_bg1_sl,                        &
    &  c_ethanol_bg1_sl,                      &
    &  c_nonsoluble_bg1_sl,                   &
    &  c_humus_1_sl,                          &
    &  c_acid_bg2_sl,                         &
    &  c_water_bg2_sl,                        &
    &  c_ethanol_bg2_sl,                      &
    &  c_nonsoluble_bg2_sl,                   &
    &  c_humus_2_sl,                          &
       ! Output fluxes
    &  c_decomp_humus_1_sum,                  &
    &  c_decomp_humus_2_sum,                  &
    &  c_into_humus_1_sum,                    &
    &  c_into_humus_2_sum,                    &
    &  soilResp_oxic_sl,                      &
    &  soilResp_anox_sl                       &
    &  )

    USE mo_carbon_constants,       ONLY: fract_green_aboveGround, fract_wood_aboveGround,    &
      &                                  LCT_HUMUS, F_RESP, F_2_HUMUS, F_F_HUMUS,            &
      &                                  M_LITTER, M_HUMUS,                                  &
      &                                  WOODLITSIZE, ROOT_DIST_BETA,                        &
      &                                  d_oxic, d_anox, d_size_leaf
    USE mo_carbon_shared,          ONLY: distribute_yasso_litter

    REAL(wp),  INTENT(in)    :: delta_time_C             !< Carbon time step [s]
    INTEGER,   INTENT(in)    :: nsoil                    !
    REAL(wp),  INTENT(in)    :: layer_depth(:)
    REAL(wp),  INTENT(in)    :: layer_thick(:)
    !
    ! Fixed input parameters for Yasso
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp),  INTENT(in)    :: LeafLit_coefV(:)
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp),  INTENT(in)    :: WoodLit_coefV(:)
    !   Other constants
    REAL(wp),  INTENT(in)    :: MiscPar(:)

    ! Carbon inputs [mol m-2 s-1]:
    REAL(wp),  INTENT(in)    :: leafLitter
    REAL(wp),  INTENT(in)    :: woodLitter
    REAL(wp),  INTENT(in)    :: rootExudates

    ! Input meteorology for Yasso
    REAL(wp),  INTENT(in)    :: temp                     !< multi-day mean temperature [K]
    REAL(wp),  INTENT(in)    :: precip                   !< multi-day mean precipitation [mm s-1]
    REAL(wp),  INTENT(in)    :: t_soil_sl(:)             !< soil layer temperatures
    REAL(wp),  INTENT(in)    :: wtr_soil_rel_sl(:)       !< relative soil moisture
    REAL(wp),  INTENT(in)    :: n_thaw_cycle             !< number of freeze-thaw cycles in previous year

    ! Input parameters from topmodel (or other wetland scheme)
    REAL(wp),  INTENT(in)    :: frac_inun                !< Fraction of tile inundated []

    ! Input soil and / or active layer depth
    REAL(wp),  INTENT(in)    :: soil_depth               !< Total soil depth until bedrock [m]
    REAL(wp),  INTENT(in)    :: soil_depth_sl(:)         !< Soil depth within each layer [m]
    REAL(wp),  INTENT(in)    :: max_depth                !< maximum available depth [m], either soil or thaw
    REAL(wp),  INTENT(in)    :: max_depth_sl(:)          !< maximum available depth [m], either soil or thaw, per layer
    REAL(wp),  INTENT(in)    :: thaw_depth               !< thaw depth [m]
    REAL(wp),  INTENT(inout) :: root_dist_sl(:)          !< root distribution, per layer
    REAL(wp),  INTENT(in)    :: permafrost_layer(:)      !< permafrost status of soil layer

    ! Yasso pools [mol(C) m-2(canopy)], input and updated on output
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp),  INTENT(inout) :: c_acid_ag1               !< Acid soluble litter pool
    REAL(wp),  INTENT(inout) :: c_water_ag1              !< Water soluble litter pool
    REAL(wp),  INTENT(inout) :: c_ethanol_ag1            !< Ethanol soluble litter pool
    REAL(wp),  INTENT(inout) :: c_nonsoluble_ag1         !< Non-soluble litter pool
    !   Belowground
    REAL(wp),  INTENT(inout) :: c_acid_bg1_sl(:)         !< Acid soluble litter pool
    REAL(wp),  INTENT(inout) :: c_water_bg1_sl(:)        !< Water soluble litter pool
    REAL(wp),  INTENT(inout) :: c_ethanol_bg1_sl(:)      !< Ethanol soluble litter pool
    REAL(wp),  INTENT(inout) :: c_nonsoluble_bg1_sl(:)   !< Non-soluble litter pool
    REAL(wp),  INTENT(inout) :: c_humus_1_sl(:)          !< Humus litter pool
    ! Size class 2; woody litter
    !   Aboveground
    REAL(wp),  INTENT(inout) :: c_acid_ag2               !< Acid soluble litter pool
    REAL(wp),  INTENT(inout) :: c_water_ag2              !< Water soluble litter pool
    REAL(wp),  INTENT(inout) :: c_ethanol_ag2            !< Ethanol soluble litter pool
    REAL(wp),  INTENT(inout) :: c_nonsoluble_ag2         !< Non-soluble litter pool
    !   Belowground
    REAL(wp),  INTENT(inout) :: c_acid_bg2_sl(:)         !< Acid soluble litter pool
    REAL(wp),  INTENT(inout) :: c_water_bg2_sl(:)        !< Water soluble litter pool
    REAL(wp),  INTENT(inout) :: c_ethanol_bg2_sl(:)      !< Ethanol soluble litter pool
    REAL(wp),  INTENT(inout) :: c_nonsoluble_bg2_sl(:)   !< Non-soluble litter pool
    REAL(wp),  INTENT(inout) :: c_humus_2_sl(:)          !< Humus litter pool

    ! Decomposition fluxes
    REAL(wp),  INTENT(inout) :: c_decomp_humus_1_sum     !< annual cflux sum of humus decomposed (leaf)
    REAL(wp),  INTENT(inout) :: c_decomp_humus_2_sum     !< annual cflux sum of humus decomposed (wood)
    REAL(wp),  INTENT(inout) :: c_into_humus_1_sum       !< annual cflux sum into humus (leaf)
    REAL(wp),  INTENT(inout) :: c_into_humus_2_sum       !< annual cflux sum into humus (wood)

    REAL(wp),  INTENT(out)   :: soilResp_oxic_sl(:)      !< soil respiration rate, oxic part
    REAL(wp),  INTENT(out)   :: soilResp_anox_sl(:)      !< soil respiration rate, anoxic part

    ! decomposition modifier terms
    REAL(wp)     :: d_temp                  !< term which accounts for the temperature influence on decomposition
    REAL(wp)     :: d_water                 !< term which accounts for the water influence on decomposition
    REAL(wp)     :: d_size_wood             !< term which accounts for the litter size influence on decomposition (woody litter)
    REAL(wp)     :: d_temp_sl(nsoil)        !< term which accounts for the temperature influence on decomposition with depth
    REAL(wp)     :: d_water_sl(nsoil)       !< term which accounts for the water influence on decomposition with depth
    REAL(wp)     :: d_depth_sl(nsoil)       !< term which accounts for the depth influence on decomposition


    REAL(wp) :: frac_inun_loc                   !< local frac_inun -- can be changed without problems

    ! transfer to YASSO subroutine
    REAL(wp), DIMENSION(5)  :: Y_in_pools       !< Yasso pools IN [mol(c)/m2]
    REAL(wp), DIMENSION(2)  :: Y_in_mod         !< Yasso rate modifiers IN
    REAL(wp), DIMENSION(5)  :: Y_out_pools      !< Yasso pools OUT [mol(c)/m2]
    REAL(wp), DIMENSION(4)  :: Y_out_fluxes     !< Yasso fluxes OUT [mol(c)/m2]

    REAL(wp) :: diff_coeff(nsoil)

    REAL(wp) :: tmp_out_humus
    INTEGER  :: is

    REAL(wp) :: soilResp_oxic                   !< ag respiration flux oxic
    REAL(wp) :: soilResp_anox                   !< ag respiration flux anoxic
    CHARACTER(len=*), PARAMETER  :: routine = modname//':calc_vertical_soil_c'


    ! ------------------------------------------------------------
    soilResp_oxic_sl(:)   = 0._wp
    soilResp_anox_sl(:)   = 0._wp

    root_dist_sl = determine_root_dist(max_depth_sl, layer_depth, MiscPar(ROOT_DIST_BETA), t_soil_sl, permafrost_layer)

    frac_inun_loc = frac_inun

    ! ------------------------------------------------------------
    ! Formerly within YASSO

    ! ----------------------------------------------------------------------------------------
    ! YASSO calls
    ! ----------------------------------------------------------------------------------------

    ! setup: initialise variables
    tmp_out_humus       = 0._wp
    soilResp_oxic       = 0._wp
    soilResp_anox       = 0._wp
    soilResp_oxic_sl(:) = 0._wp
    soilResp_anox_sl(:) = 0._wp

    ! get decomposition modifiers
    d_temp      = get_modifier_temp(temp)
    d_water     = get_modifier_precip(precip)
    d_temp_sl   = get_modifier_temp(t_soil_sl)
    d_water_sl  = get_modifier_water(wtr_soil_rel_sl)
    d_depth_sl  = get_modifier_depth(layer_depth)
    d_size_wood = get_modifier_size(MiscPar(WOODLITSIZE))

    ! -------------------------------
    ! call YASSO for dry / oxic areas
    !

    ! LEAF LITTER
    ! ag pools

    ! rate modifiers YASSO call for oxic
    Y_in_mod(M_LITTER) = d_temp * d_water * d_size_leaf * d_oxic
    Y_in_mod(M_HUMUS)  = 0._wp   ! no humus above ground

    ! C pools input
    CALL input_YASSO_pools(c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
      &                    c_nonsoluble_ag1, 0._wp, Y_in_pools)
    ! call YASSO for oxic
    CALL yasso_minimal(delta_time_C, Y_in_pools * (1._wp - frac_inun_loc), Y_in_mod, Y_out_pools, Y_out_fluxes)
    ! C pools output
    CALL output_YASSO_pools(c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
      &                     c_nonsoluble_ag1, tmp_out_humus, .false., Y_out_pools)
    ! Fluxes
    soilResp_oxic              = soilResp_oxic      + Y_out_fluxes(F_RESP) / delta_time_C
    c_into_humus_1_sum         = c_into_humus_1_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus

    ! call YASSO for inundated / anoxic
    IF (frac_inun_loc > EPSILON(1._wp)) THEN
      Y_in_mod(M_LITTER) = d_temp * d_water * d_size_leaf * d_anox

      ! call YASSO for anoxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, Y_out_pools, Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
        &                     c_nonsoluble_ag1, tmp_out_humus, .true., Y_out_pools)
      ! Fluxes
      soilResp_anox            = soilResp_anox + Y_out_fluxes(F_RESP) / delta_time_C
      c_into_humus_1_sum       = c_into_humus_1_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
    END IF

    ! -----------------
    ! YASSO bg1 calls: All soil layers, leaf litter

    DO is = 1, nsoil
      ! rate modifiers YASSO call for oxic
      Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_leaf * d_oxic
      Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_oxic

      ! C pools input
      CALL input_YASSO_pools(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is),  &
        &                    c_nonsoluble_bg1_sl(is), c_humus_1_sl(is), Y_in_pools)
      ! call YASSO for oxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * (1._wp - frac_inun_loc), Y_in_mod, Y_out_pools, Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is), &
        &                     c_nonsoluble_bg1_sl(is), c_humus_1_sl(is), .false., Y_out_pools)
      ! Fluxes
      soilResp_oxic_sl(is)     = soilResp_oxic_sl(is)  + Y_out_fluxes(F_RESP) / delta_time_C     ! respiration
      c_decomp_humus_1_sum     = c_decomp_humus_1_sum  + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
      c_into_humus_1_sum       = c_into_humus_1_sum    + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus

      ! call YASSO for inundated / anoxic
      IF (frac_inun_loc > EPSILON(1._wp)) THEN
        Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_leaf * d_anox
        Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_anox

        ! call YASSO for anoxic
        CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, Y_out_pools, Y_out_fluxes)
        ! C pools output
        CALL output_YASSO_pools(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is), &
          &                     c_nonsoluble_bg1_sl(is), c_humus_1_sl(is), .true., Y_out_pools)
        ! Fluxes
        soilResp_anox_sl(is)   = soilResp_anox_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
        c_decomp_humus_1_sum   = c_decomp_humus_1_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
        c_into_humus_1_sum     = c_into_humus_1_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
      END IF
    END DO
    ! some C goes from ag to humus -- accounted for here
    c_humus_1_sl(1) = c_humus_1_sl(1) + tmp_out_humus
    tmp_out_humus   = 0._wp

    ! -----------------
    ! WOODY LITTER
    ! ag pools

    ! rate modifiers YASSO call for oxic
    Y_in_mod(M_LITTER) = d_temp * d_water * d_size_wood * d_oxic
    Y_in_mod(M_HUMUS) = 0._wp   ! no humus above ground

    ! C pools input
    CALL input_YASSO_pools(c_acid_ag2, c_water_ag2, c_ethanol_ag2,  &
      &                    c_nonsoluble_ag2, 0._wp, Y_in_pools)
    ! call YASSO for oxic
    CALL yasso_minimal(delta_time_C, Y_in_pools * (1._wp - frac_inun_loc), Y_in_mod, Y_out_pools, Y_out_fluxes)
    ! C pools output
    CALL output_YASSO_pools(c_acid_ag2, c_water_ag2, c_ethanol_ag2, &
      &                     c_nonsoluble_ag2, tmp_out_humus, .false., Y_out_pools)
    ! Fluxes
    soilResp_oxic              = soilResp_oxic + Y_out_fluxes(F_RESP) / delta_time_C
    c_into_humus_2_sum         = c_into_humus_2_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus

    IF (frac_inun_loc > EPSILON(1._wp)) THEN
      ! call YASSO for inundated / anoxic
      Y_in_mod(M_LITTER) = d_temp * d_water * d_size_wood * d_anox

      ! call YASSO for anoxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, Y_out_pools, Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(c_acid_ag2, c_water_ag2, c_ethanol_ag2, &
        &                     c_nonsoluble_ag2, tmp_out_humus, .true., Y_out_pools)
      ! Fluxes
      soilResp_anox            = soilResp_anox + Y_out_fluxes(F_RESP) / delta_time_C
      c_into_humus_2_sum       = c_into_humus_2_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus
    END IF

    ! -----------------
    ! YASSO bg call: All layers
    DO is = 1, nsoil
      ! rate modifiers, YASSO call for oxic
      Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_wood * d_oxic
      Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_oxic

      ! C pools input
      CALL input_YASSO_pools(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is), &
        &                    c_nonsoluble_bg2_sl(is), c_humus_2_sl(is), Y_in_pools)
      ! call YASSO for oxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * (1._wp - frac_inun_loc), Y_in_mod, Y_out_pools, Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is), &
        &                     c_nonsoluble_bg2_sl(is), c_humus_2_sl(is), .false., Y_out_pools)
      ! Fluxes
      soilResp_oxic_sl(is)     = soilResp_oxic_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
      c_decomp_humus_2_sum     = c_decomp_humus_2_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
      c_into_humus_2_sum       = c_into_humus_2_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus

      IF (frac_inun_loc > EPSILON(1._wp)) THEN
        ! call YASSO for inundated / anoxic
        Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_wood * d_anox
        Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_anox

        ! call YASSO for anoxic
        CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, Y_out_pools, Y_out_fluxes)

        ! C pools output
        CALL output_YASSO_pools(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is), &
          &                     c_nonsoluble_bg2_sl(is), c_humus_2_sl(is), .true., Y_out_pools)

        ! Fluxes
        soilResp_anox_sl(is)   = soilResp_anox_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
        c_decomp_humus_2_sum   = c_decomp_humus_2_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
        c_into_humus_2_sum     = c_into_humus_2_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
      END IF
    END DO
    ! some C goes from ag to humus -- accounted for here
    c_humus_2_sl(1) = c_humus_2_sl(1) + tmp_out_humus
    tmp_out_humus   = 0._wp

    c_into_humus_1_sum   = c_into_humus_1_sum + &
      &                    LeafLit_coefV(LCT_HUMUS) * leafLitter * delta_time_C
    c_into_humus_2_sum   = c_into_humus_2_sum + &
      &                    WoodLit_coefV(LCT_HUMUS) * woodLitter * delta_time_C

    ! soil respiration: ag respiration is really part of top soil layer respiration
    soilResp_oxic_sl(1)     = soilResp_oxic_sl(1) + soilResp_oxic
    soilResp_anox_sl(1)     = soilResp_anox_sl(1) + soilResp_anox

    ! ----------------------------------------------------------------------------------------
    ! 3b. distribute litter and root exudates
    ! a) leaf litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, c_humus_1_sl(1),  &
      &                          leafLitter*delta_time_C, fract_green_aboveGround, LeafLit_coefV)

    ! bg
    DO is = 1, nsoil
      CALL distribute_yasso_litter(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is),              &
        &                          c_nonsoluble_bg1_sl(is), c_humus_1_sl(is),                                &
        &                          leafLitter*delta_time_C, (1.0_wp - fract_green_aboveGround) * root_dist_sl(is),        &
        &                          LeafLit_coefV)
      ! root exudates are carbohydrates, go directly into water_bg pools
      c_water_bg1_sl(is) = c_water_bg1_sl(is) + rootExudates * delta_time_C * root_dist_sl(is)
    END DO

    !
    ! b) woody litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, c_humus_2_sl(1),  &
      &                          woodLitter*delta_time_C, fract_wood_aboveGround, WoodLit_coefV)

    ! bg
    DO is = 1, nsoil
      CALL distribute_yasso_litter(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is),              &
        &                          c_nonsoluble_bg2_sl(is), c_humus_2_sl(is),                                &
        &                          woodLitter*delta_time_C, (1.0_wp - fract_wood_aboveGround) * root_dist_sl(is),         &
        &                          WoodLit_coefV)
    END DO

    ! ----------------------------------------------------------------------------------------
    ! vertical mixing

    CALL move_ag_2_bg(c_acid_ag1, c_acid_bg1_sl, c_acid_ag2, c_acid_bg2_sl,                         &
      &               c_water_ag1, c_water_bg1_sl, c_water_ag2, c_water_bg2_sl,                     &
      &               c_ethanol_ag1, c_ethanol_bg1_sl, c_ethanol_ag2, c_ethanol_bg2_sl,             &
      &               c_nonsoluble_ag1, c_nonsoluble_bg1_sl, c_nonsoluble_ag2, c_nonsoluble_bg2_sl, &
      &               c_humus_1_sl, c_humus_2_sl, max_depth, layer_thick, permafrost_layer)

    CALL set_diffusivity(diff_coeff, soil_depth, soil_depth_sl, thaw_depth, &
      &                 layer_depth, t_soil_sl, wtr_soil_rel_sl,            &
      &                 n_thaw_cycle, permafrost_layer)
    CALL vertical_C_transport(c_acid_bg1_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_acid_bg2_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_water_bg1_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_water_bg2_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_ethanol_bg1_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_ethanol_bg2_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_nonsoluble_bg1_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_nonsoluble_bg2_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_humus_1_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)
    CALL vertical_C_transport(c_humus_2_sl, max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_layer)

    ! remove any C that might be below bedrock boundary
    CALL remove_C_below_bedrock(c_acid_bg1_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_acid_bg2_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_water_bg1_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_water_bg2_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_ethanol_bg1_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_ethanol_bg2_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_nonsoluble_bg1_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_nonsoluble_bg2_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_humus_1_sl, soil_depth, layer_thick)
    CALL remove_C_below_bedrock(c_humus_2_sl, soil_depth, layer_thick)

  END SUBROUTINE calc_vertical_soil_c


  ! ---------------------------------------------------------------------------------------------------------

  SUBROUTINE yasso_minimal(delta_time_C, Y_in_pools, Y_in_mod, Y_out_pools, Y_out_fluxes,                   &
    &                      A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H)

    USE mo_carbon_constants, ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS, &
      &                            F_RESP, F_2_HUMUS, F_F_HUMUS, F_LITRESP,           &
      &                            M_LITTER, M_HUMUS, sec_per_year

    REAL(wp), INTENT(in)                 :: delta_time_C     !< Carbon time step [s]
    REAL(wp), DIMENSION(5),  INTENT(IN)  :: Y_in_pools       !< Yasso pools IN [mol(c)/m2]
    REAL(wp), DIMENSION(2),  INTENT(IN)  :: Y_in_mod         !< Yasso rate modifiers IN
    REAL(wp), DIMENSION(5),  INTENT(OUT) :: Y_out_pools      !< Yasso pools OUT [mol(c)/m2]
    REAL(wp), DIMENSION(4),  INTENT(OUT) :: Y_out_fluxes     !< Yasso fluxes OUT [mol(c)/m2]

    ! Possibilities to overwrite (some of) the default exchange coefficients (AWEN2H for backward compatibility)
    REAL(wp), OPTIONAL,      INTENT(in)  :: A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(wp)     :: c_acid
    REAL(wp)     :: c_water
    REAL(wp)     :: c_ethanol
    REAL(wp)     :: c_nonsoluble
    REAL(wp)     :: c_humus

    ! Decomposition rates [1/a]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus

    ! Carbon fluxes [mol(C) m-2 yr-1] - Unit indicated by names Cflx_yr (Carbon annual flux)
    REAL(wp)     :: Cflx_yr_from_acid              !< Loss flux of acid soluble pool        [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_from_water             !< Loss flux of water soluble pool       [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_from_ethanol           !< Loss flux of ethanol soluble pool     [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_from_nonsoluble        !< Loss flux of non soluble pool         [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_from_humus             !< Loss flux of humus pool               [mol(c)/m2/yr]

    REAL(wp)     :: Cflx_yr_2_acid                 !< Gain flux of acid soluble pool        [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_2_water                !< Gain flux of water soluble pool       [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_2_ethanol              !< Gain flux of ethanol soluble pool     [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_2_nonsoluble           !< Gain flux of non soluble pool         [mol(c)/m2/yr]
    REAL(wp)     :: Cflx_yr_2_humus                !< Gain flux of humus pool               [mol(c)/m2/yr]

    REAL(wp)     :: Cflx_yr_resp_from_litter
    REAL(wp)     :: Cflx_yr_respiration

    ! Fractions to split the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
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
    ! REAL(wp)   :: AWEN_2_H
    ! Alternative pathway for N pool which enters the humus pool [ ]
    REAL(wp)   :: N_2_H
    REAL(wp)   :: AWE_2_H

    ! Yasso parameters (Tuomi '11)
    !   reference decompositation rates for T = 0C and unlimited water availability in [yr-1]...
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !<  ... acid-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !<  ... water-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !<  ... ethanol-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !<  ... nonsoluble-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !<  ... humus

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates.
    REAL(wp)     :: years_per_ts, ts_per_year

    ts_per_year = REAL(sec_per_year,wp) / delta_time_C
    years_per_ts = 1.0_wp/ts_per_year

    ! Treatment of optional parameters
    !   Set exchange coefficients
    IF (PRESENT(AWEN2H)) THEN
      AWE_2_H = AWEN2H
      N_2_H   = AWEN2H
    END IF

    ! Default values in [yr-1]
    AWE_2_H  = 0.0045_wp; IF (PRESENT(AWE2H))  AWE_2_H  = AWE2H
    N_2_H    = 0.0045_wp; IF (PRESENT(  N2H))    N_2_H  =   N2H
    A_2_W    = 0.99_wp;   IF (PRESENT(  A2W))    A_2_W  =   A2W
    A_2_E    = 0.00_wp;   IF (PRESENT(  A2E))    A_2_E  =   A2E
    A_2_N    = 0.00_wp;   IF (PRESENT(  A2N))    A_2_N  =   A2N
    W_2_A    = 0.48_wp;   IF (PRESENT(  W2A))    W_2_A  =   W2A
    W_2_E    = 0.00_wp;   IF (PRESENT(  W2E))    W_2_E  =   W2E
    W_2_N    = 0.015_wp;  IF (PRESENT(  W2N))    W_2_N  =   W2N
    E_2_A    = 0.01_wp;   IF (PRESENT(  E2A))    E_2_A  =   E2A
    E_2_W    = 0.00_wp;   IF (PRESENT(  E2W))    E_2_W  =   E2W
    E_2_N    = 0.95_wp;   IF (PRESENT(  E2N))    E_2_N  =   E2N
    N_2_A    = 0.83_wp;   IF (PRESENT(  N2A))    N_2_A  =   N2A
    N_2_W    = 0.01_wp;   IF (PRESENT(  N2W))    N_2_W  =   N2W
    N_2_E    = 0.02_wp;   IF (PRESENT(  N2E))    N_2_E  =   N2E

    ! INPUTS
    ! a) pools
    c_acid       = Y_in_pools(Y_ACID)
    c_water      = Y_in_pools(Y_WATER)
    c_ethanol    = Y_in_pools(Y_ETHANOL)
    c_nonsoluble = Y_in_pools(Y_NONSOLUBLE)
    c_humus      = Y_in_pools(Y_HUMUS)

    ! b) decomposition modifiers
    ! decomposition rate is reference decomposition rate * modifier
    ! with modifier function(temperature, water, size, etc) set externally
    d_acid       = ref_decomp_rate_acid       * Y_in_mod(M_LITTER)
    d_water      = ref_decomp_rate_water      * Y_in_mod(M_LITTER)
    d_ethanol    = ref_decomp_rate_ethanol    * Y_in_mod(M_LITTER)
    d_nonsoluble = ref_decomp_rate_nonsoluble * Y_in_mod(M_LITTER)
    d_humus      = ref_decomp_rate_humus      * Y_in_mod(M_HUMUS)

    ! Calculate fluxes
    !   annual loss fluxes (negative values):
    Cflx_yr_from_acid         = c_acid       * d_acid
    Cflx_yr_from_water        = c_water      * d_water
    Cflx_yr_from_ethanol      = c_ethanol    * d_ethanol
    Cflx_yr_from_nonsoluble   = c_nonsoluble * d_nonsoluble
    Cflx_yr_from_humus        = c_humus      * d_humus

    !   annual gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_yr_2_acid            = ABS(Cflx_yr_from_water      * W_2_A  &    ! returns positive fluxes
      &                           + Cflx_yr_from_ethanol    * E_2_A  &
      &                           + Cflx_yr_from_nonsoluble * N_2_A)

    Cflx_yr_2_water           = ABS(Cflx_yr_from_acid       * A_2_W  &
      &                            + Cflx_yr_from_ethanol    * E_2_W &
      &                            + Cflx_yr_from_nonsoluble * N_2_W)

    Cflx_yr_2_ethanol         = ABS(Cflx_yr_from_acid       * A_2_E  &
      &                            + Cflx_yr_from_water      * W_2_E &
      &                            + Cflx_yr_from_nonsoluble * N_2_E)

    Cflx_yr_2_nonsoluble      = ABS(Cflx_yr_from_acid       * A_2_N  &
      &                           + Cflx_yr_from_water      * W_2_N  &
      &                           + Cflx_yr_from_ethanol    * E_2_N)

    Cflx_yr_2_humus           = ABS(Cflx_yr_from_acid                &
      &                           + Cflx_yr_from_water               &
      &                           + Cflx_yr_from_ethanol             &
      &                           ) * AWE_2_H                        &
      &                             +                                &
      &                         ABS(Cflx_yr_from_nonsoluble)         &
      &                             * N_2_H

    !   litter ag & bg respiration (needed for N cycle)
    Cflx_yr_resp_from_litter =  Cflx_yr_from_acid       * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_water      * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_ethanol    * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_nonsoluble * (1.0_wp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_respiration      =  Cflx_yr_resp_from_litter                      &
      &                       + Cflx_yr_from_humus

    ! 3. update Yasso pools

    !   3a. redistribution among yasso pools
    c_acid        = MAX(0.0_wp,                                    &
      &                 c_acid                                     & ! old pool
      &              + (Cflx_yr_from_acid                          & ! incoming flux from AWEN pools
      &               + Cflx_yr_2_acid) * years_per_ts)
    c_water       = MAX(0.0_wp,                                    & ! ...
      &                 c_water                                    &
      &              + (Cflx_yr_from_water                         &
      &               + Cflx_yr_2_water) * years_per_ts)
    c_ethanol     = MAX(0.0_wp,                                    &
      &                 c_ethanol                                  &
      &              + (Cflx_yr_from_ethanol                       &
      &               + Cflx_yr_2_ethanol) * years_per_ts)
    c_nonsoluble  = MAX(0.0_wp,                                    &
      &                 c_nonsoluble                               &
      &              + (Cflx_yr_from_nonsoluble                    &
      &               + Cflx_yr_2_nonsoluble) * years_per_ts)

    c_humus       = MAX(0.0_wp,                                    &
      &                 c_humus                                    &
      &              + (Cflx_yr_from_humus                         &
      &               + Cflx_yr_2_humus) * years_per_ts)

    ! 4.1 Write pools into the output array
    Y_out_pools(Y_ACID)       = c_acid
    Y_out_pools(Y_WATER)      = c_water
    Y_out_pools(Y_ETHANOL)    = c_ethanol
    Y_out_pools(Y_NONSOLUBLE) = c_nonsoluble
    Y_out_pools(Y_HUMUS)      = c_humus

    ! 4.2 Write respiration into the output array
    Y_out_fluxes(F_RESP)     = Cflx_yr_respiration * years_per_ts ! convert back to daily

    ! the routine might have been called to diagnose the nutrient demand, so we give out
    ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
    Y_out_fluxes(F_2_HUMUS)  =  Cflx_yr_2_humus * years_per_ts !  per timestep ! the flux is positive
    ! 2. flux for mineralisation:  this is the flux from humus to atmosphere
    Y_out_fluxes(F_F_HUMUS)  = -Cflx_yr_from_humus * years_per_ts !  per timestep ! sign CHECKED!
    ! ... as well as all fluxes from AWEN pools to atmosphere
    Y_out_fluxes(F_LITRESP)  = -Cflx_yr_resp_from_litter * years_per_ts ! per timestep ! the flux is negative
  END SUBROUTINE yasso_minimal

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> We here determine the distribution of plant roots in the soil (based on Jackson 1996), assuming an
  !> exponential decay of root material with depth.
  !>
  !> Note: The root distribution assumed here is only used to determine the vertical profile of root litter components.
  !>       It has no direct affect on the plant hydrology.
  !
  FUNCTION determine_root_dist(max_thick, layer_depth, root_dist_beta, t_soil_sl, permafrost_layer) RESULT(root_dist)
    USE mo_carbon_constants,       ONLY: permafrost_threshold
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(wp), INTENT(in)  :: max_thick(:)         !< Active layer thickness [m]
    REAL(wp), INTENT(in)  :: layer_depth(:)       !< Model layer depth [m]
    REAL(wp), INTENT(in)  :: root_dist_beta       !< PFT dependent factor for root decay with depth
    REAL(wp), INTENT(in)  :: t_soil_sl(:)         !< Soil temperature [K]
    REAL(wp), INTENT(in)  :: permafrost_layer(:)  !< permafrost status of soil layers [??]


    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':determine_root_dist'
    REAL(wp) :: root_dist(size(max_thick))     ! root distribution in each layer
    REAL(wp) :: root_dist_cum(size(max_thick)) ! cumulative root distribution for layers
    INTEGER  :: is, nsoil, nsoil_p

    nsoil = size(max_thick)
    nsoil_p = nsoil
    root_dist(:) = 0._wp
    root_dist_cum(:) = 0._wp

    DO is = 1, nsoil
      IF ((nsoil_p == nsoil) .AND. ((permafrost_layer(is) > permafrost_threshold) .OR. (t_soil_sl(is) < tmelt))) THEN
        nsoil_p = MAX(is-1, 1)
      END IF
    END DO

    IF (root_dist_beta > 0._wp) THEN
      DO is = 1, nsoil_p
        root_dist_cum(is) = 1._wp - (root_dist_beta ** (layer_depth(is) * 100._wp))
          IF (is == 1) THEN
            root_dist(is) = root_dist_cum(is)
          ELSE
            root_dist(is) = root_dist_cum(is) - root_dist_cum(is-1)
          END IF
      END DO
    END IF

    DO is = 1, nsoil
      IF (permafrost_layer(is) > permafrost_threshold) THEN
        root_dist(is) = 0._wp
      END IF
    END DO

    ! normalise root distribution -> SUM(root_dist) = 1
    IF (SUM(root_dist) > 0._wp) THEN
      root_dist(:) = root_dist(:) / SUM(root_dist(:))
    ELSE
      root_dist(:) = 0._wp
    END IF

  END FUNCTION determine_root_dist

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> The routine prepares the yasso call: The "yasso_in" buffer is filled with data from the current vertical
  !> yasso pools.
  !>
  PURE SUBROUTINE input_YASSO_pools(cpool_acid, cpool_water, cpool_ethanol,      &
    &                               cpool_nonsoluble, cpool_humus, YASSO_in)
    USE mo_carbon_constants,  ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS

    REAL(wp), INTENT(in)  :: cpool_acid
    REAL(wp), INTENT(in)  :: cpool_water
    REAL(wp), INTENT(in)  :: cpool_ethanol
    REAL(wp), INTENT(in)  :: cpool_nonsoluble
    REAL(wp), INTENT(in)  :: cpool_humus
    REAL(wp), INTENT(out) :: YASSO_in(5)

    YASSO_in(Y_ACID)       = cpool_acid
    YASSO_in(Y_WATER)      = cpool_water
    YASSO_in(Y_ETHANOL)    = cpool_ethanol
    YASSO_in(Y_NONSOLUBLE) = cpool_nonsoluble
    YASSO_in(Y_HUMUS)      = cpool_humus
  END SUBROUTINE input_YASSO_pools

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> This subroutine moves the yasso output to the respective vertical yasso pools. Depending on argument "plus"
  !> the yasso output values are added to (true) or are replacing (false) current vertical pool values.
  !>
  PURE SUBROUTINE output_YASSO_pools(cpool_acid, cpool_water, cpool_ethanol,         &
    &                                cpool_nonsoluble, cpool_humus, plus, YASSO_out)

    USE mo_carbon_constants,  ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS

    REAL(wp), INTENT(inout)  :: cpool_acid
    REAL(wp), INTENT(inout)  :: cpool_water
    REAL(wp), INTENT(inout)  :: cpool_ethanol
    REAL(wp), INTENT(inout)  :: cpool_nonsoluble
    REAL(wp), INTENT(inout)  :: cpool_humus
    LOGICAL,  INTENT(in)     :: plus
    REAL(wp), INTENT(in)     :: YASSO_out(5)

    IF (plus) THEN
      cpool_acid       = cpool_acid       + YASSO_out(Y_ACID)
      cpool_water      = cpool_water      + YASSO_out(Y_WATER)
      cpool_ethanol    = cpool_ethanol    + YASSO_out(Y_ETHANOL)
      cpool_nonsoluble = cpool_nonsoluble + YASSO_out(Y_NONSOLUBLE)
      cpool_humus      = cpool_humus      + YASSO_out(Y_HUMUS)
    ELSE
      cpool_acid       = YASSO_out(Y_ACID)
      cpool_water      = YASSO_out(Y_WATER)
      cpool_ethanol    = YASSO_out(Y_ETHANOL)
      cpool_nonsoluble = YASSO_out(Y_NONSOLUBLE)
      cpool_humus      = YASSO_out(Y_HUMUS)
    END IF
  END SUBROUTINE output_YASSO_pools

  ! ---------------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION get_modifier_temp(temp_soil, use_Q10) RESULT(d_temp)
    USE mo_carbon_constants,       ONLY: Q10, temp_p1, temp_p2
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(wp), INTENT(in)          :: temp_soil
    LOGICAL, OPTIONAL, INTENT(in) :: use_Q10

    LOGICAL :: do_Q10
    REAL(wp):: d_temp

    d_temp = 0._wp

    do_Q10 = .FALSE.
    IF (PRESENT(use_Q10)) THEN
      do_Q10 = use_Q10
    END IF

    IF (.NOT. do_Q10) THEN       ! use YASSO formulation
      d_temp  = EXP(temp_p1*(temp_soil - tmelt) + &                 ! YASSO formulation for temperature dependence
        &           temp_p2*(temp_soil- tmelt)**2.0_wp)
    ELSE                        ! use alternative Q10 formulation
      d_temp  = Q10**((temp_soil - (tmelt + 10._wp))/10._wp)        ! Temperature modifier as CLM
    END IF
  END FUNCTION get_modifier_temp

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_water(water_soil) RESULT(d_water)
    USE mo_carbon_constants,  ONLY: moist_p_a, moist_p_b, moist_p_d

    REAL(wp), INTENT(in)  :: water_soil(:)

    REAL(wp) :: d_water(size(water_soil))
    INTEGER  :: is, nsoil

    nsoil = size(water_soil)
    d_water(:) = 0._wp

    DO is = 1,nsoil
      d_water(is) = ( ( (water_soil(is) - moist_p_b) / (moist_p_a-moist_p_b) )**  &
        &             ( moist_p_d*( (moist_p_b-moist_p_a)/moist_p_a ) ) ) *       &
        &             ( (water_soil(is) / moist_p_a)**moist_p_d )
    END DO
  END FUNCTION get_modifier_water

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_depth(layer_depth) RESULT(d_depth)
    USE mo_carbon_constants,  ONLY: ref_dep

    REAL(wp), INTENT(in)  :: layer_depth(:)

    REAL(wp) :: d_depth(size(layer_depth))
    INTEGER :: is, nsoil
    nsoil = size(layer_depth)
    d_depth(:) = 0._wp

    DO is = 1, nsoil
      d_depth(is) = exp(-layer_depth(is)/ref_dep)                       ! depth modifier as CLM
    END DO

  END FUNCTION get_modifier_depth

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_size(WoodLitterSize) RESULT(d_size)
    USE mo_carbon_constants,       ONLY: size_p1, size_p2, size_p3

    REAL(wp), INTENT(in)  :: WoodLitterSize
    REAL(wp) :: d_size

    d_size = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)

  END FUNCTION get_modifier_size

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_precip(precip) RESULT(d_precip)
    USE mo_carbon_constants,       ONLY: precip_p1, sec_per_year

    REAL(wp), INTENT(in)  :: precip
    REAL(wp) :: d_precip

    d_precip = 1.0_wp - EXP(precip_p1*(precip*sec_per_year/1000._wp))      ! YASSO formulation for precipitation dependence

  END FUNCTION get_modifier_precip

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION depth_to_layers(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth), SIZE(dz))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, is

    CHARACTER(len=*), PARAMETER :: routine = modname//':depth_to_layers'

    nsoil = SIZE(dz)

    depth_l(:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO is = 1,nsoil
      bottom = bottom + dz(is)
      WHERE (depth(:) >= bottom)
        depth_l(:,is) = dz(is)
      ELSEWHERE (depth(:) > bottom_m1)
        depth_l(:,is) = depth(:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(is)
    END DO

  END FUNCTION depth_to_layers

  ! ---------------------------------------------------------------------------------------------------------
  SUBROUTINE vertical_C_transport(cpool_bg_sl, act_depth_sl, diff_coeff, delta_time_C, dz, permafrost_layer)
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    REAL(wp), INTENT(inout) :: cpool_bg_sl(:)       !< bg C pools on soil layers
    REAL(wp), INTENT(in)    :: act_depth_sl(:)      !< actual layer thickness, considering soil depth and active layer thickness
    REAL(wp), INTENT(in)    :: diff_coeff(:)        !< diffusivity
    REAL(wp), INTENT(in)    :: delta_time_C         !< Carbon time step [s]
    REAL(wp), INTENT(in)    :: dz(:)                !< layer thickness
    REAL(wp), INTENT(in)    :: permafrost_layer(:)  !< permafrost indicator

    ! layer counters, # soil layers
    INTEGER :: is, nsoil, nsoil_loc
    REAL(wp), PARAMETER :: implicit_contrib = 0.5_wp

    ! Diffusion constants averaged over adjacent layers.
    REAL(wp) :: D_below, D_above

    ! Layer-dependent weighting factors in Crank-Nicholson scheme.
    REAL(wp) :: dz_factor, W_below, W_above
    REAL(wp) :: dz_here, dz_above, dz_below
    REAL(wp) :: C_here, C_above, C_below

    ! Diagonal, left and right subdiagonals for Crank-Nicholson matrix
    REAL(wp) :: left(size(act_depth_sl))
    REAL(wp) :: diag(size(act_depth_sl))
    REAL(wp) :: right(size(act_depth_sl))

    ! Requires carbon concentrations for diffusion algorithm
    REAL(wp) :: C_conc(size(act_depth_sl))
    REAL(wp) :: src_C_conc(size(act_depth_sl))

    ! Right hand side vector for Crank-Nicholson scheme equations
    ! This is overwritten with solution by dgtsv tridiagonal solver
    REAL(wp) :: rhs(size(act_depth_sl))
    INTEGER  :: info             ! dgtsv status variable

    REAL(wp) :: C_before, C_after

    nsoil = size(act_depth_sl)
    nsoil_loc = nsoil
    src_C_conc(:) = 0._wp

    C_before = SUM(cpool_bg_sl(:))

    IF (C_before > 0._wp) THEN
      C_conc(:) = 0._wp
      DO is = 1, nsoil
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_layer(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          C_conc(is) = cpool_bg_sl(is) / act_depth_sl(is)
        ELSE
          IF (nsoil_loc == nsoil) THEN
            nsoil_loc = MAX(is - 1, 1)
          END IF
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
          dz_above = 0._wp
          C_above = 0._wp
        ELSE
          dz_above = dz(is-1)
          C_above = C_conc(is-1)
        END IF
        IF (is == nsoil_loc) THEN
          dz_below = 0._wp
          C_below = 0._wp
        ELSE
          dz_below = dz(is+1)
          C_below = C_conc(is+1)
        END IF

        ! Determine diffusion coefficients to the layers above and below the current layer.
        IF (is == 1) THEN
          ! upwards diffusivity in uppermost layer is zero
          D_above = 0._wp
        ELSE
          ! for the layers below, upwards diffusivity is the weighted mean of the layer
          ! diffusivity in layer and layer above
          D_above = (dz(is-1)*diff_coeff(is-1) + &
            &        dz(is  )*diff_coeff(is))  / &
            &        (dz(is-1) + dz(is))
        END IF
        IF (is == nsoil_loc) THEN
          ! downwards diffusivity in the bottom layer is set to zero
          ! to impose no flux boundary condition
          D_below = 0._wp
        ELSE
          ! for the layers above, downwards diffusivity is the weighted mean of the
          ! diffusivity in the layer and the layer below
          D_below = (dz(is  )*diff_coeff(is  )  + &
            &        dz(is+1)*diff_coeff(is+1)) / &
            &        (dz(is) + dz(is+1))
        END IF

        ! Calculate layer-dependent weighting factors.
        dz_factor = 0.25_wp * (dz_below + 2._wp * dz_here + dz_above)
        W_below = D_below * (delta_time_C / dz_factor) * (2._wp / (dz_below + dz_here))
        W_above = D_above * (delta_time_C / dz_factor) * (2._wp / (dz_here + dz_above))

        ! boundary conditions for surface and bottom layers: no flux through top / bottom
        IF (is == 1) THEN
          W_above = 0._wp
        END IF
        IF (is == nsoil_loc) THEN
          W_below = 0._wp
        END IF

        ! Fill in matrix diagonal and off-diagonal elements.
        diag(is) = 1.0_wp + (implicit_contrib * (W_below + W_above))
        left(is) = -implicit_contrib * W_above
        right(is) = -implicit_contrib * W_below

        ! Calculate right hand side vector values.
        rhs(is) = (1.0_wp - (1.0_wp - implicit_contrib)*(W_below + W_above)) * C_here + &
          &        (1.0_wp - implicit_contrib)*W_below*C_below +                         &
          &        (1.0_wp - implicit_contrib)*W_above*C_above +                         &
          &        src_C_conc(is)
      END DO

      ! Semi-implicit case -- whenever implicit_contrib > 0.1,
      ! fully explicit case is solced already ( = rhs)
      IF (implicit_contrib >= 0.1_wp) THEN
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
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_layer(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          cpool_bg_sl(is) = C_conc(is) * act_depth_sl(is)
        END IF
      END DO

      WHERE(cpool_bg_sl(:) < 0._wp)
        cpool_bg_sl(:) = 0._wp
      END WHERE

      C_after = SUM(cpool_bg_sl(:))

      ! re-normalise to conserve carbon
      IF (C_after > 0._wp) THEN
        cpool_bg_sl(:) = cpool_bg_sl(:)/C_after*C_before
      END IF
    END IF    ! C_before > 0

  END SUBROUTINE vertical_C_transport

  ! ---------------------------------------------------------------------------------------------------------
  SUBROUTINE set_diffusivity(diff_coeff, soil_depth, soil_depth_sl, thaw_depth, &
    &                        layer_depth, t_soil_sl, wtr_rel_sl,                &
    &                        n_thaw_cycle, permafrost_layer)
    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    REAL(wp), INTENT(out)   :: diff_coeff(:)
    REAL(wp), INTENT(in)    :: soil_depth
    REAL(wp), INTENT(in)    :: soil_depth_sl(:)
    REAL(wp), INTENT(in)    :: thaw_depth
    REAL(wp), INTENT(in)    :: layer_depth(:)
    REAL(wp), INTENT(in)    :: t_soil_sl(:)
    REAL(wp), INTENT(in)    :: wtr_rel_sl(:)
    REAL(wp), INTENT(in)    :: n_thaw_cycle
    REAL(wp), INTENT(in)    :: permafrost_layer(:)
    INTEGER :: is, isp1, nsoil

!     REAL(wp), PARAMETER :: cryo_max = 0.15_wp     ! m2/yr (?)
!     REAL(wp), PARAMETER :: bio_max  = 0.015_wp

    ! unit conversion: sec per year = 31557600. = 3e8 -->
    ! 0.15 -> 0.15/31557600 = 4.7532e-09
    !
    ! 0.15 multiplier 0.0000001 (1e-7) -> 1.5e-8
    !
    ! cm2 -> m2 --> 1m2 = 100cm*100cm = 10000cm2 ==> 1m2 = 1e5 cm2
    !
    ! 15 cm2/yr --> 15/1e5 m2/yr -> 15/1e5 m2/31557600s --> 15/1e5/31557600 m2/s --> 4.7532e-12 m2/s
    !

    REAL(wp), PARAMETER :: cryo_max = 4.7e-12_wp     ! m2/yr (?)
    REAL(wp), PARAMETER :: bio_max  = 4.7e-13_wp


    nsoil = size(diff_coeff)

    diff_coeff(:) = 0._wp

    IF (thaw_depth < soil_depth) THEN    ! active layer depth < soil depth --> permafrost
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF (  (permafrost_layer(is)   < permafrost_threshold) .AND. &
          &   (permafrost_layer(isp1) < permafrost_threshold) .AND. &
          &   (t_soil_sl(is)          > tmelt)) THEN
          diff_coeff(is) = cryo_max *                               &
            &             MIN(1._wp, wtr_rel_sl(is)) *              &
            &              MIN(1._wp, (n_thaw_cycle / 40._wp))*     &
            &              (1._wp - layer_depth(is)/thaw_depth)
        ELSE
          diff_coeff(is) = 0._wp
        END IF
      END DO
    ELSE
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF ((soil_depth_sl(is) > minimum_layer_thickness) .AND. (t_soil_sl(is) > tmelt)) THEN
          diff_coeff(is) = MAX(0._wp, bio_max *(1._wp - (layer_depth(is)/3._wp)))
        ELSE
          diff_coeff(is) = 0._wp
        END IF
      END DO
    END IF

  END SUBROUTINE set_diffusivity

  ! ---------------------------------------------------------------------------------------------------------
  ! move_ag_2_bg
  ! check column C amounts, ensure that C makes sense.
  ! Idea: peat column defines maximum C possible, if more C in 1 layer, it needs to be moved down.
    SUBROUTINE move_ag_2_bg(c_acid_ag1, c_acid_bg1_sl, c_acid_ag2, c_acid_bg2_sl,                         &
      &                     c_water_ag1, c_water_bg1_sl, c_water_ag2, c_water_bg2_sl,                     &
      &                     c_ethanol_ag1, c_ethanol_bg1_sl, c_ethanol_ag2, c_ethanol_bg2_sl,             &
      &                     c_nonsoluble_ag1, c_nonsoluble_bg1_sl, c_nonsoluble_ag2, c_nonsoluble_bg2_sl, &
      &                     c_humus_1_sl, c_humus_2_sl, max_depth, layer_thick, permafrost_layer)

    USE mo_carbon_constants, ONLY: C_2_biomass_ratio, rho_org_ag, rho_org_bg, molarMassC_kg, &
      &                            permafrost_threshold

    REAL(wp), INTENT(inout) :: c_acid_ag1
    REAL(wp), INTENT(inout) :: c_acid_bg1_sl(:)
    REAL(wp), INTENT(inout) :: c_acid_ag2
    REAL(wp), INTENT(inout) :: c_acid_bg2_sl(:)
    REAL(wp), INTENT(inout) :: c_water_ag1
    REAL(wp), INTENT(inout) :: c_water_bg1_sl(:)
    REAL(wp), INTENT(inout) :: c_water_ag2
    REAL(wp), INTENT(inout) :: c_water_bg2_sl(:)
    REAL(wp), INTENT(inout) :: c_ethanol_ag1
    REAL(wp), INTENT(inout) :: c_ethanol_bg1_sl(:)
    REAL(wp), INTENT(inout) :: c_ethanol_ag2
    REAL(wp), INTENT(inout) :: c_ethanol_bg2_sl(:)
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1_sl(:)
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2_sl(:)
    REAL(wp), INTENT(inout) :: c_humus_1_sl(:)
    REAL(wp), INTENT(inout) :: c_humus_2_sl(:)
    REAL(wp), INTENT(in)    :: max_depth
    REAL(wp), INTENT(in)    :: layer_thick(:)
    REAL(wp), INTENT(in)    :: permafrost_layer(:)

    REAL(wp),  PARAMETER :: thresh_frac_org    = 1.25_wp  !< Threshold organic fraction above which C is being moved downwards
    REAL(wp),  PARAMETER :: ag_layer_depth     = 0.13_wp
    REAL(wp) :: cpool_tot
    REAL(wp) :: depth_upper
    REAL(wp) :: rfrac
    REAL(wp) :: tmp

    REAL(wp) :: max_C_ag
    REAL(wp) :: max_C_bg(size(layer_thick))
    INTEGER  :: is, isp1, nsoil

    nsoil = SIZE(layer_thick)

    ! determine max C amount per layer
    ! C [kgC] = height [m] x A [m2] x rho [kg/m3] x C_2_bio [kgC/kg]
    ! finally converted to mol
    max_C_ag = ag_layer_depth * rho_org_ag * C_2_biomass_ratio / molarMassC_kg
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
    depth_upper = 0._wp
    DO is = 1,nsoil-1
      isp1 = MIN(is+1, nsoil)
      IF ((depth_upper < max_depth) .AND. (permafrost_layer(is) < permafrost_threshold) &
        & .AND. (permafrost_layer(isp1) < permafrost_threshold)) THEN
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

  ! --- remove_C_below_bedrock ---------------------------------------------------------------------------------------------------
  !
  SUBROUTINE remove_C_below_bedrock(c_pool_bg, soil_depth, layer_thick)
    REAL(wp), INTENT(inout)  :: c_pool_bg(:)
    REAL(wp), INTENT(in)     :: soil_depth
    REAL(wp), INTENT(in)     :: layer_thick(:)

    INTEGER  :: is, nsoil, max_layer_above
    REAL(wp) :: soil_above, C_below_bedrock, upper_bound

    nsoil = size(layer_thick)

    soil_above = 0._wp
    C_below_bedrock = 0._wp
    max_layer_above = 0
    upper_bound = 0._wp

    ! Check for C below bedrock depth
    DO is = 1, nsoil
      IF (soil_depth > upper_bound) THEN
        max_layer_above = is
        soil_above = soil_above + layer_thick(is)
      ELSE
        IF (c_pool_bg(is) > 0._wp) THEN
          C_below_bedrock = C_below_bedrock + c_pool_bg(is)
        END IF
      END IF
      upper_bound = upper_bound + layer_thick(is)
    END DO

    ! Move C to lowest layer above bedrock, the lowest layer with a soil fraction
    IF ((C_below_bedrock > 0._wp) .AND. (soil_above > 0._wp) .AND. (max_layer_above > 0)) THEN
      c_pool_bg(max_layer_above) = c_pool_bg(max_layer_above) + C_below_bedrock
      DO is = MIN(max_layer_above+1,nsoil), nsoil
        c_pool_bg(is) = 0._wp
      END DO
    END IF

  END SUBROUTINE remove_C_below_bedrock


  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  !> This routine handles the impact of fire on the carbon pools.
  !>
  !> We assume that the above ground litter carbon (yasso ag1 and ag2 pools), is completely released
  !> to the atmosphere while the carbon of the living plants (c_green, c_reserve, c_woods) is partly
  !> released to the atmosphere and partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_vertical_fire( &
    & delta_time_C,                      & ! in
    & burned_fract,                      & ! in
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
    & co2flux_fire_all_2_atm,            & ! out
    & c_acid_bg1_sl,                     & ! inout
    & c_water_bg1_sl,                    & ! inout
    & c_ethanol_bg1_sl,                  & ! inout
    & c_nonsoluble_bg1_sl,               & ! inout
    & c_acid_bg2_sl,                     & ! inout
    & c_water_bg2_sl,                    & ! inout
    & c_ethanol_bg2_sl,                  & ! inout
    & c_nonsoluble_bg2_sl,               & ! inout
    & c_humus_1_sl,                      & ! inout
    & c_humus_2_sl                       & ! inout
    & )

    USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, molarMassCO2_kg
    USE mo_carbon_shared,    ONLY: distribute_yasso_litter

    ! Arguments
    REAL(wp), INTENT(in)    :: delta_time_C                 !< Carbon time step [s]
    REAL(wp), INTENT(in)    :: burned_fract(:)              !< Fraction of the vegetated area of each tile burned till the last
                                                            !< call of this routine
    REAL(wp), INTENT(in)    :: fire_fract_wood_2_atmos      !< Fraction of above ground wood immediately emitted to atmosphere
                                                            !< in case of fire
    REAL(wp), INTENT(in)    :: fuel_moisture(:)             !
    REAL(wp), INTENT(in)    :: LeafLit_coef(:)              !< Factor to spread non woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(in)    :: WoodLit_coef(:)              !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(inout) :: c_green(:)                   !< Amount of green carbon [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods(:)                   !< Amount of woody carbon [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve(:)                 !< Amount of reserve carbon [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    !< CO2 released by fire [kg(CO2)m-2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_green_2_soil(:)   !< Carbon relocated by disturbances from non-woody biomass
                                                            !< to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   !< Carbon relocated by disturbances from woody biomass to
                                                            !< the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1: Green, originating from non-woody biomass
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag1(:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1(:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1(:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1(:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg1(:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1(:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1(:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1(:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1(:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2: Woody, originating from woody biomass
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag2(:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2(:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2(:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2(:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg2(:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2(:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2(:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2(:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2(:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    REAL(wp), INTENT(inout) :: c_acid_bg1_sl(:,:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1_sl(:,:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1_sl(:,:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1_sl(:,:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1_sl(:,:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    REAL(wp), INTENT(inout) :: c_acid_bg2_sl(:,:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2_sl(:,:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2_sl(:,:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2_sl(:,:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2_sl(:,:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Carbon released due to disturbances from the green and reserve pools, that enters the green litter and
      ! soil pools. We here add the contribution of fires.
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._wp - fract_green_aboveGround) * burned_fract(ic) / delta_time_C
      ! Carbon released due to disturbances from the woody vegetation pool, that enters the woody litter and soil
      ! pools. We here add the contribution of fires.
      cflux_dist_woods_2_soil(ic)  = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(ic) / delta_time_C

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
      fraction = 1._wp - burned_fract(ic)
      c_acid_ag1      (ic) = fraction * c_acid_ag1       (ic)
      c_water_ag1     (ic) = fraction * c_water_ag1      (ic)
      c_ethanol_ag1   (ic) = fraction * c_ethanol_ag1    (ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1 (ic)

      c_acid_ag2      (ic) = fraction * c_acid_ag2       (ic)
      c_water_ag2     (ic) = fraction * c_water_ag2      (ic)
      c_ethanol_ag2   (ic) = fraction * c_ethanol_ag2    (ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2 (ic)

      ! Add the below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._wp-fract_green_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1_sl(ic,1), c_water_bg1_sl(ic,1), c_ethanol_bg1_sl(ic,1), &
        &                          c_nonsoluble_bg1_sl(ic,1), c_humus_1_sl(ic,1), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remanings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._wp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2_sl(ic,1), c_water_bg2_sl(ic,1), c_ethanol_bg2_sl(ic,1), &
        &                          c_nonsoluble_bg2_sl(ic,1), c_humus_2_sl(ic,1), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - burned_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

      c_acid_bg1(ic)       = SUM(c_acid_bg1_sl(ic,:))
      c_water_bg1(ic)      = SUM(c_water_bg1_sl(ic,:))
      c_ethanol_bg1(ic)    = SUM(c_ethanol_bg1_sl(ic,:))
      c_nonsoluble_bg1(ic) = SUM(c_nonsoluble_bg1_sl(ic,:))
      c_humus_1(ic)        = SUM(c_humus_1_sl(ic,:))
      c_acid_bg2(ic)       = SUM(c_acid_bg2_sl(ic,:))
      c_water_bg2(ic)      = SUM(c_water_bg2_sl(ic,:))
      c_ethanol_bg2(ic)    = SUM(c_ethanol_bg2_sl(ic,:))
      c_nonsoluble_bg2(ic) = SUM(c_nonsoluble_bg2_sl(ic,:))
      c_humus_2(ic)        = SUM(c_humus_2_sl(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_vertical_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  !< In this routine we model the impact of vegetation damage (i.e. wind through) on the carbon pools:
  !< Depending on the disturbed area fraction carbon from the living vegetation pools (c_green,
  !< c_reserve and c_woods) is relocated  to the litter and soil pools.
  !
  SUBROUTINE relocate_carbon_vertical_damage(  &
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
    & c_humus_2,                      & ! inout
    & c_acid_bg1_sl,                  & ! inout
    & c_water_bg1_sl,                 & ! inout
    & c_ethanol_bg1_sl,               & ! inout
    & c_nonsoluble_bg1_sl,            & ! inout
    & c_acid_bg2_sl,                  & ! inout
    & c_water_bg2_sl,                 & ! inout
    & c_ethanol_bg2_sl,               & ! inout
    & c_nonsoluble_bg2_sl,            & ! inout
    & c_humus_1_sl,                   & ! inout
    & c_humus_2_sl                    & ! inout
    & )

    USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround
    USE mo_carbon_shared,    ONLY: distribute_yasso_litter

    ! Arguments
    REAL(wp), INTENT(in)    :: delta_time_C                !< Carbon time step [s]
    real(wp), intent(in)    :: damaged_fract(:)            !< Fraction of the vegetated area of each tile damaged since the last
                                                           !< call of this routine
    real(wp), intent(in)    :: LeafLit_coef(:)             !< Factor to spread non woody litter into yasso pools [ ]
    real(wp), intent(in)    :: WoodLit_coef(:)             !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(inout) :: c_green(:)                  !< Amount of carbon in green pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods(:)                  !< Amount carbon in wood pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve(:)                !< Amount of carbon in reserve pool [mol(C)/m^(canopy)2]
    REAL(wp), INTENT(inout) :: cflux_dist_green_2_soil(:)  !< Disturbance-induced carbon flux of non-woody living biomass
                                                           !< to the green litter and soil pools [mol(C)/m^2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)  !< Disturbance-induced carbon flux of woody living biomass to
                                                           !< the woody litter and soil pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag1(:)        !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag1(:)       !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag1(:)     !< Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag1(:)  !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg1(:)        !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg1(:)       !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg1(:)     !< Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg1(:)  !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_1(:)         !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag2(:)        !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag2(:)       !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag2(:)     !< Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag2(:)  !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg2(:)        !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg2(:)       !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg2(:)     !< Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg2(:)  !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_2(:)         !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    REAL(wp), INTENT(inout) :: c_acid_bg1_sl(:,:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1_sl(:,:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1_sl(:,:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1_sl(:,:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1_sl(:,:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    REAL(wp), INTENT(inout) :: c_acid_bg2_sl(:,:)          !< Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2_sl(:,:)         !< Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2_sl(:,:)       !< Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2_sl(:,:)    !< Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2_sl(:,:)           !< Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
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
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._wp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1_sl(ic,1), c_water_bg1_sl(ic,1), c_ethanol_bg1_sl(ic,1), &
        &                          c_nonsoluble_bg1_sl(ic,1), c_humus_1_sl(ic,1), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the carbon of damaged woody vegetation to the yasso litter and soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2_sl(ic,1), c_water_bg2_sl(ic,1), c_ethanol_bg2_sl(ic,1), &
        &                          c_nonsoluble_bg2_sl(ic,1), c_humus_2_sl(ic,1), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

      c_acid_bg1(ic)       = SUM(c_acid_bg1_sl(ic,:))
      c_water_bg1(ic)      = SUM(c_water_bg1_sl(ic,:))
      c_ethanol_bg1(ic)    = SUM(c_ethanol_bg1_sl(ic,:))
      c_nonsoluble_bg1(ic) = SUM(c_nonsoluble_bg1_sl(ic,:))
      c_humus_1(ic)        = SUM(c_humus_1_sl(ic,:))
      c_acid_bg2(ic)       = SUM(c_acid_bg2_sl(ic,:))
      c_water_bg2(ic)      = SUM(c_water_bg2_sl(ic,:))
      c_ethanol_bg2(ic)    = SUM(c_ethanol_bg2_sl(ic,:))
      c_nonsoluble_bg2(ic) = SUM(c_nonsoluble_bg2_sl(ic,:))
      c_humus_2(ic)        = SUM(c_humus_2_sl(ic,:))

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_vertical_damage

#endif
END MODULE mo_carbon_process_vertical
