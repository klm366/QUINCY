!> Contains the routines for the lateral (intertile) exchange processes
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

!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_lex_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport, calc_drainage_transport,  &
          & calc_lat_pond_soil_flux, calc_diffusion_water,   &
          & calc_diffusion_energy, calc_snow_redistribution, &
          & remove_exess_water, shift_fluxes_sl,             &
          & get_geometry, get_flow_lags, check_flow_matrix

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! Compute intertile exchange of surface runoff
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_runoff_transport(                        &
    & dtime,                                               &
    & fract_source_2_target,                               &
    & fract_source, fract_target,                          &
    & wtr_latflow_srf,                                     &
    & wtr_lex_runoff_target_acc,                           &
    & wtr_lex_runoff_source_acc,                           &
    & wtr_lex_srf_source_acc,                              &
    & wtr_lex_srf_target_acc,                              &
    & wtr_pond_target                                      &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime                          !< Time step length [s]

    REAL(wp), INTENT(in)    :: fract_source_2_target,      &  !< Fraction of outflow from source tile received by target tile [/]
      &                        fract_source,               &  !< Fraction of source tile [/]
      &                        fract_target,               &  !< Fraction of target tile [/]
      &                        wtr_latflow_srf                !< Outflow of intermediate storage on source tile [m]

    REAL(wp), INTENT(inout) :: wtr_lex_runoff_target_acc,  &  !< Inflow of intermediate storage received on target tile [kg m-2 s-1]
      &                        wtr_lex_runoff_source_acc,  &  !< Outflow of intermediate storage received by target tiles [kg m-2 s-1]
      &                        wtr_lex_srf_source_acc,     &  !< Netflux (surface) on source tile [kg m-2 s-1]
      &                        wtr_lex_srf_target_acc,     &  !< Netflux (surface) on target tile [kg m-2 s-1]
      &                        wtr_pond_target                !< Surface water storage on target tile    [m]

    !$ACC ROUTINE SEQ

    IF (fract_source < 1.0e-10_wp .OR. fract_target < 1.0e-10_wp) RETURN

    IF (fract_source_2_target > 1.0e-10_wp .AND. wtr_latflow_srf > 1.0e-10_wp) THEN

      wtr_lex_runoff_source_acc = wtr_lex_runoff_source_acc               &
        &                       + wtr_latflow_srf * fract_source_2_target &
        &                       * rhoh2o / dtime
      
      wtr_lex_srf_source_acc    = wtr_lex_srf_source_acc                  &
        &                       - wtr_latflow_srf * fract_source_2_target &
        &                       * rhoh2o / dtime

      wtr_lex_runoff_target_acc = wtr_lex_runoff_target_acc               &
        &                       + wtr_latflow_srf * fract_source_2_target &
        &                       * (fract_source / fract_target) * rhoh2o / dtime

      wtr_pond_target           = wtr_pond_target                         &
        &                       + wtr_latflow_srf * fract_source_2_target &
        &                       * (fract_source / fract_target)
 
      wtr_lex_srf_target_acc    = wtr_lex_srf_target_acc                  &
        &                       + wtr_latflow_srf * fract_source_2_target &
        &                       * (fract_source / fract_target) * rhoh2o / dtime

    END IF

  END SUBROUTINE calc_runoff_transport

  !
  ! ===============================================================================================================================
  !>
  !! Compute intertile exchange of below ground runoff
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_drainage_transport(                     &
    & dtime,                                              &
    & fract_source_2_target,                              &
    & fract_source, fract_target,                         &
    & wtr_latflow_sl,                                     &
    & wtr_lex_drain_target_acc,                           &
    & wtr_lex_drain_source_acc,                           &
    & wtr_lex_sl_source_acc,                              &
    & wtr_lex_sl_target_acc,                              &
    & wtr_soil_sl_target                                  &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime                         !< Time step length [s]

    REAL(wp), INTENT(in)    :: fract_source_2_target,     &  !< Fraction of outflow from source tile received by target tile [/]
      &                        fract_source,              &  !< Fraction of source tile [/]
      &                        fract_target,              &  !< Fraction of target tile [/]
      &                        wtr_latflow_sl                !< Outflow of intermediate storage on source tile [m]


    REAL(wp), INTENT(inout) :: wtr_lex_drain_target_acc,  &  !< Inflow of intermediate storage received on target tile [kg m-2 s-1]
      &                        wtr_lex_drain_source_acc,  &  !< Outflow of intermediate storage received by target tiles [kg m-2 s-1]
      &                        wtr_lex_sl_source_acc,     &  !< Netflux (below ground) on source tile [kg m-2 s-1]
      &                        wtr_lex_sl_target_acc,     &  !< Netflux (below ground) on target tile [kg m-2 s-1]
      &                        wtr_soil_sl_target            !< Soil moisture on target tile [m]

    !$ACC ROUTINE SEQ

    IF (fract_source < 1.0e-10_wp .OR. fract_target < 1.0e-10_wp) RETURN

    IF (fract_source_2_target > 1.0e-10_wp .AND. wtr_latflow_sl > 1.0e-10_wp) THEN

      wtr_lex_drain_source_acc = wtr_lex_drain_source_acc               &
        &                      + wtr_latflow_sl * fract_source_2_target &
        &                      * rhoh2o / dtime
      
      wtr_lex_sl_source_acc    = wtr_lex_sl_source_acc                     &
        &                      - wtr_latflow_sl * fract_source_2_target    &
        &                      * rhoh2o / dtime

      wtr_lex_drain_target_acc = wtr_lex_drain_target_acc               &
        &                      + wtr_latflow_sl * fract_source_2_target &
        &                      * (fract_source / fract_target) * rhoh2o / dtime

      wtr_lex_sl_target_acc    = wtr_lex_sl_target_acc                     &
        &                      + wtr_latflow_sl * fract_source_2_target    &
        &                      * (fract_source / fract_target) * rhoh2o / dtime

      wtr_soil_sl_target       = wtr_soil_sl_target                        &
        &                      + wtr_latflow_sl * fract_source_2_target    &
        &                      * (fract_source / fract_target)
    END IF

  END SUBROUTINE calc_drainage_transport
  !
  ! ===============================================================================================================================
  !>
  !! Compute intertile lateral ingress from pond to adjacent soil
  !!
  SUBROUTINE calc_lat_pond_soil_flux(             &
      &      nsoil, dtime,                        &
      &      h_offset, length_scale, len_inter,   &
      &      lvt_fract, dpr_fract, dsoil_lvt,     &
      &      field_cap_lvt,  k_lvt, mpot_lvt,     &
      &      ingress_lvt,                         &
      &      egress_dpr,                          &
      &      wtr_latflow_dpr, wtr_latflow_lvt,    &
      &      wtr_pond_dpr, wtr_soil_lvt           &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER,  INTENT(in)    :: nsoil                 !< Number of soil ayers

    REAL(wp), INTENT(in)    :: dtime,             &  !< Timestep length [s]
      &                        h_offset,          &  !< Differece in surface elevation between elevation and depression [m]
      &                        length_scale,      &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,         &  !< Length of interface between elevation and depression [m m-2]
      &                        lvt_fract,         &  !< Area fraction of elevation tile [/]
      &                        dpr_fract             !< Area fraction of depression tile [/]

    REAL(wp), INTENT(in)    :: dsoil_lvt(:),      &  !< Soil depth on layers of elevation tile [m]
      &                        field_cap_lvt(:),  &  !< Soil field capacity of elevation tile [m]
      &                        k_lvt(:),          &  !< Soil layer hyd. conductivity on source tile [m s-1]
      &                        mpot_lvt(:)           !< Matric potential on source tile [m]


    REAL(wp), INTENT(inout) :: ingress_lvt,       &  !< Diffusion flux on elevation tile [kg m-2 s-1]
      &                        egress_dpr            !< Diffusion flux on depression tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_soil_lvt(:),   &  !< Soil water storage on elevation tile [m]
      &                        wtr_latflow_lvt(:)    !< Netflux (below ground) on elevation tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_pond_dpr,      &  !< Pond storage on depression tile [m]
      &                        wtr_latflow_dpr       !< Netflux (surface) on depression tile [kg m-2 s-1]


    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,             &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, i_lvt              !< Index variables for layers on helper grid (i) and in elevation (i_lvt)  

    REAL(wp)                :: fr_lvt_rel,        &  !< Relative fraction of elevation tile [/]
      &                        fr_dpr_rel,        &  !< Relative fraction of depression tile [/]
      &                        wtr_pond_dpr_ini,  &  !< Initial surface water storage of depression tile [m]
      &                        free_space,        &  !< Free pore space on elevation tile on a given layer [m]
      &                        wtr_avail,         &  !< Available water on depression tile on a given layer [m]
      &                        wtr_flx_lat,       &  !< Potential flux (based on k and mpot) between depression and elevation [m/s]
      &                        wtr_lat               !< Actual flux between depression and elevation [m]

    REAL(wp)                :: diffusion(nsoil)      !< Diffusive flux between depression and elevation on soil layers [m]


    REAL(wp)                :: z_lvt,             &  !< Layer bottom of current layer on elevation tile [m]
      &                        h,                 &  !< Height of current layer [m]
      &                        a_inter,           &  !< Interface area on current layer [m2 m-2]
      &                        z_rgrd_top,        &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,        &  !< Bottom of current layer on helper grid [m] 
      &                        z_rgrd_bot_new,    &  !< New bottom of current layer on helper grid [m]
      &                        z_bot_lvt,         &  !< Bedrock boundary on elevation tile [m]
      &                        wtr_soil_lvt_rgrd, &  !< Soil water content on elevation tile on helper grid layer [m]
      &                        field_cap_lvt_rgrd,&  !< Field capacity on elevation tile on helper grid layer [m]
      &                        mpot_lvt_rgrd,     &  !< Matric potential  on elevation tile on helper grid layer [m]
      &                        k_lvt_rgrd            !< Hydraulic conductivity on elevation tile on helper grid layer [[m s-1]

    nrgrd = 2*nsoil

    ! Initilize fields for accumulating regridded output
    diffusion(:)        = 0._wp

    IF (lvt_fract < 1.0e-10_wp .OR. dpr_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    wtr_pond_dpr_ini = wtr_pond_dpr


    IF (h_offset < 1.0e-10_wp .OR. wtr_pond_dpr < 1.0e-10_wp) RETURN

    ! Start setting up the interpolated helper grid
    i_lvt = 1
    z_bot_lvt = SUM(dsoil_lvt(:))
    z_lvt     = dsoil_lvt(i_lvt)
    z_rgrd_top = 0._wp

    IF (h_offset > wtr_pond_dpr_ini) THEN
      z_rgrd_bot = MIN(z_lvt, h_offset - wtr_pond_dpr_ini)
    ELSE
      z_rgrd_bot = MIN(z_lvt, h_offset)
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO i = 1, nrgrd
      IF (z_rgrd_bot <= z_bot_lvt) THEN ! Within soil of elevation
        h                     = z_rgrd_bot - z_rgrd_top
        a_inter               = len_inter * h
        wtr_soil_lvt_rgrd     = wtr_soil_lvt(i_lvt)  * h / dsoil_lvt(i_lvt)
        wtr_lat               = 0._wp
        IF (z_rgrd_bot <= h_offset) THEN ! Above depression surface
          field_cap_lvt_rgrd  = field_cap_lvt(i_lvt) * h / dsoil_lvt(i_lvt)
          k_lvt_rgrd          = k_lvt(i_lvt)
          mpot_lvt_rgrd       = mpot_lvt(i_lvt)
          IF (z_rgrd_bot > h_offset - wtr_pond_dpr_ini) THEN ! Within the ponded surface water in the depression
            ! Water flows from the pond reservoir of the depression tile into the soil of the elevation tile
            wtr_flx_lat       = abs(k_lvt_rgrd  * mpot_lvt_rgrd) / length_scale *  a_inter
            free_space        = MAX(0._wp,(field_cap_lvt_rgrd - wtr_soil_lvt_rgrd) * fr_lvt_rel)
            wtr_avail         = MIN(h, wtr_pond_dpr) * fr_dpr_rel
            wtr_lat           = MIN(wtr_flx_lat * dtime, free_space, wtr_avail)
          END IF
        END IF
        diffusion(i_lvt)      = diffusion(i_lvt) + wtr_lat
      END IF  
      ! Layer bottom was on soil layer
      IF (z_rgrd_bot == z_lvt) THEN
        i_lvt = i_lvt + 1
        IF (i_lvt <= nsoil) THEN
          IF (dsoil_lvt(i_lvt) > 0._wp) THEN
            z_lvt = z_lvt + dsoil_lvt(i_lvt)
          ELSE
            z_lvt = 10000._wp
            i_lvt = nsoil + 1
          END IF
        ELSE
          z_lvt = 10000._wp
          i_lvt = nsoil + 1
        END IF
      END IF

      ! Finding new layer Bottom
      IF (z_rgrd_bot < h_offset - wtr_pond_dpr_ini) THEN ! Still above water level in depression
        z_rgrd_bot_new = MIN(h_offset - wtr_pond_dpr_ini, z_lvt)
      ELSE IF (z_rgrd_bot < h_offset ) THEN              ! Within the ponded surface water in the depression 
        z_rgrd_bot_new = MIN(h_offset, z_lvt)
      ELSE
        z_rgrd_bot_new = z_lvt
      END IF       

      z_rgrd_top=z_rgrd_bot
      z_rgrd_bot=z_rgrd_bot_new      
    END DO
    !$ACC END PARALLEL LOOP

    wtr_soil_lvt     = wtr_soil_lvt    +     diffusion  / fr_lvt_rel
    wtr_pond_dpr     = wtr_pond_dpr    - SUM(diffusion) / fr_dpr_rel

    wtr_latflow_lvt  = wtr_latflow_lvt +     diffusion  / fr_lvt_rel * rhoh2o / dtime
    wtr_latflow_dpr  = wtr_latflow_dpr - SUM(diffusion) / fr_dpr_rel * rhoh2o / dtime

    ingress_lvt      = ingress_lvt     + SUM(diffusion) / fr_lvt_rel * rhoh2o / dtime
    egress_dpr       = egress_dpr      + SUM(diffusion) / fr_dpr_rel * rhoh2o / dtime

  END SUBROUTINE calc_lat_pond_soil_flux

  ! ===============================================================================================================================
  !>
  !! Compute intertile lateral diffusion (water)
  !!
  SUBROUTINE calc_diffusion_water(                  &
      &      nsoil, dtime,                          &
      &      h_offset, length_scale, len_inter,     &
      &      lvt_fract, dpr_fract,                  &
      &      dsoil_lvt, dsoil_dpr,                  &
      &      field_cap_lvt, field_cap_dpr,          &
      &      wsat_lvt, wsat_dpr,                    &
      &      res_lvt, res_dpr,                      &
      &      k_sat_lvt, k_sat_dpr,                  &
      &      k_lvt, k_dpr,                          &
      &      mpot_lvt, mpot_dpr,                    &
      &      diffus_in_lvt, diffus_in_dpr,          &
      &      diffus_out_lvt, diffus_out_dpr,        &
      &      wtr_latflow_lvt, wtr_latflow_dpr,      &
      &      wtr_soil_lvt, wtr_soil_dpr             &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER, INTENT(in)     :: nsoil                   !< Number of soil layers

    REAL(wp), INTENT(in)    :: dtime,               &  !< Timestep length [s]
      &                        h_offset,            &  !< Differece in surface elevation between elevation and depression [m]
      &                        length_scale,        &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,           &  !< Length of interface between elevation and depression [m m-2]
      &                        lvt_fract,           &  !< Area fraction of elevation tile [/]
      &                        dpr_fract               !< Area fraction of depression tile [/]

    REAL(wp), INTENT(in)    :: dsoil_lvt(:),        &  !< Soil depth on layers of elevation tile [m]
      &                        dsoil_dpr(:),        &  !< Soil depth on layers of depression tile [m]
      &                        field_cap_lvt(:),    &  !< Soil field capacity of elevation tile [m]
      &                        field_cap_dpr(:),    &  !< Soil field capacity of depression tile [m]
      &                        wsat_lvt(:),         &  !< Soil moisture at saturation depression tile [m]
      &                        wsat_dpr(:),         &  !< Soil moisture at saturation elevation tile [m]
      &                        res_lvt(:),          &  !< Residual soil moisture elevation tile [m]
      &                        res_dpr(:),          &  !< Residual soil moisture depression tile [m]
      &                        k_sat_lvt(:),        &  !< Saturation hydraulic conductivity elevation tile [m s-1]
      &                        k_sat_dpr(:),        &  !< Saturation hydraulic conductivity depression tile [m s-1]
      &                        k_lvt(:),            &  !< Soil layer hyd. conductivity on source tile [m s-1]
      &                        k_dpr(:),            &  !< Soil layer hyd. conductivity on target tile [m s-1]
      &                        mpot_lvt(:),         &  !< Soil layer matric potential on source tile [m]
      &                        mpot_dpr(:)             !< Soil layer matric potential on target tile [m]

    REAL(wp), INTENT(inout) :: diffus_in_lvt(:),    &  !< Diffusion influx at elevation tile [kg m-2 s-1]
      &                        diffus_in_dpr(:),    &  !< Diffusion influx at depression tile [kg m-2 s-1]
      &                        diffus_out_lvt(:),   &  !< Diffusion influx at elevation tile [kg m-2 s-1]
      &                        diffus_out_dpr(:)       !< Diffusion influx at depression tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_latflow_lvt(:),  &  !< Netflux (below ground) on elevation tile [kg m-2 s-1]
      &                        wtr_latflow_dpr(:)      !< Netflux (below ground) on depression tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_soil_lvt(:),     &  !< Soil water storage on elevation tile [m]
      &                        wtr_soil_dpr(:)         !< Soil water storage on depression tile [m]



    ! Local variables ------------------------------------------------------------------------------------

    INTEGER                 :: nrgrd,               &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, i_lvt, i_dpr         !< Index variables for layers on helper grid (i), in elevation (i_lvt) and depression tile (i_dpr)  

    REAL(wp)                :: fr_lvt_rel,          &  !< Relative fraction of elevation tile [/]
      &                        fr_dpr_rel,          &  !< Relative fraction of depression tile [/]
      &                        h,                   &  !< Height of current layer [m]
      &                        a_inter,             &  !< Interface area on current layer [m2 m-2]
      &                        k_inter,             &  !< Hydraulic conductivity at interface [m s-1]
      &                        d_mpot,              &  !< Difference in matric potential [m]
      &                        wtr_rel_qlbrm,       &  !< Relative saturation at equilibrium between connected tiles [/]
      &                        free_space,          &  !< Free pore space on elevation tile on a given layer [m]  
      &                        wtr_avail,           &  !< Available water on depression tile on a given layer [m]
      &                        wtr_2_qlbrm,         &  !< Amount of water moved in order to reach equilibrium state [m]
      &                        wtr_flx_lat,         &  !< potential flux (based on k and mpot) between depression and elevation [m/s]
      &                        wtr_lat,             &  !< actual flux between depression and elevation [m]
      &                        z_rgrd_top,          &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,          &  !< Bottom of current layer on helper grid [m] 
      &                        z_lvt,               &  !< Layer bottom of current layer on elevation tile [m]
      &                        z_dpr,               &  !< Layer bottom of current layer on depression tile [m]
      &                        z_bot_lvt,           &  !< Bedrock boundary on elevation tile [m]
      &                        z_bot_dpr,           &  !< Bedrock boundary on depression tile [m]
      &                        wtr_soil_lvt_rgrd,   &  !< Soil water content on elevation tile on helper grid layer [m]
      &                        wsat_lvt_rgrd,       &  !< Soil moisture at saturation on elevation tile on helper grid layer [m] 
      &                        field_cap_lvt_rgrd,  &  !< Field capacity on elevation tile on helper grid layer [m]
      &                        res_lvt_rgrd,        &  !< Residual soil moisture on elevation tile on helper grid layer [m]
      &                        mpot_lvt_rgrd,       &  !< Matric potential on elevation tile on helper grid layer [m]
      &                        k_lvt_rgrd,          &  !< Hydraulic conductivity on elevation tile on helper grid layer [[m s-1]
      &                        wtr_soil_dpr_rgrd,   &  !< Soil water content on depression tile on helper grid layer [m]
      &                        wsat_dpr_rgrd,       &  !< Soil moisture at saturation on depression tile on helper grid layer [m]
      &                        field_cap_dpr_rgrd,  &  !< Field capacity on depression tile on helper grid layer [m]
      &                        res_dpr_rgrd,        &  !< Residual soil moisture on depression tile on helper grid layer [m]
      &                        mpot_dpr_rgrd,       &  !< Matric potential on deperssion tile on helper grid layer [m]
      &                        k_dpr_rgrd              !< Hydraulic conductivity on depression tile on helper grid layer [[m s-1]

    REAL(wp)                :: diffus_lvt(nsoil),   &  !< Net diffusive flux on elevation on soil layers [m]
      &                        diffus_dpr(nsoil)       !< Net diffusive flux on depression on soil layers [m]

    nrgrd         = nsoil * 2 + 3

    ! Initilize fields for accumulating regridded output
    diffus_lvt(:)        = 0._wp
    diffus_dpr(:)        = 0._wp

    IF (h_offset > SUM(dsoil_lvt(:))) RETURN
    IF (lvt_fract < 1.0e-10_wp .OR. dpr_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    ! Start setting up the interpolated helper grid
    ! Initial values for elevation tile
    i_lvt = 1
    z_bot_lvt = SUM(dsoil_lvt(:))
    z_lvt = dsoil_lvt(i_lvt)

    ! Initial values for depression tile 
    IF (h_offset > 0._wp) THEN
      i_dpr = 0
      z_dpr = h_offset
      z_bot_dpr = SUM(dsoil_dpr(:)) + h_offset
    ELSE
      i_dpr = 1
      z_dpr = dsoil_dpr(i_dpr)
      z_bot_dpr = SUM(dsoil_dpr(:))
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, z_dpr)    
    
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO i = 1, nrgrd     
      IF (z_rgrd_bot <= z_bot_lvt .AND. z_rgrd_bot <= z_bot_dpr .AND. & ! Above Bedrock in both tiles
          i_lvt      <= nsoil     .AND. i_dpr      <= nsoil) THEN           ! within soil column
        wtr_lat              = 0._wp
        h                    = z_rgrd_bot - z_rgrd_top
        a_inter              = len_inter * h
        wtr_soil_lvt_rgrd    = wtr_soil_lvt(i_lvt)   * h / dsoil_lvt(i_lvt)
        wsat_lvt_rgrd        = wsat_lvt(i_lvt)       * h / dsoil_lvt(i_lvt)
        field_cap_lvt_rgrd   = field_cap_lvt(i_lvt)  * h / dsoil_lvt(i_lvt)
        res_lvt_rgrd         = res_lvt(i_lvt)        * h / dsoil_lvt(i_lvt)
        mpot_lvt_rgrd        = mpot_lvt(i_lvt)
        k_lvt_rgrd           = k_lvt(i_lvt)

        IF (z_rgrd_bot > h_offset ) THEN ! Within the soil column of both tiles
          ! Flux is determined by the difference in matric potential between the two tiles
          wtr_soil_dpr_rgrd  = wtr_soil_dpr(i_dpr) * h / dsoil_dpr(i_dpr)
          wsat_dpr_rgrd      = wsat_dpr(i_dpr)  * h / dsoil_dpr(i_dpr)
          field_cap_dpr_rgrd = field_cap_dpr(i_dpr)  * h / dsoil_dpr(i_dpr)
          res_dpr_rgrd       = res_dpr(i_dpr)   * h / dsoil_dpr(i_dpr)
          mpot_dpr_rgrd      = mpot_dpr(i_dpr)
          k_dpr_rgrd         = k_dpr(i_dpr)
          k_inter            =(k_dpr_rgrd * k_lvt_rgrd)**0.5_wp
          d_mpot             =  abs(mpot_lvt_rgrd) - abs(mpot_dpr_rgrd)
          wtr_flx_lat        = - k_inter  * d_mpot / length_scale *  a_inter

          wtr_rel_qlbrm      = MIN(1._wp, (wtr_soil_lvt_rgrd * fr_lvt_rel       &
            &                +  wtr_soil_dpr_rgrd * fr_dpr_rel)      &
            &                / (wsat_lvt_rgrd * fr_lvt_rel + wsat_dpr_rgrd * fr_dpr_rel ))
          IF (wtr_flx_lat > 0._wp) THEN ! Flux is going from elevation to depression

            wtr_2_qlbrm      = MAX(0._wp, &
              &               (wtr_rel_qlbrm * wsat_dpr_rgrd - wtr_soil_dpr_rgrd) * fr_dpr_rel)
            free_space       = MAX(0._wp,(field_cap_dpr_rgrd - wtr_soil_dpr_rgrd) * fr_dpr_rel)
            wtr_avail        = MAX(0._wp,(wtr_soil_lvt_rgrd - res_lvt_rgrd) * fr_lvt_rel)
            wtr_lat          = MIN(wtr_flx_lat * dtime, free_space, wtr_avail, wtr_2_qlbrm)
            wtr_lat          = MAX(0._wp, wtr_lat)

          ELSE ! Flux is going from depression to elevation

            wtr_2_qlbrm      = MAX(0._wp, &
              &               (wtr_rel_qlbrm  * wsat_lvt_rgrd - wtr_soil_lvt_rgrd) * fr_lvt_rel)
            free_space       = MAX(0._wp, (field_cap_lvt_rgrd - wtr_soil_lvt_rgrd) * fr_lvt_rel)
            wtr_avail        = MAX(0._wp, (wtr_soil_dpr_rgrd - res_dpr_rgrd) * fr_dpr_rel)
            wtr_lat          = MAX(wtr_flx_lat * dtime, -free_space, -wtr_avail, -wtr_2_qlbrm)
            wtr_lat          = MIN(0._wp, wtr_lat)
          END IF
        END IF

        diffus_lvt(i_lvt)   = diffus_lvt(i_lvt) - wtr_lat / fr_lvt_rel

        IF (i_dpr > 0) THEN
          diffus_dpr(i_dpr) = diffus_dpr(i_dpr) + wtr_lat / fr_dpr_rel
        ENDIF

        IF (z_rgrd_bot == z_lvt) THEN
          i_lvt = i_lvt + 1
          IF (i_lvt <= nsoil) THEN
            IF (dsoil_lvt(i_lvt) > 0._wp) THEN
              z_lvt = z_lvt + dsoil_lvt(i_lvt)
            ELSE
              z_lvt = 10000._wp
              i_lvt = nsoil + 10
            END IF
          ELSE
            z_lvt = 10000._wp
            i_lvt = nsoil + 10
          END IF
        END IF

        IF (z_rgrd_bot == z_dpr) THEN
          i_dpr = i_dpr + 1
          IF (i_dpr <= nsoil) THEN
            IF (dsoil_dpr(i_dpr) > 0._wp) THEN
              z_dpr = z_dpr + dsoil_dpr(i_dpr)
            ELSE
              z_dpr = 10000._wp
              i_dpr = nsoil + 10
            END IF
          ELSE
            z_dpr = 10000._wp
            i_dpr = nsoil + 10
          END IF
        END IF      

        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=MIN(z_lvt, z_dpr)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! Adjusting the water contents and calculating the net- in- and outgoing fluxes for the tiles
    wtr_soil_lvt(:)   =  wtr_soil_lvt(:)             + diffus_lvt(:)
    wtr_soil_dpr(:)   =  wtr_soil_dpr(:)             + diffus_dpr(:)

    wtr_latflow_lvt   = wtr_latflow_lvt              + diffus_lvt(:)  * rhoh2o / dtime
    wtr_latflow_dpr   = wtr_latflow_dpr              + diffus_dpr(:)  * rhoh2o / dtime

    diffus_in_lvt(:)  = diffus_in_lvt(:)  + MAX(0._wp, diffus_lvt(:)) * rhoh2o / dtime
    diffus_in_dpr(:)  = diffus_in_dpr(:)  + MAX(0._wp, diffus_dpr(:)) * rhoh2o / dtime
    diffus_out_lvt(:) = diffus_out_lvt(:) - MIN(0._wp, diffus_lvt(:)) * rhoh2o / dtime
    diffus_out_dpr(:) = diffus_out_dpr(:) - MIN(0._wp, diffus_dpr(:)) * rhoh2o / dtime

  END SUBROUTINE calc_diffusion_water

! ===============================================================================================================================
  !>
  !! Compute intertile lateral diffusion (energy)
  !!
  SUBROUTINE calc_diffusion_energy(                    &
      &      nsoil, dtime,                             &
      &      dz, h_offset,                             &
      &      length_scale, len_inter,                  &
      &      lvt_fract, dpr_fract,                     &
      &      t_soil_lvt, t_soil_dpr,                   &
      &      heat_cond_soil_lvt, heat_cond_soil_dpr,   &
      &      heat_cap_soil_lvt, heat_cap_soil_dpr,     &
      &      diffus_in_lvt, diffus_in_dpr,             &
      &      diffus_out_lvt, diffus_out_dpr,           &
      &      heat_flx_lat_lvt, heat_flx_lat_dpr        &
      &      )

    INTEGER, INTENT(in)     :: nsoil                      !< Number of oil layers

    REAL(wp), INTENT(in)    :: dtime,                  &  !< Timestep length [s]
      &                        h_offset,               &  !< Differece in surface elevation between elevation and depression [m]
      &                        length_scale,           &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,              &  !< Length of interface between elevation and depression [m m-2]
      &                        lvt_fract,              &  !< Area fraction of elevation tile [/]
      &                        dpr_fract                  !< Area fraction of depression tile [/]

    REAL(wp), INTENT(in)    :: dz(:),                  &  !< Layer thickness [m]
      &                        t_soil_lvt(:),          &  !< Soil tempertures of elevation tile [K]
      &                        t_soil_dpr(:),          &  !< Soil tempertures of depression tile [K]
      &                        heat_cond_soil_lvt(:),  &  !< Soil heat conductivity of elevation tile [J m-1 s-1 K-1]
      &                        heat_cond_soil_dpr(:),  &  !< Soil heat conductivity of depression tile [J m-1 s-1 K-1]
      &                        heat_cap_soil_lvt(:),   &  !< Soil heat capacity of elevation tile [J m-3 K-1]
      &                        heat_cap_soil_dpr(:)       !< Soil heat capacity of depression tile [J m-3 K-1]

    REAL(wp), INTENT(inout) :: diffus_in_lvt(:),       &  !< Diffusion influx at elevation tile [J m-2 s-1]
      &                        diffus_in_dpr(:),       &  !< Diffusion outflux at depression tile [J m-2 s-1]
      &                        diffus_out_lvt(:),      &  !< Diffusion influx at elevation tile [J m-2 s-1]
      &                        diffus_out_dpr(:)          !< Diffusion outflux at depression tile [J m-2 s-1]

    REAL(wp), INTENT(inout) :: heat_flx_lat_lvt(:),    &  !< Net Heat flux at elevation tile [J m-2 s-1]
      &                        heat_flx_lat_dpr(:)        !< Net Heat flux at depression tile [J m-2 s-1]



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,                  &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, i_lvt, i_dpr            !< Index variables for layers on helper grid (i), in elevation (i_lvt) and depression tile (i_dpr)  

    REAL(wp)                :: fr_lvt_rel,             &  !< Relative fraction of elevation tile [/]
      &                        fr_dpr_rel,             &  !< Relative fraction of depression tile [/]
      &                        z_lvt,                  &  !< Layer bottom of current layer on elevation tile [m]
      &                        z_dpr,                  &  !< Layer bottom of current layer on depression tile [m]
      &                        z_rgrd_top,             &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,             &  !< Bottom of current layer on helper grid [m]
      &                        h,                      &  !< Height of current layer [m]
      &                        a_inter,                &  !< Interface area on current layer [m2 m-2]
      &                        t_lvt_rgrd,             &  !< Soil tempertures of elevation tile on current layer on helper grid [K]
      &                        t_dpr_rgrd,             &  !< Soil tempertures of depression tile on current layer on helper grid [K]
      &                        h_cond_soil_lvt_rgrd,   &  !< Soil heat conductivity of elevation tile on current layer on helper grid [J m-1 s-1 K-1]
      &                        h_cond_soil_dpr_rgrd,   &  !< Soil heat conductivity of depression tile on current layer on helper grid [J m-1 s-1 K-1]
      &                        h_cap_soil_lvt_rgrd,    &  !< Soil heat capacity of elevation tile on current layer on helper grid [J m-3 K-1]
      &                        h_cap_soil_dpr_rgrd,    &  !< Soil heat capacity of depression tile on current layer on helper grid [J m-3 K-1]
      &                        h_cond_inter,           &  !< Soil heat conductivity at interface on current layer on helper grid [J m-1 s-1 K-1]
      &                        d_t,                    &  !< Temperature difference on current layer on helper grid [K]
      &                        heat_flx_lat,           &  !< Energy flux (based on h_cond_inter and d_t) between depression and elevation [J m-2 s-1]
      &                        t_eq,                   &  !< Equilibrium temperature [K]
      &                        heat_flx_eq_lat            !< Equilirium energy flux between depression and elevation [J m-2 s-1]

    nrgrd         = nsoil * 2 + 3

    IF (h_offset > SUM(dz(:))) RETURN
    IF (lvt_fract < 1.0e-10_wp .OR. dpr_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    ! Start setting up the interpolated helper grid
    ! As ponded water / ice do not have a heat content we only account for energy exchange in the soil:
    i_lvt = 1
    z_lvt = dz(i_lvt)
    i_dpr = 1
    IF (h_offset > 0._wp) THEN
      i_dpr = 0
      z_dpr = h_offset
    ELSE
      z_dpr = dz(i_dpr)
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, z_dpr)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO i = 1, nrgrd
      IF (i_lvt <= nsoil .AND. i_dpr <= nsoil) THEN ! Within the soil column of both tiles

        heat_flx_lat           = 0._wp

        IF (z_rgrd_bot > h_offset ) THEN ! Within the soil column of both tiles
          h                    = z_rgrd_bot - z_rgrd_top
          a_inter              = len_inter * h

          t_lvt_rgrd           = t_soil_lvt(i_lvt)
          h_cond_soil_lvt_rgrd = heat_cond_soil_lvt(i_lvt)
          h_cap_soil_lvt_rgrd  = heat_cap_soil_lvt(i_lvt)

          t_dpr_rgrd           = t_soil_dpr(i_dpr)
          h_cond_soil_dpr_rgrd = heat_cond_soil_dpr(i_dpr)
          h_cap_soil_dpr_rgrd  = heat_cap_soil_dpr(i_dpr)

          ! Flux is determined by the difference in temperature between the two tiles
          h_cond_inter         =(h_cond_soil_lvt_rgrd * h_cond_soil_dpr_rgrd)**0.5_wp
          d_t                  = t_lvt_rgrd - t_dpr_rgrd
          heat_flx_lat         = h_cond_inter  * d_t / length_scale *  a_inter

          ! Limiting heat flux to avoid instablilities for long time steps
          ! Determine flux required to equilibriate Temps within timestep
          ! Equilibrium temperature
          t_eq = (t_lvt_rgrd * h_cap_soil_lvt_rgrd * fr_lvt_rel  &
            &  +  t_dpr_rgrd * h_cap_soil_dpr_rgrd * fr_dpr_rel) &
            &  / (h_cap_soil_lvt_rgrd * fr_lvt_rel + h_cap_soil_dpr_rgrd * fr_dpr_rel)
          ! Heat flux to reach equilibrium temp
          heat_flx_eq_lat = MIN( &
            &               ABS((t_lvt_rgrd - t_eq)  * h_cap_soil_lvt_rgrd * fr_lvt_rel), &
            &               ABS((t_dpr_rgrd - t_eq)  * h_cap_soil_dpr_rgrd * fr_dpr_rel)) &
            &               *  h / dtime
          ! Limiting heat flux to equilibrium flux
          IF (heat_flx_lat > 0._wp .AND. heat_flx_lat > ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = ABS(heat_flx_eq_lat)
          ELSE IF (heat_flx_lat < 0._wp .AND. heat_flx_lat < -ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = -ABS(heat_flx_eq_lat)
          END IF

          ! Diagnosing incoming-, outgoing- and net fluxes on tiles
          diffus_in_lvt(i_lvt)    = diffus_in_lvt(i_lvt)  + MAX(0._wp, heat_flx_lat / fr_lvt_rel)
          diffus_in_dpr(i_dpr)    = diffus_in_dpr(i_dpr)  + MAX(0._wp, heat_flx_lat / fr_dpr_rel)
          diffus_out_lvt(i_lvt)   = diffus_out_lvt(i_lvt) - MIN(0._wp, heat_flx_lat / fr_lvt_rel)
          diffus_out_dpr(i_dpr)   = diffus_out_dpr(i_dpr) - MIN(0._wp, heat_flx_lat / fr_dpr_rel)

          heat_flx_lat_lvt(i_lvt) = heat_flx_lat_lvt(i_lvt)          - heat_flx_lat / fr_lvt_rel
          heat_flx_lat_dpr(i_dpr) = heat_flx_lat_dpr(i_dpr)          + heat_flx_lat / fr_dpr_rel
        END IF

        ! Determening the top and bottem of the next layer of the regridded
        IF (z_rgrd_bot == z_lvt) THEN
          i_lvt = i_lvt + 1
          IF (i_lvt <= nsoil) THEN
            z_lvt = z_lvt + dz(i_lvt)
          END IF
        END IF

        IF (z_rgrd_bot == z_dpr) THEN
          i_dpr = i_dpr + 1
          IF (i_dpr <= nsoil) THEN
            z_dpr = z_dpr + dz(i_dpr)
          END IF
        END IF
        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=MIN(z_lvt, z_dpr)
      END IF
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE calc_diffusion_energy



  ! ===============================================================================================================================
  !>
  !! Compute intertile snow redistribution
  !!
  SUBROUTINE calc_snow_redistribution(              &
      &      l_dynsnow,                             &
      &      dtime,                                 &
      &      h_offset,                              &
      &      lvt_fract, dpr_fract,                  &
      &      weq_pond_lvt, weq_pond_dpr,            &
      &      snow_rdst_in_lvt, snow_rdst_in_dpr,    &
      &      snow_rdst_out_lvt, snow_rdst_out_dpr,  &
      &      wtr_latflow_lvt, wtr_latflow_dpr,      &
      &      snow_soil_lvt, snow_soil_dpr,          &
      &      snow_new_lvt, snow_new_dpr,            &
      &      snow_lvt, snow_dpr,                    &
      &      snow_dens_lvt, snow_dens_dpr           &
      &      )


    USE mo_jsb_physical_constants, ONLY: rhoh2o, dens_snow_min, dens_snow_max, dens_snow

    LOGICAL,  INTENT(in)    :: l_dynsnow               !< Dynamic snow density is being used

    REAL(wp), INTENT(in)    :: dtime,               &  !< Timestep length [s]
      &                        h_offset,            &  !< Differece in surface elevation between elevation and depression [m]
      &                        lvt_fract,           &  !< Area fraction of elevation tile [/]
      &                        dpr_fract               !< Area fraction of depression tile [/]

    REAL(wp), INTENT(in)    :: weq_pond_lvt,        &  !< Pond storage on elevation tile [m]
      &                        weq_pond_dpr            !< Pond storage on depression tile [m]

    REAL(wp), INTENT(inout) :: snow_rdst_in_lvt,    &  !< Snow redistribution influx at elevation tile [kg m-2 s-1]
      &                        snow_rdst_in_dpr,    &  !< Snow redistribution influx at depression tile [kg m-2 s-1]
      &                        snow_rdst_out_lvt,   &  !< Snow redistribution outflux at elevation tile [kg m-2 s-1]
      &                        snow_rdst_out_dpr       !< Snow redistribution outflux at depression tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: &
      &                        wtr_latflow_lvt,     &  !< Netflux (surface) on elevation tile [kg m-2 s-1]  
      &                        wtr_latflow_dpr,     &  !< Netflux (surface) on depression tile [kg m-2 s-1]
      &                        snow_soil_lvt,       &  !< Snow soil (WEQ) on elevation tile [m]
      &                        snow_soil_dpr,       &  !< Snow soil (WEQ) on depression tile [m]
      &                        snow_new_lvt,        &  !< New snow soil (WEQ) on elevation tile [m]
      &                        snow_new_dpr,        &  !< New snow soil (WEQ) on depression tile [m]
      &                        snow_lvt,            &  !< Snow (WEQ; can and soil) on elevation tile [m]
      &                        snow_dpr,            &  !< Snow (WEQ; can and soil) on depression tile [m]
      &                        snow_dens_lvt,       &  !< Snow density on elevation tile [kg m-3]
      &                        snow_dens_dpr           !< Snow density on depression tile [kg m-3]

    ! Local variables ------------------------------------------------------------------------------------

    REAL(wp)                :: fr_lvt_rel,          &  !< Relative fraction of elevation tile [/]
      &                        fr_dpr_rel,          &  !< Relative fraction of depression tile [/]
      &                        h_offset_snow_fix,   &  !< Differece in surface elevation including immobile snow layer between elevation and depression [m]
      &                        snow_fix_lvt ,       &  !< Immobile snow on elevation [m]
      &                        snow_fix_dpr ,       &  !< Immobile snow on depression [m]
      &                        snow_fix_rel_lvt,    &  !< Fraction of immobile snow on elevation [/]
      &                        snow_fix_rel_dpr,    &  !< Fraction of immobile snow on depression [/]
      &                        snow_depth_mob_avg,  &  !< Average depth of mobile snow [m]
      &                        snow_depth_mob_lvt,  &  !< Depth of mobile snow on elevation [m]
      &                        snow_depth_mob_dpr,  &  !< Depth of mobile snow on depression [m]
      &                        snow_mob_lvt,        &  !< Mobile snow (WEQ) on elevation [m]
      &                        snow_mob_dpr,        &  !< Mobile snow (WEQ)  on depression [m]
      &                        snow_mob_upd_lvt,    &  !< Updated mobile snow (WEQ) on elevation [m]
      &                        snow_mob_upd_dpr,    &  !< Updated mobile snow (WEQ) on depression [m]
      &                        snow_dens_max,       &  !< Maximum snow density [kg m-3]
      &                        snow_dens_min,       &  !< Minimum snow density [kg m-3]
      &                        snow_dens_fix_lvt,   &  !< Density of immobile snow on elevation [kg m-3]
      &                        snow_dens_fix_dpr,   &  !< Density of immobile snow on depression [kg m-3]
      &                        snow_rdst_lvt,       &  !< Snow (WEQ) redistribution flux on elevation [kg m-2 s-1]
      &                        snow_rdst_dpr           !< Snow (WEQ) redistribution flux on depression [kg m-2 s-1]

    IF (lvt_fract < 1.0e-10_wp .OR. dpr_fract < 1.0e-10_wp) RETURN

    ! Setting initial values to 0.0
    snow_rdst_lvt = 0._wp
    snow_rdst_dpr = 0._wp

    ! Snow density paraeters depending on snow scheme
    IF (l_dynsnow) THEN
      snow_dens_max = dens_snow_max
      snow_dens_min = dens_snow_min
    ELSE
      snow_dens_max = dens_snow
      snow_dens_min = dens_snow
    END IF

    ! Determine relative fractions from absolute fractions
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    ! Determine amount of fixed and mobile fractions of the snow cover
    ! (assuming that in the previous time step the old snow has melted before the new snow)
    snow_fix_lvt =  snow_soil_lvt - snow_new_lvt
    snow_fix_dpr =  snow_soil_dpr - snow_new_dpr
    snow_mob_lvt =  snow_new_lvt
    snow_mob_dpr =  snow_new_dpr

    IF (snow_fix_lvt < 0._wp) THEN
      snow_mob_lvt = snow_new_lvt + snow_fix_lvt
      snow_fix_lvt = 0._wp
    END IF
    IF (snow_fix_dpr < 0._wp) THEN
      snow_mob_dpr = snow_new_dpr + snow_fix_dpr
      snow_fix_dpr = 0._wp
    END IF

    ! Determine snow density of fixed snow
    IF (snow_fix_lvt > 0._wp) THEN
      snow_fix_rel_lvt  = snow_fix_lvt / snow_soil_lvt
      snow_dens_fix_lvt = (snow_dens_lvt - (snow_dens_min * (1._wp - snow_fix_rel_lvt))) &
        &               / snow_fix_rel_lvt
    ELSE
      snow_dens_fix_lvt = snow_dens_lvt
    END IF
    IF (snow_fix_dpr > 0._wp) THEN
      snow_fix_rel_dpr  = snow_fix_dpr / snow_soil_dpr
      snow_dens_fix_dpr = (snow_dens_dpr - (snow_dens_min * (1._wp - snow_fix_rel_dpr))) &
        &               / snow_fix_rel_dpr
    ELSE
      snow_dens_fix_dpr = snow_dens_dpr
    END IF

    ! Determine height differences of the top of fixed snow
    h_offset_snow_fix = h_offset                                                   &
      &               + weq_pond_lvt + snow_fix_lvt * (rhoh2o / snow_dens_fix_lvt) &
      &               - weq_pond_dpr - snow_fix_dpr * (rhoh2o / snow_dens_fix_dpr)

    ! Calculate snow redistribution for two cases (1) from elevation to depression or (2) from depression to elevation
    ! Transport only applies if in the source tile there is fixed snow
    IF  ((h_offset_snow_fix >  1.0e-10_wp .AND. snow_fix_lvt > 0.01_wp) .OR. &
         (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_dpr > 0.01_wp))THEN

      snow_depth_mob_avg = (snow_mob_lvt * fr_lvt_rel + snow_mob_dpr * fr_dpr_rel) &
        &                * (rhoh2o / snow_dens_min)

      ! 1. (Fixed) snow surface is lower in derpression and there is fixed and mobile snow on elevation
      IF  (h_offset_snow_fix > 1.0e-10_wp .AND. snow_fix_lvt > 0.01_wp)THEN
        IF (snow_depth_mob_avg / fr_dpr_rel <= h_offset_snow_fix ) THEN ! New snow fits into depression
          snow_depth_mob_dpr = snow_depth_mob_avg / fr_dpr_rel
          snow_depth_mob_lvt = 0._wp
        ELSE
          snow_depth_mob_dpr = h_offset_snow_fix  &
            &                + (snow_depth_mob_avg - h_offset_snow_fix * fr_dpr_rel)
          snow_depth_mob_lvt =  snow_depth_mob_avg - h_offset_snow_fix * fr_dpr_rel
        END IF
      ! 2. (Fixed) snow surface is lower on elevation and there is fixed snow in depression
      ELSE IF (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_dpr > 0.01_wp) THEN
        h_offset_snow_fix = - h_offset_snow_fix ! Inversing sign of height offset to simplify calculations
        IF (snow_depth_mob_avg / fr_lvt_rel <= h_offset_snow_fix ) THEN  ! New snow fits into elevation
          snow_depth_mob_lvt = snow_depth_mob_avg / fr_lvt_rel
          snow_depth_mob_dpr = 0._wp
        ELSE
          snow_depth_mob_lvt = h_offset_snow_fix &
                           & + (snow_depth_mob_avg - h_offset_snow_fix * fr_lvt_rel)
          snow_depth_mob_dpr =  snow_depth_mob_avg - h_offset_snow_fix * fr_lvt_rel
        END IF
      END IF

      ! Recalculation of snow weqs
      snow_mob_upd_dpr = snow_depth_mob_dpr / (rhoh2o / snow_dens_min)
      snow_mob_upd_lvt = snow_depth_mob_lvt / (rhoh2o / snow_dens_min)

      snow_rdst_dpr    = snow_mob_upd_dpr - snow_mob_dpr
      snow_rdst_lvt    = snow_mob_upd_lvt - snow_mob_lvt

      snow_soil_dpr    = snow_soil_dpr    + snow_rdst_dpr
      snow_soil_lvt    = snow_soil_lvt    + snow_rdst_lvt

      snow_dpr         = snow_dpr         + snow_rdst_dpr
      snow_lvt         = snow_lvt         + snow_rdst_lvt

      snow_new_dpr     = snow_new_dpr     + snow_rdst_dpr
      snow_new_lvt     = snow_new_lvt     + snow_rdst_lvt

      ! Converting redistribution flux from m/step -> kg m-2 s-1
      snow_rdst_dpr     = snow_rdst_dpr * rhoh2o / dtime
      snow_rdst_lvt     = snow_rdst_lvt * rhoh2o / dtime

      wtr_latflow_lvt   = wtr_latflow_lvt   + snow_rdst_lvt
      wtr_latflow_dpr   = wtr_latflow_dpr   + snow_rdst_dpr

      snow_rdst_in_lvt  = snow_rdst_in_lvt  + MAX(0._wp, snow_rdst_lvt)
      snow_rdst_in_dpr  = snow_rdst_in_dpr  + MAX(0._wp, snow_rdst_dpr)
      snow_rdst_out_lvt = snow_rdst_out_lvt - MIN(0._wp, snow_rdst_lvt)
      snow_rdst_out_dpr = snow_rdst_out_dpr - MIN(0._wp, snow_rdst_dpr)

      ! Updating snow densities based on redistributed snow
      IF (snow_rdst_dpr .NE. 0_wp .AND. snow_soil_dpr > 0._wp) THEN
        snow_dens_dpr   = (snow_new_dpr * snow_dens_min &
          &             + snow_fix_dpr * min(snow_dens_fix_dpr, snow_dens_max)) / snow_soil_dpr
      END IF

      IF (snow_rdst_lvt .NE. 0_wp .AND. snow_soil_lvt > 0._wp) THEN
        snow_dens_lvt   = (snow_new_lvt * snow_dens_min &
          &             + snow_fix_lvt * min(snow_dens_fix_lvt, snow_dens_max)) / snow_soil_lvt
      END IF
    END IF

  END SUBROUTINE calc_snow_redistribution

  !  
  ! ===============================================================================================================================
  !>
  !! Limit soil water and pond storage content
  !!
  SUBROUTINE remove_exess_water(              &    
    & weq_pond_max,                           &
    & ice_pond,                               &
    & wtr_pond,                               &
    & wtr_latflow_res_srf,                    &
    & nsoil,                                  &
    & wtr_soil_sat_sl,                        &
    & wtr_soil_sl                             &
    )
     
    REAL(wp), INTENT(in)    :: weq_pond_max,  &              !< Maximum pond water storage [m]
      &                        ice_pond                      !< Ice in pond water storage [m]

    REAL(wp), INTENT(inout) :: wtr_pond,      &              !< Liquid water in pond water storage [m]
      &                        wtr_latflow_res_srf           !< Surface water retention storage [m]


    INTEGER , OPTIONAL, INTENT(in)    :: nsoil               !< Number of soil ayers
    REAL(wp), OPTIONAL, INTENT(in)    :: wtr_soil_sat_sl(:)  !< Saturation soil moisture on soil layers [m]
    REAL(wp), OPTIONAL, INTENT(inout) :: wtr_soil_sl(:)      !< Water content on soil layers [m]]


    INTEGER                 :: i                             !< Counter for loop over soil layers
    REAL(wp)                :: wtr_pond_max                  !< Maximum liquid water in pond storage [m]
    
    wtr_pond_max = MAX(0.0_wp, weq_pond_max - ice_pond)

    IF (PRESENT(nsoil)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO i=nsoil,2, -1
        IF (wtr_soil_sl(i) > wtr_soil_sat_sl(i)) THEN
          wtr_soil_sl(i-1) = wtr_soil_sl(i-1) &
            &              + wtr_soil_sl(i) - wtr_soil_sat_sl(i)
          wtr_soil_sl(i)   = wtr_soil_sat_sl(i)
        END IF
      END DO
      !$ACC END PARALLEL

      IF (wtr_soil_sl(1) > wtr_soil_sat_sl(1)) THEN
        wtr_pond = wtr_pond &
          &      + wtr_soil_sl(1) - wtr_soil_sat_sl(1)
        wtr_soil_sl(1) = wtr_soil_sat_sl(1)
      END IF
    END IF

    IF (wtr_pond > wtr_pond_max) THEN
      wtr_latflow_res_srf = wtr_latflow_res_srf &
        &                 + wtr_pond - wtr_pond_max
      wtr_pond = wtr_pond_max
    END IF

  END SUBROUTINE remove_exess_water

  !  
  ! ===============================================================================================================================
  !>
  !! Adjust lateral fluxes according to differences in surface elevation and bedrock boundary
  !!
  SUBROUTINE shift_fluxes_sl(                 &
    & nsoil,                                  &
    & h_offset,                               &
    & dsoil_lvt, dsoil_dpr,                   &   
    & flx_lat_sl                              &
    )
     

    INTEGER , INTENT(in)    :: nsoil             !< Number of soil ayers


    REAL(wp), INTENT(in)    :: h_offset          !< Differece in surface elevation between elevation and depression [m]

    REAL(wp), INTENT(in)    :: dsoil_lvt(:),  &  !< Soil depth on layers of elevation tile [m]
      &                        dsoil_dpr(:)      !< Soil depth on layers of depression tile [m]

    REAL(wp), INTENT(inout) :: flx_lat_sl(:)     !< Lateralal flux to be adjusted [?]


    INTEGER                 :: nrgrd,         &  !< Number of layers on helper grid determined by layer overlaps
      &                        i,             &  !< Counter variable
      &                        i_lvt, i_dpr      !< Index variables for layers on in elevation (i_lvt) and in depression (i_dpr)

    REAL(wp)                :: z_lvt,         &  !< Layer bottom of current layer on elevation tile [m]
      &                        z_dpr,         &  !< Layer bottom of current layer on depression tile [m]
      &                        z_rgrd_top,    &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,    &  !< Bottom of current layer on helper grid [m]
      &                        h                 !< Height of current layer [m]


    REAL(wp)                :: flx_lat_sl_upd(nsoil) !< Updated lateral flux [?]

    flx_lat_sl_upd(:) = 0._wp

    nrgrd       = nsoil * 2 + 1    

    i_lvt = 0
    z_lvt = dsoil_lvt(1)

    i_dpr = 0
    z_dpr = dsoil_dpr(1) + h_offset

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, z_dpr)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO i  = 1, nrgrd
      h                     = z_rgrd_bot-z_rgrd_top
      flx_lat_sl_upd(i_dpr) = flx_lat_sl_upd(i_dpr) + flx_lat_sl(i_lvt) * h / dsoil_lvt(i_lvt)

      IF (z_rgrd_bot == z_lvt .AND. z_lvt < SUM(dsoil_lvt(:))) THEN
        i_lvt = i_lvt + 1
        IF (i_lvt <= nsoil) THEN
          z_lvt = z_lvt + dsoil_lvt(i_lvt)
        END IF
      END IF

      IF (z_rgrd_bot == z_dpr .AND. z_dpr < SUM(dsoil_dpr(:)) + h_offset) THEN
        i_dpr = i_dpr + 1
        IF (i_dpr <= nsoil) THEN
          z_dpr = z_dpr + dsoil_dpr(i_dpr)
        END IF
      ELSE IF (z_rgrd_bot >= SUM(dsoil_dpr(:)) + h_offset) THEN
        z_dpr = z_lvt
      END IF
      
      z_rgrd_top=z_rgrd_bot
      z_rgrd_bot=MIN(z_lvt, z_dpr)

    END DO
    !$ACC END PARALLEL   

    flx_lat_sl(:) = flx_lat_sl_upd(:)

  END SUBROUTINE shift_fluxes_sl


  !
  ! ===============================================================================================================================
  !>
  !! Compute geometric relationships for tile and connected sibling
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE get_geometry(                      &
    & microcon_slp_tile,                        &
    & microcon_dfsn_tile,                       &
    & microcon_dfsn_sibl,                       &
    & macro_connect,                            &
    & fract, fract_sibl,                        &
    & fract_nclsd, fract_nclsd_sibl,            &
    & area, area_sibl,                          &
    & area_nclsd, area_nclsd_sibl,              &
    & h_surf, h_surf_sibl,                      &
    & h_offset,                                 &
    & length_scale,                             &
    & length_iface,                             &
    & micro_slope                               & 
    )

    USE mo_jsb_math_constants,     ONLY: pi

    REAL(wp), INTENT(in)    :: microcon_slp_tile,   &  !< Micro-scale connectivity between tile and sibling  along slope[/]
      &                        microcon_dfsn_tile,  &  !< Micro-scale connectivity between tile and sibling (if tile is encompased by sibl) [/]
      &                        microcon_dfsn_sibl,  &  !< Micro-scale connectivity between sibling and tile (if sibl is encompased by tile) [/]
      &                        macro_connect           !< Macro-scale connectivity between tile and sibling [/]

    REAL(wp), INTENT(in)    :: fract,               &  !< Cover fraction of tile [/]
      &                        fract_sibl,          &  !< Cover fraction of sibling [/]
      &                        fract_nclsd,         &  !< Cover fraction enclosed by tile [/]
      &                        fract_nclsd_sibl,    &  !< Cover fraction enclosed by sibling [/]
      &                        area,                &  !< Characteristic area of tile [m2]
      &                        area_sibl,           &  !< Characteristic area of sibling [m2]
      &                        area_nclsd,          &  !< Characteristic area enclosed by tile [m2]
      &                        area_nclsd_sibl,     &  !< Characteristic area enclosed by sibling [m2]
      &                        h_surf,              &  !< Surface elevation of tile [/]
      &                        h_surf_sibl             !< Surface elevation of sibling [/]


    REAL(wp), INTENT(out)   :: h_offset,            &  !< Differece in surface elevation between tile and siblig [m]
      &                        length_scale,        &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        length_iface,        &  !< Length of interface between tile and sibling [m m-2]
      &                        micro_slope             !< Micro-scale slope between tile and sibling [/]

    REAL(wp)                :: fract_rel,           &  !< Relative fraction of tile [/]
      &                        fract_rel_sibl,      &  !< Relative fraction of sibling [/]
      &                        rad_in,              &  !< Radius of inner circle [m]
      &                        rad_out,             &  !< Radius outer circle [m]
      &                        area_cmbnd,          &  !< Combined area of nested tiles [m2]
      &                        fract_rel_cmbnd,     &  !< Relative fraction of nested tiles [/]
      &                        area_fact               !< Conversion factor to get from area of inner to that of outer tile

    !$ACC ROUTINE SEQ

    ! Calculate difference in surface elevation between tile and sibling
    h_offset = h_surf - h_surf_sibl

    ! Calculate relative fractions
    IF ((fract + fract_sibl) == 0._wp) THEN
      fract_rel      = 0._wp
      fract_rel_sibl = 0._wp
    ELSE
      fract_rel      = fract      / (fract + fract_sibl)
      fract_rel_sibl = fract_sibl / (fract + fract_sibl)
    END IF      

    ! Calculate characteristic length scale and interface length between tile and sibling:
    !
    ! Assumption is that the smaller tile constitutes a circle (C1) which allows the calulation of the 
    ! interface length as the minor arc whose fraction of the circle's diameter is equal to the connected fraction.
    !
    ! The smaller circle is encompassed by a larger circle (C2) which contains the areas of both tiles and whose radius is used 
    ! to approximate the characteristic length scale. 
    ! The latter is assumed to be equal to half the distance between circumference of C1 and C2 added to the radius of C1.
    ! The area of C2 is determined so that a segment with a central angle (X/360) equal to the connected fraction contains an area
    ! equal to that of the larger tile and the area of a segment of C1 with the same central angle. 
    IF ((microcon_dfsn_tile == 0._wp .AND. microcon_dfsn_sibl == 0._wp) .OR.  &
      &  fract_rel          == 0._wp .OR.  fract_rel_sibl == 0._wp      .OR.  &
      &  area               == 0._wp .OR.  area_sibl      == 0._wp) THEN
      length_scale = 0._wp
      length_iface = 0._wp

    ! Areas of outer tiles are not used directly, but approximated by: 
    ! outer_area    = inner_area * ( 1. / fract_rel_inner - 1.)
    ! combined_area = inner_area *   1. / fract_rel_inner
    ! this is done to be (at least) consistent with the cover fractions and 
    ! the assumption of having basic geometry of nested circles.
    !
    ! Interface length is given [as m per m2] relative to the area containing both tiles
    ! i.e. scaled with fract_rel / area or fract_rel_sibl / area_sibl
    !
    ELSE IF (microcon_dfsn_tile > microcon_dfsn_sibl) THEN ! Tile is encompassed by sibling
      IF (area_nclsd > 0._wp .AND. fract_nclsd > 1.0e-10_wp) THEN ! Tile encloses another tile

        area_cmbnd      = area_nclsd &
          &             * (fract_nclsd + fract) / fract_nclsd     
        fract_rel_cmbnd = (fract_nclsd + fract) &
          &             / (fract_nclsd + fract + fract_sibl)

        area_fact = 1._wp / fract_rel_cmbnd - 1._wp

        length_iface = ((4._wp * pi * area_cmbnd) ** 0.5_wp) &
          &          * microcon_dfsn_tile * fract_rel / area

        rad_in  = (area_nclsd    / pi)**0.5_wp

        rad_out = (area_cmbnd * (1._wp + area_fact / microcon_dfsn_tile) / pi)**0.5_wp

        length_scale = (rad_in + rad_out) / 2._wp

      ELSE  ! Tile does not enclose other tile

        area_fact = 1._wp / fract_rel - 1._wp

        length_iface = ((4._wp * pi * area) ** 0.5_wp) & 
          &          * microcon_dfsn_tile * fract_rel / area 

        rad_in  = (area / pi)**0.5_wp

        rad_out = (area *(1._wp + area_fact / microcon_dfsn_tile) / pi)**0.5_wp
               
        length_scale = rad_in + (rad_out - rad_in) / 2._wp
      END IF        
    ELSE  ! Sibling is encompassed by tile
      IF (area_nclsd_sibl > 0._wp .AND. fract_nclsd_sibl > 1.0e-10_wp) THEN  ! Sibling encloses another tile

        area_cmbnd      = area_nclsd_sibl &
          &             * (fract_nclsd_sibl + fract_sibl) / fract_nclsd_sibl

        fract_rel_cmbnd = (fract_nclsd_sibl + fract_sibl) &
          &             / (fract_nclsd_sibl + fract + fract_sibl)

        area_fact = 1._wp / fract_rel_cmbnd - 1._wp

        length_iface = ((4._wp * pi * area_cmbnd) ** 0.5_wp) & 
          &          * microcon_dfsn_sibl * fract_rel_sibl / area_sibl 

        rad_in  = (area_nclsd_sibl / pi)**0.5_wp

        rad_out = (area_cmbnd *(1._wp + area_fact / microcon_dfsn_sibl) / pi)**0.5_wp

        length_scale = (rad_in + rad_out) / 2._wp

      ELSE   ! Sibling does not enclose other tile

        area_fact = 1._wp / fract_rel_sibl - 1._wp
 
        length_iface = ((4._wp * pi * area_sibl) ** 0.5_wp) & 
          &          * microcon_dfsn_sibl * fract_rel_sibl / area_sibl 

        rad_in  = (area_sibl / pi)**0.5_wp

        rad_out = (area_sibl *(1._wp + area_fact / microcon_dfsn_sibl) / pi)**0.5_wp

        length_scale = rad_in + (rad_out - rad_in) / 2._wp

      END IF         
    END IF

    ! calculate micro-scale slope based on difference in surfce elevation and characteristic length scale    
    IF (microcon_slp_tile > 0._wp .AND. length_scale > 0._wp) THEN
      micro_slope = MIN(0.5_wp, h_offset / length_scale)
    ELSE
      micro_slope = 0._wp
    END IF
  END SUBROUTINE get_geometry

  !
  ! ===============================================================================================================================
  !>
  !! Compute lag factors for lateral exchange
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE get_flow_lags(                          &
    & lex_idt, lex_ids,                              &
    & velocity_srf, velocity_blg,                    &
    & ret_micro_uncon_srf, ret_micro_uncon_blg,      &
    & ret_macro_srf, ret_macro_blg,                  &
    & length_scale,                                  &
    & micro_connect, macro_connect,                  &
    & fract, fract_sibl,                             &
    & fract_microcon, fract_macrocon,                &
    & fact_lag_micro_srf, fact_lag_micro_blg,        &
    & fact_lag_macro_srf, fact_lag_macro_blg         &
    )

    USE mo_jsb_math_constants,     ONLY: pi

    INTEGER,  INTENT(in)    :: lex_idt,              &  !< lex id of tile
      &                        lex_ids                  !< lex id of sibling

    REAL(wp), INTENT(in)    ::                       &
      &                        velocity_srf,         &  !< Potential velocity of lateral flow through connected flow pathways [m s-1]
      &                        velocity_blg,         &  !< Potential velocity of lateral flow at the surface [m s-1]
      &                        ret_micro_uncon_srf,  &  !< Retention time for micro-scale fluxes at the surfcace, when no connections exist[s]
      &                        ret_micro_uncon_blg,  &  !< Retention time for micro-scale fluxes below the surfcace, when no connections exist [s]
      &                        ret_macro_srf,        &  !< Retention time for macro-scale fluxes at the surfcace [s]
      &                        ret_macro_blg            !< Retention time for macro-scale fluxes below the surfcace [s]

    REAL(wp), INTENT(in)    :: length_scale,         &  !< Assumed distance between centers of micro-scale connected tiles [m]
      &                        micro_connect,        &  !< Micro-scale connectivity between tile and sibling [/]
      &                        macro_connect,        &  !< Macro-scale connectivity between tile and sibling [/]
      &                        fract,                &  !< cover fraction of tile [/]
      &                        fract_sibl,           &  !< cover fraction of sibling [/]
      &                        fract_microcon,       &  !< (Weighted) Sum of micro-scale connection fractions on tile [/]
      &                        fract_macrocon           !< (Weighted) Sum of macro-scale connection fractions on tile [/]


    REAL(wp), INTENT(out)   :: fact_lag_micro_srf,   &  !< lag factor for micro-scale surface fluxes between tile and sibling [s-1]  
                               fact_lag_micro_blg,   &  !< lag factor for micro-scale below-ground fluxes between tile and sibling [s-1] 
                               fact_lag_macro_srf,   &  !< lag factor for macro-scale surface fluxes between tile and sibling [s-1]
                               fact_lag_macro_blg       !< lag factor for macro-scale below-ground fluxes between tile and sibling [s-1] 

   ! Local variables
   !
   REAL(wp)                 ::                       &
      &                        ret_micro_srf,        &  !< Retention time for micro-scale fluxes at the surfcace [s]
      &                        ret_micro_blg            !< Retention time for micro-scale fluxes below the surfcace [s]

    !$ACC ROUTINE SEQ

    ! Micro scale lag factors are determined by the (inverse) retention time, the micro-connectivity between tile and sibling, 
    ! and the cover fraction of sibling relative to the total connected farction of tile

    ret_micro_srf = MAX(length_scale / velocity_srf, 3600._wp)
    ret_micro_blg = MAX(length_scale / velocity_blg, 3600._wp)

    IF (micro_connect > 0._wp) THEN
      fact_lag_micro_srf = 1._wp / ret_micro_srf &
        &                * micro_connect * fract_sibl / fract_microcon
      fact_lag_micro_blg = 1._wp / ret_micro_blg &
        &                * micro_connect * fract_sibl / fract_microcon
    ELSE
      fact_lag_micro_srf = 0.0_wp
      fact_lag_micro_blg = 0.0_wp   
    END IF      

    ! Micro scale lag factors are determined by the (inverse) retention time, the micro-connectivity between tile and sibling, 
    ! and the cover fraction of sibling relative to the total connected farction of tile.
    ! Additionally, the lag foctor is decreased accoring to the micro-scale connectivity, assuming that the fractions of the tile 
    ! that have a micro-scale connectivity do not additionally feature macro-scale connections. 
    IF ((macro_connect > 0._wp) .AND. (fract > 0._wp)) THEN

      fact_lag_macro_srf = 1._wp / ret_macro_srf                       &
        &                * macro_connect * fract_sibl / fract_macrocon &
        &                * MAX(0._wp, (fract - fract_microcon) / fract)


      fact_lag_macro_blg = 1._wp / ret_macro_blg                       &
        &                * macro_connect * fract_sibl / fract_macrocon &
        &                * MAX(0._wp, (fract - fract_microcon) / fract)


    ELSE
      fact_lag_macro_srf = 0.0_wp
      fact_lag_macro_blg = 0.0_wp   
    END IF  

    ! For outflow on tiles that have no downstream tile, lag factors are calculated based on the time-scales of the macro-scale retention times.
    ! The respective fluxes are given to surface runoff and drainage.
    IF ((lex_idt == lex_ids)       .AND. &
      & (fract_microcon == 0.0_wp) .AND. (fract_macrocon == 0.0_wp)) THEN
      fact_lag_micro_srf = 1._wp / ret_micro_uncon_srf
      fact_lag_micro_blg = 1._wp / ret_micro_uncon_blg
      fact_lag_macro_srf = 0._wp
      fact_lag_macro_blg = 0._wp     
    END IF

  END SUBROUTINE get_flow_lags

  !  
  ! ===============================================================================================================================
  !>
  !! Check wether matrix contains closed loop
  !!
  SUBROUTINE check_flow_matrix(               &    
    & nconti,                                 &
    & matrix,                                 &
    & has_loop                                &
    )
     
    INTEGER,  INTENT(in)    :: nconti              !< number of connected tiles

    REAL(wp), INTENT(inout) :: matrix(:,:)         !< Flow matrix

    LOGICAL,  INTENT(out)   :: has_loop            !< .TRUE. if matrix contains closed loop


    INTEGER                 :: it, is, cntr,  &    !< Counter variables for loop over soil layers
      &                        src_id, trgt_id     !< Id of current source and target tile

    has_loop = .FALSE.

    DO it = 1, nconti
      DO is = 1, nconti
        IF (matrix(is,it) > 0._wp) THEN ! sibling connects to tile
          IF (SUM(matrix(:,is)) > 0._wp) THEN ! sibling has tiles connected to it
            trgt_id=is ! sibling is currenttarget          
            src_id=MAXLOC(matrix(:,is), 1)  ! getting the largest incomming connection
            cntr=0  ! Starting counter for closed loop check
            DO WHILE (cntr < nconti)
              IF (SUM(matrix(:,src_id)) == 0._wp) THEN ! (incomming) connected tile is end of flow net branch
                matrix(src_id,trgt_id)=0._wp ! Setting the current source -> target to 0, to not get on the branch again 
                IF (sum(matrix(:,is)) > 0._wp) THEN ! Check if connected sibling still has incomming connections
                  trgt_id=is ! Set target back to siblig
                  src_id=MAXLOC(matrix(:,is), 1) ! and find new source 
                  cntr=0 ! Goind down new branch so counter is reset 
                ELSE
                  matrix(is,it) = 0._wp ! set sibling -> tile conection to 0 (sibling branch has no closed loop)
                  cntr=nconti-1 ! nothing more to be checked 
                END IF
              ELSE
                trgt_id=src_id
                src_id=MAXLOC(matrix(:,src_id), 1)
              END IF
              cntr=cntr+1
              IF ((cntr == nconti) .AND. (SUM(matrix(:,src_id)) /= 0._wp)) THEN
                has_loop = .TRUE.
              END IF
            END DO
          ELSE
            matrix(is,it) = 0._wp
          END IF 
        END IF 
      END DO
    END DO


  END SUBROUTINE check_flow_matrix


#endif
END MODULE mo_lex_process
