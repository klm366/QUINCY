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

  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_control, ONLY: acc_stream

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport, calc_drainage_transport,                &
          & calc_pond_soil_flux_water, calc_pond_soil_flux_heat,           &
          & calc_diffusion_water, calc_diffusion_heat,                     &
          & calc_convection_heat, calc_snow_redistribution,                &
          & remove_excess_water, shift_fluxes_sl,                          &
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
  !! Compute ingress of water from pond to adjacent soil
  !!
  SUBROUTINE calc_pond_soil_flux_water(           &
      &      nc, nsoil, dtime, connect,           &
      &      h_offset, length_scale,              &
      &      length_iface_ct, fract_pond_max_dpr, &
      &      lvt_fract, dpr_fract, dsoil_lvt,     &
      &      field_cap_lvt,  k_lvt, mpot_lvt,     &
      &      egress_dpr, ingress_lvt,             &
      &      wtr_latflow_dpr, wtr_latflow_lvt,    &
      &      wtr_pond_dpr, wtr_soil_lvt           &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER,  INTENT(in)    :: nc,                   &  !< Number of cells
      &                        nsoil                    !< Number of soil ayers

    REAL(wp), INTENT(in)    :: dtime,                &  !< Timestep length [s]
      &                        connect(:),           &  !< Connectivity array for lateral exchange
      &                        h_offset(:),          &  !< Differece in surface elevation between elevation and depression [m]
      &                        length_scale(:),      &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        length_iface_ct(:),   &
      &                        fract_pond_max_dpr(:),&  !< Maximum fraction of ponded water on depression tile [/]
      &                        lvt_fract(:),         &  !< Area fraction of elevation tile [/]
      &                        dpr_fract(:)             !< Area fraction of depression tile [/]

    REAL(wp), INTENT(in)    :: dsoil_lvt(:,:),       &  !< Soil depth on layers of elevation tile [m]
      &                        field_cap_lvt(:,:),   &  !< Soil field capacity of elevation tile [m]
      &                        k_lvt(:,:),           &  !< Soil layer hyd. conductivity on source tile [m s-1]
      &                        mpot_lvt(:,:)            !< Matric potential on source tile [m]


    REAL(wp), INTENT(inout) :: egress_dpr(:),        &  !< Diffusion flux on depression tile [kg m-2 s-1]
      &                        wtr_pond_dpr(:),      &  !< Pond storage on depression tile [m]
      &                        wtr_latflow_dpr(:)       !< Netflux (surface) on depression tile [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: ingress_lvt(:,:),     &  !< Diffusion flux on elevation tile [kg m-2 s-1]
      &                        wtr_soil_lvt(:,:),    &  !< Soil water storage on elevation tile [m]
      &                        wtr_latflow_lvt(:,:)     !< Netflux (below ground) on elevation tile [kg m-2 s-1]

    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,             &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, ic, is, i_lvt      !< Index variables for layers on helper grid (i), grid cell, soil level, and in elevation (i_lvt)

    LOGICAL                 :: l_compute(nc)         !< Flag to indicate whether fluxes should be computed for a given cell

    REAL(wp)                :: wtr_pond_dpr_ini,  &  !< Initial surface water storage of depression tile [m]
      &                        free_space,        &  !< Free pore space on elevation tile on a given layer [m]
      &                        wtr_avail,         &  !< Available water on depression tile on a given layer [m]
      &                        wtr_flx_lat,       &  !< Potential flux (based on k and mpot) between depression and elevation [m/s]
      &                        wtr_lat               !< Actual flux between depression and elevation [m]

    REAL(wp)                :: len_inter(nc),     &  !< Length of interface between elevation and depression [m m-2]
      &                        fr_lvt_rel(nc),    &  !< Relative fraction of elevation tile [/]
      &                        fr_dpr_rel(nc),    &  !< Relative fraction of depression tile [/]
      &                        z_bot_lvt(nc),     &  !< Bedrock boundary on elevation tile [m]
      &                        diffusion_sum(nc)     !< Sum of diffusive fluxes on all soil layers [m]

    REAL(wp)                :: diffusion(nc, nsoil)  !< Diffusive flux between depression and elevation on soil layers [m]

    REAL(wp)                :: z_lvt,             &  !< Layer bottom of current layer on elevation tile [m]
      &                        h,                 &  !< Height of current layer [m]
      &                        a_inter,           &  !< Interface area on current layer [m2 m-2]
      &                        z_rgrd_top,        &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,        &  !< Bottom of current layer on helper grid [m]
      &                        wtr_soil_lvt_rgrd, &  !< Soil water content on elevation tile on helper grid layer [m]
      &                        field_cap_lvt_rgrd,&  !< Field capacity on elevation tile on helper grid layer [m]
      &                        mpot_lvt_rgrd,     &  !< Matric potential  on elevation tile on helper grid layer [m]
      &                        k_lvt_rgrd            !< Hydraulic conductivity on elevation tile on helper grid layer [[m s-1]

    ! Initialize 1D fields
    DO ic = 1, nc
      IF (lvt_fract(ic) >= 1.0e-10_wp .AND. dpr_fract(ic)    >= 1.0e-10_wp .AND. &
        & h_offset(ic)  >= 1.0e-10_wp .AND. wtr_pond_dpr(ic) >= 1.0e-10_wp .AND. &
        & connect(ic) > 0._wp) THEN
        ! Conditions for lateral transport are met
        l_compute(ic)  = .TRUE.
        ! Compute length of interface between soil on elevation and ponds on depression tile
        len_inter(ic)  = length_iface_ct(ic) * fract_pond_max_dpr(ic)
        ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
        fr_lvt_rel(ic) = lvt_fract(ic) / (lvt_fract(ic) + dpr_fract(ic))
        fr_dpr_rel(ic) = dpr_fract(ic) / (lvt_fract(ic) + dpr_fract(ic))
      ELSE
        ! No lateral transport possible
        l_compute(ic) = .FALSE.
      END IF
      diffusion_sum(ic) = 0._wp
      z_bot_lvt(ic)     = 0._wp
    END DO

    ! Initilize 2D fields
    DO is = 1, nsoil
      DO  ic = 1, nc
        diffusion(ic, is) = 0._wp
        z_bot_lvt(ic)     = z_bot_lvt(ic) + dsoil_lvt(ic,is)
      END DO
    END DO

    nrgrd = 2 * nsoil
    DO i = 1, nrgrd
      DO ic = 1, nc
        IF (l_compute(ic)) THEN
          ! Initialize first
          IF (i == 1) THEN
            wtr_pond_dpr_ini = wtr_pond_dpr(ic)
            i_lvt      = 1
            z_lvt      = dsoil_lvt(ic,1)
            z_rgrd_top = 0._wp
            IF (h_offset(ic) > wtr_pond_dpr_ini) THEN
              z_rgrd_bot = MIN(z_lvt, h_offset(ic) - wtr_pond_dpr_ini)
            ELSE
              z_rgrd_bot = MIN(z_lvt, h_offset(ic))
            END IF
          END IF

          IF (z_rgrd_bot <= z_bot_lvt(ic)) THEN ! Within soil of elevation
            h                     = z_rgrd_bot - z_rgrd_top
            a_inter               = len_inter(ic) * h
            wtr_soil_lvt_rgrd     = wtr_soil_lvt(ic,i_lvt)  * h / dsoil_lvt(ic,i_lvt)
            wtr_lat               = 0._wp
            IF (z_rgrd_bot <= h_offset(ic)) THEN ! Above depression surface
              field_cap_lvt_rgrd  = field_cap_lvt(ic,i_lvt) * h / dsoil_lvt(ic,i_lvt)
              k_lvt_rgrd          = k_lvt(ic,i_lvt)
              mpot_lvt_rgrd       = mpot_lvt(ic,i_lvt)
              IF (z_rgrd_bot > h_offset(ic) - wtr_pond_dpr_ini) THEN ! Within the ponded surface water in the depression
                ! Water flows from the pond reservoir of the depression tile into the soil of the elevation tile
                wtr_flx_lat       = abs(k_lvt_rgrd  * mpot_lvt_rgrd) / length_scale(ic) *  a_inter
                free_space        = MAX(0._wp,(field_cap_lvt_rgrd - wtr_soil_lvt_rgrd) * fr_lvt_rel(ic))
                wtr_avail         = MIN(h, wtr_pond_dpr(ic)) * fr_dpr_rel(ic)
                wtr_lat           = MIN(wtr_flx_lat * dtime, free_space, wtr_avail)
              END IF
            END IF
            diffusion(ic,i_lvt)   = diffusion(ic,i_lvt) + wtr_lat
          END IF
          ! Layer bottom was on soil layer
          IF (z_rgrd_bot == z_lvt) THEN
            i_lvt = i_lvt + 1
            IF (i_lvt <= nsoil) THEN
              IF (dsoil_lvt(ic,i_lvt) > 0._wp) THEN
                z_lvt = z_lvt + dsoil_lvt(ic,i_lvt)
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
          z_rgrd_top = z_rgrd_bot
          IF (z_rgrd_bot < h_offset(ic) - wtr_pond_dpr_ini) THEN ! Still above water level in depression
            z_rgrd_bot = MIN(h_offset(ic) - wtr_pond_dpr_ini, z_lvt)
          ELSE IF (z_rgrd_bot < h_offset(ic) ) THEN              ! Within the ponded surface water in the depression
            z_rgrd_bot = MIN(h_offset(ic), z_lvt)
          ELSE
            z_rgrd_bot = z_lvt
          END IF
        END IF
      END DO
    END DO

    DO is = 1, nsoil
      DO  ic = 1, nc
        diffusion_sum(ic) = diffusion_sum(ic) + diffusion(ic,is)
      END DO
    END DO

    DO is = 1, nsoil
      DO ic = 1, nc
        IF (l_compute(ic)) THEN
          wtr_soil_lvt(ic,is)    = wtr_soil_lvt(ic,is)    + diffusion(ic,is) / fr_lvt_rel(ic)
          wtr_latflow_lvt(ic,is) = wtr_latflow_lvt(ic,is) + diffusion(ic,is) / fr_lvt_rel(ic) * rhoh2o / dtime
          ingress_lvt(ic,is)     = ingress_lvt(ic,is)     + diffusion(ic,is) / fr_lvt_rel(ic) * rhoh2o / dtime
        END IF
      END DO
    END DO

    DO ic = 1, nc
      IF (l_compute(ic)) THEN
        wtr_pond_dpr(ic)    = wtr_pond_dpr(ic)    - diffusion_sum(ic) / fr_dpr_rel(ic)
        wtr_latflow_dpr(ic) = wtr_latflow_dpr(ic) - diffusion_sum(ic) / fr_dpr_rel(ic) * rhoh2o / dtime
        egress_dpr(ic)      = egress_dpr(ic)      + diffusion_sum(ic) / fr_dpr_rel(ic) * rhoh2o / dtime
      END IF
    END DO

  END SUBROUTINE calc_pond_soil_flux_water


! ===============================================================================================================================
  !>
  !! Compute diffusive heat exchange between pond and adjacent soil
  !!
  SUBROUTINE calc_pond_soil_flux_heat(                 &
      &      nsoil, dtime,                             &
      &      dz, h_offset,                             &
      &      length_scale, len_inter,                  &
      &      lvt_fract, dpr_fract,                     &
      &      t_dpr_top, heat_cap_top_dpr,              &
      &      wtr_pond_dpr, ice_pond_dpr,               &
      &      t_soil_lvt,                               &
      &      heat_cond_soil_lvt, heat_cap_soil_lvt,    &
      &      diffus_in_dpr_top,  diffus_out_dpr_top,   &
      &      heat_flx_lat_dpr_top,                     &
      &      diffus_in_lvt, diffus_out_lvt,            &
      &      heat_flx_lat_lvt                          &
      &      )

    !$ACC ROUTINE SEQ

    USE mo_sse_constants, ONLY: hcond_water, hcond_ice

    INTEGER, INTENT(in)     :: nsoil                      !< Number of oil layers

    REAL(wp), INTENT(in)    :: dtime,                  &  !< Timestep length [s]
      &                        h_offset,               &  !< Differece in surface elevation between elevation and depression [m]
      &                        length_scale,           &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,              &  !< Length of interface between elevation and depression [m m-2]
      &                        lvt_fract,              &  !< Area fraction of elevation [/]
      &                        dpr_fract,              &  !< Area fraction of depression [/]
      &                        t_dpr_top,              &  !< Top layer soil temperture of depression [K]
      &                        heat_cap_top_dpr,       &  !< Top layer soil heat capacity of sibing [J m-3 K-1]
      &                        wtr_pond_dpr,           &  !< Water content of surface reservoir on depression [m]
      &                        ice_pond_dpr               !< Ice content of surface reservoir on depression [m (water equivalent)]

    REAL(wp), INTENT(in)    :: dz(:),                  &  !< Layer thickness [m]
      &                        t_soil_lvt(:),          &  !< Soil tempertures of elevation [K]
      &                        heat_cond_soil_lvt(:),  &  !< Soil heat conductivity of elevation [J m-1 s-1 K-1]
      &                        heat_cap_soil_lvt(:)       !< Soil heat capacity of elevation [J m-3 K-1]

    REAL(wp), INTENT(inout) :: diffus_in_dpr_top,      &  !< Diffusion influx at depression [J m-2 s-1]
      &                        diffus_out_dpr_top,     &  !< Diffusion outflux at depression [J m-2 s-1]
      &                        heat_flx_lat_dpr_top       !< Net Heat flux at depression [J m-2 s-1]

    REAL(wp), INTENT(inout) :: diffus_in_lvt(:),       &  !< Diffusion influx at elevation [J m-2 s-1]
      &                        diffus_out_lvt(:),      &  !< Diffusion outflux at elevation [J m-2 s-1]
      &                        heat_flx_lat_lvt(:)        !< Net Heat flux at elevation [J m-2 s-1]




    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: i, i_lvt                   !< Index variables for layers on grid (i) on in elevation (i_lvt)

    REAL(wp)                :: fr_lvt_rel,             &  !< Relative fraction of elevation [/]
      &                        fr_dpr_rel,             &  !< Relative fraction of depression [/]
      &                        z_lvt,                  &  !< Layer bottom of current layer on elevation [m]
      &                        z_rgrd_top,             &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,             &  !< Bottom of current layer on helper grid [m]
      &                        h,                      &  !< Height of current layer [m]
      &                        a_inter,                &  !< Interface area on current layer [m2 m-2]
      &                        weq_pond_dpr,           &  !< Water and ice content of surface reservoir on depression [m]
      &                        h_cond_pond_dpr,        &  !< Heat conductivity of surface reservoir on depression [J m-1 s-1 K-1]
      &                        t_lvt_rgrd,             &  !< Soil tempertures of elevation on current layer on helper grid [K]
      &                        h_cond_soil_lvt_rgrd,   &  !< Soil heat conductivity of elevation on current layer on helper grid [J m-1 s-1 K-1]
      &                        h_cap_soil_lvt_rgrd,    &  !< Soil heat capacity of elevation on current layer on helper grid [J m-3 K-1]
      &                        h_cond_inter,           &  !< Soil heat conductivity at interface on current layer on helper grid [J m-1 s-1 K-1]
      &                        d_t,                    &  !< Temperature difference on current layer on helper grid [K]
      &                        heat_flx_lat,           &  !< Energy flux (based on h_cond_inter and d_t) between depression and elevation [J m-2 s-1]
      &                        t_eq,                   &  !< Equilibrium temperature [K]
      &                        heat_flx_eq_lat            !< Equilirium energy flux between depression and elevation [J m-2 s-1]

    weq_pond_dpr    = wtr_pond_dpr + ice_pond_dpr

    IF (h_offset < 1.0e-10_wp .OR. weq_pond_dpr < 1.0e-10_wp) RETURN
    IF (lvt_fract < 1.0e-10_wp .OR. dpr_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    ! Start setting up the interpolated helper grid
    ! As ponded water / ice do not have a heat content we only account for energy exchange in the soil:

    i_lvt     = 1
    z_lvt     = dz(i_lvt)
    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, h_offset)

    h_cond_pond_dpr = (wtr_pond_dpr * hcond_water + ice_pond_dpr * hcond_ice) / weq_pond_dpr

    DO i = 1, nsoil
      IF (i_lvt <= nsoil .AND.                        &
        & z_rgrd_bot >  h_offset - weq_pond_dpr .AND. &
        & z_rgrd_bot <= h_offset) THEN ! Within the ponded surface water in the depression
        h                    = z_rgrd_bot - z_rgrd_top
        a_inter              = len_inter * h

        t_lvt_rgrd           = t_soil_lvt(i_lvt)
        h_cond_soil_lvt_rgrd = heat_cond_soil_lvt(i_lvt)
        h_cap_soil_lvt_rgrd  = heat_cap_soil_lvt(i_lvt)

        ! Flux is determined by the difference in temperature between the two tiles
        h_cond_inter         =(h_cond_soil_lvt_rgrd * h_cond_pond_dpr)**0.5_wp
        d_t                  = t_lvt_rgrd - t_dpr_top
        heat_flx_lat         = h_cond_inter  * d_t / length_scale *  a_inter

        ! Limiting heat flux to avoid instablilities for long time steps
        ! Determine flux required to equilibriate Temps within timestep
        ! Equilibrium temperature
        t_eq = (t_lvt_rgrd * h_cap_soil_lvt_rgrd * fr_lvt_rel  &
          &  +  t_dpr_top  * heat_cap_top_dpr    * fr_dpr_rel) &
          &  / (h_cap_soil_lvt_rgrd * fr_lvt_rel + heat_cap_top_dpr * fr_dpr_rel)
        ! Heat flux to reach equilibrium temp
        heat_flx_eq_lat = MIN( &
          &               ABS((t_lvt_rgrd - t_eq)  * h_cap_soil_lvt_rgrd * fr_lvt_rel), &
          &               ABS((t_dpr_top  - t_eq)  * heat_cap_top_dpr    * fr_dpr_rel)) &
          &               *  h / dtime
        ! Limiting heat flux to equilibrium flux
        IF (heat_flx_lat > 0._wp .AND. heat_flx_lat > ABS(heat_flx_eq_lat)) THEN
          heat_flx_lat = ABS(heat_flx_eq_lat)
        ELSE IF (heat_flx_lat < 0._wp .AND. heat_flx_lat < -ABS(heat_flx_eq_lat)) THEN
          heat_flx_lat = -ABS(heat_flx_eq_lat)
        END IF

        ! Diagnosing incoming-, outgoing- and net fluxes on tiles
        diffus_in_lvt(i_lvt)    = diffus_in_lvt(i_lvt)  + MAX(0._wp, -heat_flx_lat / fr_lvt_rel)
        diffus_in_dpr_top       = diffus_in_dpr_top     + MAX(0._wp,  heat_flx_lat / fr_dpr_rel)
        diffus_out_lvt(i_lvt)   = diffus_out_lvt(i_lvt) - MIN(0._wp, -heat_flx_lat / fr_lvt_rel)
        diffus_out_dpr_top      = diffus_out_dpr_top    - MIN(0._wp,  heat_flx_lat / fr_dpr_rel)
        heat_flx_lat_lvt(i_lvt) = heat_flx_lat_lvt(i_lvt)           - heat_flx_lat / fr_lvt_rel
        heat_flx_lat_dpr_top    = heat_flx_lat_dpr_top              + heat_flx_lat / fr_dpr_rel
      END IF

      ! Determening the top and bottem of the next layer of the regridded
      IF (z_rgrd_bot == h_offset) THEN
        i_lvt = nsoil + 1
      ELSE
        i_lvt = i_lvt + 1
        IF (i_lvt <= nsoil) THEN
          z_lvt = z_lvt + dz(i_lvt)
        END IF
      END IF

      z_rgrd_top=z_rgrd_bot
      z_rgrd_bot=MIN(z_lvt, h_offset)
    END DO

  END SUBROUTINE calc_pond_soil_flux_heat

  ! ===============================================================================================================================
  !>
  !! Compute lateral diffusion between two adjacent soil columns (water)
  !!
  SUBROUTINE calc_diffusion_water(                  &
      &      nsoil, dtime,                          &
      &      h_offset, length_scale, len_inter,     &
      &      tle_fract, sbl_fract,                  &
      &      dsoil_tle, dsoil_sbl,                  &
      &      field_cap_tle, field_cap_sbl,          &
      &      wsat_tle, wsat_sbl,                    &
      &      res_tle, res_sbl,                      &
      &      k_sat_tle, k_sat_sbl,                  &
      &      k_tle, k_sbl,                          &
      &      mpot_tle, mpot_sbl,                    &
      &      diffus_in_tle, diffus_in_sbl,          &
      &      diffus_out_tle, diffus_out_sbl,        &
      &      wtr_latflow_tle, wtr_latflow_sbl,      &
      &      wtr_soil_tle, wtr_soil_sbl             &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER, INTENT(in)     :: nsoil                   !< Number of soil layers

    REAL(wp), INTENT(in)    :: dtime,               &  !< Timestep length [s]
      &                        h_offset,            &  !< Differece in surface elevation between tile and sibling [m]
      &                        length_scale,        &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,           &  !< Length of interface between tile and sibling [m m-2]
      &                        tle_fract,           &  !< Area fraction of tile [/]
      &                        sbl_fract               !< Area fraction of sibling [/]

    REAL(wp), INTENT(in)    :: dsoil_tle(:),        &  !< Soil depth on layers of tile [m]
      &                        dsoil_sbl(:),        &  !< Soil depth on layers of sibling [m]
      &                        field_cap_tle(:),    &  !< Soil field capacity of tile [m]
      &                        field_cap_sbl(:),    &  !< Soil field capacity of sibling [m]
      &                        wsat_tle(:),         &  !< Soil moisture at saturation sibling [m]
      &                        wsat_sbl(:),         &  !< Soil moisture at saturation tile [m]
      &                        res_tle(:),          &  !< Residual soil moisture tile [m]
      &                        res_sbl(:),          &  !< Residual soil moisture sibling [m]
      &                        k_sat_tle(:),        &  !< Saturation hydraulic conductivity tile [m s-1]
      &                        k_sat_sbl(:),        &  !< Saturation hydraulic conductivity sibling [m s-1]
      &                        k_tle(:),            &  !< Soil layer hyd. conductivity on source tile [m s-1]
      &                        k_sbl(:),            &  !< Soil layer hyd. conductivity on target tile [m s-1]
      &                        mpot_tle(:),         &  !< Soil layer matric potential on source tile [m]
      &                        mpot_sbl(:)             !< Soil layer matric potential on target tile [m]

    REAL(wp), INTENT(inout) :: diffus_in_tle(:),    &  !< Diffusion influx at tile [kg m-2 s-1]
      &                        diffus_in_sbl(:),    &  !< Diffusion influx at sibling [kg m-2 s-1]
      &                        diffus_out_tle(:),   &  !< Diffusion influx at tile [kg m-2 s-1]
      &                        diffus_out_sbl(:)       !< Diffusion influx at sibling [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_latflow_tle(:),  &  !< Netflux (below ground) on tile [kg m-2 s-1]
      &                        wtr_latflow_sbl(:)      !< Netflux (below ground) on sibling [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: wtr_soil_tle(:),     &  !< Soil water storage on tile [m]
      &                        wtr_soil_sbl(:)         !< Soil water storage on sibling [m]



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,               &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, i_tle, i_sbl         !< Index variables for layers on helper grid (i), in tile (i_tle) and sibling (i_sbl)

    REAL(wp)                :: fr_tle_rel,          &  !< Relative fraction of tile [/]
      &                        fr_sbl_rel,          &  !< Relative fraction of sibling [/]
      &                        h,                   &  !< Height of current layer [m]
      &                        a_inter,             &  !< Interface area on current layer [m2 m-2]
      &                        k_inter,             &  !< Hydraulic conductivity at interface [m s-1]
      &                        d_mpot,              &  !< Difference in matric potential [m]
      &                        wtr_rel_qlbrm,       &  !< Relative saturation at equilibrium between connected tiles [/]
      &                        free_space,          &  !< Free pore space on tile on a given layer [m]
      &                        wtr_avail,           &  !< Available water on sibling on a given layer [m]
      &                        wtr_2_qlbrm,         &  !< Amount of water moved in order to reach equilibrium state [m]
      &                        wtr_flx_lat,         &  !< potential flux (based on k and mpot) between sibling and tile [m/s]
      &                        wtr_lat,             &  !< actual flux between sibling and tile [m]
      &                        z_rgrd_top,          &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,          &  !< Bottom of current layer on helper grid [m]
      &                        z_tle,               &  !< Layer bottom of current layer on tile [m]
      &                        z_sbl,               &  !< Layer bottom of current layer on sibling [m]
      &                        z_bot_tle,           &  !< Bedrock boundary on tile [m]
      &                        z_bot_sbl,           &  !< Bedrock boundary on sibling [m]
      &                        wtr_soil_tle_rgrd,   &  !< Soil water content on tile on helper grid layer [m]
      &                        wsat_tle_rgrd,       &  !< Soil moisture at saturation on tile on helper grid layer [m]
      &                        field_cap_tle_rgrd,  &  !< Field capacity on tile on helper grid layer [m]
      &                        res_tle_rgrd,        &  !< Residual soil moisture on tile on helper grid layer [m]
      &                        mpot_tle_rgrd,       &  !< Matric potential on tile on helper grid layer [m]
      &                        k_tle_rgrd,          &  !< Hydraulic conductivity on tile on helper grid layer [[m s-1]
      &                        wtr_soil_sbl_rgrd,   &  !< Soil water content on sibling on helper grid layer [m]
      &                        wsat_sbl_rgrd,       &  !< Soil moisture at saturation on sibling on helper grid layer [m]
      &                        field_cap_sbl_rgrd,  &  !< Field capacity on sibling on helper grid layer [m]
      &                        res_sbl_rgrd,        &  !< Residual soil moisture on sibling on helper grid layer [m]
      &                        mpot_sbl_rgrd,       &  !< Matric potential on deperssion tile on helper grid layer [m]
      &                        k_sbl_rgrd              !< Hydraulic conductivity on sibling on helper grid layer [[m s-1]

    REAL(wp)                :: diffus_tle(nsoil),   &  !< Net diffusive flux on tile on soil layers [m]
      &                        diffus_sbl(nsoil)       !< Net diffusive flux on sibling on soil layers [m]

    nrgrd         = nsoil * 2 + 3

    ! Initilize fields for accumulating regridded output
    diffus_tle(:)        = 0._wp
    diffus_sbl(:)        = 0._wp

    IF (h_offset > 0._wp .AND. h_offset > SUM(dsoil_tle(:))) RETURN
    IF (h_offset < 0._wp .AND. abs(h_offset) > SUM(dsoil_sbl(:))) RETURN
    IF (tle_fract < 1.0e-10_wp .OR. sbl_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_tle_rel = tle_fract / (tle_fract+sbl_fract)
    fr_sbl_rel = sbl_fract / (tle_fract+sbl_fract)

    ! Start setting up the interpolated helper grid
    IF (h_offset > 0._wp) THEN
      i_tle     = 1
      i_sbl     = 0
      z_tle     = dsoil_tle(i_tle)
      z_sbl     = abs(h_offset)
      z_bot_tle = SUM(dsoil_tle(:))
      z_bot_sbl = SUM(dsoil_sbl(:)) + abs(h_offset)
    ELSE IF (h_offset < 0._wp) THEN
      i_tle     = 0
      i_sbl     = 1
      z_tle     = abs(h_offset)
      z_sbl     = dsoil_sbl(i_sbl)
      z_bot_tle = SUM(dsoil_tle(:)) + abs(h_offset)
      z_bot_sbl = SUM(dsoil_sbl(:))
    ELSE
      i_tle     = 1
      i_sbl     = 1
      z_tle     = dsoil_tle(i_tle)
      z_sbl     = dsoil_sbl(i_sbl)
      z_bot_tle = SUM(dsoil_tle(:))
      z_bot_sbl = SUM(dsoil_sbl(:))
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_tle, z_sbl)

    DO i = 1, nrgrd
      IF (z_rgrd_bot <= z_bot_tle .AND. z_rgrd_bot <= z_bot_sbl .AND. & ! Above Bedrock in both tiles
          i_tle      <= nsoil     .AND. i_sbl      <= nsoil     ) THEN  ! within soil column
        IF (wsat_tle(i_tle) > 0._wp .AND. wsat_sbl(i_sbl) > 0._wp) THEN ! Space for liquid water exists in both tiles

          IF (z_rgrd_bot > abs(h_offset) ) THEN ! Within the soil column of both tiles
            ! Flux is determined by the difference in matric potential between the two tiles

            wtr_lat            = 0._wp
            h                  = z_rgrd_bot - z_rgrd_top
            a_inter            = len_inter * h

            wtr_soil_tle_rgrd  = wtr_soil_tle(i_tle)   * h / dsoil_tle(i_tle)
            wsat_tle_rgrd      = wsat_tle(i_tle)       * h / dsoil_tle(i_tle)
            field_cap_tle_rgrd = field_cap_tle(i_tle)  * h / dsoil_tle(i_tle)
            res_tle_rgrd       = res_tle(i_tle)        * h / dsoil_tle(i_tle)
            mpot_tle_rgrd      = mpot_tle(i_tle)
            k_tle_rgrd         = k_tle(i_tle)

            wtr_soil_sbl_rgrd  = wtr_soil_sbl(i_sbl)   * h / dsoil_sbl(i_sbl)
            wsat_sbl_rgrd      = wsat_sbl(i_sbl)       * h / dsoil_sbl(i_sbl)
            field_cap_sbl_rgrd = field_cap_sbl(i_sbl)  * h / dsoil_sbl(i_sbl)
            res_sbl_rgrd       = res_sbl(i_sbl)        * h / dsoil_sbl(i_sbl)
            mpot_sbl_rgrd      = mpot_sbl(i_sbl)
            k_sbl_rgrd         = k_sbl(i_sbl)
            k_inter            = (k_sbl_rgrd * k_tle_rgrd)**0.5_wp
            d_mpot             =  abs(mpot_tle_rgrd) - abs(mpot_sbl_rgrd)
            wtr_flx_lat        = - k_inter  * d_mpot / length_scale *  a_inter
            wtr_rel_qlbrm      = MIN(1._wp, (wtr_soil_tle_rgrd * fr_tle_rel  &
              &                +  wtr_soil_sbl_rgrd * fr_sbl_rel)            &
              &                / (wsat_tle_rgrd * fr_tle_rel + wsat_sbl_rgrd * fr_sbl_rel ))
            IF (wtr_flx_lat > 0._wp) THEN ! Flux is going from tile to sibling

              wtr_2_qlbrm      = MAX(0._wp, &
                &               (wtr_rel_qlbrm * wsat_sbl_rgrd - wtr_soil_sbl_rgrd) * fr_sbl_rel)
              free_space       = MAX(0._wp,(field_cap_sbl_rgrd - wtr_soil_sbl_rgrd) * fr_sbl_rel)
              wtr_avail        = MAX(0._wp,(wtr_soil_tle_rgrd - res_tle_rgrd) * fr_tle_rel)
              wtr_lat          = MIN(wtr_flx_lat * dtime, free_space, wtr_avail, wtr_2_qlbrm)
              wtr_lat          = MAX(0._wp, wtr_lat)

            ELSE ! Flux is going from sibling to tile

              wtr_2_qlbrm      = MAX(0._wp, &
                &               (wtr_rel_qlbrm  * wsat_tle_rgrd - wtr_soil_tle_rgrd) * fr_tle_rel)
              free_space       = MAX(0._wp, (field_cap_tle_rgrd - wtr_soil_tle_rgrd) * fr_tle_rel)
              wtr_avail        = MAX(0._wp, (wtr_soil_sbl_rgrd - res_sbl_rgrd) * fr_sbl_rel)
              wtr_lat          = MAX(wtr_flx_lat * dtime, -free_space, -wtr_avail, -wtr_2_qlbrm)
              wtr_lat          = MIN(0._wp, wtr_lat)

            END IF
            diffus_tle(i_tle)  = diffus_tle(i_tle) - wtr_lat / fr_tle_rel
            diffus_sbl(i_sbl)  = diffus_sbl(i_sbl) + wtr_lat / fr_sbl_rel
          END IF

          IF (z_rgrd_bot == z_tle) THEN
            i_tle = i_tle + 1
            IF (i_tle <= nsoil) THEN
              IF (dsoil_tle(i_tle) > 0._wp) THEN
                z_tle = z_tle + dsoil_tle(i_tle)
              ELSE
                i_tle = nsoil + 1
              END IF
            END IF
          END IF

          IF (z_rgrd_bot == z_sbl) THEN
            i_sbl = i_sbl + 1
            IF (i_sbl <= nsoil) THEN
              IF (dsoil_sbl(i_sbl) > 0._wp) THEN
                z_sbl = z_sbl + dsoil_sbl(i_sbl)
              ELSE
                i_sbl = nsoil + 1
              END IF
            END IF
          END IF

          z_rgrd_top=z_rgrd_bot
          z_rgrd_bot=MIN(z_tle, z_sbl)
        END IF
      END IF
    END DO

    ! Adjusting the water contents and calculating the net- in- and outgoing fluxes for the tiles
    wtr_soil_tle(:)   = wtr_soil_tle(:)              + diffus_tle(:)
    wtr_soil_sbl(:)   = wtr_soil_sbl(:)              + diffus_sbl(:)

    wtr_latflow_tle   = wtr_latflow_tle              + diffus_tle(:)  * rhoh2o / dtime
    wtr_latflow_sbl   = wtr_latflow_sbl              + diffus_sbl(:)  * rhoh2o / dtime

    diffus_in_tle(:)  = diffus_in_tle(:)  + MAX(0._wp, diffus_tle(:)) * rhoh2o / dtime
    diffus_in_sbl(:)  = diffus_in_sbl(:)  + MAX(0._wp, diffus_sbl(:)) * rhoh2o / dtime
    diffus_out_tle(:) = diffus_out_tle(:) - MIN(0._wp, diffus_tle(:)) * rhoh2o / dtime
    diffus_out_sbl(:) = diffus_out_sbl(:) - MIN(0._wp, diffus_sbl(:)) * rhoh2o / dtime

  END SUBROUTINE calc_diffusion_water

! ===============================================================================================================================
  !>
  !! Compute lateral diffusion between two adjacent soil columns (energy)
  !!
  SUBROUTINE calc_diffusion_heat(                      &
      &      nsoil, dtime,                             &
      &      dz, h_offset,                             &
      &      length_scale, len_inter,                  &
      &      tle_fract, sbl_fract,                     &
      &      t_soil_tle, t_soil_sbl,                   &
      &      heat_cond_soil_tle, heat_cond_soil_sbl,   &
      &      heat_cap_soil_tle, heat_cap_soil_sbl,     &
      &      diffus_in_tle, diffus_in_sbl,             &
      &      diffus_out_tle, diffus_out_sbl,           &
      &      heat_flx_lat_tle, heat_flx_lat_sbl        &
      &      )

    !$ACC ROUTINE SEQ

    INTEGER, INTENT(in)     :: nsoil                      !< Number of oil layers

    REAL(wp), INTENT(in)    :: dtime,                  &  !< Timestep length [s]
      &                        h_offset,               &  !< Differece in surface elevation between tile and sibling [m]
      &                        length_scale,           &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      &                        len_inter,              &  !< Length of interface between tile and sibling [m m-2]
      &                        tle_fract,              &  !< Area fraction of tile [/]
      &                        sbl_fract                  !< Area fraction of sibling [/]

    REAL(wp), INTENT(in)    :: dz(:),                  &  !< Layer thickness [m]
      &                        t_soil_tle(:),          &  !< Soil tempertures of tile [K]
      &                        t_soil_sbl(:),          &  !< Soil tempertures of sibling [K]
      &                        heat_cond_soil_tle(:),  &  !< Soil heat conductivity of tile [J m-1 s-1 K-1]
      &                        heat_cond_soil_sbl(:),  &  !< Soil heat conductivity of sibling [J m-1 s-1 K-1]
      &                        heat_cap_soil_tle(:),   &  !< Soil heat capacity of tile [J m-3 K-1]
      &                        heat_cap_soil_sbl(:)       !< Soil heat capacity of sibling [J m-3 K-1]

    REAL(wp), INTENT(inout) :: diffus_in_tle(:),       &  !< Diffusion influx at tile [J m-2 s-1]
      &                        diffus_in_sbl(:),       &  !< Diffusion outflux at sibling [J m-2 s-1]
      &                        diffus_out_tle(:),      &  !< Diffusion influx at tile [J m-2 s-1]
      &                        diffus_out_sbl(:)          !< Diffusion outflux at sibling [J m-2 s-1]

    REAL(wp), INTENT(inout) :: heat_flx_lat_tle(:),    &  !< Net Heat flux at tile [J m-2 s-1]
      &                        heat_flx_lat_sbl(:)        !< Net Heat flux at sibling [J m-2 s-1]



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,                  &  !< Number of layers on helper grid determined by layer overlaps
      &                        i, i_tle, i_sbl            !< Index variables for layers on helper grid (i), in tile (i_tle) and sibling (i_sbl)

    REAL(wp)                :: fr_tle_rel,             &  !< Relative fraction of tile [/]
      &                        fr_sbl_rel,             &  !< Relative fraction of sibling [/]
      &                        z_tle,                  &  !< Layer bottom of current layer on tile [m]
      &                        z_sbl,                  &  !< Layer bottom of current layer on sibling [m]
      &                        z_rgrd_top,             &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,             &  !< Bottom of current layer on helper grid [m]
      &                        h,                      &  !< Height of current layer [m]
      &                        a_inter,                &  !< Interface area on current layer [m2 m-2]
      &                        t_tle_rgrd,             &  !< Soil tempertures of tile on current layer on helper grid [K]
      &                        t_sbl_rgrd,             &  !< Soil tempertures of sibling on current layer on helper grid [K]
      &                        h_cond_soil_tle_rgrd,   &  !< Soil heat conductivity of tile on current layer on helper grid [J m-1 s-1 K-1]
      &                        h_cond_soil_sbl_rgrd,   &  !< Soil heat conductivity of sibling on current layer on helper grid [J m-1 s-1 K-1]
      &                        h_cap_soil_tle_rgrd,    &  !< Soil heat capacity of tile on current layer on helper grid [J m-3 K-1]
      &                        h_cap_soil_sbl_rgrd,    &  !< Soil heat capacity of sibling on current layer on helper grid [J m-3 K-1]
      &                        h_cond_inter,           &  !< Soil heat conductivity at interface on current layer on helper grid [J m-1 s-1 K-1]
      &                        d_t,                    &  !< Temperature difference on current layer on helper grid [K]
      &                        heat_flx_lat,           &  !< Energy flux (based on h_cond_inter and d_t) between sibling and tile [J m-2 s-1]
      &                        t_eq,                   &  !< Equilibrium temperature [K]
      &                        heat_flx_eq_lat            !< Equilirium energy flux between sibling and tile [J m-2 s-1]

    nrgrd         = nsoil * 2 + 3

    IF (abs(h_offset) > SUM(dz(:))) RETURN
    IF (tle_fract < 1.0e-10_wp .OR. sbl_fract < 1.0e-10_wp) RETURN

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_tle_rel = tle_fract / (tle_fract+sbl_fract)
    fr_sbl_rel = sbl_fract / (tle_fract+sbl_fract)

    ! Start setting up the interpolated helper grid
    ! As ponded water / ice do not have a heat content we only account for energy exchange in the soil:

    IF (h_offset > 0._wp) THEN
      i_tle    = 1
      i_sbl    = 0
      z_tle    = dz(i_tle)
      z_sbl    = abs(h_offset)
    ELSE IF (h_offset < 0._wp) THEN
      i_tle    = 0
      i_sbl    = 1
      z_tle    = abs(h_offset)
      z_sbl    = dz(i_sbl)
    ELSE
      i_tle    = 1
      i_sbl    = 1
      z_tle    = dz(i_tle)
      z_sbl    = dz(i_sbl)
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_tle, z_sbl)

    DO i = 1, nrgrd
      IF (i_tle <= nsoil .AND. i_sbl <= nsoil) THEN ! Within the soil column of both tiles

        IF (z_rgrd_bot > abs(h_offset)) THEN ! Within the soil column of both tiles
          h                    = z_rgrd_bot - z_rgrd_top
          a_inter              = len_inter * h

          t_tle_rgrd           = t_soil_tle(i_tle)
          h_cond_soil_tle_rgrd = heat_cond_soil_tle(i_tle)
          h_cap_soil_tle_rgrd  = heat_cap_soil_tle(i_tle)

          t_sbl_rgrd           = t_soil_sbl(i_sbl)
          h_cond_soil_sbl_rgrd = heat_cond_soil_sbl(i_sbl)
          h_cap_soil_sbl_rgrd  = heat_cap_soil_sbl(i_sbl)

          ! Flux is determined by the difference in temperature between the two tiles
          h_cond_inter         =(h_cond_soil_tle_rgrd * h_cond_soil_sbl_rgrd)**0.5_wp
          d_t                  = t_tle_rgrd - t_sbl_rgrd
          heat_flx_lat         = h_cond_inter  * d_t / length_scale *  a_inter

          ! Limiting heat flux to avoid instablilities for long time steps
          ! Determine flux required to equilibriate Temps within timestep
          ! Equilibrium temperature
          t_eq = (t_tle_rgrd * h_cap_soil_tle_rgrd * fr_tle_rel  &
            &  +  t_sbl_rgrd * h_cap_soil_sbl_rgrd * fr_sbl_rel) &
            &  / (h_cap_soil_tle_rgrd * fr_tle_rel + h_cap_soil_sbl_rgrd * fr_sbl_rel)
          ! Heat flux to reach equilibrium temp
          heat_flx_eq_lat = MIN( &
            &               ABS((t_tle_rgrd - t_eq)  * h_cap_soil_tle_rgrd * fr_tle_rel), &
            &               ABS((t_sbl_rgrd - t_eq)  * h_cap_soil_sbl_rgrd * fr_sbl_rel)) &
            &               *  h / dtime
          ! Limiting heat flux to equilibrium flux
          IF (heat_flx_lat > 0._wp .AND. heat_flx_lat > ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = ABS(heat_flx_eq_lat)
          ELSE IF (heat_flx_lat < 0._wp .AND. heat_flx_lat < -ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = -ABS(heat_flx_eq_lat)
          END IF

          ! Diagnosing incoming-, outgoing- and net fluxes on tiles
          diffus_in_tle(i_tle)    = diffus_in_tle(i_tle)  + MAX(0._wp, -heat_flx_lat / fr_tle_rel)
          diffus_in_sbl(i_sbl)    = diffus_in_sbl(i_sbl)  + MAX(0._wp,  heat_flx_lat / fr_sbl_rel)
          diffus_out_tle(i_tle)   = diffus_out_tle(i_tle) - MIN(0._wp, -heat_flx_lat / fr_tle_rel)
          diffus_out_sbl(i_sbl)   = diffus_out_sbl(i_sbl) - MIN(0._wp,  heat_flx_lat / fr_sbl_rel)

          heat_flx_lat_tle(i_tle) = heat_flx_lat_tle(i_tle)           - heat_flx_lat / fr_tle_rel
          heat_flx_lat_sbl(i_sbl) = heat_flx_lat_sbl(i_sbl)           + heat_flx_lat / fr_sbl_rel
        END IF

        ! Determening the top and bottem of the next layer of the regridded
        IF (z_rgrd_bot == z_tle) THEN
          i_tle = i_tle + 1
          IF (i_tle <= nsoil) THEN
            z_tle = z_tle + dz(i_tle)
          END IF
        END IF

        IF (z_rgrd_bot == z_sbl) THEN
          i_sbl = i_sbl + 1
          IF (i_sbl <= nsoil) THEN
            z_sbl = z_sbl + dz(i_sbl)
          END IF
        END IF

        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=MIN(z_tle, z_sbl)
      END IF
    END DO

  END SUBROUTINE calc_diffusion_heat

! ===============================================================================================================================
  !>
  !! Compute intertile heat exchenge due to lateral movement of water
  !!
  SUBROUTINE calc_convection_heat(                         &
      &      nc, nsoil, dtime,                             &
      &      dz, h_offset,                                 &
      &      tle_fract, sbl_fract,                         &
      &      wtr_inflx_lat_srf_tle, wtr_inflx_lat_srf_sbl, &
      &      t_soil_tle, t_soil_sbl,                       &
      &      wtr_inflx_lat_sl_tle, wtr_inflx_lat_sl_sbl,   &
      &      heat_cap_soil_tle, heat_cap_soil_sbl,         &
      &      cnvctv_in_tle, cnvctv_in_sbl,                 &
      &      cnvctv_out_tle, cnvctv_out_sbl,               &
      &      heat_flx_lat_tle, heat_flx_lat_sbl            &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_sse_constants, ONLY: vol_hcap_water

    INTEGER, INTENT(in)     :: nc,                        &  !< Number of grid cells
      &                        nsoil                         !< Number of oil layers

    REAL(wp), INTENT(in)    :: dtime                         !< Timestep length [s]

    REAL(wp), INTENT(in)    :: h_offset(:),               &  !< Difference in surface elevation between tile and sibling [m]
      &                        tle_fract(:),              &  !< Area fraction of tile [/]
      &                        sbl_fract(:),              &  !< Area fraction of sibling [/]
      &                        wtr_inflx_lat_srf_tle(:),  &  !< Lateral flux (surface) of water from sibling to tile [kg m-2 s-1]
      &                        wtr_inflx_lat_srf_sbl(:),  &  !< Lateral flux (surface) of water from tile to sibling [kg m-2 s-1]
      &                        dz(:)                         !< Layer thickness [m]

    REAL(wp), INTENT(in)    :: t_soil_tle(:,:),           &  !< Soil tempertures of tile [K]
      &                        t_soil_sbl(:,:),           &  !< Soil tempertures of sibling [K]
      &                        wtr_inflx_lat_sl_tle(:,:), &  !< Lateral flux (soil) of water from sibling to tile [kg m-2 s-1]
      &                        wtr_inflx_lat_sl_sbl(:,:), &  !< Lateral flux (soil) of water from tile to sibling [kg m-2 s-1]
      &                        heat_cap_soil_tle(:,:),    &  !< Soil heat capacity of tile [J m-3 K-1]
      &                        heat_cap_soil_sbl(:,:)        !< Soil heat capacity of sibling [J m-3 K-1]

    REAL(wp), INTENT(inout) :: cnvctv_in_tle(:,:),        &  !< Convective influx at tile [J m-2 s-1]
      &                        cnvctv_in_sbl(:,:),        &  !< Convective outflux at sibling [J m-2 s-1]
      &                        cnvctv_out_tle(:,:),       &  !< Convective influx at tile [J m-2 s-1]
      &                        cnvctv_out_sbl(:,:)           !< Convective outflux at sibling [J m-2 s-1]

    REAL(wp), INTENT(inout) :: heat_flx_lat_tle(:,:),     &  !< Net Heat flux at tile [J m-2 s-1]
      &                        heat_flx_lat_sbl(:,:)         !< Net Heat flux at sibling [J m-2 s-1]



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd,                  &  !< Number of layers on helper grid determined by layer overlaps
      &                        ic,                     &  !< Index variable for grid cells
      &                        is,                     &  !< Index variable for soil layers
      &                        i, i_tle, i_sbl            !< Index variables for layers on helper grid (i), on tile (i_tle) and sibling tile (i_sbl)

    LOGICAL                 :: l_compute(nc)              !< Flag to indicate whether fluxes should be computed for a given cell

    REAL(wp)                :: dz_sum,                 &  !< Total soil column depth [m]
      &                        z_tle,                  &  !< Layer bottom of current layer on tile [m]
      &                        z_sbl,                  &  !< Layer bottom of current layer on sibling [m]
      &                        z_rgrd_top,             &  !< Top of current layer on helper grid [m]
      &                        z_rgrd_bot,             &  !< Bottom of current layer on helper grid [m]
      &                        h,                      &  !< Height of current layer [m]
      &                        t_tle_rgrd,             &  !< Soil tempertures of tile on current layer on helper grid [K]
      &                        t_sbl_rgrd,             &  !< Soil tempertures of sibling on current layer on helper grid [K]
      &                        d_t_wtr_soil_tle,       &  !< Temperture difference between incomming water and soil on tile on current layer on helper grid [K]
      &                        d_t_wtr_soil_sbl,       &  !< Temperture difference between incomming water and soil on sibling on current layer on helper grid [K]
      &                        h_cap_soil_tle_rgrd,    &  !< Soil heat capacity of tile on current layer on helper grid [J m-3 K-1]
      &                        h_cap_soil_sbl_rgrd,    &  !< Soil heat capacity of sibling on current layer on helper grid [J m-3 K-1]
      &                        wtr_flx_in_tle_regrid,  &  !< Incomming water on tile (from sibling) on current layer on helper grid [m3 m-2 s]
      &                        wtr_flx_in_sbl_regrid,  &  !< Incomming water on sibling (from tile) on current layer on helper grid [m3 m-2 s]
      &                        heat_flx_tle,           &  !< Incomming heat flux on tile on current layer on helper grid [J m-3]
      &                        heat_flx_sbl               !< Incomming heat flux on sibling on current layer on helper grid [J m-3]


    ! Initialize 1D fields
    dz_sum = 0._wp
    DO is = 1, nsoil
      dz_sum = dz_sum + dz(is)
    END DO

    ! Idenify grid cells where fluxes should be computed
    DO ic = 1, nc
      IF (abs(h_offset(ic)) > dz_sum .OR. &
        & tle_fract(ic) < 1.0e-10_wp .OR. sbl_fract(ic) < 1.0e-10_wp) THEN
        l_compute(ic) = .FALSE.
      ELSE
        l_compute(ic) = .TRUE.
      END IF
    END DO

    ! Determine number of layers on helper grid
    nrgrd = nsoil * 2 + 3

    ! Start setting up the interpolated helper grid
    ! As ponded water / ice do not have a heat content we only account for energy exchange in the soil:
    DO i = 1, nrgrd
      DO ic = 1, nc
        IF (l_compute(ic)) THEN
          ! Initialize first
          IF (i == 1) THEN
            IF (h_offset(ic) < 0._wp) THEN
              i_tle = 0
              i_sbl = 1
              z_tle = abs(h_offset(ic))
              z_sbl = dz(i_sbl)
            ELSE IF  (h_offset(ic) > 0._wp) THEN
              i_tle = 1
              i_sbl = 0
              z_tle = dz(i_tle)
              z_sbl = abs(h_offset(ic))
            ELSE
              i_tle = 1
              i_sbl = 1
              z_tle = dz(i_tle)
              z_sbl = dz(i_sbl)
            END IF

            z_rgrd_top = 0._wp
            z_rgrd_bot = MIN(z_tle, z_sbl)
          END IF

          ! First determine the lateral heat fluxes for water movements within the soil
          IF (i_tle <= nsoil .AND. i_sbl <= nsoil) THEN ! Within the soil column of both tiles

            h = z_rgrd_bot - z_rgrd_top

            IF (i_tle > 0) THEN
              t_tle_rgrd            = t_soil_tle(ic,i_tle)
              h_cap_soil_tle_rgrd   = heat_cap_soil_tle(ic,i_tle)
              wtr_flx_in_tle_regrid = wtr_inflx_lat_sl_tle(ic,i_tle) * h / dz(i_tle) / rhoh2o
            ELSE
              t_tle_rgrd            = t_soil_tle(ic,1)
              h_cap_soil_tle_rgrd   = heat_cap_soil_tle(ic,1)
              wtr_flx_in_tle_regrid = 0._wp
            END IF

            IF (i_sbl > 0) THEN
              t_sbl_rgrd            = t_soil_sbl(ic,i_sbl)
              h_cap_soil_sbl_rgrd   = heat_cap_soil_sbl(ic,i_sbl)
              wtr_flx_in_sbl_regrid = wtr_inflx_lat_sl_sbl(ic,i_sbl) * h / dz(i_sbl) / rhoh2o
            ELSE
              t_sbl_rgrd            = t_soil_sbl(ic,1)
              h_cap_soil_sbl_rgrd   = heat_cap_soil_sbl(ic,1)
              wtr_flx_in_sbl_regrid = 0._wp
            END IF

            ! Effective heat flux is determined so that:
            ! t_sl_new = t_wtr + (t_sl-t_wtr) * (h - h_wtr) / h * vol_heat_cap_sl_old / vol_heat_cap_sl_new.
            ! Since there is no way of knowing what vol_heat_cap_sl_new will be, the scheme assumes that
            ! the water content will change according to wtr_flx * dtime:
            ! vol_heat_cap_new = (vol_heat_cap_sl_old * (h-h_wtr) + vol_heat_cap_h2o * h_wtr) / h
            ! Incomming water is assumed to have the same temperature as the soil layer it originates from.

            d_t_wtr_soil_tle = t_sbl_rgrd - t_tle_rgrd
            d_t_wtr_soil_sbl = t_tle_rgrd - t_sbl_rgrd

            heat_flx_tle     = d_t_wtr_soil_tle * vol_hcap_water * wtr_flx_in_tle_regrid
            heat_flx_sbl     = d_t_wtr_soil_sbl * vol_hcap_water * wtr_flx_in_sbl_regrid

            ! Soil may not warm above or cool below the temperature of the incomming water
            ! This could happen in the above formulation if wtr_flx_in_tle_regrid * dtime > h
            ! Here, we use the current soil heat capacities instead of assuming that
            ! the water content will change according to wtr_flx * dtime.
            ! This makes the limits more restrictive, given that these are anyway very large margins.
            IF (     heat_flx_tle < 0._wp) THEN
              heat_flx_tle   = MAX(heat_flx_tle, d_t_wtr_soil_tle * h_cap_soil_tle_rgrd * h / dtime)
            ELSE IF (heat_flx_tle > 0._wp) THEN
              heat_flx_tle   = MIN(heat_flx_tle, d_t_wtr_soil_tle * h_cap_soil_tle_rgrd * h / dtime)
            END IF
            IF (     heat_flx_sbl < 0._wp) THEN
              heat_flx_sbl   = MAX(heat_flx_sbl, d_t_wtr_soil_sbl * h_cap_soil_sbl_rgrd * h / dtime)
            ELSE IF (heat_flx_sbl > 0._wp) THEN
              heat_flx_sbl   = MIN(heat_flx_sbl, d_t_wtr_soil_sbl * h_cap_soil_sbl_rgrd * h / dtime)
            END IF

            ! Diagnosing incoming-, outgoing- and net fluxes on tiles
            ! Note that the the outgoing fluxes are not included in the net fluxes
            ! since these do not change the temperature (only the heat capacity)
            IF (i_tle > 0) THEN
              cnvctv_in_tle(ic,i_tle)    = cnvctv_in_tle(ic,i_tle)    + heat_flx_tle
              cnvctv_out_tle(ic,i_tle)   = cnvctv_out_tle(ic,i_tle)   + heat_flx_sbl * (sbl_fract(ic) / tle_fract(ic))
              heat_flx_lat_tle(ic,i_tle) = heat_flx_lat_tle(ic,i_tle) + heat_flx_tle
            END IF

            IF (i_sbl > 0) THEN
              cnvctv_in_sbl(ic,i_sbl)    = cnvctv_in_sbl(ic,i_sbl)    + heat_flx_sbl
              cnvctv_out_sbl(ic,i_sbl)   = cnvctv_out_sbl(ic,i_sbl)   + heat_flx_tle * (tle_fract(ic) / sbl_fract(ic))
              heat_flx_lat_sbl(ic,i_sbl) = heat_flx_lat_sbl(ic,i_sbl) + heat_flx_sbl
            END IF

            ! Determening the top and bottem of the next layer of the regridded
            IF (z_rgrd_bot == z_tle) THEN
              i_tle = i_tle + 1
              IF (i_tle <= nsoil) THEN
                z_tle = z_tle + dz(i_tle)
              END IF
            END IF

            IF (z_rgrd_bot == z_sbl) THEN
              i_sbl = i_sbl + 1
              IF (i_sbl <= nsoil) THEN
                z_sbl = z_sbl + dz(i_sbl)
              END IF
            END IF
            z_rgrd_top = z_rgrd_bot
            z_rgrd_bot = MIN(z_tle, z_sbl)
          END IF
        END IF
      END DO
    END DO

    ! Add the heat fluxes due to surface runoff / runon
    DO ic = 1, nc
      IF (l_compute(ic)) THEN
        d_t_wtr_soil_tle = t_soil_sbl(ic,1) - t_soil_tle(ic,1)
        d_t_wtr_soil_sbl = t_soil_tle(ic,1) - t_soil_sbl(ic,1)
        heat_flx_tle = d_t_wtr_soil_tle * vol_hcap_water * wtr_inflx_lat_srf_tle(ic) / rhoh2o
        heat_flx_sbl = d_t_wtr_soil_sbl * vol_hcap_water * wtr_inflx_lat_srf_sbl(ic) / rhoh2o

        IF (     heat_flx_tle < 0._wp) THEN
          heat_flx_tle   = MAX(heat_flx_tle, d_t_wtr_soil_tle * heat_cap_soil_tle(ic,1) * dz(1) / dtime)
        ELSE IF (heat_flx_tle > 0._wp) THEN
          heat_flx_tle   = MIN(heat_flx_tle, d_t_wtr_soil_tle * heat_cap_soil_tle(ic,1) * dz(1) / dtime)
        END IF

        IF (     heat_flx_sbl < 0._wp) THEN
          heat_flx_sbl   = MAX(heat_flx_sbl, d_t_wtr_soil_sbl * heat_cap_soil_sbl(ic,1) * dz(1) / dtime)
        ELSE IF (heat_flx_sbl > 0._wp) THEN
          heat_flx_sbl   = MIN(heat_flx_sbl, d_t_wtr_soil_sbl * heat_cap_soil_sbl(ic,1) * dz(1) / dtime)
        END IF

        cnvctv_in_tle(ic,1)    = cnvctv_in_tle(ic,1)  + heat_flx_tle
        cnvctv_in_sbl(ic,1)    = cnvctv_in_sbl(ic,1)  + heat_flx_sbl
        cnvctv_out_tle(ic,1)   = cnvctv_out_tle(ic,1) + heat_flx_sbl * (sbl_fract(ic) / tle_fract(ic))
        cnvctv_out_sbl(ic,1)   = cnvctv_out_sbl(ic,1) + heat_flx_tle * (tle_fract(ic) / sbl_fract(ic))

        heat_flx_lat_tle(ic,1) = heat_flx_lat_tle(ic,1) + heat_flx_tle
        heat_flx_lat_sbl(ic,1) = heat_flx_lat_sbl(ic,1) + heat_flx_sbl
      END IF
    END DO

  END SUBROUTINE calc_convection_heat


  ! ===============================================================================================================================
  !>
  !! Compute intertile snow redistribution
  !!
  SUBROUTINE calc_snow_redistribution(              &
      &      l_dynsnow,                             &
      &      dtime,                                 &
      &      cnct_rel,                              &
      &      h_offset,                              &
      &      tle_fract, sbl_fract,                  &
      &      weq_pond_tle, weq_pond_sbl,            &
      &      snow_soil_tle_tm1, snow_new_tle_tm1,   &
      &      snow_dens_tle_tm1,                     &
      &      snow_rdst_in_tle, snow_rdst_in_sbl,    &
      &      snow_rdst_out_tle, snow_rdst_out_sbl,  &
      &      wtr_latflow_tle, wtr_latflow_sbl,      &
      &      snow_soil_tle, snow_soil_sbl,          &
      &      snow_new_tle, snow_new_sbl,            &
      &      snow_tle, snow_sbl,                    &
      &      snow_dens_tle, snow_dens_sbl           &
      &      )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: rhoh2o, dens_snow_min, dens_snow_max, dens_snow

    LOGICAL,  INTENT(in)    :: l_dynsnow               !< Dynamic snow density is being used

    REAL(wp), INTENT(in)    :: dtime,               &  !< Timestep length [s]
      &                        cnct_rel,            &  !< connectivity of tile and sibing relative to the sum of connectivities of tile [/]
      &                        h_offset,            &  !< Differece in surface elevation between tile and sibling [m]
      &                        tle_fract,           &  !< Area fraction of tile [/]
      &                        sbl_fract               !< Area fraction of sibling [/]

    REAL(wp), INTENT(in)    :: weq_pond_tle,        &  !< Pond storage on tile [m]
      &                        weq_pond_sbl            !< Pond storage on sibling [m]

    REAL(wp), INTENT(in)    :: snow_soil_tle_tm1,   &  !< Snow soil (WEQ) on tile at end of previous timestep [m]
      &                        snow_new_tle_tm1,    &  !< New snow soil (WEQ) on tile at end of previous timestep [m]
      &                        snow_dens_tle_tm1       !< Snow density on tile at end of previous timestep [kg m-3]


    REAL(wp), INTENT(inout) :: snow_rdst_in_tle,    &  !< Snow redistribution influx at tile [kg m-2 s-1]
      &                        snow_rdst_in_sbl,    &  !< Snow redistribution influx at sibling [kg m-2 s-1]
      &                        snow_rdst_out_tle,   &  !< Snow redistribution outflux at tile [kg m-2 s-1]
      &                        snow_rdst_out_sbl       !< Snow redistribution outflux at sibling [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: &
      &                        wtr_latflow_tle,     &  !< Netflux (surface) on tile [kg m-2 s-1]
      &                        wtr_latflow_sbl,     &  !< Netflux (surface) on sibling [kg m-2 s-1]
      &                        snow_soil_tle,       &  !< Snow soil (WEQ) on tile [m]
      &                        snow_soil_sbl,       &  !< Snow soil (WEQ) on sibling [m]
      &                        snow_new_tle,        &  !< New snow soil (WEQ) on tile [m]
      &                        snow_new_sbl,        &  !< New snow soil (WEQ) on sibling [m]
      &                        snow_tle,            &  !< Snow (WEQ; can and soil) on tile [m]
      &                        snow_sbl,            &  !< Snow (WEQ; can and soil) on sibling [m]
      &                        snow_dens_tle,       &  !< Snow density on tile [kg m-3]
      &                        snow_dens_sbl           !< Snow density on sibling [kg m-3]

    ! Local variables ------------------------------------------------------------------------------------

    REAL(wp)                :: fr_tle_rel,          &  !< Relative fraction of tile [/]
      &                        fr_sbl_rel,          &  !< Relative fraction of sibling [/]
      &                        h_offset_snow_fix,   &  !< Differece in surface elevation including immobile snow layer between tile and sibling [m]
      &                        snow_fix_tle ,       &  !< Immobile snow on tile [m]
      &                        snow_fix_sbl ,       &  !< Immobile snow on sibling [m]
      &                        snow_fix_rel_tle,    &  !< Fraction of immobile snow on tile [/]
      &                        snow_fix_rel_sbl,    &  !< Fraction of immobile snow on sibling [/]
      &                        snow_depth_mob_avg,  &  !< Average depth of mobile snow [m]
      &                        snow_depth_mob_tle,  &  !< Depth of mobile snow on tile [m]
      &                        snow_depth_mob_sbl,  &  !< Depth of mobile snow on sibling [m]
      &                        snow_mob_tle,        &  !< Mobile snow (WEQ) on tile [m]
      &                        snow_mob_sbl,        &  !< Mobile snow (WEQ)  on sibling [m]
      &                        snow_mob_upd_tle,    &  !< Updated mobile snow (WEQ) on tile [m]
      &                        snow_mob_upd_sbl,    &  !< Updated mobile snow (WEQ) on sibling [m]
      &                        snow_dens_max,       &  !< Maximum snow density [kg m-3]
      &                        snow_dens_min,       &  !< Minimum snow density [kg m-3]
      &                        snow_dens_fix_tle,   &  !< Density of immobile snow on tile [kg m-3]
      &                        snow_dens_fix_sbl,   &  !< Density of immobile snow on sibling [kg m-3]
      &                        snow_rdst_tle,       &  !< Snow (WEQ) redistribution flux on tile [kg m-2 s-1]
      &                        snow_rdst_sbl           !< Snow (WEQ) redistribution flux on sibling [kg m-2 s-1]

    IF (tle_fract < 1.0e-10_wp .OR. sbl_fract < 1.0e-10_wp) RETURN

    ! Setting initial values to 0.0
    snow_rdst_tle = 0._wp
    snow_rdst_sbl = 0._wp

    ! Snow density paraeters depending on snow scheme
    IF (l_dynsnow) THEN
      snow_dens_max = dens_snow_max
      snow_dens_min = dens_snow_min
    ELSE
      snow_dens_max = dens_snow
      snow_dens_min = dens_snow
    END IF

    ! Determine relative fractions from absolute fractions
    fr_tle_rel = tle_fract / (tle_fract+sbl_fract)
    fr_sbl_rel = sbl_fract / (tle_fract+sbl_fract)

    ! Determine amount of fixed and mobile fractions of the snow cover
    ! (assuming that in the previous time step the old snow has melted before the new snow)
    snow_fix_tle =  snow_soil_tle_tm1 - snow_new_tle_tm1
    snow_fix_sbl =  snow_soil_sbl - snow_new_sbl
    snow_mob_tle =  snow_new_tle_tm1
    snow_mob_sbl =  snow_new_sbl

    IF (snow_fix_tle < 0._wp) THEN
      snow_mob_tle = snow_new_tle_tm1 + snow_fix_tle
      snow_fix_tle = 0._wp
    END IF
    IF (snow_fix_sbl < 0._wp) THEN
      snow_mob_sbl = snow_new_sbl + snow_fix_sbl
      snow_fix_sbl = 0._wp
    END IF

    ! Determine snow density of fixed snow
    IF (snow_fix_tle > 0._wp) THEN
      snow_fix_rel_tle  = snow_fix_tle / snow_soil_tle_tm1
      snow_dens_fix_tle = (snow_dens_tle_tm1 - (snow_dens_min * (1._wp - snow_fix_rel_tle))) &
        &               / snow_fix_rel_tle
    ELSE
      snow_dens_fix_tle = snow_dens_tle_tm1
    END IF
    IF (snow_fix_sbl > 0._wp) THEN
      snow_fix_rel_sbl  = snow_fix_sbl / snow_soil_sbl
      snow_dens_fix_sbl = (snow_dens_sbl - (snow_dens_min * (1._wp - snow_fix_rel_sbl))) &
        &               / snow_fix_rel_sbl
    ELSE
      snow_dens_fix_sbl = snow_dens_sbl
    END IF

    ! Determine height differences of the top of fixed snow
    h_offset_snow_fix = h_offset                                                   &
      &               + weq_pond_tle + snow_fix_tle * (rhoh2o / snow_dens_fix_tle) &
      &               - weq_pond_sbl - snow_fix_sbl * (rhoh2o / snow_dens_fix_sbl)

    ! Calculate snow redistribution for two cases (1) from tile to sibling or (2) from sibling to tile
    ! Transport only applies if in the source tile there is fixed snow
    IF  ((h_offset_snow_fix >  1.0e-10_wp .AND. snow_fix_tle > 0.01_wp) .OR. &
         (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_sbl > 0.01_wp))THEN

      snow_depth_mob_avg = (snow_mob_tle * fr_tle_rel + snow_mob_sbl * fr_sbl_rel) &
        &                * (rhoh2o / snow_dens_min)

      ! 1. (Fixed) snow surface is lower in derpression and there is fixed and mobile snow on tile
      IF  (h_offset_snow_fix > 1.0e-10_wp .AND. snow_fix_tle > 0.01_wp)THEN
        IF (snow_depth_mob_avg / fr_sbl_rel <= h_offset_snow_fix ) THEN ! New snow fits into sibling
          snow_depth_mob_sbl = snow_depth_mob_avg / fr_sbl_rel
          snow_depth_mob_tle = 0._wp
        ELSE
          snow_depth_mob_sbl = h_offset_snow_fix  &
            &                + (snow_depth_mob_avg - h_offset_snow_fix * fr_sbl_rel)
          snow_depth_mob_tle =  snow_depth_mob_avg - h_offset_snow_fix * fr_sbl_rel
        END IF
      ! 2. (Fixed) snow surface is lower on tile and there is fixed snow in sibling
      ELSE IF (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_sbl > 0.01_wp) THEN
        h_offset_snow_fix = - h_offset_snow_fix ! Inversing sign of height offset to simplify calculations
        IF (snow_depth_mob_avg / fr_tle_rel <= h_offset_snow_fix ) THEN  ! New snow fits into tile
          snow_depth_mob_tle = snow_depth_mob_avg / fr_tle_rel
          snow_depth_mob_sbl = 0._wp
        ELSE
          snow_depth_mob_tle = h_offset_snow_fix &
                           & + (snow_depth_mob_avg - h_offset_snow_fix * fr_tle_rel)
          snow_depth_mob_sbl =  snow_depth_mob_avg - h_offset_snow_fix * fr_tle_rel
        END IF
      END IF

      ! Recalculation of snow weqs
      snow_mob_upd_sbl = snow_depth_mob_sbl / (rhoh2o / snow_dens_min)
      snow_mob_upd_tle = snow_depth_mob_tle / (rhoh2o / snow_dens_min)

      snow_rdst_sbl    = (snow_mob_upd_sbl - snow_mob_sbl) * cnct_rel
      snow_rdst_tle    = (snow_mob_upd_tle - snow_mob_tle) * cnct_rel

      snow_soil_sbl    = snow_soil_sbl    + snow_rdst_sbl
      snow_soil_tle    = snow_soil_tle    + snow_rdst_tle

      snow_sbl         = snow_sbl         + snow_rdst_sbl
      snow_tle         = snow_tle         + snow_rdst_tle

      snow_new_sbl     = snow_new_sbl     + snow_rdst_sbl
      snow_new_tle     = snow_new_tle     + snow_rdst_tle

      ! Converting redistribution flux from m/step -> kg m-2 s-1
      snow_rdst_sbl     = snow_rdst_sbl * rhoh2o / dtime
      snow_rdst_tle     = snow_rdst_tle * rhoh2o / dtime

      wtr_latflow_tle   = wtr_latflow_tle   + snow_rdst_tle
      wtr_latflow_sbl   = wtr_latflow_sbl   + snow_rdst_sbl

      snow_rdst_in_tle  = snow_rdst_in_tle  + MAX(0._wp, snow_rdst_tle)
      snow_rdst_in_sbl  = snow_rdst_in_sbl  + MAX(0._wp, snow_rdst_sbl)
      snow_rdst_out_tle = snow_rdst_out_tle - MIN(0._wp, snow_rdst_tle)
      snow_rdst_out_sbl = snow_rdst_out_sbl - MIN(0._wp, snow_rdst_sbl)

      ! Updating snow densities based on redistributed snow
      IF (snow_rdst_sbl .NE. 0_wp .AND. snow_soil_sbl > 0._wp) THEN
        snow_dens_sbl   = (snow_new_sbl * snow_dens_min &
          &             + snow_fix_sbl * min(snow_dens_fix_sbl, snow_dens_max)) / snow_soil_sbl
      END IF

      IF (snow_rdst_tle .NE. 0_wp .AND. snow_soil_tle > 0._wp) THEN
        snow_dens_tle   = (snow_new_tle * snow_dens_min &
          &             + snow_fix_tle * min(snow_dens_fix_tle, snow_dens_max)) / snow_soil_tle
      END IF
    END IF

  END SUBROUTINE calc_snow_redistribution

  !
  ! ===============================================================================================================================
  !>
  !! Limit soil water and pond storage content
  !!
  SUBROUTINE remove_excess_water(             &
    & nc,                                     &
    & nsoil,                                  &
    & weq_pond_max,                           &
    & ice_pond,                               &
    & wtr_pond,                               &
    & wtr_latflow_res_srf,                    &
    & wtr_soil_sat_sl,                        &
    & wtr_soil_sl                             &
    )

    INTEGER,  INTENT(in)    :: nc                             !< Number of cells
    INTEGER , INTENT(in)    :: nsoil                          !< Number of soil ayers
    REAL(wp), INTENT(in)    :: weq_pond_max(:),  &            !< Maximum pond water storage [m]
      &                        ice_pond(:)                    !< Ice in pond water storage [m]

    REAL(wp), INTENT(inout) :: wtr_pond(:),      &            !< Liquid water in pond water storage [m]
      &                        wtr_latflow_res_srf(:)         !< Surface water retention storage [m]

    REAL(wp), OPTIONAL, INTENT(in)    :: wtr_soil_sat_sl(:,:) !< Saturation soil moisture on soil layers [m]
    REAL(wp), OPTIONAL, INTENT(inout) :: wtr_soil_sl(:,:)     !< Water content on soil layers [m]]

    INTEGER                 :: ic, is                         !< Counters for loops over cells and soil layers
    REAL(wp)                :: wtr_pond_max                   !< Maximum liquid water in pond storage [m]

    IF (nsoil > 0) THEN
      DO is = nsoil, 2, -1
        DO ic = 1, nc
          IF (wtr_soil_sl(ic,is) > wtr_soil_sat_sl(ic,is)) THEN
            wtr_soil_sl(ic,is-1) = wtr_soil_sl(ic,is-1) &
              &                  + wtr_soil_sl(ic,is) - wtr_soil_sat_sl(ic,is)
            wtr_soil_sl(ic,is)   = wtr_soil_sat_sl(ic,is)
          END IF
        END DO
      END DO

      DO ic = 1, nc
        IF (wtr_soil_sl(ic,1) > wtr_soil_sat_sl(ic,1)) THEN
          wtr_pond(ic)      = wtr_pond(ic) &
            &               + wtr_soil_sl(ic,1) - wtr_soil_sat_sl(ic,1)
          wtr_soil_sl(ic,1) = wtr_soil_sat_sl(ic,1)
        END IF
      END DO
    END IF

    DO ic = 1, nc
      wtr_pond_max = MAX(0.0_wp, weq_pond_max(ic) - ice_pond(ic))
      IF (wtr_pond(ic) > wtr_pond_max) THEN
        wtr_latflow_res_srf(ic) = wtr_latflow_res_srf(ic) &
          &                     + wtr_pond(ic) - wtr_pond_max
        wtr_pond(ic)            = wtr_pond_max
      END IF
    END DO

  END SUBROUTINE remove_excess_water

  !
  ! ===============================================================================================================================
  !>
  !! Adjust lateral fluxes according to differences in surface elevation and bedrock boundary
  !!
  SUBROUTINE shift_fluxes_sl( &
    & nc,                     &
    & nsoil,                  &
    & hyd_flow_path,          &
    & h_offset,               &
    & dsoil_lvt, dsoil_dpr,   &
    & flx_lat_sl              &
    )

    INTEGER , INTENT(in)    :: nc,               &  !< Number of cells
      &                        nsoil                !< Number of soil ayers

    REAL(wp), INTENT(in)    :: hyd_flow_path(:), &  !< Connectivity factor for below-ground flow paths [/]
      &                        h_offset(:)          !< Difference in surface hights between elevation and depression [m]

    REAL(wp), INTENT(in)    :: dsoil_lvt(:,:),   &  !< Soil depth on layers of elevation tile [m]
      &                        dsoil_dpr(:,:)       !< Soil depth on layers of depression tile [m]

    REAL(wp), INTENT(inout) :: flx_lat_sl(:,:)     !< Lateralal flux to be adjusted [?]

    INTEGER  :: nrgrd,                &  !< Number of layers on helper grid determined by layer overlaps
      &         i, ic, is,            &  !< Counter variables
      &         i_lvt, i_dpr             !< Index variables for layers on in elevation (i_lvt) and in depression (i_dpr)

    REAL(wp) :: z_lvt,                &  !< Layer bottom of current layer on elevation tile [m]
      &         z_dpr,                &  !< Layer bottom of current layer on depression tile [m]
      &         z_rgrd_top,           &  !< Top of current layer on helper grid [m]
      &         z_rgrd_bot,           &  !< Bottom of current layer on helper grid [m]
      &         h                        !< Height of current layer [m]
    REAL(wp) :: dsoil_lvt_sum(nc),    &  !< Acc. soil depth on layers of elevation tile [m]
      &         dsoil_dpr_sum(nc)        !< Acc. soil depth on layers of depression tile + tile height diff [m]

    LOGICAL  :: l_regrid(nc)             !< Flag to indicate if regridding is needed

    REAL(wp) :: flx_lat_sl_upd(nc,nsoil) !< Updated lateral flux [?]


    ! Check where regridding is needed and initialize non-layered variables
    DO ic = 1, nc
      IF (hyd_flow_path(ic) > 0._wp .AND. h_offset(ic) > 0._wp) THEN
        l_regrid(ic) = .TRUE.
        dsoil_lvt_sum(ic) = 0._wp
        dsoil_dpr_sum(ic) = h_offset(ic)
      ELSE
        l_regrid(ic) = .FALSE.
      END IF
    END DO

    ! Initialize layered variables
    DO is = 1, nsoil
      DO ic = 1, nc
        IF (l_regrid(ic)) THEN
          flx_lat_sl_upd(ic,is) = 0._wp
          dsoil_lvt_sum(ic)     = dsoil_lvt_sum(ic) + dsoil_lvt(ic,is)
          dsoil_dpr_sum(ic)     = dsoil_dpr_sum(ic) + dsoil_dpr(ic,is)
        ELSE
          flx_lat_sl_upd(ic,is) = flx_lat_sl(ic,is)
        END IF
      END DO
    END DO

    ! For cells without flow paths or height differences, move lateral fluxes upwards in case
    ! they happen in layers below the total soil depth
    DO is = nsoil, 2, -1
      DO ic = 1, nc
        IF (.NOT. l_regrid(ic) .AND. dsoil_dpr(ic,is) <= 0._wp) THEN
          flx_lat_sl_upd(ic,is-1) = flx_lat_sl_upd(ic,is-1) + flx_lat_sl_upd(ic,is)
          flx_lat_sl_upd(ic,is)   = 0._wp
        END IF
      END DO
    END DO

    ! Compute cells with flow paths and height differences
    nrgrd = nsoil * 2 + 1
    DO i = 1, nrgrd
      DO ic = 1, nc
        IF (l_regrid(ic)) THEN
          ! Initialize first
          IF (i == 1) THEN
            i_lvt      = 1
            z_lvt      = dsoil_lvt(ic,1)
            i_dpr      = 1
            z_dpr      = dsoil_dpr(ic,1) + h_offset(ic)
            z_rgrd_top = 0._wp
            z_rgrd_bot = MIN(z_lvt, z_dpr)
          END IF
          ! Compute layer thickness and update lateral fluxes depending on layer fraction
          h = z_rgrd_bot - z_rgrd_top
          flx_lat_sl_upd(ic,i_dpr) = flx_lat_sl_upd(ic,i_dpr) + flx_lat_sl(ic,i_lvt) * h / dsoil_lvt(ic,i_lvt)

          ! Update layer information
          IF (z_rgrd_bot == z_lvt .AND. z_lvt < dsoil_lvt_sum(ic)) THEN
            i_lvt = i_lvt + 1
            IF (i_lvt <= nsoil) THEN
              z_lvt = z_lvt + dsoil_lvt(ic,i_lvt)
            END IF
          END IF
          IF (z_rgrd_bot == z_dpr .AND. z_dpr < dsoil_dpr_sum(ic)) THEN
            i_dpr = i_dpr + 1
            IF (i_dpr <= nsoil) THEN
              z_dpr = z_dpr + dsoil_dpr(ic,i_dpr)
            END IF
          ELSE IF (z_rgrd_bot >= dsoil_dpr_sum(ic)) THEN
            z_dpr = z_lvt
          END IF

          z_rgrd_top = z_rgrd_bot
          z_rgrd_bot = MIN(z_lvt, z_dpr)

        END IF
      END DO
    END DO

    ! Write updated lateral fluxes back to the original array
    DO is = 1, nsoil
      DO ic = 1, nc
        flx_lat_sl(ic,is) = flx_lat_sl_upd(ic,is)
      END DO
    END DO

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
    & microcon_ctc_tile,                        &
    & microcon_ctc_sibl,                        &
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

    !$ACC ROUTINE SEQ

    USE mo_jsb_math_constants,     ONLY: pi

    REAL(wp), INTENT(in)    :: microcon_ctc_tile,   &  !< Micro-scale connectivity between tile and sibling (if tile is encompased by sibl) [/]
      &                        microcon_ctc_sibl       !< Micro-scale connectivity between sibling and tile (if sibl is encompased by tile) [/]

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

    ! Calculate relative fractions
    IF ((fract + fract_sibl) == 0._wp) THEN
      fract_rel      = 0._wp
      fract_rel_sibl = 0._wp
    ELSE
      fract_rel      = fract      / (fract + fract_sibl)
      fract_rel_sibl = fract_sibl / (fract + fract_sibl)
    END IF

    ! Calculate difference in surface elevation between tile and sibling (only valid for any micro-connections)
    IF ((microcon_ctc_tile  == 0._wp .AND. microcon_ctc_sibl  == 0._wp) .OR.  &
      &  fract_rel          == 0._wp .OR.  fract_rel_sibl     == 0._wp) THEN
      h_offset = 0._wp
    ELSE
      h_offset = h_surf - h_surf_sibl
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
    IF ((microcon_ctc_tile == 0._wp .AND. microcon_ctc_sibl == 0._wp) .OR.  &
      &  fract_rel         == 0._wp .OR.  fract_rel_sibl    == 0._wp  .OR.  &
      &  area              == 0._wp .OR.  area_sibl         == 0._wp) THEN
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
    ELSE IF (microcon_ctc_tile > microcon_ctc_sibl) THEN ! Tile is encompassed by sibling
      IF (area_nclsd > 0._wp .AND. fract_nclsd > 1.0e-10_wp) THEN ! Tile encloses another tile

        area_cmbnd      = area_nclsd &
          &             * (fract_nclsd + fract) / fract_nclsd
        fract_rel_cmbnd = (fract_nclsd + fract) &
          &             / (fract_nclsd + fract + fract_sibl)

        area_fact = 1._wp / fract_rel_cmbnd - 1._wp

        length_iface = ((4._wp * pi * area_cmbnd) ** 0.5_wp) &
          &          * microcon_ctc_tile * fract_rel / area

        rad_in  = (area_nclsd    / pi)**0.5_wp

        rad_out = (area_cmbnd * (1._wp + area_fact / microcon_ctc_tile) / pi)**0.5_wp

        length_scale = (rad_in + rad_out) / 2._wp

      ELSE  ! Tile does not enclose other tile

        area_fact = 1._wp / fract_rel - 1._wp

        length_iface = ((4._wp * pi * area) ** 0.5_wp) &
          &          * microcon_ctc_tile * fract_rel / area

        rad_in  = (area / pi)**0.5_wp

        rad_out = (area *(1._wp + area_fact / microcon_ctc_tile) / pi)**0.5_wp

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
          &          * microcon_ctc_sibl * fract_rel_sibl / area_sibl

        rad_in  = (area_nclsd_sibl / pi)**0.5_wp

        rad_out = (area_cmbnd *(1._wp + area_fact / microcon_ctc_sibl) / pi)**0.5_wp

        length_scale = (rad_in + rad_out) / 2._wp

      ELSE   ! Sibling does not enclose other tile

        area_fact = 1._wp / fract_rel_sibl - 1._wp

        length_iface = ((4._wp * pi * area_sibl) ** 0.5_wp) &
          &          * microcon_ctc_sibl * fract_rel_sibl / area_sibl

        rad_in  = (area_sibl / pi)**0.5_wp

        rad_out = (area_sibl *(1._wp + area_fact / microcon_ctc_sibl) / pi)**0.5_wp

        length_scale = rad_in + (rad_out - rad_in) / 2._wp

      END IF
    END IF

    ! calculate micro-scale slope based on difference in surfce elevation and characteristic length scale
    ! ... also in case that microcon_ctc_sibl > since the microcon_slp will be reversed later on ...
    IF ((microcon_ctc_tile > 0._wp .OR. microcon_ctc_sibl > 0._wp) .AND. &
      &  h_offset > 0._wp .AND. length_scale > 0._wp) THEN
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
    & ret_macro_srf, ret_macro_blg,                  &
    & length_scale,                                  &
    & fract, fract_sibl,                             &
    & micro_connect_srf, micro_connect_blg,          &
    & macro_connect_srf, macro_connect_blg,          &
    & fract_microcon_srf, fract_microcon_blg,        &
    & fract_macrocon_srf, fract_macrocon_blg,        &
    & fact_lag_micro_srf, fact_lag_micro_blg,        &
    & fact_lag_macro_srf, fact_lag_macro_blg         &
    )

    USE mo_jsb_math_constants,     ONLY: pi

    INTEGER,  INTENT(in)    :: lex_idt,              &  !< lex id of tile
      &                        lex_ids                  !< lex id of sibling

    REAL(wp), INTENT(in)    ::                       &
      &                        velocity_srf,         &  !< Potential velocity of lateral flow through connected flow pathways [m s-1]
      &                        velocity_blg,         &  !< Potential velocity of lateral flow at the surface [m s-1]
      &                        ret_macro_srf,        &  !< Retention time for macro-scale fluxes at the surfcace [s]
      &                        ret_macro_blg            !< Retention time for macro-scale fluxes below the surfcace [s]

    REAL(wp), INTENT(in)    :: length_scale,         &  !< Assumed distance between centers of micro-scale connected tiles [m]
      &                        fract,                &  !< cover fraction of tile [/]
      &                        fract_sibl,           &  !< cover fraction of sibling [/]
      &                        micro_connect_srf,    &  !< Micro-scale connectivity (surface) between tile and sibling [/]
      &                        micro_connect_blg,    &  !< Micro-scale connectivity (below-ground) between tile and sibling [/]
      &                        macro_connect_srf,    &  !< Macro-scale connectivity (surface) between tile and sibling [/]
      &                        macro_connect_blg,    &  !< Macro-scale connectivity (below-ground) between tile and sibling [/]
      &                        fract_microcon_srf,   &  !< (Weighted) Sum of micro-scale connect. fractions (surface) on tile [/]
      &                        fract_microcon_blg,   &  !< (Weighted) Sum of micro-scale connect. fractions (below-ground) on tile [/]
      &                        fract_macrocon_srf,   &  !< (Weighted) Sum of macro-scale connect. fractions (surface) on tile [/]
      &                        fract_macrocon_blg       !< (Weighted) Sum of macro-scale connect. fractions (below-ground) on tile [/]

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

    IF (micro_connect_srf > 0._wp) THEN
      fact_lag_micro_srf = 1._wp / ret_micro_srf &
        &                * micro_connect_srf * fract_sibl / fract_microcon_srf
    ELSE
      fact_lag_micro_srf = 0.0_wp
    END IF

    IF (micro_connect_blg > 0._wp) THEN
      fact_lag_micro_blg = 1._wp / ret_micro_blg &
        &                * micro_connect_blg * fract_sibl / fract_microcon_blg
    ELSE
      fact_lag_micro_blg = 0.0_wp
    END IF

    ! Micro scale lag factors are determined by the (inverse) retention time, the micro-connectivity between tile and sibling,
    ! and the cover fraction of sibling relative to the total connected farction of tile.
    ! Additionally, the lag foctor is decreased accoring to the micro-scale connectivity, assuming that the fractions of the tile
    ! that have a micro-scale connectivity do not additionally feature macro-scale connections.
    IF ((macro_connect_srf > 0._wp) .AND. (fract > 0._wp)) THEN

      fact_lag_macro_srf = 1._wp / ret_macro_srf                       &
        &                * macro_connect_srf * fract_sibl / fract_macrocon_srf &
        &                * MAX(0._wp, (fract - fract_microcon_srf) / fract)
    ELSE
      fact_lag_macro_srf = 0.0_wp
    END IF

    IF ((macro_connect_blg > 0._wp) .AND. (fract > 0._wp)) THEN
      fact_lag_macro_blg = 1._wp / ret_macro_blg                       &
        &                * macro_connect_blg * fract_sibl / fract_macrocon_blg &
        &                * MAX(0._wp, (fract - fract_microcon_blg) / fract)
    ELSE
      fact_lag_macro_blg = 0.0_wp
    END IF

    ! For outflow on tiles that have no downstream tile, lag factors are calculated based on the time-scales of the macro-scale retention times.
    ! The respective fluxes are given to surface runoff and drainage.
    IF ((lex_idt == lex_ids)       .AND. &
      & (fract_microcon_srf == 0.0_wp) .AND. (fract_macrocon_srf == 0.0_wp)) THEN
      fact_lag_micro_srf = 0._wp
      fact_lag_macro_srf = 1._wp / ret_macro_srf
    END IF

    IF ((lex_idt == lex_ids)       .AND. &
      & (fract_microcon_blg == 0.0_wp) .AND. (fract_macrocon_blg == 0.0_wp)) THEN
      fact_lag_micro_blg = 0._wp
      fact_lag_macro_blg = 1._wp / ret_macro_blg
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
            trgt_id=is ! sibling is current target
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
