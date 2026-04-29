!> Contains the routines for the lateral exchange processes
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

MODULE mo_lex_process
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_control, ONLY: acc_stream

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport, calc_drainage_transport,                &
          & calc_lat_fluxes_pond, calc_lat_fluxes_soil,                    &
          & calc_convection_heat, calc_snow_redistribution,                &
          & remove_excess_water, shift_fluxes_sl,                          &
          & get_geometry, get_flow_lags, check_flow_matrix

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_process'

CONTAINS

  ! =================================================================================================================
  !>
  !>#### Setup intersected grid for the soil layers of two adjacent tiles
  !>
  !> As connected tiles can be vertically offset with respect to each other, the soil layers boundaries of both tiles
  !> may differ. This routine computes a new intersected grid that maps the soil layers of both tiles
  !> onto a common grid. The routine returns the mapping indices from the intersected grid to the soil layers of both
  !> tiles taking into account the vertical offset as well as the pond storage in the tile with lower elevation.
  !> It also returns information about the connection type (None, Pond, Soil, Rock) between the two tiles
  !> for each layer.
  !>
  SUBROUTINE get_intersected_grid(dz, nsoil, ngrid, offset, d_bedrock, d_pond, &
    &  index_tile, index_sibl, layer_thickness, status, layer_top, layer_bot)

    USE mo_util_sort,   ONLY: quicksort

    REAL(wp), INTENT(in)  :: dz(:)                !< Standard soil layer bottom depths in both tiles [m]
    INTEGER,  INTENT(in)  :: nsoil,            &  !< Number of standard soil layers per tile
      &                      ngrid                !< Maximum number of layers in intersected grid
    REAL(wp), INTENT(in)  :: offset,           &  !< Vertical offset between tile and sibling soil layers [m]
      &                      d_bedrock,        &  !< Depth to bedrock in both tiles [m]
      &                      d_pond               !< Water content of surface ponds [m]
    INTEGER,  INTENT(out) :: index_tile(:),    &  !< Index mapping from intersected grid to tile soil layers
      &                      index_sibl(:)        !< Index mapping from intersected grid to sibling soil layers
    REAL(wp), INTENT(out) :: layer_thickness(:)   !< Soil layer thicknesses on intersected grid [m]
    REAL(wp), INTENT(out), OPTIONAL ::         &
      & layer_top(:), &  !< Upper depth bound of soil layer on intersected grid [m]
      & layer_bot(:)     !< Lower depth bound of soil layer on intersected grid [m]
    CHARACTER(LEN=4), INTENT(out) :: status(:)    !< Connection type (None, Pond, Soil, Rock) between tile and sibling
    INTEGER               :: i,                &  !< Index variable for layers on any vertical grid
      &                      i_tile,           &  !< Current index position in tile soil grid
      &                      i_sibl               !< Current index position in sibling soil grid
    REAL(wp)              :: top, bot,         &  !< Depth boundaries of current layer on intersected grid [m]
      &                      br_tile,          &  !< Bedrock depth on tile [m]
      &                      br_sibl,          &  !< Bedrock depth on sibling [m]
      &                      dz_pond              !< Upper elevation of pond wrt tile [m]
    REAL(wp)              :: dz_tile(nsoil),   &  !< Bottom depth of soil layers at tile [m]
      &                      dz_sibl(nsoil)       !< Bottom depth of soil layers at sibling [m]
    REAL(wp)              :: new_layers(ngrid)    !< Intersecting grids of both soil layer [m]

    ! Initialize new grid variables
    DO i = 1, ngrid
      index_tile(i)      = -1
      index_sibl(i)      = -1
      layer_thickness(i) = 0._wp
      status(i)          = 'None'
    END DO

    ! Setup soil profiles using the tile as reference
    dz_tile(:) = dz(:)
    dz_sibl(:) = dz(:) + offset
    br_tile    = d_bedrock
    br_sibl    = d_bedrock + offset

    ! Setup top level elevation and indices
    IF (offset < 0._wp) THEN
      ! Sibling tile is above tile
      top = offset
      IF (d_pond + offset < 0) THEN
        ! Depression not completely filled with ponded water
        i_tile = -1
      ELSE
        ! Depression completely filled with ponded water
        i_tile = 0
      END IF
      i_sibl = 1
      dz_pond = -d_pond
    ELSE IF (offset > 0._wp) THEN
      ! Sibling tile is below tile
      top = 0._wp
      i_tile = 1
      IF (offset - d_pond > 0) THEN
        ! Depression not completely filled with ponded water
        i_sibl = -1
      ELSE
        ! Depression completely filled with ponded water
        i_sibl = 0
      END IF
      dz_pond = offset - d_pond
    ELSE
      ! Both tiles are at the same level
      top = 0._wp
      i_tile = 1
      i_sibl = 1
      dz_pond = 0._wp
    END IF

    ! Merge soil depth profiles into one common grid
    DO i = 1, nsoil
      new_layers(i)         = dz_tile(i)
      new_layers(nsoil + i) = dz_sibl(i)
    END DO
    ! Add pond depth to layers
    new_layers(2*nsoil + 1) = dz_pond
    ! Add lower soil surface to layers
    new_layers(2*nsoil + 2) = MAX(0._wp, offset)
    ! Add bedrock depth to layers
    new_layers(2*nsoil + 3) = br_tile
    new_layers(2*nsoil + 4) = br_sibl

    ! Order new layers
    CALL quicksort(new_layers)

    ! Loop over all layers in intersected grid
    DO i = 1, ngrid
      bot = new_layers(i)
      IF (PRESENT(layer_top)) layer_top(i) = top
      IF (PRESENT(layer_bot)) layer_bot(i) = bot
      IF (bot > top .AND. i_tile <= nsoil .AND. i_sibl <= nsoil) THEN
        ! There is a layer here --> store thickness and status and update indices
        layer_thickness(i) = bot - top
        index_tile(i) = i_tile
        index_sibl(i) = i_sibl
        ! Determine connection status between tile and sibling soil layers
        IF ((i_tile == 0 .AND. i_sibl > 0 .AND. bot <= br_sibl) .OR. &
          & (i_sibl == 0 .AND. i_tile > 0 .AND. bot <= br_tile)) THEN
          status(i) = 'Pond'
        ELSE IF (i_tile > 0 .AND. i_sibl > 0 .AND. bot <= br_tile .AND. bot <= br_sibl) THEN
          status(i) = 'Soil'
        ELSE IF (i_tile > 0 .AND. i_sibl > 0) THEN
          status(i) = 'Rock'
        END IF
        ! Update indices for next layer for tile
        IF (i_tile < 0 .AND. bot >= dz_pond .AND. d_pond > 0._wp) THEN
          ! Entering pond on tile --> set index to 0
          i_tile = 0
        ELSE IF (i_tile < 1 .AND. bot >= 0._wp) THEN
          ! Depression on tile but now entering soil column --> set index to 1
          i_tile = 1
        ELSE IF (i_tile > 0 .AND. bot >= dz_tile(MAX(i_tile,1))) THEN
          ! Entering next soil layer on tile
          i_tile = i_tile + 1
        END IF
        ! Update indices for next layer for sibling tile
        IF (i_sibl < 0 .AND. bot >= dz_pond .AND. d_pond > 0._wp) THEN
          ! Entering pond layer on sibling tile --> set index to 0
          i_sibl = 0
        ELSE IF (i_sibl < 1 .AND. bot >= offset) THEN
          ! Depression on sibling but now entering soil column --> set index to 1
          i_sibl = 1
        ELSE IF (i_sibl > 0 .AND. bot >= dz_sibl(MAX(i_sibl,1))) THEN
          ! Entering next soil layer on tile
          i_sibl = i_sibl + 1
        END IF
      END IF
      ! Set new layer top for next iteration
      top = MAX(bot, top)
    END DO

  END SUBROUTINE get_intersected_grid

  ! =================================================================================================================
  !>
  !>#### Compute lateral transport of surface runoff between connected tiles
  !>
  !> Surface runoff stored in the intermediate storage on the source tile is transported to the target tile based on
  !> a predefined flow fraction values (representing a flow lag). The transported water increases the pond storage
  !> on the target tile.
  !>
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

    REAL(wp), INTENT(in)    :: dtime    !< Time step length [s]

    REAL(wp), INTENT(in)    ::       &
      &  fract_source_2_target,      &  !< Fraction of outflow from source tile received by target tile [/]
      &  fract_source,               &  !< Fraction of source tile [/]
      &  fract_target,               &  !< Fraction of target tile [/]
      &  wtr_latflow_srf                !< Outflow of intermediate storage on source tile [m]

    REAL(wp), INTENT(inout) ::       &
      &  wtr_lex_runoff_target_acc,  &  !< Inflow of intermediate storage received on target tile [kg m-2 s-1]
      &  wtr_lex_runoff_source_acc,  &  !< Outflow of intermediate storage received by target tiles [kg m-2 s-1]
      &  wtr_lex_srf_source_acc,     &  !< Netflux (surface) on source tile [kg m-2 s-1]
      &  wtr_lex_srf_target_acc,     &  !< Netflux (surface) on target tile [kg m-2 s-1]
      &  wtr_pond_target                !< Surface water storage on target tile    [m]

    ! Skip if either source or target tile fraction is negligible
    IF (fract_source < 1.0e-10_wp .OR. fract_target < 1.0e-10_wp) RETURN

    ! Move water from source tile intermediate storage to target tile pond storage
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

  ! =================================================================================================================
  !>
  !>#### Compute lateral transport of subsurface drainage between connected tiles
  !>
  !> Subsurface drainage stored in the intermediate storage on the source tile is transported to the target tile
  !> based on a predefined flow fraction values (representing a flow lag). Note that the water in the layers of
  !> the intermediate storage was already rescaled taking into account the vertical offset between the two tiles,
  !> and the corresponding offset between the soil layers. For this reason, the rescaled transported water can arrive
  !> in the same layer on the target tile as it left the source tile, and still can be assumed to flow horizontally
  !> through the soil.
  !>
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

    REAL(wp), INTENT(in)    :: dtime  !< Time step length [s]

    REAL(wp), INTENT(in)    ::     &
      & fract_source_2_target,     &  !< Fraction of outflow from source tile received by target tile [/]
      & fract_source,              &  !< Fraction of source tile [/]
      & fract_target,              &  !< Fraction of target tile [/]
      & wtr_latflow_sl                !< Outflow of intermediate storage on source tile [m]


    REAL(wp), INTENT(inout) ::     &
      & wtr_lex_drain_target_acc,  &  !< Inflow of intermediate storage received on target tile [kg m-2 s-1]
      & wtr_lex_drain_source_acc,  &  !< Outflow of intermediate storage received by target tiles [kg m-2 s-1]
      & wtr_lex_sl_source_acc,     &  !< Netflux (below ground) on source tile [kg m-2 s-1]
      & wtr_lex_sl_target_acc,     &  !< Netflux (below ground) on target tile [kg m-2 s-1]
      & wtr_soil_sl_target            !< Soil moisture on target tile [m]

    ! Skip if either source or target tile fraction is negligible
    IF (fract_source < 1.0e-10_wp .OR. fract_target < 1.0e-10_wp) RETURN

    ! Move water from the source tile intermediate storage to the same layer of the target tile
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

  ! =================================================================================================================
  !>
  !>#### Compute lateral transport of water and associated heat from ponds into neighbouring soil layers
  !>
  !> If the tiles are offset against each other and the tile with lower elevation has ponded water,
  !> this water will infiltrate laterally into the soil layers of the tile with higher elevation. Furthermore,
  !> the heat contained in the ponded water is also transported into the soil layers of the adjacent tile.
  !>
  SUBROUTINE calc_lat_fluxes_pond(                       &
      &      nc, nsoil, nrgrd, dtime,                    &
      &      l_diff_water, l_diff_heat,                  &
      &      connect,                                    &
      &      h_offset, length_scale, length_iface,       &
      &      lvt_fract, dpr_fract, dpr_fract_pond_max,   &
      &      dsoil_lvt, field_cap_lvt,  k_lvt, mpot_lvt, &
      &      dz, t_soil_lvt, t_soil_dpr,                 &
      &      heat_cond_soil_lvt, heat_cap_soil_lvt,      &
      &      heat_cap_soil_dpr, ice_pond_dpr,            &
      &      egress_dpr,                                 &
      &      ingress_lvt, wtr_lat_in_sl_lvt,             &
      &      wtr_latflow_dpr, wtr_latflow_lvt,           &
      &      wtr_pond_dpr, wtr_soil_lvt,                 &
      &      diffus_in_dpr_top,  diffus_out_dpr_top,     &
      &      heat_flx_lat_dpr_top,                       &
      &      diffus_in_lvt, diffus_out_lvt,              &
      &      heat_flx_lat_lvt                            &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_sse_constants,          ONLY: hcond_water, hcond_ice

    INTEGER,  INTENT(in)    ::    &
      &  nc,                      &  !< Number of grid cells
      &  nsoil,                   &  !< Number of soil layers
      &  nrgrd                       !< Number of layers on helper grid determined by layer overlaps
    REAL(wp), INTENT(in)    :: dtime !< Timestep length [s]
    LOGICAL,  INTENT(in)    ::    &
      & l_diff_water,             &  !< T: compute diffusive water fluxes
      & l_diff_heat                  !< T: compute diffusive heat fluxes
    REAL(wp), INTENT(in)    ::    &
      & connect(:),               &  !< Connectivity array for lateral exchange
      & h_offset(:),              &  !< Difference in surface elevation between elevation and depression [m]
      & length_scale(:),          &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      & length_iface(:),          &  !< Length of interface between elevation and depression [m m-2]
      & lvt_fract(:),             &  !< Area fraction of elevation tile [/]
      & dpr_fract(:),             &  !< Area fraction of depression tile [/]
      & dpr_fract_pond_max(:)        !< Maximum area fraction of ponded water on depression tile [/]
    REAL(wp), INTENT(in)    ::    &
      & dsoil_lvt(:,:),           &  !< Soil depth on layers of elevation tile [m]
      & field_cap_lvt(:,:),       &  !< Soil field capacity of elevation tile [m]
      & k_lvt(:,:),               &  !< Soil layer hyd. conductivity on source tile [m s-1]
      & mpot_lvt(:,:)                !< Matric potential on source tile [m]
    REAL(wp), INTENT(in)    :: dz(:) !< Lower depth bounds for soil layers [m]
    REAL(wp), INTENT(in)    ::    &
      & t_soil_lvt(:,:),          &  !< Soil temperature of tile [K]
      & t_soil_dpr(:,:),          &  !< Soil temperature of sibling tile [K]
      & heat_cond_soil_lvt(:,:),  &  !< Soil heat conductivity on tile [J m-1 s-1 K-1]
      & heat_cap_soil_lvt(:,:),   &  !< Soil heat capacity on tile [J m-3 K-1]
      & heat_cap_soil_dpr(:,:)       !< Soil heat capacity on sibling [J m-3 K-1]
    REAL(wp), INTENT(in)    :: ice_pond_dpr(:) !< Ice content of ponded water on depression tile [m]
    REAL(wp), INTENT(inout) :: egress_dpr(:)   !< Diffusion water flux from depression tile [kg m-2 s-1]
    REAL(wp), INTENT(inout) ::    &
      & ingress_lvt(:,:),         &  !< Diffusion water flux to elevation tile [kg m-2 s-1]
      & wtr_soil_lvt(:,:),        &  !< Soil water storage on elevation tile [m]
      & wtr_latflow_lvt(:,:),     &  !< Netflux (below ground) on elevation tile [kg m-2 s-1]
      & wtr_lat_in_sl_lvt(:,:)       !< Lateral water influx in soil layers at tile [kg m-2 s-1]
    REAL(wp), INTENT(inout) ::    &
      & wtr_pond_dpr(:),          &  !< Pond storage on depression tile [m]
      & wtr_latflow_dpr(:)           !< Netflux (surface) on depression tile [kg m-2 s-1]
    REAL(wp), INTENT(inout) ::    &
      & diffus_in_dpr_top(:),     &  !< Diffusion influx at depression [J m-2 s-1]
      & diffus_out_dpr_top(:),    &  !< Diffusion outflux at depression [J m-2 s-1]
      & heat_flx_lat_dpr_top(:)      !< Net Heat flux at depression [J m-2 s-1]
    REAL(wp), INTENT(inout) ::    &
      & diffus_in_lvt(:,:),       &  !< Diffusion influx at elevation [J m-2 s-1]
      & diffus_out_lvt(:,:),      &  !< Diffusion outflux at elevation [J m-2 s-1]
      & heat_flx_lat_lvt(:,:)        !< Net Heat flux at elevation [J m-2 s-1]

    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 ::    &
      & i,                        &  !< Loop index for layers on intersected grid
      & ic,                       &  !< Loop index for grid cells
      & i_lvt                        !< Loop index soil layer in in elevated tile
    INTEGER                 ::    &
      & index_tile_wtr(nc,nrgrd), &  !< Mapping index from helper grid to tile layers for pond water
      & index_tile_weq(nc,nrgrd), &  !< Mapping index from helper grid to tile layers for pond water equivalent
      & index_sibl(nc,nrgrd)         !< Mapping index from helper grid to sibling layers
    REAL(wp)                ::    &
      & fr_lvt_rel(nc),           &  !< Relative fraction of elevation tile [/]
      & fr_dpr_rel(nc),           &  !< Relative fraction of depression tile [/]
      & weq_pond_dpr(nc),         &  !< Water and ice content of surface reservoir on depression [m]
      & h_cond_pond_dpr(nc),      &  !< Heat conductivity of surface reservoir on depression [J m-1 s-1 K-1]
      & sum_soil_tle(nc)             !< Total soil depth on tile [m]
    REAL(wp)                ::    &
      & diffusion(nc, nsoil)         !< Diffusive flux between depression and elevation on soil layers [m]
                                     ! TODO: check whether 'diffusion' is really the correct name for this flux
    REAL(wp)                ::    &
      & layer_thickness_wtr(nc,nrgrd), &  !< Thickness of layers on helper grid for pond water [m]
      & layer_thickness_weq(nc,nrgrd)     !< Thickness of layers on helper grid for pond water equivalent [m]
    CHARACTER(LEN=4)        ::    &
      & contype_wtr(nc,nrgrd),    &  !< Type of connection between tile and sibling for pond water
      & contype_weq(nc,nrgrd)        !< Type of connection between tile and sibling for pond water equivalent
    REAL(wp)                :: &
      & h,                     &  !< Height of current layer [m]
      & a_inter,               &  !< Interface area on current layer [m2 m-2]
      & wtr_soil_lvt_rgrd,     &  !< Soil water content on elevation tile on helper grid layer [m]
      & field_cap_lvt_rgrd,    &  !< Field capacity on elevation tile on helper grid layer [m]
      & mpot_lvt_rgrd,         &  !< Matric potential  on elevation tile on helper grid layer [m]
      & k_lvt_rgrd,            &  !< Hydraulic conductivity on elevation tile on helper grid layer [[m s-1]
      & free_space,            &  !< Free pore space on elevation tile on a given layer [m]
      & wtr_avail,             &  !< Available water on depression tile on a given layer [m]
      & wtr_flx_lat,           &  !< Potential flux (based on k and mpot) between depression and elevation [m/s]
      & wtr_lat,               &  !< Actual flux between depression and elevation [m]
      & t_lvt_rgrd,            &  !< Soil tempertures of elevation on current layer on helper grid [K]
      & h_cond_soil_lvt_rgrd,  &  !< Soil heat conductivity of elevation on current layer on helper grid [J m-1 s-1 K-1]
      & h_cap_soil_lvt_rgrd,   &  !< Soil heat capacity of elevation on current layer on helper grid [J m-3 K-1]
      & h_cond_inter,          &  !< Soil heat conductivity at interface on current layer on helper grid [J m-1 s-1 K-1]
      & d_t,                   &  !< Temperature difference on current layer on helper grid [K]
      & heat_flx_lat,          &  !< Heat flux between depression and elevation [J m-2 s-1]
      & t_eq,                  &  !< Equilibrium temperature [K]
      & heat_flx_eq_lat           !< Equilirium energy flux between depression and elevation [J m-2 s-1]
    LOGICAL                 :: &
      & l_compute_water(nc),   &  !< T: water fluxes should be computed for a given cell
      & l_compute_heat(nc)        !< T: heat fluxes should be computed for a given cell

    ! Initilize fields for accumulating regridded output
    diffusion(:,:)  = 0._wp

    weq_pond_dpr(:) = wtr_pond_dpr(:) + ice_pond_dpr(:)
    sum_soil_tle(:) = SUM(dsoil_lvt(:,:), DIM=2)

    ! Loop over grid cells and check general conditions to compute lateral fluxes
    DO ic = 1, nc
      IF (lvt_fract(ic) < 1.0e-10_wp .OR. dpr_fract(ic) < 1.0e-10_wp .OR. &
        & .NOT. connect(ic) > 0._wp .OR. h_offset(ic) < 1.0e-10_wp) THEN
        ! No lateral fluxes possible under these conditions
        l_compute_water(ic) = .FALSE.
        l_compute_heat(ic)  = .FALSE.
      ELSE
        ! Basic conditions for lateral transport are met
        l_compute_water(ic) = .TRUE.
        l_compute_heat(ic)  = .TRUE.
        ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
        fr_lvt_rel(ic) = lvt_fract(ic) / (lvt_fract(ic) + dpr_fract(ic))
        fr_dpr_rel(ic) = dpr_fract(ic) / (lvt_fract(ic) + dpr_fract(ic))
      END IF

      ! Check special conditions to not compute water or heat fluxes
      IF ((l_compute_water(ic) .AND. wtr_pond_dpr(ic) < 1.0e-10_wp) .OR. .NOT. l_diff_water) THEN
        l_compute_water(ic) = .FALSE.
      END IF
      IF ((l_compute_heat(ic) .AND. weq_pond_dpr(ic) < 1.0e-10_wp) .OR. .NOT. l_diff_heat) THEN
        l_compute_heat(ic) = .FALSE.
      END IF

      ! Get intersected grids for water and heat flux computations
      IF (l_compute_water(ic)) THEN
        wtr_lat_in_sl_lvt(ic,:) = wtr_lat_in_sl_lvt(ic,:) - ingress_lvt(ic,:)
        CALL get_intersected_grid(dz(:), nsoil, nrgrd, h_offset(ic), sum_soil_tle(ic), wtr_pond_dpr(ic), &
          & index_tile_wtr(ic,:), index_sibl(ic,:), layer_thickness_wtr(ic,:), contype_wtr(ic,:))
      END IF
      IF (l_compute_heat(ic)) THEN
        h_cond_pond_dpr(ic) = (wtr_pond_dpr(ic) * hcond_water + ice_pond_dpr(ic) * hcond_ice) / weq_pond_dpr(ic)
        CALL get_intersected_grid(dz(:), nsoil, nrgrd, h_offset(ic), sum_soil_tle(ic), weq_pond_dpr(ic), &
          & index_tile_weq(ic,:), index_sibl(ic,:), layer_thickness_weq(ic,:), contype_weq(ic,:))
      END IF
    END DO

    ! Loop over intersected grid layers and grid cells to compute lateral fluxes
    DO i = 1, nrgrd
      DO ic = 1, nc

        IF (l_compute_water(ic)) THEN
          ! Get layer index and thickness on elevation tile
          i_lvt = index_tile_wtr(ic,i)
          h = layer_thickness_wtr(ic,i)
          ! Only compute lateral water fluxes if ponded water is present on depression tile
          IF (contype_wtr(ic,i) == 'Pond') THEN
            IF (dsoil_lvt(ic,i_lvt) > 0._wp) THEN
              a_inter = length_iface(ic) * dpr_fract_pond_max(ic) * h
              ! Computes water states and properties
              wtr_soil_lvt_rgrd   = wtr_soil_lvt(ic,i_lvt)  * h / dsoil_lvt(ic,i_lvt)
              wtr_lat             = 0._wp
              field_cap_lvt_rgrd  = field_cap_lvt(ic,i_lvt) * h / dsoil_lvt(ic,i_lvt)
              k_lvt_rgrd          = k_lvt(ic,i_lvt)
              mpot_lvt_rgrd       = mpot_lvt(ic,i_lvt)
              ! Water flows from the pond reservoir of the depression tile into the soil of the elevated tile
              wtr_flx_lat         = abs(k_lvt_rgrd  * mpot_lvt_rgrd) / length_scale(ic) * a_inter
              free_space          = MAX(0._wp,(field_cap_lvt_rgrd - wtr_soil_lvt_rgrd) * fr_lvt_rel(ic))
              wtr_avail           = MIN(h, wtr_pond_dpr(ic)) * fr_dpr_rel(ic)
              wtr_lat             = MIN(wtr_flx_lat * dtime, free_space, wtr_avail)
              diffusion(ic,i_lvt) = diffusion(ic,i_lvt) + wtr_lat
            END IF
          END IF
        END IF

        IF (l_compute_heat(ic)) THEN
          ! Get layer index and thickness on elevation tile
          i_lvt = index_tile_weq(ic,i)
          h = layer_thickness_weq(ic,i)
          ! Only compute lateral water fluxes if ponded water is present on depression tile
          IF (contype_weq(ic,i) == 'Pond') THEN
            a_inter = length_iface(ic) * dpr_fract_pond_max(ic) * h
            ! Computes temperature states and thermal properties
            t_lvt_rgrd           = t_soil_lvt(ic,i_lvt)
            h_cond_soil_lvt_rgrd = heat_cond_soil_lvt(ic,i_lvt)
            h_cap_soil_lvt_rgrd  = heat_cap_soil_lvt(ic,i_lvt)
            ! Heat flux is determined by the difference in temperature between the two tiles
            h_cond_inter         =(h_cond_soil_lvt_rgrd * h_cond_pond_dpr(ic))**0.5_wp
            d_t                  = t_lvt_rgrd - t_soil_dpr(ic,1)
            heat_flx_lat         = h_cond_inter  * d_t / length_scale(ic) * a_inter

            ! Limiting heat flux to avoid instablilities for long time steps
            ! Determine flux required to equilibriate temperature within a timestep
            ! Equilibrium temperature
            t_eq = (t_lvt_rgrd       * h_cap_soil_lvt_rgrd     * fr_lvt_rel(ic)  &
              &  +  t_soil_dpr(ic,1) * heat_cap_soil_dpr(ic,1) * fr_dpr_rel(ic)) &
              &  / (h_cap_soil_lvt_rgrd * fr_lvt_rel(ic) + heat_cap_soil_dpr(ic,1) * fr_dpr_rel(ic))
            ! Heat flux to reach equilibrium temperature within a timestep
            heat_flx_eq_lat = MIN( &
              &               ABS((t_lvt_rgrd        - t_eq) * h_cap_soil_lvt_rgrd     * fr_lvt_rel(ic)), &
              &               ABS(( t_soil_dpr(ic,1) - t_eq) * heat_cap_soil_dpr(ic,1) * fr_dpr_rel(ic))) &
              &               *  h / dtime
            ! Limiting heat flux to equilibrium flux
            IF (heat_flx_lat > 0._wp .AND. heat_flx_lat > ABS(heat_flx_eq_lat)) THEN
              heat_flx_lat = ABS(heat_flx_eq_lat)
            ELSE IF (heat_flx_lat < 0._wp .AND. heat_flx_lat < -ABS(heat_flx_eq_lat)) THEN
              heat_flx_lat = -ABS(heat_flx_eq_lat)
            END IF

            ! Diagnosing incoming-, outgoing- and net heat fluxes on tiles
            diffus_in_lvt(ic,i_lvt)    = diffus_in_lvt(ic,i_lvt)  + MAX(0._wp, -heat_flx_lat / fr_lvt_rel(ic))
            diffus_in_dpr_top(ic)      = diffus_in_dpr_top(ic)    + MAX(0._wp,  heat_flx_lat / fr_dpr_rel(ic))
            diffus_out_lvt(ic,i_lvt)   = diffus_out_lvt(ic,i_lvt) - MIN(0._wp, -heat_flx_lat / fr_lvt_rel(ic))
            diffus_out_dpr_top(ic)     = diffus_out_dpr_top(ic)   - MIN(0._wp,  heat_flx_lat / fr_dpr_rel(ic))
            heat_flx_lat_lvt(ic,i_lvt) = heat_flx_lat_lvt(ic,i_lvt)           - heat_flx_lat / fr_lvt_rel(ic)
            heat_flx_lat_dpr_top(ic)   = heat_flx_lat_dpr_top(ic)             + heat_flx_lat / fr_dpr_rel(ic)
          END IF
        END IF
      END DO
    END DO

    DO ic = 1, nc
      IF (l_compute_water(ic)) THEN
        ! Diagnosing incoming-, outgoing- and net water fluxes on tiles and update storages
        wtr_soil_lvt(ic,:)      = wtr_soil_lvt(ic,:)    +   diffusion(ic,:)  / fr_lvt_rel(ic)
        wtr_pond_dpr(ic)        = wtr_pond_dpr(ic)    - SUM(diffusion(ic,:)) / fr_dpr_rel(ic)
        wtr_latflow_lvt(ic,:)   = wtr_latflow_lvt(ic,:) +   diffusion(ic,:)  / fr_lvt_rel(ic) * rhoh2o / dtime
        wtr_latflow_dpr(ic)     = wtr_latflow_dpr(ic) - SUM(diffusion(ic,:)) / fr_dpr_rel(ic) * rhoh2o / dtime
        ingress_lvt(ic,:)       = ingress_lvt(ic,:)     +   diffusion(ic,:)  / fr_lvt_rel(ic) * rhoh2o / dtime
        egress_dpr(ic)          = egress_dpr(ic)      + SUM(diffusion(ic,:)) / fr_dpr_rel(ic) * rhoh2o / dtime
        wtr_lat_in_sl_lvt(ic,:) = wtr_lat_in_sl_lvt(ic,:) + ingress_lvt(ic,:)
      END IF
    END DO

  END SUBROUTINE calc_lat_fluxes_pond

  ! =================================================================================================================
  !>
  !>#### Compute lateral fluxes of water and heat between two adjacent soil columns
  !>
  !> This routine compute the lateral diffusion of water and heat between the soil columns of
  !> two adjacent tiles. This needs to take into account that the soil layers of the two tiles
  !> may be offset against each other, and that the soil layers may not align. In this case, an
  !> intersected grid is created to compute the lateral fluxes between the overlapping parts of
  !> the two soil column layers.
  !>
  SUBROUTINE calc_lat_fluxes_soil(                       &
      &      nc, nsoil, nrgrd, dtime,                    &
      &      l_diff_water, l_diff_heat,                  &
      &      connect, h_offset, length_scale, len_inter, &
      &      tle_fract, sbl_fract,                       &
      &      dsoil_tle, dsoil_sbl,                       &
      &      field_cap_tle, field_cap_sbl,               &
      &      wsat_tle, wsat_sbl,                         &
      &      res_tle, res_sbl,                           &
      &      k_sat_tle, k_sat_sbl,                       &
      &      k_tle, k_sbl,                               &
      &      mpot_tle, mpot_sbl,                         &
      &      dz,                                         &
      &      t_soil_tle, t_soil_sbl,                     &
      &      heat_cond_soil_tle, heat_cond_soil_sbl,     &
      &      heat_cap_soil_tle, heat_cap_soil_sbl,       &
      &      wtr_diff_in_tle, wtr_diff_in_sbl,           &
      &      wtr_diff_out_tle, wtr_diff_out_sbl,         &
      &      wtr_lat_in_sl_tle,  wtr_lat_in_sl_sbl,      &
      &      wtr_latflow_tle, wtr_latflow_sbl,           &
      &      wtr_soil_tle, wtr_soil_sbl,                 &
      &      nrg_diff_in_tle, nrg_diff_in_sbl,           &
      &      nrg_diff_out_tle, nrg_diff_out_sbl,         &
      &      heat_flx_lat_tle, heat_flx_lat_sbl          &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER, INTENT(in)     ::     &
      & nc,                        &  !< Number of grid cells
      & nsoil,                     &  !< Number of soil layers
      & nrgrd                         !< Number of layers on helper grid determined by layer overlaps
    REAL(wp), INTENT(in)    :: dtime  !< Timestep length [s]
    LOGICAL,  INTENT(in)    ::     &
      & l_diff_water,              &  !< T: compute diffusive water fluxes
      & l_diff_heat                   !< T: compute diffusive heat fluxes
    REAL(wp), INTENT(in)    ::     &
      & connect(:),                &  !< Connectivity array for lateral exchange
      & h_offset(:),               &  !< Differece in surface elevation between tile and sibling [m]
      & length_scale(:),           &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      & len_inter(:),              &  !< Length of interface between tile and sibling [m m-2]
      & tle_fract(:),              &  !< Area fraction of tile [/]
      & sbl_fract(:)                  !< Area fraction of sibling [/]
    REAL(wp), INTENT(in)    ::     &
      & dsoil_tle(:,:),            &  !< Soil depth on layers of tile [m]
      & dsoil_sbl(:,:),            &  !< Soil depth on layers of sibling [m]
      & field_cap_tle(:,:),        &  !< Soil field capacity of tile [m]
      & field_cap_sbl(:,:),        &  !< Soil field capacity of sibling [m]
      & wsat_tle(:,:),             &  !< Soil moisture at saturation sibling [m]
      & wsat_sbl(:,:),             &  !< Soil moisture at saturation tile [m]
      & res_tle(:,:),              &  !< Residual soil moisture tile [m]
      & res_sbl(:,:),              &  !< Residual soil moisture sibling [m]
      & k_sat_tle(:,:),            &  !< Saturation hydraulic conductivity tile [m s-1]
      & k_sat_sbl(:,:),            &  !< Saturation hydraulic conductivity sibling [m s-1]
      & k_tle(:,:),                &  !< Soil layer hyd. conductivity on source tile [m s-1]
      & k_sbl(:,:),                &  !< Soil layer hyd. conductivity on target tile [m s-1]
      & mpot_tle(:,:),             &  !< Soil layer matric potential on source tile [m]
      & mpot_sbl(:,:)                 !< Soil layer matric potential on target tile [m]
    REAL(wp), INTENT(in)    :: dz(:)  !< Lower depth bounds for soil layers [m]
    REAL(wp), INTENT(in)    ::     &
      & t_soil_tle(:,:),           &  !< Soil tempertures of tile [K]
      & t_soil_sbl(:,:),           &  !< Soil tempertures of sibling [K]
      & heat_cond_soil_tle(:,:),   &  !< Soil heat conductivity of tile [J m-1 s-1 K-1]
      & heat_cond_soil_sbl(:,:),   &  !< Soil heat conductivity of sibling [J m-1 s-1 K-1]
      & heat_cap_soil_tle(:,:),    &  !< Soil heat capacity of tile [J m-3 K-1]
      & heat_cap_soil_sbl(:,:)        !< Soil heat capacity of sibling [J m-3 K-1]
    REAL(wp), INTENT(inout) ::     &
      & wtr_diff_in_tle(:,:),      &  !< Diffusive water influx at tile [kg m-2 s-1]
      & wtr_diff_in_sbl(:,:),      &  !< Diffusive water influx at sibling [kg m-2 s-1]
      & wtr_diff_out_tle(:,:),     &  !< Diffusive water outflux at tile [kg m-2 s-1]
      & wtr_diff_out_sbl(:,:),     &  !< Diffusive water outflux at sibling [kg m-2 s-1]
      & wtr_lat_in_sl_tle(:,:),    &  !< Lateral water influx in soil layers at tile [kg m-2 s-1]
      & wtr_lat_in_sl_sbl(:,:),    &  !< Lateral water influx in soil layers at sibling [kg m-2 s-1]
      & wtr_latflow_tle(:,:),      &  !< Netflux (below ground) on tile [kg m-2 s-1]
      & wtr_latflow_sbl(:,:),      &  !< Netflux (below ground) on sibling [kg m-2 s-1]
      & wtr_soil_tle(:,:),         &  !< Soil water storage on tile [m]
      & wtr_soil_sbl(:,:)             !< Soil water storage on sibling [m]
    REAL(wp), INTENT(inout) ::     &
      & nrg_diff_in_tle(:,:),      &  !< Diffusive energy influx at tile [J m-2 s-1]
      & nrg_diff_in_sbl(:,:),      &  !< Diffusive energy influx at sibling [J m-2 s-1]
      & nrg_diff_out_tle(:,:),     &  !< Diffusive energy outflux at tile [J m-2 s-1]
      & nrg_diff_out_sbl(:,:),     &  !< Diffusive energy outflux at sibling [J m-2 s-1]
      & heat_flx_lat_tle(:,:),     &  !< Net lateral heat flux (below ground) on tile [J m-2 s-1]
      & heat_flx_lat_sbl(:,:)         !< Net lateral heat flux (below ground) on sibling [J m-2 s-1]

    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 ::     &
      & i,                         &  !< Index variable for layers on helper grid
      & ic,                        &  !< Index variable for grid cells
      & i_tle,                     &  !< Index variable for layers on tile
      & i_sbl                         !< Index variable for layers sibling
    INTEGER                 ::     &
      & index_tile(nc,nrgrd),      &  !< Mapping index from helper grid to tile layers
      & index_sibl(nc,nrgrd)          !< Mapping index from helper grid to sibling layers

    REAL(wp)                :: layer_thickness(nc,nrgrd)  !< Thickness of layers on helper grid [m]
    CHARACTER(LEN=4)        :: contype(nc,nrgrd)          !< Type of connection between tile and sibling
    REAL(wp)                ::     &
      & fr_tle_rel(nc),            &  !< Relative fraction of tile [/]
      & fr_sbl_rel(nc),            &  !< Relative fraction of sibling [/]
      & sum_soil_tle(nc),          &  !< Total soil depth on tile [m]
      & sum_soil_sbl(nc),          &  !< Total soil depth on sibling [m]
      & h,                         &  !< Height of current layer [m]
      & a_inter,                   &  !< Interface area on current layer [m2 m-2]
      & k_inter,                   &  !< Hydraulic conductivity at interface [m s-1]
      & d_mpot,                    &  !< Difference in matric potential [m]
      & wtr_rel_qlbrm,             &  !< Relative saturation at equilibrium between connected tiles [/]
      & free_space,                &  !< Free pore space on tile on a given layer [m]
      & wtr_avail,                 &  !< Available water on sibling on a given layer [m]
      & wtr_2_qlbrm,               &  !< Amount of water moved in order to reach equilibrium state [m]
      & wtr_flx_lat,               &  !< potential flux (based on k and mpot) between sibling and tile [m/s]
      & wtr_lat,                   &  !< actual flux between sibling and tile [m]
      & wtr_soil_tle_rgrd,         &  !< Soil water content on tile on helper grid layer [m]
      & wsat_tle_rgrd,             &  !< Soil moisture at saturation on tile on helper grid layer [m]
      & field_cap_tle_rgrd,        &  !< Field capacity on tile on helper grid layer [m]
      & res_tle_rgrd,              &  !< Residual soil moisture on tile on helper grid layer [m]
      & mpot_tle_rgrd,             &  !< Matric potential on tile on helper grid layer [m]
      & k_tle_rgrd,                &  !< Hydraulic conductivity on tile on helper grid layer [[m s-1]
      & wtr_soil_sbl_rgrd,         &  !< Soil water content on sibling on helper grid layer [m]
      & wsat_sbl_rgrd,             &  !< Soil moisture at saturation on sibling on helper grid layer [m]
      & field_cap_sbl_rgrd,        &  !< Field capacity on sibling on helper grid layer [m]
      & res_sbl_rgrd,              &  !< Residual soil moisture on sibling on helper grid layer [m]
      & mpot_sbl_rgrd,             &  !< Matric potential on deperssion tile on helper grid layer [m]
      & k_sbl_rgrd,                &  !< Hydraulic conductivity on sibling on helper grid layer [[m s-1]
      & h_cond_inter,              &  !< Heat conductivity at interface on current layer on helper grid [J m-1 s-1 K-1]
      & d_t,                       &  !< Temperature difference on current layer on helper grid [K]
      & heat_flx_lat,              &  !< Heat flux between sibling and tile [J m-2 s-1]
      & t_eq,                      &  !< Equilibrium temperature [K]
      & heat_flx_eq_lat               !< Equilirium energy flux between sibling and tile [J m-2 s-1]
    LOGICAL                 ::     &
      & l_compute_water(nc),       &  !< T: water fluxes should be computed for a given cell
      & l_compute_heat(nc)            !< T: heat fluxes should be computed for a given cell
    REAL(wp)                ::     &
      & wtr_diff_tle(nc,nsoil),    &  !< Net diffusive flux on tile on soil layers [m]
      & wtr_diff_sbl(nc,nsoil)        !< Net diffusive flux on sibling on soil layers [m]

    ! Reset lateral flux accumulators
    wtr_lat_in_sl_tle(:,:) = wtr_lat_in_sl_tle(:,:) - wtr_diff_in_tle(:,:)
    wtr_lat_in_sl_sbl(:,:) = wtr_lat_in_sl_sbl(:,:) - wtr_diff_in_sbl(:,:)

    ! Initilize fields for accumulating regridded output
    wtr_diff_tle(:,:) = 0._wp
    wtr_diff_sbl(:,:) = 0._wp

    ! Compute bedrock depth
    sum_soil_tle(:) = SUM(dsoil_tle(:,:), DIM=2)
    sum_soil_sbl(:) = SUM(dsoil_sbl(:,:), DIM=2)

    ! Loop over chells and check basic conditions for lateral fluxes
    DO ic = 1, nc
      IF (tle_fract(ic) < 1.0e-10_wp .OR. sbl_fract(ic) < 1.0e-10_wp .OR. .NOT. connect(ic) > 0._wp) THEN
        ! No lateral fluxes possible under these conditions
        l_compute_water(ic) = .FALSE.
        l_compute_heat(ic)  = .FALSE.
      ELSE
        ! Basic conditions for lateral transport are met
        l_compute_water(ic) = .TRUE.
        l_compute_heat(ic)  = .TRUE.
        ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
        fr_tle_rel(ic) = tle_fract(ic) / (tle_fract(ic) + sbl_fract(ic))
        fr_sbl_rel(ic) = sbl_fract(ic) / (tle_fract(ic) + sbl_fract(ic))
      END IF

      ! Test special conditions for water and heat fluxes
      IF ((h_offset(ic) > 0._wp .AND. h_offset(ic) >  sum_soil_tle(ic))  .OR. &
        & (h_offset(ic) < 0._wp .AND. h_offset(ic) < -sum_soil_sbl(ic))  .OR. &
        & .NOT. l_diff_water) THEN
        ! No lateral fluxes computed under these conditions
        l_compute_water(ic) = .FALSE.
      END IF
      IF ((abs(h_offset(ic)) > dz(nsoil))  .OR. &
        & .NOT. l_diff_heat) THEN
        ! No lateral fluxes computed under these conditions
        l_compute_heat(ic) = .FALSE.
      END IF

      ! TODO: to run on GPUs this routine later needs to use arrays instead of scalars
      IF (l_compute_water(ic) .OR. l_compute_heat(ic)) THEN
        ! Get vertical grid information for intersecting soil layers
        ! Currently assumes soil depth and layering is equal on all tiles
        ! Pond depth is set to zero, as no pond interaction is needed for soil diffusion
        CALL get_intersected_grid(dz(:), nsoil, nrgrd, h_offset(ic), sum_soil_tle(ic), 0._wp, &
          & index_tile(ic,:), index_sibl(ic,:), layer_thickness(ic,:), contype(ic,:))
      END IF
    END DO

    ! Loop over intersected grid layers and grid cells, and compute water and energy
    ! diffusion between tiles on regridded soil layers
    DO i = 1, nrgrd
      DO ic = 1, nc

        ! Get layer index and thickness on tile and sibling
        i_tle = index_tile(ic,i)
        i_sbl = index_sibl(ic,i)
        h = layer_thickness(ic,i)

        ! Compute lateral diffusion of water within bedrock layers
        IF (l_compute_water(ic) .AND. h > 0._wp .AND. contype(ic,i) == "Soil") THEN
          ! Space for liquid water must exists in both tiles
          IF (wsat_tle(ic,i_tle) > 0._wp .AND. wsat_sbl(ic,i_sbl) > 0._wp) THEN

            ! Flux is determined by the difference in matric potential between the two tiles
            wtr_lat            = 0._wp
            a_inter            = len_inter(ic) * h

            wtr_soil_tle_rgrd  = wtr_soil_tle(ic,i_tle)   * h / dsoil_tle(ic,i_tle)
            wsat_tle_rgrd      = wsat_tle(ic,i_tle)       * h / dsoil_tle(ic,i_tle)
            field_cap_tle_rgrd = field_cap_tle(ic,i_tle)  * h / dsoil_tle(ic,i_tle)
            res_tle_rgrd       = res_tle(ic,i_tle)        * h / dsoil_tle(ic,i_tle)
            mpot_tle_rgrd      = mpot_tle(ic,i_tle)
            k_tle_rgrd         = k_tle(ic,i_tle)

            wtr_soil_sbl_rgrd  = wtr_soil_sbl(ic,i_sbl)   * h / dsoil_sbl(ic,i_sbl)
            wsat_sbl_rgrd      = wsat_sbl(ic,i_sbl)       * h / dsoil_sbl(ic,i_sbl)
            field_cap_sbl_rgrd = field_cap_sbl(ic,i_sbl)  * h / dsoil_sbl(ic,i_sbl)
            res_sbl_rgrd       = res_sbl(ic,i_sbl)        * h / dsoil_sbl(ic,i_sbl)
            mpot_sbl_rgrd      = mpot_sbl(ic,i_sbl)
            k_sbl_rgrd         = k_sbl(ic,i_sbl)

            k_inter            = (k_sbl_rgrd * k_tle_rgrd)**0.5_wp
            d_mpot             =  abs(mpot_tle_rgrd) - abs(mpot_sbl_rgrd)
            wtr_flx_lat        = - k_inter  * d_mpot / length_scale(ic) *  a_inter
            wtr_rel_qlbrm      = MIN(1._wp, (wtr_soil_tle_rgrd * fr_tle_rel(ic)  &
              &                +  wtr_soil_sbl_rgrd * fr_sbl_rel(ic))            &
              &                / (wsat_tle_rgrd * fr_tle_rel(ic) + wsat_sbl_rgrd * fr_sbl_rel(ic) ))
            IF (wtr_flx_lat > 0._wp) THEN
              ! Flux is going from tile to sibling
              wtr_2_qlbrm      = MAX(0._wp, &
                &               (wtr_rel_qlbrm * wsat_sbl_rgrd - wtr_soil_sbl_rgrd) * fr_sbl_rel(ic))
              free_space       = MAX(0._wp,(field_cap_sbl_rgrd - wtr_soil_sbl_rgrd) * fr_sbl_rel(ic))
              wtr_avail        = MAX(0._wp,(wtr_soil_tle_rgrd - res_tle_rgrd) * fr_tle_rel(ic))
              wtr_lat          = MIN(wtr_flx_lat * dtime, free_space, wtr_avail, wtr_2_qlbrm)
              wtr_lat          = MAX(0._wp, wtr_lat)

            ELSE
              ! Flux is going from sibling to tile
              wtr_2_qlbrm      = MAX(0._wp, &
                &               (wtr_rel_qlbrm  * wsat_tle_rgrd - wtr_soil_tle_rgrd) * fr_tle_rel(ic))
              free_space       = MAX(0._wp, (field_cap_tle_rgrd - wtr_soil_tle_rgrd) * fr_tle_rel(ic))
              wtr_avail        = MAX(0._wp, (wtr_soil_sbl_rgrd - res_sbl_rgrd) * fr_sbl_rel(ic))
              wtr_lat          = MAX(wtr_flx_lat * dtime, -free_space, -wtr_avail, -wtr_2_qlbrm)
              wtr_lat          = MIN(0._wp, wtr_lat)

            END IF
            wtr_diff_tle(ic,i_tle)  = wtr_diff_tle(ic,i_tle) - wtr_lat / fr_tle_rel(ic)
            wtr_diff_sbl(ic,i_sbl)  = wtr_diff_sbl(ic,i_sbl) + wtr_lat / fr_sbl_rel(ic)
          END IF
        END IF

        ! Compute lateral diffusion of heat within soil column
        IF (l_compute_heat(ic) .AND. h > 0._wp .AND. (contype(ic,i) == "Soil" .OR. contype(ic,i) == "Rock")) THEN
          a_inter              = len_inter(ic) * h

          ! Flux is determined by the difference in temperature between the two tiles
          h_cond_inter         =(heat_cond_soil_tle(ic,i_tle) * heat_cond_soil_sbl(ic,i_sbl))**0.5_wp
          d_t                  = t_soil_tle(ic,i_tle) - t_soil_sbl(ic,i_sbl)
          heat_flx_lat         = h_cond_inter * d_t / length_scale(ic) * a_inter

          ! Limiting heat flux to avoid instablilities for long time steps
          ! Determine flux required to equilibriate temperature within a timestep
          ! Equilibrium temperature
          t_eq = (t_soil_tle(ic,i_tle) * heat_cap_soil_tle(ic,i_tle) * fr_tle_rel(ic)  &
            &  +  t_soil_sbl(ic,i_sbl) * heat_cap_soil_sbl(ic,i_sbl) * fr_sbl_rel(ic)) &
            &  / (heat_cap_soil_tle(ic,i_tle) * fr_tle_rel(ic) + heat_cap_soil_sbl(ic,i_sbl) * fr_sbl_rel(ic))
          ! Heat flux to reach equilibrium temp
          heat_flx_eq_lat = MIN( &
            &               ABS((t_soil_tle(ic,i_tle) - t_eq)  * heat_cap_soil_tle(ic,i_tle) * fr_tle_rel(ic)), &
            &               ABS((t_soil_sbl(ic,i_sbl) - t_eq)  * heat_cap_soil_sbl(ic,i_sbl) * fr_sbl_rel(ic))) &
            &               *  h / dtime
          ! Limiting heat flux to equilibrium flux
          IF (heat_flx_lat > 0._wp .AND. heat_flx_lat > ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = ABS(heat_flx_eq_lat)
          ELSE IF (heat_flx_lat < 0._wp .AND. heat_flx_lat < -ABS(heat_flx_eq_lat)) THEN
            heat_flx_lat = -ABS(heat_flx_eq_lat)
          END IF

          ! Diagnosing incoming-, outgoing- and net heat fluxes on tiles
          nrg_diff_in_tle(ic,i_tle)  = nrg_diff_in_tle(ic,i_tle)  + MAX(0._wp, -heat_flx_lat / fr_tle_rel(ic))
          nrg_diff_in_sbl(ic,i_sbl)  = nrg_diff_in_sbl(ic,i_sbl)  + MAX(0._wp,  heat_flx_lat / fr_sbl_rel(ic))
          nrg_diff_out_tle(ic,i_tle) = nrg_diff_out_tle(ic,i_tle) - MIN(0._wp, -heat_flx_lat / fr_tle_rel(ic))
          nrg_diff_out_sbl(ic,i_sbl) = nrg_diff_out_sbl(ic,i_sbl) - MIN(0._wp,  heat_flx_lat / fr_sbl_rel(ic))

          heat_flx_lat_tle(ic,i_tle)   = heat_flx_lat_tle(ic,i_tle)   - heat_flx_lat / fr_tle_rel(ic)
          heat_flx_lat_sbl(ic,i_sbl)   = heat_flx_lat_sbl(ic,i_sbl)   + heat_flx_lat / fr_sbl_rel(ic)
        END IF
      END DO
    END DO

    ! Adjusting the water contents and calculating the net- in- and outgoing water fluxes for the tiles
    wtr_soil_tle(:,:)     = wtr_soil_tle(:,:)     + wtr_diff_tle(:,:)
    wtr_soil_sbl(:,:)     = wtr_soil_sbl(:,:)     + wtr_diff_sbl(:,:)
    wtr_latflow_tle(:,:)  = wtr_latflow_tle(:,:)  + wtr_diff_tle(:,:)  * rhoh2o / dtime
    wtr_latflow_sbl(:,:)  = wtr_latflow_sbl(:,:)  + wtr_diff_sbl(:,:)  * rhoh2o / dtime
    wtr_diff_in_tle(:,:)  = wtr_diff_in_tle(:,:)  + MAX(0._wp, wtr_diff_tle(:,:)) * rhoh2o / dtime
    wtr_diff_in_sbl(:,:)  = wtr_diff_in_sbl(:,:)  + MAX(0._wp, wtr_diff_sbl(:,:)) * rhoh2o / dtime
    wtr_diff_out_tle(:,:) = wtr_diff_out_tle(:,:) - MIN(0._wp, wtr_diff_tle(:,:)) * rhoh2o / dtime
    wtr_diff_out_sbl(:,:) = wtr_diff_out_sbl(:,:) - MIN(0._wp, wtr_diff_sbl(:,:)) * rhoh2o / dtime

    ! Adding updated water diffusion fluxes to lateral fluxes
    wtr_lat_in_sl_tle(:,:)   = wtr_lat_in_sl_tle(:,:) + wtr_diff_in_tle(:,:)
    wtr_lat_in_sl_sbl(:,:)   = wtr_lat_in_sl_sbl(:,:) + wtr_diff_in_sbl(:,:)

  END SUBROUTINE calc_lat_fluxes_soil

  ! =================================================================================================================
  !>
  !>#### Compute lateral heat advection due to surface and subsurface water flows
  !>
  !> This routine compute the lateral advection of heat due to water exchange between the soil
  !> columns of two adjacent tiles. This needs to take into account that the soil layers of the
  !> two tiles may be offset against each other, and that the soil layers may not align. In this
  !> case, an intersected grid is created to compute the lateral fluxes between the overlapping
  !> parts of the two soil column layers.
  !>
  SUBROUTINE calc_convection_heat(                         &
      &      nc, nsoil, nrgrd, dtime,                      &
      &      dz, h_offset,                                 &
      &      tle_fract, sbl_fract,                         &
      &      wtr_inflx_lat_srf_tle, wtr_inflx_lat_srf_sbl, &
      &      t_soil_tle, t_soil_sbl,                       &
      &      wtr_pond_tle, wtr_pond_sbl,                   &
      &      wtr_inflx_lat_sl_tle, wtr_inflx_lat_sl_sbl,   &
      &      dsoil_tle, dsoil_sbl,                         &
      &      heat_cap_soil_tle, heat_cap_soil_sbl,         &
      &      cnvctv_in_tle, cnvctv_in_sbl,                 &
      &      cnvctv_out_tle, cnvctv_out_sbl,               &
      &      heat_flx_lat_tle, heat_flx_lat_sbl            &
      &      )

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o
    USE mo_sse_constants, ONLY: vol_hcap_water

    INTEGER, INTENT(in)     ::     &
      & nc,                        &  !< Number of grid cells
      & nsoil,                     &  !< Number of oil layers
      & nrgrd                         !< Number of layers on helper grid determined by layer overlaps
    REAL(wp), INTENT(in)    :: dtime  !< Timestep length [s]
    REAL(wp), INTENT(in)    ::     &
      & h_offset(:),               &  !< Difference in surface elevation between tile and sibling [m]
      & tle_fract(:),              &  !< Area fraction of tile [/]
      & sbl_fract(:),              &  !< Area fraction of sibling [/]
      & wtr_inflx_lat_srf_tle(:),  &  !< Lateral flux (surface) of water from sibling to tile [kg m-2 s-1]
      & wtr_inflx_lat_srf_sbl(:),  &  !< Lateral flux (surface) of water from tile to sibling [kg m-2 s-1]
      & wtr_pond_tle(:),           &  !< Ponded water on tile [m]
      & wtr_pond_sbl(:),           &  !< Ponded water on sibling [m]
      & dz(:)                         !< Layer thickness [m]
    REAL(wp), INTENT(in)    ::     &
      & t_soil_tle(:,:),           &  !< Soil tempertures of tile [K]
      & t_soil_sbl(:,:),           &  !< Soil tempertures of sibling [K]
      & wtr_inflx_lat_sl_tle(:,:), &  !< Lateral flux (soil) of water from sibling to tile [kg m-2 s-1]
      & wtr_inflx_lat_sl_sbl(:,:), &  !< Lateral flux (soil) of water from tile to sibling [kg m-2 s-1]
      & dsoil_tle(:,:),            &  !< Soil depth on layers of tile [m]
      & dsoil_sbl(:,:),            &  !< Soil depth on layers of sibling [m]
      & heat_cap_soil_tle(:,:),    &  !< Soil heat capacity of tile [J m-3 K-1]
      & heat_cap_soil_sbl(:,:)        !< Soil heat capacity of sibling [J m-3 K-1]

    REAL(wp), INTENT(inout) ::     &
      & cnvctv_in_tle(:,:),        &  !< Advective influx at tile [J m-2 s-1]
      & cnvctv_in_sbl(:,:),        &  !< Advective outflux at sibling [J m-2 s-1]
      & cnvctv_out_tle(:,:),       &  !< Advective influx at tile [J m-2 s-1]
      & cnvctv_out_sbl(:,:)           !< Advective outflux at sibling [J m-2 s-1]

    REAL(wp), INTENT(inout) ::     &
      & heat_flx_lat_tle(:,:),     &  !< Net Heat flux at tile [J m-2 s-1]
      & heat_flx_lat_sbl(:,:)         !< Net Heat flux at sibling [J m-2 s-1]

    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 ::     &
      & ic,                        &  !< Index variable for grid cells
      & is,                        &  !< Index variable for soil layers
      & i,                         &  !< Index variable for layers on helper grid
      & i_tle,                     &  !< Index variables for layers on tile
      & i_sbl                         !< Index variables for layers on sibling tile
    INTEGER                 ::     &
      & index_tile(nc,nrgrd),      &  !< Mapping index from helper grid to tile layers
      & index_sibl(nc,nrgrd)          !< Mapping index from helper grid to sibling layers
    REAL(wp)                :: layer_thickness(nc,nrgrd)  !< Thickness of layers on helper grid [m]
    CHARACTER(LEN=4)        :: contype(nc,nrgrd)          !< Type of connection between tile and sibling
    LOGICAL                 ::     &
      & l_compute(nc)                 !< Flag to indicate whether fluxes should be computed for a given cell
    REAL(wp)                ::     &
      & dz_bedrock(nc),            &  !< Depth until bedrock (assumed to be equal on tile and sibling) [m]
      & dz_pond(nc)                   !< Pond depth [m]
    REAL(wp)                :: dz_lb(nsoil)               !< Lower depth bounds of soil layers [m]
    REAL(wp)                ::     &
      & dz_sum,                    &  !< Total soil column depth [m]
      & z_tle,                     &  !< Layer bottom of current layer on tile [m]
      & z_sbl,                     &  !< Layer bottom of current layer on sibling [m]
      & z_rgrd_top,                &  !< Top of current layer on helper grid [m]
      & z_rgrd_bot,                &  !< Bottom of current layer on helper grid [m]
      & h,                         &  !< Height of current layer [m]
      & t_tle_rgrd,                &  !< Soil tempertures of tile on current layer on helper grid [K]
      & t_sbl_rgrd,                &  !< Soil tempertures of sibling on current layer on helper grid [K]
      & d_t_wtr_soil_tle,          &  !< Temperature difference between incomming water and soil on tile on current
                                      !< layer on helper grid [K]
      & d_t_wtr_soil_sbl,          &  !< Temperature difference between incomming water and soil on sibling on current
                                      !< layer on helper grid [K]
      & h_cap_soil_tle_rgrd,       &  !< Soil heat capacity of tile on current layer on helper grid [J m-3 K-1]
      & h_cap_soil_sbl_rgrd,       &  !< Soil heat capacity of sibling on current layer on helper grid [J m-3 K-1]
      & wtr_flx_in_tle_regrid,     &  !< Incomming water on tile (from sibling) on current layer on helper grid [m3 m-2 s]
      & wtr_flx_in_sbl_regrid,     &  !< Incomming water on sibling (from tile) on current layer on helper grid [m3 m-2 s]
      & heat_flx_tle,              &  !< Incomming heat flux on tile on current layer on helper grid [J m-3]
      & heat_flx_sbl                  !< Incomming heat flux on sibling on current layer on helper grid [J m-3]

    ! Initialize 1D fields
    dz_sum = 0._wp
    dz_lb(:) = 0._wp
    DO is = 1, nsoil
      dz_sum = dz_sum + dz(is)
      IF (is == 1) THEN
        dz_lb(is) = dz(is)
      ELSE
        dz_lb(is) = dz_lb(is-1) + dz(is)
      END IF
    END DO

    ! Compute bedrock depth (assumed equal on tile and sibling)
    dz_bedrock(:) = SUM(dsoil_tle(:,:), DIM=2)

    ! Loop over grid cells and idenify grid cells where fluxes should be computed
    DO ic = 1, nc
      IF (abs(h_offset(ic)) > dz_sum .OR. &
        & tle_fract(ic) < 1.0e-10_wp .OR. sbl_fract(ic) < 1.0e-10_wp) THEN
        l_compute(ic) = .FALSE.
        dz_pond(ic) = 0._wp
      ELSE
        l_compute(ic) = .TRUE.
        dz_pond(ic) = MAX(wtr_pond_tle(ic), wtr_pond_sbl(ic))
        ! Get vertical grid information for intersecting soil layers
        ! Currently assumes soil depth and layering is equal on all tiles
        CALL get_intersected_grid(dz_lb(:), nsoil, nrgrd, h_offset(ic), dz_bedrock(ic), dz_pond(ic), &
          & index_tile(ic,:), index_sibl(ic,:), layer_thickness(ic,:), contype(ic,:))
      END IF
    END DO

    ! Loop over intersected grid layers and grid cells. As ponded water / ice currently does not have a heat content
    ! we only account for heat exchange in the soil:
    DO i = 1, nrgrd
      DO ic = 1, nc
        IF (l_compute(ic)) THEN

          ! Get indices and layer thickness for actual cell and layer
          i_tle = index_tile(ic,i)
          i_sbl = index_sibl(ic,i)
          h     = layer_thickness(ic,i)

          ! First determine the lateral heat fluxes for water movements within the soil
          IF (contype(ic,i) /= "None") THEN ! Within the pond or soil column of both tiles

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
            ! Incoming water is assumed to have the same temperature as the soil layer it originates from
            ! or a minimum of tmelt to be liquid.

            d_t_wtr_soil_tle = MAX(t_sbl_rgrd, tmelt) - t_tle_rgrd
            d_t_wtr_soil_sbl = MAX(t_tle_rgrd, tmelt) - t_sbl_rgrd

            heat_flx_tle     = d_t_wtr_soil_tle * vol_hcap_water * wtr_flx_in_tle_regrid
            heat_flx_sbl     = d_t_wtr_soil_sbl * vol_hcap_water * wtr_flx_in_sbl_regrid

            ! Soil may not warm above or cool below the temperature of the incomming water
            ! This could happen in the above formulation if wtr_flx_in_tle_regrid * dtime > h
            ! Here, we use the current soil heat capacities instead of assuming that
            ! the water content will change according to wtr_flx * dtime.
            ! This makes the limits more restrictive, given that these are anyway very large margins.
            IF (     heat_flx_tle < 0._wp) THEN
              heat_flx_tle   = MAX(heat_flx_tle, MAX(d_t_wtr_soil_tle, -10._wp) * h_cap_soil_tle_rgrd * h / dtime)
            ELSE IF (heat_flx_tle > 0._wp) THEN
              heat_flx_tle   = MIN(heat_flx_tle, MIN(d_t_wtr_soil_tle,  10._wp) * h_cap_soil_tle_rgrd * h / dtime)
            END IF
            IF (     heat_flx_sbl < 0._wp) THEN
              heat_flx_sbl   = MAX(heat_flx_sbl, MAX(d_t_wtr_soil_sbl, -10._wp) * h_cap_soil_sbl_rgrd * h / dtime)
            ELSE IF (heat_flx_sbl > 0._wp) THEN
              heat_flx_sbl   = MIN(heat_flx_sbl, MIN(d_t_wtr_soil_sbl,  10._wp) * h_cap_soil_sbl_rgrd * h / dtime)
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

          END IF
        END IF
      END DO
    END DO

    ! Add the heat fluxes due to surface runoff / runon
    DO ic = 1, nc
      IF (l_compute(ic)) THEN
        d_t_wtr_soil_tle = MAX(t_soil_sbl(ic,1), tmelt) - t_soil_tle(ic,1)
        d_t_wtr_soil_sbl = MAX(t_soil_tle(ic,1), tmelt) - t_soil_sbl(ic,1)
        heat_flx_tle = d_t_wtr_soil_tle * vol_hcap_water * wtr_inflx_lat_srf_tle(ic) / rhoh2o
        heat_flx_sbl = d_t_wtr_soil_sbl * vol_hcap_water * wtr_inflx_lat_srf_sbl(ic) / rhoh2o

        IF (     heat_flx_tle < 0._wp) THEN
          heat_flx_tle   = MAX(heat_flx_tle, MAX(d_t_wtr_soil_tle, -10._wp) * heat_cap_soil_tle(ic,1) * dz(1) / dtime)
        ELSE IF (heat_flx_tle > 0._wp) THEN
          heat_flx_tle   = MIN(heat_flx_tle, MIN(d_t_wtr_soil_tle,  10._wp) * heat_cap_soil_tle(ic,1) * dz(1) / dtime)
        END IF

        IF (     heat_flx_sbl < 0._wp) THEN
          heat_flx_sbl   = MAX(heat_flx_sbl, MAX(d_t_wtr_soil_sbl, -10._wp) * heat_cap_soil_sbl(ic,1) * dz(1) / dtime)
        ELSE IF (heat_flx_sbl > 0._wp) THEN
          heat_flx_sbl   = MIN(heat_flx_sbl, MIN(d_t_wtr_soil_sbl,  10._wp) * heat_cap_soil_sbl(ic,1) * dz(1) / dtime)
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

  ! =================================================================================================================
  !>
  !>#### Compute lateral resdistribution of snow between adjacent tiles
  !>
  !> This routine compute the lateral transport of snow from elevated to depression tiles. Here, only
  !> this simple redistribution is considered, while more complex processes such as wind-driven, vegetation hampered
  !> transport or complex topography effects are not included.
  !>
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

    USE mo_jsb_physical_constants, ONLY: rhoh2o, dens_snow_min, dens_snow
    USE mo_hydro_constants,        ONLY: dens_snow_max

    LOGICAL,  INTENT(in)    :: &
      & l_dynsnow                 !< Dynamic snow density is being used
    REAL(wp), INTENT(in)    :: &
      & dtime,                 &  !< Timestep length [s]
      & cnct_rel,              &  !< connectivity of tile and sibing relative to the sum of connectivities of tile [/]
      & h_offset,              &  !< Differece in surface elevation between tile and sibling [m]
      & tle_fract,             &  !< Area fraction of tile [/]
      & sbl_fract                 !< Area fraction of sibling [/]
    REAL(wp), INTENT(in)    :: &
      & weq_pond_tle,          &  !< Pond storage on tile [m]
      & weq_pond_sbl              !< Pond storage on sibling [m]
    REAL(wp), INTENT(in)    :: &
      & snow_soil_tle_tm1,     &  !< Snow soil (WEQ) on tile at end of previous timestep [m]
      & snow_new_tle_tm1,      &  !< New snow soil (WEQ) on tile at end of previous timestep [m]
      & snow_dens_tle_tm1         !< Snow density on tile at end of previous timestep [kg m-3]
    REAL(wp), INTENT(inout) :: &
      & snow_rdst_in_tle,      &  !< Snow redistribution influx at tile [kg m-2 s-1]
      & snow_rdst_in_sbl,      &  !< Snow redistribution influx at sibling [kg m-2 s-1]
      & snow_rdst_out_tle,     &  !< Snow redistribution outflux at tile [kg m-2 s-1]
      & snow_rdst_out_sbl         !< Snow redistribution outflux at sibling [kg m-2 s-1]
    REAL(wp), INTENT(inout) :: &
      & wtr_latflow_tle,       &  !< Netflux (surface) on tile [kg m-2 s-1]
      & wtr_latflow_sbl,       &  !< Netflux (surface) on sibling [kg m-2 s-1]
      & snow_soil_tle,         &  !< Snow soil (WEQ) on tile [m]
      & snow_soil_sbl,         &  !< Snow soil (WEQ) on sibling [m]
      & snow_new_tle,          &  !< New snow soil (WEQ) on tile [m]
      & snow_new_sbl,          &  !< New snow soil (WEQ) on sibling [m]
      & snow_tle,              &  !< Snow (WEQ; can and soil) on tile [m]
      & snow_sbl,              &  !< Snow (WEQ; can and soil) on sibling [m]
      & snow_dens_tle,         &  !< Snow density on tile [kg m-3]
      & snow_dens_sbl             !< Snow density on sibling [kg m-3]

    ! Local variables ------------------------------------------------------------------------------------
    REAL(wp)                :: &
      & fr_tle_rel,          &  !< Relative fraction of tile [/]
      & fr_sbl_rel,          &  !< Relative fraction of sibling [/]
      & h_offset_snow_fix,   &  !< Difference in elevation incl. immobile snow layer between tile and sibling [m]
      & snow_fix_tle ,       &  !< Immobile snow on tile [m]
      & snow_fix_sbl ,       &  !< Immobile snow on sibling [m]
      & snow_fix_rel_tle,    &  !< Fraction of immobile snow on tile [/]
      & snow_fix_rel_sbl,    &  !< Fraction of immobile snow on sibling [/]
      & snow_depth_mob_avg,  &  !< Average depth of mobile snow [m]
      & snow_depth_mob_tle,  &  !< Depth of mobile snow on tile [m]
      & snow_depth_mob_sbl,  &  !< Depth of mobile snow on sibling [m]
      & snow_mob_tle,        &  !< Mobile snow (WEQ) on tile [m]
      & snow_mob_sbl,        &  !< Mobile snow (WEQ)  on sibling [m]
      & snow_mob_upd_tle,    &  !< Updated mobile snow (WEQ) on tile [m]
      & snow_mob_upd_sbl,    &  !< Updated mobile snow (WEQ) on sibling [m]
      & snow_dens_max,       &  !< Maximum snow density [kg m-3]
      & snow_dens_min,       &  !< Minimum snow density [kg m-3]
      & snow_dens_fix_tle,   &  !< Density of immobile snow on tile [kg m-3]
      & snow_dens_fix_sbl,   &  !< Density of immobile snow on sibling [kg m-3]
      & snow_rdst_tle,       &  !< Snow (WEQ) redistribution flux on tile [kg m-2 s-1]
      & snow_rdst_sbl           !< Snow (WEQ) redistribution flux on sibling [kg m-2 s-1]

    ! Return if area fractions are too small
    IF (tle_fract < 1.0e-10_wp .OR. sbl_fract < 1.0e-10_wp) RETURN

    ! Setting initial values to 0.0
    snow_rdst_tle = 0._wp
    snow_rdst_sbl = 0._wp

    ! Snow density parameters depending on snow scheme
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
      ! TODO: Needs to be updated according to Mittals new snow density scheme
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

  ! =================================================================================================================
  !>
  !>#### Remove excess water from soil water and pond storages
  !>
  !> This routine make sure that neither soil moisture nor the pond storage exceed their maximum values
  !> after lateral fluxes are applied.
  !>
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
    INTEGER , INTENT(in)    :: nsoil                          !< Number of soil layers
    REAL(wp), INTENT(in)    :: weq_pond_max(:),  &            !< Maximum pond water storage [m]
      &                        ice_pond(:)                    !< Ice in pond water storage [m]

    REAL(wp), INTENT(inout) :: wtr_pond(:),      &            !< Liquid water in pond water storage [m]
      &                        wtr_latflow_res_srf(:)         !< Surface water retention storage [m]

    REAL(wp), OPTIONAL, INTENT(in)    :: wtr_soil_sat_sl(:,:) !< Saturation soil moisture on soil layers [m]
    REAL(wp), OPTIONAL, INTENT(inout) :: wtr_soil_sl(:,:)     !< Water content on soil layers [m]]

    INTEGER                 :: ic,               &            !< Loop index for cells
      &                        is                             !< Loop index for soil layers
    REAL(wp)                :: wtr_pond_max                   !< Maximum liquid water in pond storage [m]

    ! Loop over soil layers from bottom to top
    IF (nsoil > 0) THEN
      ! Push excess water from lower to upper soil layers
      DO is = nsoil, 2, -1
        DO ic = 1, nc
          IF (wtr_soil_sl(ic,is) > wtr_soil_sat_sl(ic,is)) THEN
            wtr_soil_sl(ic,is-1) = wtr_soil_sl(ic,is-1) &
              &                  + wtr_soil_sl(ic,is) - wtr_soil_sat_sl(ic,is)
            wtr_soil_sl(ic,is)   = wtr_soil_sat_sl(ic,is)
          END IF
        END DO
      END DO

      ! Excess water in top soil layer goes into pond storage
      DO ic = 1, nc
        IF (wtr_soil_sl(ic,1) > wtr_soil_sat_sl(ic,1)) THEN
          wtr_pond(ic)      = wtr_pond(ic) &
            &               + wtr_soil_sl(ic,1) - wtr_soil_sat_sl(ic,1)
          wtr_soil_sl(ic,1) = wtr_soil_sat_sl(ic,1)
        END IF
      END DO
    END IF

    DO ic = 1, nc
      ! Excess water in pond storage goes into retention storage
      wtr_pond_max = MAX(0.0_wp, weq_pond_max(ic) - ice_pond(ic))
      IF (wtr_pond(ic) > wtr_pond_max) THEN
        wtr_latflow_res_srf(ic) = wtr_latflow_res_srf(ic) &
          &                     + wtr_pond(ic) - wtr_pond_max
        wtr_pond(ic)            = wtr_pond_max
      END IF
    END DO

  END SUBROUTINE remove_excess_water

  ! =================================================================================================================
  !>
  !>#### Adjust lateral fluxes to compensate for elevation differences between tiles
  !>
  !> Subsurface drainage is assumed to flow horizontally from the elevation tile into the depression tile
  !> without following the topography. Therefore, water from a lower layer on the elevation tile may end up on
  !> a heighter layer in the depression tile if there is a positive elevation offset between the two tiles.
  !> Futhermore, the offset does not corresponds to the layer thickness on the soil grid, therefore interpolation
  !> between layers is necessary.
  !>
  SUBROUTINE shift_fluxes_sl(                 &
    & nc, nsoil, nrgrd,                       &
    & connect, h_offset, dz,                  &
    & dsoil_lvt, dsoil_dpr,                   &
    & flx_lat_sl                              &
    )

    INTEGER , INTENT(in)    ::     &
      & nc,                        &  !< Number of cells
      & nsoil,                     &  !< Number of soil ayers
      & nrgrd                         !< Number of layers on helper grid determined by layer overlaps
    REAL(wp), INTENT(in)    ::     &
      & connect(:),                &  !< Connectivity array for lateral exchange
      & h_offset(:)                   !< Difference in surface elevation between elevation and depression [m]
    REAL(wp), INTENT(in)    ::     &
      & dz(:)                         !< Standard soil layer bottom depths in both tiles [m]
    REAL(wp), INTENT(in)    ::     &
      & dsoil_lvt(:,:),            &  !< Soil depth on layers of elevation tile [m]
      & dsoil_dpr(:,:)                !< Soil depth on layers of depression tile [m]
    REAL(wp), INTENT(inout) ::     &
      & flx_lat_sl(:,:)               !< Lateralal flux to be adjusted [?]
    INTEGER                 ::     &
      & i,                         &  !< Index variable for intersected grid layers
      & ic,                        &  !< Index variable for grid cells
      & is,                        &  !< Index variable for soil layers
      & i_lvt,                     &  !< Index variable for layers on elevation tile
      & i_dpr                         !< Index variable for layers on depression tile
    REAL(wp)                ::     &
      & dz_bedrock(nc)                !< Depth of bedrock boundary [m]
    REAL(wp)                ::     &
      & flx_lat_sl_upd(nc,nsoil)      !< Updated lateral flux [?]
    LOGICAL                 ::     &
      & l_compute(nc)                 !< Flag to indicate whether adjusting the fluxes is necessary
    INTEGER                 ::     &
      & index_tile(nc,nrgrd),      &  !< Mapping index from helper grid to tile layers
      & index_sibl(nc,nrgrd)          !< Mapping index from helper grid to sibling layers
    REAL(wp)                ::     &
      & layer_thickness(nc,nrgrd), &  !< Thickness of layers on helper grid [m]
      & layer_bot(nc,nrgrd)           !< Bottom depth of soil layers on helper grid [m]
    CHARACTER(LEN=4)        ::     &
      & contype(nc,nrgrd)             !< Type of connection between tile and sibling

    ! Initialize fields
    dz_bedrock(:) = SUM(dsoil_dpr(:,:), DIM=2)

    ! Determine whether adjustment of fluxes is necessary
    DO ic = 1, nc
      IF (connect(ic) > 0._wp .AND. h_offset(ic) > 0._wp)THEN
        l_compute(ic) = .TRUE.
        flx_lat_sl_upd(ic,:) = 0._wp
        ! Get intersected grid between elevation and depression tile. Pond depth is set to zero here
        ! since fluxes are only shifted within the soil
        ! TODO: Any reason, why lateral drainage should not discharge into ponds? At least from those
        !       layers, that are above the pond water level?
        CALL get_intersected_grid(dz(:), nsoil, nrgrd, h_offset(ic), dz_bedrock(ic), 0._wp, &
          & index_tile(ic,:), index_sibl(ic,:), layer_thickness(ic,:), contype(ic,:),       &
          & layer_bot=layer_bot(ic,:))
      ELSE
        ! No adjustment necessary, either because there is no positive offset, or no connectivity
        l_compute(ic) = .FALSE.
        flx_lat_sl_upd(ic,:) = flx_lat_sl(ic,:)
      END IF
    END DO

    ! Remap lateral drainage onto depression tile
    DO i = 1, nrgrd
      DO ic = 1, nc
        IF (l_compute(ic)) THEN
          i_lvt = index_tile(ic,i)
          i_dpr = MAX(index_sibl(ic,i), 1)
          ! Map fluxes on soil grid of depression tile while soil depth is above bedrock
          IF (layer_bot(ic, i) <= dz_bedrock(ic) .AND. i_lvt > 0) THEN
            flx_lat_sl_upd(ic,i_dpr) = flx_lat_sl_upd(ic,i_dpr) &
              & + flx_lat_sl(ic,i_lvt) * layer_thickness(ic,i) / dsoil_lvt(ic,i_lvt)
          END IF
        END IF
      END DO
    END DO

    ! Without offset or connectivity, just check no incoming flux is below the bedrock boundary
    DO is = nsoil, 2, -1
      DO ic = 1, nc
        IF (dsoil_dpr(ic,is) <= 0._wp .AND. .NOT. l_compute(ic)) THEN
          flx_lat_sl_upd(ic,is-1) = flx_lat_sl_upd(ic,is-1) + flx_lat_sl_upd(ic,is)
          flx_lat_sl_upd(ic,is)   = 0._wp
        END IF
      END DO
    END DO

    ! Update lateral fluxes with mapped values
    flx_lat_sl(:,:) = flx_lat_sl_upd(:,:)

  END SUBROUTINE shift_fluxes_sl

  ! =================================================================================================================
  !>
  !>#### Compute geometric relationships for tile and connected sibling
  !>
  !> Lateral transport and diffusion require assumptions about the average distance between two adjacent tiles
  !> as well as the length of their interface. This routine computes these geometric relationships.
  !>
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

    REAL(wp), INTENT(in)  :: &
      & microcon_ctc_tile,   &  !< Micro-scale connectivity between tile and sibling (if tile is encompased by sibl) []
      & microcon_ctc_sibl       !< Micro-scale connectivity between sibling and tile (if sibl is encompased by tile) []
    REAL(wp), INTENT(in)  :: &
      & fract,               &  !< Cover fraction of tile []
      & fract_sibl,          &  !< Cover fraction of sibling []
      & fract_nclsd,         &  !< Cover fraction enclosed by tile []
      & fract_nclsd_sibl,    &  !< Cover fraction enclosed by sibling []
      & area,                &  !< Characteristic area of tile [m2]
      & area_sibl,           &  !< Characteristic area of sibling [m2]
      & area_nclsd,          &  !< Characteristic area enclosed by tile [m2]
      & area_nclsd_sibl,     &  !< Characteristic area enclosed by sibling [m2]
      & h_surf,              &  !< Surface elevation of tile []
      & h_surf_sibl             !< Surface elevation of sibling []
    REAL(wp), INTENT(out) :: &
      & h_offset,            &  !< Difference in surface elevation between tile and sibling [m]
      & length_scale,        &  !< Characteristic length scale of gradient between two adjacent tiles [m]
      & length_iface,        &  !< Length of interface between tile and sibling [m m-2]
      & micro_slope             !< Micro-scale slope between tile and sibling []

    ! Local variables ------------------------------------------------------------------------------------
    REAL(wp)              :: &
      & fract_rel,           &  !< Relative fraction of tile []
      & fract_rel_sibl,      &  !< Relative fraction of sibling []
      & rad_in,              &  !< Radius of inner circle [m]
      & rad_out,             &  !< Radius outer circle [m]
      & area_cmbnd,          &  !< Combined area of nested tiles [m2]
      & fract_rel_cmbnd,     &  !< Relative fraction of nested tiles []
      & area_fact               !< Conversion factor to get from area of inner to that of outer tile

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
    ! The smaller circle is encompassed by a larger circle (C2) which contains the areas of both tiles and whose radius
    ! is used to approximate the characteristic length scale. The latter is assumed to be equal to half the distance
    ! between circumference of C1 and C2 added to the radius of C1. The area of C2 is determined so that a segment
    ! with a central angle (X/360) equal to the connected fraction contains an area equal to that of the larger tile
    ! and the area of a segment of C1 with the same central angle.
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
    ELSE IF (microcon_ctc_tile > microcon_ctc_sibl) THEN
      ! Tile is encompassed by sibling
      IF (area_nclsd > 0._wp .AND. fract_nclsd > 1.0e-10_wp) THEN
        ! Tile encloses another tile
        area_cmbnd      = area_nclsd &
          &             * (fract_nclsd + fract) / fract_nclsd
        fract_rel_cmbnd = (fract_nclsd + fract) &
          &             / (fract_nclsd + fract + fract_sibl)
        area_fact       = 1._wp / fract_rel_cmbnd - 1._wp
        length_iface    = ((4._wp * pi * area_cmbnd) ** 0.5_wp) &
          &             * microcon_ctc_tile * fract_rel / area
        rad_in          = (area_nclsd    / pi)**0.5_wp
        rad_out         = (area_cmbnd * (1._wp + area_fact / microcon_ctc_tile) / pi)**0.5_wp
        length_scale    = (rad_in + rad_out) / 2._wp
      ELSE
        ! Tile does not enclose other tile
        area_fact       = 1._wp / fract_rel - 1._wp
        length_iface    = ((4._wp * pi * area) ** 0.5_wp) &
          &             * microcon_ctc_tile * fract_rel / area
        rad_in          = (area / pi)**0.5_wp
        rad_out         = (area *(1._wp + area_fact / microcon_ctc_tile) / pi)**0.5_wp
        length_scale    = rad_in + (rad_out - rad_in) / 2._wp
      END IF

    ELSE
      ! Sibling is encompassed by tile
      IF (area_nclsd_sibl > 0._wp .AND. fract_nclsd_sibl > 1.0e-10_wp) THEN
        ! Sibling encloses another tile
        area_cmbnd      = area_nclsd_sibl &
          &             * (fract_nclsd_sibl + fract_sibl) / fract_nclsd_sibl
        fract_rel_cmbnd = (fract_nclsd_sibl + fract_sibl) &
          &             / (fract_nclsd_sibl + fract + fract_sibl)
        area_fact       = 1._wp / fract_rel_cmbnd - 1._wp
        length_iface    = ((4._wp * pi * area_cmbnd) ** 0.5_wp) &
          &             * microcon_ctc_sibl * fract_rel_sibl / area_sibl
        rad_in          = (area_nclsd_sibl / pi)**0.5_wp
        rad_out         = (area_cmbnd *(1._wp + area_fact / microcon_ctc_sibl) / pi)**0.5_wp
        length_scale    = (rad_in + rad_out) / 2._wp
      ELSE
        ! Sibling does not enclose other tile
        area_fact       = 1._wp / fract_rel_sibl - 1._wp
        length_iface    = ((4._wp * pi * area_sibl) ** 0.5_wp) &
          &             * microcon_ctc_sibl * fract_rel_sibl / area_sibl
        rad_in          = (area_sibl / pi)**0.5_wp
        rad_out         = (area_sibl *(1._wp + area_fact / microcon_ctc_sibl) / pi)**0.5_wp
        length_scale    = rad_in + (rad_out - rad_in) / 2._wp
      END IF
    END IF

    ! calculate micro-scale slope based on difference in surface elevation and characteristic length scale
    ! ... also in case that microcon_ctc_sibl > since the microcon_slp will be reversed later on ...
    IF ((microcon_ctc_tile > 0._wp .OR. microcon_ctc_sibl > 0._wp) .AND. &
      &  h_offset > 0._wp .AND. length_scale > 0._wp) THEN
      micro_slope = MIN(0.5_wp, h_offset / length_scale)
    ELSE
      micro_slope = 0._wp
    END IF
  END SUBROUTINE get_geometry

  ! =================================================================================================================
  !>
  !>#### Compute lag factors for lateral exchange
  !>
  !> Macro- and micro-scale lag factors are determined by the (inverse) retention time,
  !> the connectivity between tile and sibling, and the cover fraction of a sibling relative to
  !> the total connected fraction of tile. Additionally, the lag foctor is decreased according
  !> to the micro-scale connectivity, assuming that the fractions of the tile that have a micro-scale
  !> connectivity do not additionally feature macro-scale connections.
  ! TODO: this routine needs to be double-checked with Philipp. Difference between micro and macro lag factors
  !       are rather unclear and the comments only talk about micro-scale.
  !>
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

    INTEGER,  INTENT(in)    :: &
      & lex_idt,               &  !< lex id of tile
      & lex_ids                   !< lex id of sibling
    REAL(wp), INTENT(in)    :: &
      & velocity_srf,          &  !< Potential velocity of lateral surface flow along pathways [m s-1]
      & velocity_blg,          &  !< Potential velocity of lateral below-ground flow along flowpaths [m s-1]
      & ret_macro_srf,         &  !< Retention time for macro-scale fluxes at the surface [s]
      & ret_macro_blg             !< Retention time for macro-scale fluxes below the surface [s]
    REAL(wp), INTENT(in)    :: &
      & length_scale,          &  !< Assumed distance between centers of micro-scale connected tiles [m]
      & fract,                 &  !< cover fraction of tile []
      & fract_sibl,            &  !< cover fraction of sibling []
      & micro_connect_srf,     &  !< Micro-scale connectivity (surface) between tile and sibling []
      & micro_connect_blg,     &  !< Micro-scale connectivity (below-ground) between tile and sibling []
      & macro_connect_srf,     &  !< Macro-scale connectivity (surface) between tile and sibling []
      & macro_connect_blg,     &  !< Macro-scale connectivity (below-ground) between tile and sibling []
      & fract_microcon_srf,    &  !< (Weighted) Sum of micro-scale connect. fractions (surface) on tile []
      & fract_microcon_blg,    &  !< (Weighted) Sum of micro-scale connect. fractions (below-ground) on tile []
      & fract_macrocon_srf,    &  !< (Weighted) Sum of macro-scale connect. fractions (surface) on tile []
      & fract_macrocon_blg        !< (Weighted) Sum of macro-scale connect. fractions (below-ground) on tile []
    REAL(wp), INTENT(out)   :: &
      & fact_lag_micro_srf,    &  !< lag factor for micro-scale surface fluxes between tile and sibling [s-1]
      & fact_lag_micro_blg,    &  !< lag factor for micro-scale below-ground fluxes between tile and sibling [s-1]
      & fact_lag_macro_srf,    &  !< lag factor for macro-scale surface fluxes between tile and sibling [s-1]
      & fact_lag_macro_blg        !< lag factor for macro-scale below-ground fluxes between tile and sibling [s-1]

    ! Local variables ------------------------------------------------------------------------------------
    REAL(wp)                :: &
       & ret_micro_srf,        &  !< Retention time for micro-scale fluxes at the surfcace [s]
       & ret_micro_blg            !< Retention time for micro-scale fluxes below the surfcace [s]

    !$ACC ROUTINE SEQ

    ! Micro scale lag factors are determined by the (inverse) retention time, the micro-connectivity
    ! between tile and sibling, and the cover fraction of sibling relative to the total connected
    ! fraction of the tile.

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

    ! Micro scale lag factors are determined by the (inverse) retention time, the micro-connectivity
    ! between tile and sibling, and the cover fraction of sibling relative to the total connected
    ! fraction of the tile.
    ! Additionally, the lag factor is decreased according to the micro-scale connectivity, assuming
    ! that the fractions of the tile that have a micro-scale connectivity do not additionally feature
    ! macro-scale connections.
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

    ! For outflow on tiles that have no downstream tile, lag factors are calculated based on the
    ! time-scales of the macro-scale retention times. The respective fluxes are given to surface
    ! runoff and drainage.
    IF ((lex_idt == lex_ids) .AND. &
      & (fract_microcon_srf == 0.0_wp) .AND. (fract_macrocon_srf == 0.0_wp)) THEN
      fact_lag_micro_srf = 0._wp
      fact_lag_macro_srf = 1._wp / ret_macro_srf
    END IF

    IF ((lex_idt == lex_ids) .AND. &
      & (fract_microcon_blg == 0.0_wp) .AND. (fract_macrocon_blg == 0.0_wp)) THEN
      fact_lag_micro_blg = 0._wp
      fact_lag_macro_blg = 1._wp / ret_macro_blg
    END IF

  END SUBROUTINE get_flow_lags

  ! =================================================================================================================
  !>
  !>#### Check flow matrix for closed loops
  !>
  !> This routine checks whether any of the tiles in a flow matrix points back to itself
  !> through a series of connections.
  !>
  SUBROUTINE check_flow_matrix(               &
    & nconti,                                 &
    & matrix,                                 &
    & has_loop                                &
    )

    INTEGER,  INTENT(in)    :: nconti              !< number of connected tiles
    REAL(wp), INTENT(inout) :: matrix(:,:)         !< Flow matrix
    LOGICAL,  INTENT(out)   :: has_loop            !< .TRUE. if matrix contains closed loop
    INTEGER                 :: it,            &    !< First index variable for loop over connections
      &                        is,            &    !< Seconds index variable for loop over connections
      &                        cntr,          &    !< Index variable for closed loop check
      &                        src_id,        &    !< Id of current tile
      &                        trgt_id             !< Id of current sibling

    ! Initialize
    has_loop = .FALSE.

    ! Loop over flow matrix
    DO it = 1, nconti
      DO is = 1, nconti
        IF (matrix(is,it) > 0._wp) THEN
          ! sibling connects to tile
          IF (SUM(matrix(:,is)) > 0._wp) THEN
            ! sibling has tiles connected to it
            trgt_id = is                         ! sibling is current target
            src_id  = MAXLOC(matrix(:,is), 1)    ! getting the largest incomming connection
            cntr    = 0                          ! Starting counter for closed loop check
            DO WHILE (cntr < nconti)
              IF (SUM(matrix(:,src_id)) == 0._wp) THEN
                ! (incomming) connected tile is the end of flow network branch
                matrix(src_id,trgt_id) = 0._wp   ! Set the current source -> target to 0, to not get on the branch again
                IF (sum(matrix(:,is)) > 0._wp) THEN
                  ! Check if connected sibling still has incomming connections
                  trgt_id = is                   ! Set target back to sibling
                  src_id  = MAXLOC(matrix(:,is), 1) ! and find new source
                  cntr    = 0                    ! Going down new branch so counter is reset
                ELSE
                  matrix(is,it) = 0._wp          ! set sibling -> tile connection to 0 (no closed loop here)
                  cntr          = nconti - 1     ! nothing more to be checked
                END IF
              ELSE
                trgt_id = src_id
                src_id  = MAXLOC(matrix(:,src_id), 1)
              END IF
              cntr = cntr + 1
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
