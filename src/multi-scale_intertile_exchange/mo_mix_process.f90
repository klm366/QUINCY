!>
!! @brief Contains the routines for the mix processes
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-08-24)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_mix_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_runoff_transport, calc_drainage_transport, calc_intertile_diffusion

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of surface runoff
  !!
  PURE ELEMENTAL SUBROUTINE calc_runoff_transport(            &
    & dtime, fract_to_retention, retention_time,                 &
    & source_wtr_pond_max, source_fract, target_fract, runoff_source,       &
    & wtr_latflow_source, wtr_latflow_target, wtr_pond_source, &
    & runoff_retention_target,  wtr_pond_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        fract_to_retention, &
      &                        retention_time, &
      &                        source_wtr_pond_max,    & !< maximum surface water storage on source tile
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract          !< area fraction of water target tile
    REAL(wp), INTENT(inout) :: runoff_source,      & !< surface runoff on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        wtr_pond_source,    & !< surface water storage on source tile
      &                        runoff_retention_target, &
      &                        wtr_pond_target       !< surface water storage on target tile

    REAL(wp)                :: overland_flx,       &
      &                        wtr_to_retention,  &
      &                        wtr_from_retention

    !> Update target ponds with surface runoff from source tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Surface runoff is scaled with the tile fractions and moved into the target tile's
    !!   pond storage
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    overland_flx    = MAX(0._wp, wtr_pond_source + runoff_source * dtime / rhoh2o - source_wtr_pond_max)
    wtr_pond_source = MIN(source_wtr_pond_max, wtr_pond_source + runoff_source * dtime / rhoh2o)

    IF (target_fract > 1.0e-10_wp) THEN
      wtr_latflow_source = wtr_latflow_source - overland_flx * rhoh2o / dtime
      runoff_source      = 0._wp
      wtr_latflow_target = wtr_latflow_target + overland_flx * (source_fract / target_fract) * rhoh2o / dtime
      wtr_to_retention   = fract_to_retention * overland_flx * (source_fract / target_fract)
      wtr_pond_target    = wtr_pond_target    + overland_flx * (source_fract / target_fract) - wtr_to_retention
      runoff_retention_target = runoff_retention_target + wtr_to_retention
      wtr_from_retention = runoff_retention_target * dtime / retention_time
      wtr_pond_target    = wtr_pond_target + wtr_from_retention
      runoff_retention_target = runoff_retention_target - wtr_from_retention
    END IF

  END SUBROUTINE calc_runoff_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of lateral drainage
  !!
  PURE ELEMENTAL SUBROUTINE calc_drainage_transport(        &
    & dtime, fract_to_retention, retention_time,            &
    & source_fract, target_fract, field_cap_source,  field_cap_target,         &
    & latdrain_source, drainage_source, wtr_latflow_source, &
    & wtr_latflow_target, latdrain_retention_target,        &
    & wtr_soil_source, wtr_soil_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        fract_to_retention, &
      &                        retention_time,     &
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract,       & !< area fraction of water target tile
      &                        field_cap_source,   & !< soil field capacity of source tile
      &                        field_cap_target      !< soil field capacity of target tile
    REAL(wp), INTENT(inout) :: latdrain_source,    & !< lateral drainage on source tile
      &                        drainage_source,    & !< bottom drainage on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        latdrain_retention_target, &
      &                        wtr_soil_source,    & !< soil water storage on source tile
      &                        wtr_soil_target       !< soil water storage on target tile
    ! Local variables
    REAL(wp)                :: src2trgt,           &
      &                        trgt2src,           &
      &                        latflow_act,        &
      &                        wtr_to_retention,   &
      &                        wtr_from_retention, &
      &                        wtr_remoist

    !> Update target soil with lateral drainage from source tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Lateral drainage is scaled with the tile fractions and moved into the target tile's
    !!   soil storage. It flows into the same layer as it originates from on the source tile.
    !!   Excess drainage is added to the total drainage of the source tile.
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    IF (target_fract > 1.0e-10_wp) THEN
      src2trgt           = source_fract / target_fract
      ! Substract diagnosed lateral drainage from total drainage
      ! and from latflow of srouce tile
      drainage_source    = drainage_source    - latdrain_source
      wtr_latflow_source = wtr_latflow_source - latdrain_source
      ! Scale and add all lateral drainage to target tile soil layer
      latflow_act        = latdrain_source * src2trgt
      wtr_latflow_target = wtr_latflow_target + latflow_act
      latdrain_source    = 0._wp

      ! converting fluxes into volume
      latflow_act        = latflow_act * dtime / rhoh2o
      ! and add to the retention resersoir of the target tile
      ! or to the tile directly if retention reservoir is not taken into considderation (micro scale)
      wtr_to_retention   = fract_to_retention * latflow_act
      wtr_soil_target    = wtr_soil_target + latflow_act - wtr_to_retention
      latdrain_retention_target = latdrain_retention_target + wtr_to_retention
      ! calculate fluxes from retention reservoir to soil column in target tile
      ! and adding water to soil and subtracting it from reservoir
      wtr_from_retention = latdrain_retention_target * dtime / retention_time
      latdrain_retention_target = latdrain_retention_target - wtr_from_retention
      wtr_soil_target           = wtr_soil_target           + wtr_from_retention
      ! water that has not left the reservoir can remoist the soil in the source tile
      ! if water is below the field capacity
      IF (source_fract > 1.0e-10_wp .AND. &
          latdrain_retention_target > 1.0e-10_wp .AND. &
          wtr_soil_source < field_cap_source) THEN
        trgt2src    = target_fract / source_fract
        ! calculating the maximum amount of water that re-enters the source tile
        wtr_remoist = MIN(latdrain_retention_target*trgt2src, &
                          field_cap_source - wtr_soil_source)
        ! add it to the soil in the tile
        wtr_soil_source = wtr_soil_source + wtr_remoist
        ! adjust the diagnosed lateral flux from source to sink
        latflow_act        = wtr_remoist * rhoh2o / dtime
        wtr_latflow_source = wtr_latflow_source + latflow_act
        !latdrain_source    = latdrain_source + latflow_act
        ! update the retention reservoir on the target tile
        latdrain_retention_target = latdrain_retention_target - wtr_remoist * src2trgt
        ! and adjust the lateral flux to the target tile
        latflow_act        = wtr_remoist * src2trgt * rhoh2o / dtime
        wtr_latflow_target = wtr_latflow_target - latflow_act
      END IF
    END IF

  END SUBROUTINE calc_drainage_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles with different surface elevation (elevation/depression)
  !!
  SUBROUTINE calc_intertile_diffusion(nc, nsoil, soilhydmodel, dtime,   &
      &      h_offset, lvt_fract, dpr_fract, dsoil_lvt,dsoil_dpr,          &
      &      field_cap_lvt, field_cap_dpr, wsat_lvt, wsat_dpr,   &
      &      pwp_lvt, pwp_dpr, mpot_sat_lvt, mpot_sat_dpr,       &
      &      k_sat_lvt, k_sat_dpr, bclapp_lvt, bclapp_dpr, ps_index_lvt, ps_index_dpr,             &
      &      ice_soil_lvt, ice_soil_dpr, wtr_scool_lvt, wtr_scool_dpr, wtr_pond_lvt, wtr_pond_dpr, &
      &      wtr_latflow_lvt, wtr_latflow_dpr, wtr_soil_lvt, wtr_soil_dpr, diffus_lvt, diffus_dpr  &
      &      )

    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: get_soilhyd_cond_properties

    INTEGER, INTENT(in)     :: nc,                  & !< number of grid cells
      &                        nsoil,               & !< number of soil layers
                               soilhydmodel           !< Choice of soil hydrological model

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
      &                        h_offset,            &
      &                        lvt_fract,           & !< area fraction of water elevation tile
      &                        dpr_fract              !< area fraction of water depression tile

    REAL(wp), INTENT(in)    :: dsoil_lvt(:),        &
      &                        dsoil_dpr(:),        &
      &                        field_cap_lvt(:),    & !< soil field capacity of elevation tile
      &                        field_cap_dpr(:),    & !< soil field capacity of depression tile
      &                        wsat_lvt(:),         &
      &                        wsat_dpr(:),         &
      &                        pwp_lvt(:),          &
      &                        pwp_dpr(:),          &
      &                        mpot_sat_lvt(:),     &
      &                        mpot_sat_dpr(:),     &
      &                        k_sat_lvt(:),        &
      &                        k_sat_dpr(:),        &
      &                        bclapp_lvt(:),       &
      &                        bclapp_dpr(:),       &
      &                        ps_index_lvt(:),     &
      &                        ps_index_dpr(:),     &
      &                        ice_soil_lvt(:),     &
      &                        ice_soil_dpr(:),     &
      &                        wtr_scool_lvt(:),    &
      &                        wtr_scool_dpr(:)


    REAL(wp), INTENT(inout) :: wtr_pond_lvt,    & !< lateral water flow from elevation tile
      &                        wtr_pond_dpr,    & !< lateral water flow towards depression tile
      &                        wtr_latflow_lvt, & !< lateral water flow from elevation tile
      &                        wtr_latflow_dpr    !< lateral water flow towards depression tile

    REAL(wp), INTENT(inout) :: wtr_soil_lvt(:),     & !< soil water storage on elevation tile
      &                        wtr_soil_dpr(:)        !< soil water storage on depression tile

    REAL(wp), INTENT(out)   :: diffus_lvt(:),       & !< Diffusion flux at elevation tile [kg m-2 s-1]
      &                        diffus_dpr(:)          !< Diffusion flux at depression tile [kg m-2 s-1]



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: nrgrd, i, i_lvt, i_dpr

    REAL(wp)                :: fr_lvt_rel,       &
      &                        fr_dpr_rel,       &
      &                        pond_flx_lvt,     &
      &                        pond_flx_dpr,     &
      &                        wtr_pond_total,   &
      &                        z_lvt,            &
      &                        z_dpr,            &
      &                        len_inter,        &
      &                        h,                &
      &                        a_inter,          &
      &                        k_inter,          &
      &                        d_mpot,           &
      &                        wtr_rel_qlbrm,    &
      &                        free_space,       &
      &                        wtr_avail,        &
      &                        wtr_2_qlbrm

    REAL(wp)                ::                           &
      &                        mpot_lvt(nsoil),          &
      &                        mpot_dpr(nsoil),          &
      &                        ice_impede_lvt(nsoil),    &
      &                        ice_impede_dpr(nsoil),    &
      &                        k_lvt(nsoil),             &
      &                        k_dpr(nsoil),             &
      &                        wtr_soil_lvt_updt(nsoil), &
      &                        wtr_soil_dpr_updt(nsoil), &
      &                        field_cap_wi_lvt(nsoil),    & !< soil field capacity incl. ice of elevation tile
      &                        field_cap_wi_dpr(nsoil),    & !< soil field capacity incl. ice of depression tile
      &                        wsat_wi_lvt(nsoil),         &
      &                        wsat_wi_dpr(nsoil),         &
      &                        pwp_wi_lvt(nsoil),          &
      &                        pwp_wi_dpr(nsoil)


    REAL(wp)                :: z_rgrd_top,        &
      &                        z_rgrd_bot,        &
      &                        wtr_soil_lvt_rgrd, &
      &                        ice_soil_lvt_rgrd, &
      &                        wsat_lvt_rgrd,     &
      &                        field_cap_lvt_rgrd,&
      &                        pwp_lvt_rgrd,      &
      &                        mpot_lvt_rgrd,     &
      &                        k_lvt_rgrd,        &
      &                        k_sat_lvt_rgrd,    &
      &                        wtr_soil_dpr_rgrd, &
      &                        ice_soil_dpr_rgrd, &
      &                        wsat_dpr_rgrd,     &
      &                        field_cap_dpr_rgrd,&
      &                        pwp_dpr_rgrd,      &
      &                        mpot_dpr_rgrd,     &
      &                        k_dpr_rgrd,        &
      &                        wtr_flx_lat,       &
      &                        wtr_lat

    REAL(wp)                :: ice_impede_lvt_sl(nc,nsoil), &
      &                        ice_impede_dpr_sl(nc,nsoil), &
      &                        k_lvt_sl(nc,nsoil),          &
      &                        k_dpr_sl(nc,nsoil),          &
      &                        mpot_lvt_sl(nc,nsoil),       &
      &                        mpot_dpr_sl(nc,nsoil)



!    REAL(wp), PARAMETER ::     length_scale = 1.0_wp,        & ! assumed diameter of depressions (could be input field / namelist para)
!      &                        h_offset     = 0.1_wp,        & ! hight difference between ridge and depression (should later be input field)
!      &                        ck           = 8.0_wp           ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    REAL(wp) ::     length_scale,    & ! assumed diameter of depressions (could be input field / namelist para)
!      &             h_offset,        & ! hight difference between ridge and depression (should later be input field)
      &             ck                 ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    nrgrd = nsoil * 2 + 3

    ! Initilize fields for accumulating regridded output
    wtr_soil_lvt_updt(:) = 0._wp
    wtr_soil_dpr_updt(:) = 0._wp

    pond_flx_lvt         = 0._wp
    pond_flx_dpr         = 0._wp

    diffus_lvt(:)        = 0._wp
    diffus_dpr(:)        = 0._wp

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    !Assumption of circular depression with A = 1 m2 with in a circular elevation
    !length_scale = ((1._wp/pi)**0.5_wp+((1._wp/fr_dpr_rel)/pi)**0.5_wp)/2._wp ! From center depression to half of outer rim
    length_scale = ((1._wp/fr_dpr_rel)/pi)**0.5_wp                      ! From center to outer rim

    ! Fist calculate the standing water at the surface of both tiles ... with the assumption that runoff is instant
    wtr_pond_total = wtr_pond_lvt * fr_lvt_rel + wtr_pond_dpr * fr_dpr_rel

    IF (wtr_pond_total / fr_dpr_rel <= h_offset ) THEN
      pond_flx_dpr =   wtr_pond_total / fr_dpr_rel - wtr_pond_dpr
      pond_flx_lvt = - wtr_pond_lvt
      wtr_pond_dpr =   wtr_pond_total / fr_dpr_rel
      wtr_pond_lvt =   0._wp
    ELSE
      pond_flx_dpr =   h_offset  + (wtr_pond_total - h_offset * fr_dpr_rel) - wtr_pond_dpr
      pond_flx_lvt =   wtr_pond_total - h_offset * fr_dpr_rel - wtr_pond_lvt
      wtr_pond_dpr =   h_offset  + (wtr_pond_total - h_offset * fr_dpr_rel)
      wtr_pond_lvt =   wtr_pond_total - h_offset * fr_dpr_rel
    END IF

    ! Determine soil parameters on the native grids

    !!!!!! Getting back from soil charctierstics excluding the ice to general characteristics
    field_cap_wi_lvt(:) = field_cap_lvt(:) + ice_soil_lvt(:)
    field_cap_wi_dpr(:) = field_cap_dpr(:) + ice_soil_dpr(:)
    wsat_wi_lvt(:)      = wsat_lvt(:)      + ice_soil_lvt(:)
    wsat_wi_dpr(:)      = wsat_dpr(:)      + ice_soil_dpr(:)
    WHERE(wsat_lvt(:) > 0._wp)
      pwp_wi_lvt(:)       = pwp_lvt(:) * wsat_wi_lvt(:) / wsat_lvt(:)
    ELSEWHERE(wsat_wi_lvt(:) > 0._wp)
      pwp_wi_lvt(:)       = ice_soil_lvt(:)
    ELSEWHERE
      pwp_wi_lvt(:)       = 0._wp
    END WHERE

    WHERE(wsat_dpr(:) > 0._wp)
      pwp_wi_dpr(:)       = pwp_lvt(:) * wsat_wi_dpr(:) / wsat_dpr(:)
    ELSEWHERE(wsat_wi_dpr(:) > 0._wp)
      pwp_wi_dpr(:)       = ice_soil_dpr(:)
    ELSEWHERE
      pwp_wi_dpr(:)       = 0._wp
    END WHERE

    call get_soilhyd_cond_properties(soilhydmodel, 1, -1._wp, nc, nsoil,        &
      & SPREAD(1,                   DIM=1, ncopies=nc), & ! active soil layer set to one
      & SPREAD(dsoil_lvt(:),        DIM=1, ncopies=nc), & ! depth of acutal soil layer
      & SPREAD(wtr_soil_lvt(:),     DIM=1, ncopies=nc), &
      & SPREAD(ice_soil_lvt(:),     DIM=1, ncopies=nc), &
      & SPREAD(wtr_scool_lvt(:),    DIM=1, ncopies=nc), &
      & SPREAD(wsat_wi_lvt(:),      DIM=1, ncopies=nc), &
      & SPREAD(field_cap_wi_lvt(:), DIM=1, ncopies=nc), &
      & SPREAD(pwp_wi_lvt(:),       DIM=1, ncopies=nc), &
      & SPREAD(k_sat_lvt(:),        DIM=1, ncopies=nc), &
      & SPREAD(mpot_sat_lvt(:),     DIM=1, ncopies=nc), &
      & SPREAD(bclapp_lvt(:),       DIM=1, ncopies=nc), &
      & SPREAD(ps_index_lvt(:),     DIM=1, ncopies=nc), &
      & ice_impedance=ice_impede_lvt_sl(:,:),           &
      & K=k_lvt_sl(:,:),                                &
      & mpot_act=mpot_lvt_sl(:,:)                       &
      )
      ice_impede_lvt(:) = ice_impede_lvt_sl(1,:)
      k_lvt(:)          = k_lvt_sl(1,:)
      mpot_lvt(:)       = mpot_lvt_sl(1,:)

    call get_soilhyd_cond_properties(soilhydmodel, 1, -1._wp, nc, nsoil,        &
      & SPREAD(1,                   DIM=1, ncopies=nc), & ! active soil layer set to one
      & SPREAD(dsoil_dpr(:),        DIM=1, ncopies=nc), & ! depth of acutal soil layer
      & SPREAD(wtr_soil_dpr(:),     DIM=1, ncopies=nc), &
      & SPREAD(ice_soil_dpr(:),     DIM=1, ncopies=nc), &
      & SPREAD(wtr_scool_dpr(:),    DIM=1, ncopies=nc), &
      & SPREAD(wsat_wi_dpr(:),      DIM=1, ncopies=nc), &
      & SPREAD(field_cap_wi_dpr(:), DIM=1, ncopies=nc), &
      & SPREAD(pwp_wi_dpr(:),       DIM=1, ncopies=nc), &
      & SPREAD(k_sat_dpr(:),        DIM=1, ncopies=nc), &
      & SPREAD(mpot_sat_dpr(:),     DIM=1, ncopies=nc), &
      & SPREAD(bclapp_dpr(:),       DIM=1, ncopies=nc), &
      & SPREAD(ps_index_dpr(:),     DIM=1, ncopies=nc), &
      & ice_impedance=ice_impede_lvt_sl(:,:),           &
      & K=k_dpr_sl(:,:),                                &
      & mpot_act=mpot_dpr_sl(:,:)                       &
      )
      ice_impede_dpr(:) = ice_impede_dpr_sl(1,:)
      k_dpr(:)          = k_dpr_sl(1,:)
      mpot_dpr(:)       = mpot_dpr_sl(1,:)

    ! Start setting up the interpolated grid
    i_lvt = 1
    z_lvt = dsoil_lvt(i_lvt)
    i_dpr = 1
    IF (h_offset > 0._wp) THEN
      i_dpr = 0
      IF (h_offset > wtr_pond_dpr) THEN
        z_dpr = h_offset - wtr_pond_dpr
      ELSE
        z_dpr = h_offset
      END IF
    ELSE
      z_dpr = dsoil_dpr(i_dpr)
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, z_dpr)

    len_inter = ((4._wp * pi) ** 0.5_wp) * fr_dpr_rel ! base assumption is the circumference of circle with area 1 m2 within elevation area x

    DO i = 1, nrgrd
      IF (i_lvt <= nsoil .AND. i_dpr <= nsoil .AND. dsoil_lvt(i_lvt) > 0._wp .AND. dsoil_dpr(MAX(i_dpr,1)) > 0._wp) THEN ! Above Bedrock in both tiles
        h                   = z_rgrd_bot - z_rgrd_top
        a_inter             = len_inter * h

        wtr_soil_lvt_rgrd   = wtr_soil_lvt(i_lvt)      * h / dsoil_lvt(i_lvt)
        ice_soil_lvt_rgrd   = ice_soil_lvt(i_lvt)      * h / dsoil_lvt(i_lvt)
        wsat_lvt_rgrd       = wsat_wi_lvt(i_lvt)       * h / dsoil_lvt(i_lvt)
        field_cap_lvt_rgrd  = field_cap_wi_lvt(i_lvt)  * h / dsoil_lvt(i_lvt)
        pwp_lvt_rgrd        = pwp_wi_lvt(i_lvt)        * h / dsoil_lvt(i_lvt)
        mpot_lvt_rgrd       = mpot_lvt(i_lvt)
        k_sat_lvt_rgrd      = k_sat_lvt(i_lvt)
        k_lvt_rgrd          = k_lvt(i_lvt)

        IF (z_rgrd_bot <= h_offset - wtr_pond_dpr) THEN ! Still above water level in depression
          ! Water flows out of the soil in the elevation tile into the pond reservoir of the depression tile
          wtr_flx_lat       = k_lvt_rgrd * a_inter
          wtr_avail         = MAX(0._wp,(wtr_soil_lvt_rgrd - pwp_lvt_rgrd) * fr_lvt_rel)
          wtr_lat           = MIN(wtr_flx_lat * dtime, wtr_avail)
          wtr_flx_lat       = wtr_lat / dtime
          wtr_pond_dpr      = wtr_pond_dpr + wtr_lat / fr_dpr_rel
          wtr_soil_lvt_rgrd = wtr_soil_lvt_rgrd - wtr_lat / fr_lvt_rel
          wtr_soil_dpr_rgrd = 0.0_wp

        ELSE IF (z_rgrd_bot <= h_offset ) THEN ! Within the ponded surface water in the depression
          ! Water flows from the pond reservoir of the depression tile into the soil of the elevation tile
!          wtr_flx_lat       = k_lvt_rgrd * a_inter
          wtr_flx_lat       = ((k_lvt_rgrd * k_sat_lvt_rgrd)**0.5_wp)*a_inter
!          free_space        = MAX(0._wp,(wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          free_space        = MAX(0._wp,(field_cap_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          wtr_avail         = MIN(h, wtr_pond_dpr) * fr_dpr_rel
          wtr_lat           = MAX(-wtr_flx_lat * dtime, -free_space, -wtr_avail)
          wtr_flx_lat       = wtr_lat / dtime
          wtr_pond_dpr      = wtr_pond_dpr + wtr_lat / fr_dpr_rel
          wtr_soil_lvt_rgrd = wtr_soil_lvt_rgrd - wtr_lat / fr_lvt_rel
          wtr_soil_dpr_rgrd = 0.0_wp

        ELSE ! Within the soil column of both tiles
          ! Flux is determined by the difference in matric potential between the two tiles
          wtr_soil_dpr_rgrd = wtr_soil_dpr(i_dpr) * h / dsoil_dpr(i_dpr)
          ice_soil_dpr_rgrd = ice_soil_dpr(i_dpr) * h / dsoil_dpr(i_dpr)

          wsat_dpr_rgrd     = wsat_wi_dpr(i_dpr)  * h / dsoil_dpr(i_dpr)
          field_cap_dpr_rgrd= field_cap_wi_dpr(i_dpr)  * h / dsoil_dpr(i_dpr)
          pwp_dpr_rgrd      = pwp_wi_dpr(i_dpr)   * h / dsoil_dpr(i_dpr)
          mpot_dpr_rgrd     = mpot_dpr(i_dpr)
          k_dpr_rgrd        = k_dpr(i_dpr)

          !k_inter  =    k_dpr_rgrd * fr_dpr_rel + k_lvt_rgrd * fr_lvt_rel
          !k_inter           = MIN(k_dpr_rgrd, k_lvt_rgrd)
          k_inter           =(k_dpr_rgrd * k_lvt_rgrd)**0.5_wp
          d_mpot            = mpot_dpr_rgrd - mpot_lvt_rgrd
          wtr_flx_lat       = - k_inter  * d_mpot / length_scale *  a_inter

          wtr_rel_qlbrm     = MIN(1._wp, ((wtr_soil_lvt_rgrd+ice_soil_lvt_rgrd) * fr_lvt_rel       &
                          & +  (wtr_soil_dpr_rgrd+ice_soil_dpr_rgrd) * fr_dpr_rel)      &
                          & / (wsat_lvt_rgrd * fr_lvt_rel +wsat_dpr_rgrd * fr_dpr_rel ))

          IF (wtr_flx_lat > 0._wp) THEN ! Flux is going from elevation to depression

            wtr_2_qlbrm     = MAX(0._wp, &
                          &  (wtr_rel_qlbrm * wsat_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
!            free_space      = MAX(0._wp,(wsat_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
            free_space      = MAX(0._wp,(field_cap_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
            wtr_avail       = MAX(0._wp,(wtr_soil_lvt_rgrd - pwp_lvt_rgrd) * fr_lvt_rel)
            wtr_lat         = MIN(wtr_flx_lat * dtime, free_space, wtr_avail, wtr_2_qlbrm)
            wtr_lat         = MAX(0._wp, wtr_lat)
            wtr_flx_lat     = wtr_lat / dtime

          ELSE ! Flux is going from depression to elevation

            wtr_2_qlbrm     = MAX(0._wp, &
                           & (wtr_rel_qlbrm  * wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
!            free_space      = MAX(0._wp, (wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
            free_space      = MAX(0._wp, (field_cap_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
            wtr_avail       = MAX(0._wp, (wtr_soil_dpr_rgrd - pwp_dpr_rgrd) * fr_dpr_rel)
            wtr_lat         = MAX(wtr_flx_lat * dtime, -free_space, -wtr_avail, -wtr_2_qlbrm)
            wtr_lat         = MIN(0._wp, wtr_lat)
            wtr_flx_lat     = wtr_lat / dtime

          END IF

          wtr_soil_lvt_rgrd = wtr_soil_lvt_rgrd - wtr_lat / fr_lvt_rel
          wtr_soil_dpr_rgrd = wtr_soil_dpr_rgrd + wtr_lat / fr_dpr_rel

        END IF

        wtr_soil_lvt_updt(i_lvt)  = wtr_soil_lvt_updt(i_lvt) + wtr_soil_lvt_rgrd
        diffus_lvt(i_lvt)         = diffus_lvt(i_lvt) - wtr_lat / fr_lvt_rel

        IF (i_dpr > 0) THEN
          wtr_soil_dpr_updt(i_dpr)  = wtr_soil_dpr_updt(i_dpr) + wtr_soil_dpr_rgrd
          diffus_dpr(i_dpr)         = diffus_dpr(i_dpr) + wtr_lat / fr_dpr_rel
        ELSE
          pond_flx_dpr              = pond_flx_dpr + wtr_lat / fr_dpr_rel
        END IF

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

        IF (z_rgrd_bot == z_dpr) THEN
          IF (z_dpr < h_offset) THEN  ! Are we still above the surface ? Then set z_dpr to surface
            z_dpr = h_offset
          ELSE
            i_dpr = i_dpr + 1
            IF (i_dpr <= nsoil) THEN
              IF (dsoil_dpr(i_dpr) > 0._wp) THEN
                z_dpr = z_dpr + dsoil_dpr(i_dpr)
              ELSE
                z_dpr = 10000._wp
                i_dpr = nsoil + 1
              END IF
            ELSE
              z_dpr = 10000._wp
              i_dpr = nsoil + 1
            END IF
          END IF
        END IF

        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=MIN(z_lvt, z_dpr)


      ELSE IF (i_lvt <= nsoil .AND. dsoil_lvt(i_lvt) > 0._wp ) THEN  ! Above Bedrock in lvt tile only (should only be the case when using different soil depths)
        ! hence there are no fluxes and wtr is just mapped to wtr_updt
        h                        = z_rgrd_bot - z_rgrd_top
        wtr_soil_lvt_rgrd        = wtr_soil_lvt(i_lvt) * h / dsoil_lvt(i_lvt)
        wtr_soil_lvt_updt(i_lvt) = wtr_soil_lvt_updt(i_lvt) + wtr_soil_lvt_rgrd

        IF (z_rgrd_bot == z_lvt) THEN
          i_lvt=i_lvt+1
          IF (i_lvt <= nsoil) THEN
            IF (dsoil_dpr(i_lvt) > 0._wp) THEN
              z_lvt = z_lvt + dsoil_lvt(i_lvt)
            END IF
          END IF
        END IF
        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=z_lvt

      ELSE IF (i_dpr > 0 .AND. i_dpr <= nsoil .AND. dsoil_dpr(i_dpr) > 0._wp ) THEN   ! Above Bedrock in dpr tile only
        ! hence there are no fluxes and wtr is just mapped to wtr_updt
        h                        = z_rgrd_bot - z_rgrd_top
        wtr_soil_dpr_rgrd        = wtr_soil_dpr(i_dpr) * h / dsoil_dpr(i_dpr)
        wtr_soil_dpr_updt(i_dpr) = wtr_soil_dpr_updt(i_dpr) + wtr_soil_dpr_rgrd

        IF (z_rgrd_bot == z_dpr) THEN
          i_dpr=i_dpr+1
          IF (i_dpr <= nsoil) THEN
            IF (dsoil_dpr(i_dpr) > 0._wp) THEN
              z_dpr = z_dpr + dsoil_dpr(i_dpr)
            END IF
          END IF
        END IF

        z_rgrd_top=z_rgrd_bot
        z_rgrd_bot=z_dpr

      END IF
    END DO

    wtr_soil_lvt(:)  = wtr_soil_lvt_updt(:)
    wtr_soil_dpr(:)  = wtr_soil_dpr_updt(:)

    wtr_latflow_lvt = wtr_latflow_lvt + (SUM(diffus_lvt(:)) + pond_flx_lvt) * rhoh2o / dtime
    wtr_latflow_dpr = wtr_latflow_dpr + (SUM(diffus_dpr(:)) + pond_flx_dpr) * rhoh2o / dtime

    diffus_lvt(:) = diffus_lvt(:) * rhoh2o / dtime
    diffus_dpr(:) = diffus_dpr(:) * rhoh2o / dtime

  END SUBROUTINE calc_intertile_diffusion

#endif
END MODULE mo_mix_process
