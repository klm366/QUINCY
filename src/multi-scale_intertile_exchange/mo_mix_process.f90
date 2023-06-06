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

  PUBLIC :: calc_runoff_transport, calc_drainage_transport, calc_intertile_diffusion, calc_snow_redistribution

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of surface runoff
  !!
  PURE ELEMENTAL SUBROUTINE calc_runoff_transport(             &
    & l_lag, dtime, retention_time, wtr_pond_max_source,       &
    & source_fract, target_fract, runoff_source,               &
    & wtr_latflow_source, wtr_latflow_target, wtr_pond_source, &
    & wtr_pond_target, runoff_retention_source                 &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    LOGICAL,  INTENT(in)    :: l_lag                      !< T: apply lag to overland flow
    REAL(wp), INTENT(in)    :: dtime,                   & !< timestep length                        [s]
      &                        retention_time,          & !< retention time for lateral transport   [s]
      &                        wtr_pond_max_source,     & !< maximum surface water storage on source tile [m]
      &                        source_fract,            & !< area fraction of water source tile
      &                        target_fract               !< area fraction of water target tile
    REAL(wp), INTENT(inout) :: runoff_source,           & !< surface runoff on source tile          [kg m-2 s-1]
      &                        wtr_latflow_source,      & !< lateral water flow from source tile    [kg m-2 s-1]
      &                        wtr_latflow_target,      & !< lateral water flow towards target tile [kg m-2 s-1]
      &                        wtr_pond_source,         & !< surface water storage on source tile   [m]
      &                        wtr_pond_target,         & !< surface water storage on target tile   [m]
      &                        runoff_retention_source    !< intermediate storage for overland flow [m]

    REAL(wp)                :: overland_flx_source,     & !< surface flow leaving source tile       [m]
      &                        overland_flx_target        !< surface flow arriving at target tile   [m]

    !> Update target ponds with surface runoff/pond overflow from source tile
    !!   Surface runoff is scaled with the tile fractions and moved into the target tile's
    !!   pond storage.
    !!   In case of long transport routes (upland to lowland) a retention storage is used.
    !! @todo the retention store lag time needs to be physics based, not fixed to some arbitrary value

    IF (target_fract > 1.0e-10_wp) THEN
      ! Compute pond overflow on source tile and update source tile pond storage
      overland_flx_source = MAX(0._wp, wtr_pond_source + runoff_source * dtime / rhoh2o - wtr_pond_max_source)
      wtr_pond_source     = MIN(wtr_pond_max_source, wtr_pond_source + runoff_source * dtime / rhoh2o)
      ! Set source tile runoff to zero as it goes down the hill
      runoff_source       = 0._wp

      ! Compute flows for and update runoff retention storage
      IF (l_lag) THEN
        runoff_retention_source = runoff_retention_source + overland_flx_source
        overland_flx_source     = runoff_retention_source * MIN(1.0_wp, dtime / retention_time)
        runoff_retention_source = MAX(0._wp, runoff_retention_source - overland_flx_source)
      END IF

      ! Scale overland flow to account for tile fraction differences
      overland_flx_target = overland_flx_source * (source_fract / target_fract)
      ! Update target pond storage
      wtr_pond_target     = wtr_pond_target + overland_flx_target
      ! Update lateral flows
      wtr_latflow_source  = wtr_latflow_source - overland_flx_source * rhoh2o / dtime
      wtr_latflow_target  = wtr_latflow_target + overland_flx_target * rhoh2o / dtime

    END IF

  END SUBROUTINE calc_runoff_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of lateral drainage
  !!
  PURE ELEMENTAL SUBROUTINE calc_drainage_transport(        &
    & l_lag, l_latdrain, dtime, retention_time,             &
    & field_cap_source, source_fract, target_fract,         &
    & drainlat_source, drainage_source, wtr_latflow_source, &
    & wtr_latflow_target, wtr_soil_source, wtr_soil_target, &
    & drainlat_retention_source                             &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    LOGICAL,  INTENT(in)    :: l_lag                       !< T: apply lag to overland flow
    LOGICAL,  INTENT(in)    :: l_latdrain                  !< T: allow lateral drainage
    REAL(wp), INTENT(in)    :: dtime,                    & !< timestep length                        [s]
      &                        retention_time,           & !< retention time for lateral transport   [s]
      &                        field_cap_source,         & !< soil field capacity of source tile     [m]
      &                        source_fract,             & !< area fraction of water source tile     []
      &                        target_fract                !< area fraction of water target tile     []
    REAL(wp), INTENT(inout) :: drainlat_source,          & !< lateral drainage on source tile        [kg m-2 s-1]
      &                        drainage_source,          & !< bottom drainage on source tile         [kg m-2 s-1]
      &                        wtr_latflow_source,       & !< lateral water flow from source tile    [kg m-2 s-1]
      &                        wtr_latflow_target,       & !< lateral water flow towards target tile [kg m-2 s-1]
      &                        wtr_soil_source,          & !< soil water storage on source tile      [m]
      &                        wtr_soil_target,          & !< soil water storage on target tile      [m]
      &                        drainlat_retention_source   !< lateral drainage retention storage     [m]

    ! Local variables
    REAL(wp)                :: drainlat_flx_source,      & !< lateral drainage flux leaving source tile                        [m]
      &                        drainlat_flx_target,      & !< lateral drainage flux arriving at target tile                    [m]
      &                        wtr_remoist_source          !< water from retention storage re-entering the soil at source tile [m]

    !> Update target soil with lateral drainage from source tile
    !!   Lateral drainage is scaled with the tile fractions and moved into the target tile's
    !!   soil storage. It flows into the same layer as it originates from on the source tile.
    !!   Excess drainage is added to the total drainage of the source tile.
    !!   In case of long transport routes (upland to lowland) a retention storage is used.
    !!   This reservoir can remoist the soil water storage in the source tile in later time steps.
    !! @todo the retention store lag time needs to be physics based, not fixed to some arbitrary value
    !! @todo where does the retention storage reside physically? On the hillslope? Doesn't this concept
    !!       bypass the upper limit for soil moisture storage

    IF (l_latdrain) THEN
      ! Water within the retention storage is allowed to remoist the soil layer water storage
      IF (drainlat_retention_source > 1.0e-10_wp .AND. wtr_soil_source < field_cap_source) THEN
        wtr_remoist_source        = MIN(drainlat_retention_source, field_cap_source - wtr_soil_source)
        wtr_soil_source           = wtr_soil_source + wtr_remoist_source
        drainlat_retention_source = MAX(0._wp, drainlat_retention_source - wtr_remoist_source)
      END IF

      ! Compute lateral flow
      IF (target_fract > 1.0e-10_wp) THEN
        ! Compute lateral drainage on source tile and update total drainage and total lateral water flow
        drainlat_flx_source = drainlat_source    * dtime / rhoh2o
        drainage_source     = drainage_source    - drainlat_source
        ! Set source tile lateral drainage to zero as it left the tile
        drainlat_source     = 0._wp

        ! Compute flows for and update drainage retention storage (only macro scale flow)
        IF (l_lag) THEN
          drainlat_retention_source = drainlat_retention_source + drainlat_flx_source
          drainlat_flx_source       = drainlat_retention_source * MIN(1.0_wp, dtime / retention_time)
          drainlat_retention_source = MAX(0._wp, drainlat_retention_source - drainlat_flx_source)
        END IF

        ! Scale lateral drainage flow to account for tile fraction differences
        drainlat_flx_target = drainlat_flx_source * (source_fract / target_fract)
        ! Update target soil layer water content
        wtr_soil_target    = wtr_soil_target + drainlat_flx_target
        ! Update lateral flows
        wtr_latflow_source  = wtr_latflow_source - drainlat_flx_source * rhoh2o / dtime
        wtr_latflow_target  = wtr_latflow_target + drainlat_flx_target * rhoh2o / dtime

      END IF
    END IF

  END SUBROUTINE calc_drainage_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles
  !!
  PURE ELEMENTAL SUBROUTINE calc_intertile_diffusion_old(           &
    & dtime, source_fract, target_fract, field_cap_source,      &
    & field_cap_target, wtr_latflow_source, wtr_latflow_target, &
    & wtr_soil_source, wtr_soil_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
      &                        source_fract,        & !< area fraction of water source tile
      &                        target_fract,        & !< area fraction of water target tile
      &                        field_cap_source,    & !< soil field capacity of source tile
      &                        field_cap_target       !< soil field capacity of target tile
    REAL(wp), INTENT(inout) :: wtr_latflow_source,  & !< lateral water flow from source tile
      &                        wtr_latflow_target,  & !< lateral water flow towards target tile
      &                        wtr_soil_source,     & !< soil water storage on source tile
      &                        wtr_soil_target        !< soil water storage on target tile

    ! Local variables
    REAL(wp)                :: diffus_source,       & !< Diffusion flux at source tile per time step
      &                        diffus_target,       & !< Diffusion flux at target tile per time step
      &                        wtr_soil_rel_source, & !< Relative soil water content in source tile
      &                        wtr_soil_rel_target, & !< Relative soil water content in target tile
      &                        wtr_soil_rel_mean      !< Mean relative soil water content of both tiles

    !> Compute lateral diffusion flux between soil water storages of source and target tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Lateral diffusion is implemented such that is compensates the the differences
    !!   in relative soil water content between both tiles. Thus, the flow is not necessarily
    !!   from the source to the target tile but follows the moisture gradient.
    !! @todo This dummy implementation sets the lateral diffusion flux to 10% of what would
    !!   be needed to fully compensate the soil moisture difference. Later on, soil hydraulic
    !!   conductivity should be used.

    IF (source_fract > 1.0e-10_wp .AND. target_fract > 1.0e-10_wp .AND. &
        field_cap_source > 1.0e-10_wp .AND. field_cap_target > 1.0e-10_wp) THEN
      ! Compute relative soil water content as well as the fraction-weighted relative soil water content mean between both tiles
      wtr_soil_rel_source = wtr_soil_source / field_cap_source
      wtr_soil_rel_target = wtr_soil_target / field_cap_target
      wtr_soil_rel_mean  = 1.0 / (source_fract + target_fract) * &
        & (wtr_soil_rel_source * source_fract + wtr_soil_rel_target * target_fract)
      ! Compute diffusive flow between both tiles
      diffus_target = (wtr_soil_rel_mean - wtr_soil_rel_target) * field_cap_target * 0.1
      diffus_source = (wtr_soil_rel_mean - wtr_soil_rel_source) * field_cap_source * 0.1
      ! Add diffusive flow to soil storages
      wtr_soil_target = wtr_soil_target + diffus_target
      wtr_soil_source = wtr_soil_source + diffus_source
      ! Add diffusive flow to overal lateral flow
      wtr_latflow_target = wtr_latflow_target + diffus_target * rhoh2o / dtime
      wtr_latflow_source = wtr_latflow_source + diffus_source * rhoh2o / dtime
    END IF

  END SUBROUTINE calc_intertile_diffusion_old


  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles with different surface elevation (elevation/depression)
  !!
  SUBROUTINE calc_intertile_diffusion(nsoil, soilhydmodel, dtime,   &
      &      h_offset, lvt_fract, dpr_fract, dsoil_lvt,dsoil_dpr,          &
      &      field_cap_lvt, field_cap_dpr, wsat_lvt, wsat_dpr,   &
      &      pwp_lvt, pwp_dpr, mpot_sat_lvt, mpot_sat_dpr,       &
      &      k_sat_lvt, k_sat_dpr, bclapp_lvt, bclapp_dpr, ps_index_lvt, ps_index_dpr,             &
      &      ice_impede_lvt, ice_impede_dpr, k_lvt, k_dpr, mpot_lvt, mpot_dpr,                     &
      &      ice_soil_lvt, ice_soil_dpr, wtr_scool_lvt, wtr_scool_dpr, wtr_pond_lvt, wtr_pond_dpr, &
      &      wtr_latflow_lvt, wtr_latflow_dpr, wtr_soil_lvt, wtr_soil_dpr, diffus_lvt, diffus_dpr  &
      &      )

    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    INTEGER, INTENT(in)     :: nsoil,               & !< number of soil layers
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
      &                        ice_impede_lvt(:),   & !< Soil layer ice impedance factor on source tile []
      &                        ice_impede_dpr(:),   & !< Soil layer ice impedance factor on target tile []
      &                        k_lvt(:),            & !< Soil layer hyd. conductivity on source tile [m s-1]
      &                        k_dpr(:),            & !< Soil layer hyd. conductivity on target tile [m s-1]
      &                        mpot_lvt(:),         & !< Soil layer matric potential on source tile [?]
      &                        mpot_dpr(:),         & !< Soil layer matric potential on target tile [?]
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



!    REAL(wp), PARAMETER ::     length_scale = 1.0_wp,        & ! assumed diameter of depressions (could be input field / namelist para)
!      &                        h_offset     = 0.1_wp,        & ! hight difference between ridge and depression (should later be input field)
!      &                        ck           = 8.0_wp           ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    REAL(wp) ::     length_scale,    & ! assumed diameter of depressions (could be input field / namelist para)
!      &             h_offset,        & ! hight difference between ridge and depression (should later be input field)
      &             ck                 ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    nrgrd         = nsoil * 2 + 3

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

    !Assumption of circular depression with A = 1m2 with in a circular elevation
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

    len_inter = ((4._wp * pi) ** 0.5_wp) * fr_dpr_rel ! base assumption is the circumference of circle with area 1m2 within elevation area x

    DO i = 1, nrgrd
      IF (i_lvt <= nsoil .AND. i_dpr <= nsoil .AND. dsoil_lvt(MIN(i_lvt,nsoil)) > 0._wp &
        & .AND. dsoil_dpr(MIN(MAX(i_dpr,1),nsoil)) > 0._wp) THEN ! Above Bedrock in both tiles
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


      ELSE IF (i_lvt <= nsoil .AND. dsoil_lvt(MIN(i_lvt,nsoil)) > 0._wp ) THEN ! Above Bedrock in lvt tile only (should only be the case when using different soil depths)
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

      ELSE IF (i_dpr > 0 .AND. i_dpr <= nsoil .AND. dsoil_dpr(MIN(i_dpr,nsoil)) > 0._wp ) THEN   ! Above Bedrock in dpr tile only
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

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile snow redistribution between tiles with different surface elevation (elevation/depression)
  !!
  SUBROUTINE calc_snow_redistribution(l_dynsnow,   &
      &      dtime,                              &
      &      h_offset,                           &
      &      lvt_fract, dpr_fract,               &
      &      weq_pond_lvt, weq_pond_dpr,         &
      &      snow_soil_lvt, snow_soil_dpr,       &
      &      snow_new_lvt, snow_new_dpr,         &
      &      snow_lvt, snow_dpr,                 &
      &      snow_dens_lvt, snow_dens_dpr,       &
      &      snow_rdst_lvt,snow_rdst_dpr         &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o, dens_snow_min, dens_snow_max, dens_snow

    LOGICAL,  INTENT(in)    :: l_dynsnow              !< Dynamic snow density

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
      &                        h_offset,            & !< Height difference between tiles
      &                        lvt_fract,           & !< area fraction of water elevation tile
      &                        dpr_fract              !< area fraction of water depression tile

    REAL(wp), INTENT(in)    :: weq_pond_lvt,        & !< water & ice lvl in pond reservoir on lvt
      &                        weq_pond_dpr           !< water & ice lvl in pond reservoir in dpr

    REAL(wp), INTENT(inout) :: snow_soil_lvt,       & !< snow soil (WEQ) on elevation tile
      &                        snow_soil_dpr,       & !< snow soil (WEQ) on depression tile
      &                        snow_new_lvt,        & !< new snow soil (WEQ) on elevation tile
      &                        snow_new_dpr,        & !< new snow soil (WEQ) on depression tile
      &                        snow_lvt,            & !< snow (WEQ; can and soil) on elevation tile
      &                        snow_dpr,            & !< snow (WEQ; can and soil) on depression tile
      &                        snow_dens_lvt,       & !< snow density on elevation tile
      &                        snow_dens_dpr          !< snow density on depression tile

    REAL(wp), INTENT(out)   :: snow_rdst_lvt,       & !< Snow redistribution flux at elevation tile [kg m-2 s-1]
      &                        snow_rdst_dpr          !< Snow redistribution flux at depression tile [kg m-2 s-1]



    ! Local variables ------------------------------------------------------------------------------------

    REAL(wp)                :: fr_lvt_rel,          &
      &                        fr_dpr_rel,          &
      &                        h_offset_snow_fix,   &
      &                        snow_fix_lvt ,       &
      &                        snow_fix_dpr ,       &
      &                        snow_fix_rel_lvt,    &
      &                        snow_fix_rel_dpr,    &
      &                        snow_depth_mob_avg,  &
      &                        snow_depth_mob_lvt,  &
      &                        snow_depth_mob_dpr,  &
      &                        snow_mob_lvt,        &
      &                        snow_mob_dpr,        &
      &                        snow_mob_upd_lvt,    &
      &                        snow_mob_upd_dpr,    &
      &                        snow_dens_max,       &
      &                        snow_dens_min,       &
      &                        snow_dens_fix_lvt,   &
      &                        snow_dens_fix_dpr

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
                      & / snow_fix_rel_lvt
    ELSE
      snow_dens_fix_lvt = snow_dens_lvt
    END IF
    IF (snow_fix_dpr > 0._wp) THEN
      snow_fix_rel_dpr = snow_fix_dpr / snow_soil_dpr
      snow_dens_fix_dpr = (snow_dens_dpr - (snow_dens_min * (1._wp - snow_fix_rel_dpr))) &
                      & / snow_fix_rel_dpr
    ELSE
      snow_dens_fix_dpr = snow_dens_dpr
    END IF

    !! Determine height differences of the top of fixed snow
    h_offset_snow_fix = h_offset &
                    & + weq_pond_lvt + snow_fix_lvt * (rhoh2o / snow_dens_fix_lvt) &
                    & - weq_pond_dpr - snow_fix_dpr * (rhoh2o / snow_dens_fix_dpr)

    ! Calculate snow redistribution for two cases (1) from elevation to depression or (2) from depression to elevation
    ! Transport only applies if in the source tile there is fixed and mobile (new) snow
    IF  ((h_offset_snow_fix > 1.0e-10_wp .AND. snow_fix_lvt > 0.01_wp .AND. snow_mob_lvt > 0._wp) .OR. &
         (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_dpr > 0.01_wp .AND. snow_mob_dpr > 0._wp))THEN

      snow_depth_mob_avg = (snow_mob_lvt * fr_lvt_rel + snow_mob_dpr * fr_dpr_rel) * (rhoh2o / snow_dens_min)

      !! 1. (Fixed) snow surface is lower in derpression and there is fixed and mobile snow on elevation    
      IF  (h_offset_snow_fix > 1.0e-10_wp .AND. snow_fix_lvt > 0.01_wp .AND. snow_mob_lvt > 0._wp)THEN
        IF (snow_depth_mob_avg / fr_dpr_rel <= h_offset_snow_fix ) THEN ! New snow fits into depression
          snow_depth_mob_dpr = snow_depth_mob_avg / fr_dpr_rel
          snow_depth_mob_lvt = 0._wp
        ELSE
          snow_depth_mob_dpr = h_offset_snow_fix  &
                           & + (snow_depth_mob_avg - h_offset_snow_fix * fr_dpr_rel)
          snow_depth_mob_lvt =  snow_depth_mob_avg - h_offset_snow_fix * fr_dpr_rel
        END IF
      !! 2. (Fixed) snow surface is lower on elevation and there is fixed and mobile snow in depression
      ELSE IF (h_offset_snow_fix < -1.0e-10_wp .AND. snow_fix_dpr > 0.01_wp .AND. snow_mob_dpr > 0._wp) THEN
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

      ! Recalculation of up snow weqs 
      snow_mob_upd_dpr = snow_depth_mob_dpr / (rhoh2o / snow_dens_min)
      snow_mob_upd_lvt = snow_depth_mob_lvt / (rhoh2o / snow_dens_min)

      snow_rdst_dpr    = snow_mob_upd_dpr - snow_mob_dpr
      snow_rdst_lvt    = snow_mob_upd_lvt - snow_mob_lvt  
     
      !snow_soil_dpr    = snow_fix_dpr + snow_mob_upd_dpr  
      !snow_soil_lvt    = snow_fix_lvt + snow_mob_upd_lvt

      !snow_dpr         = snow_dpr + snow_rdst_dpr
      !snow_lvt         = snow_lvt + snow_rdst_lvt 
      
      !snow_new_dpr     = snow_mob_upd_dpr
      !snow_new_lvt     = snow_mob_upd_lvt 

      snow_soil_dpr    = snow_soil_dpr + snow_rdst_dpr
      snow_soil_lvt    = snow_soil_lvt + snow_rdst_lvt

      snow_dpr         = snow_dpr + snow_rdst_dpr
      snow_lvt         = snow_lvt + snow_rdst_lvt 
      
      snow_new_dpr     = snow_new_dpr + snow_rdst_dpr
      snow_new_lvt     = snow_new_lvt + snow_rdst_lvt
      
      ! Converting redistribution flux from m/step -> kg m-2 s-1 
      snow_rdst_dpr    = snow_rdst_dpr * rhoh2o / dtime
      snow_rdst_lvt    = snow_rdst_lvt * rhoh2o / dtime

      ! Updating snow densities based on redistributed snow
      IF (snow_rdst_dpr .NE. 0_wp .AND. snow_soil_dpr > 0._wp) THEN         
        snow_dens_dpr = (snow_new_dpr * snow_dens_min &
                    & + snow_fix_dpr * min(snow_dens_fix_dpr, snow_dens_max)) / snow_soil_dpr               
      END IF
    
      IF (snow_rdst_lvt .NE. 0_wp .AND. snow_soil_lvt > 0._wp) THEN          
        snow_dens_lvt = (snow_new_lvt * snow_dens_min &
                    & + snow_fix_lvt * min(snow_dens_fix_lvt, snow_dens_max)) / snow_soil_lvt 
      END IF
    END IF

  END SUBROUTINE calc_snow_redistribution

#endif
END MODULE mo_mix_process
