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

  PUBLIC :: calc_runoff_transport, calc_drainage_transport, calc_intertile_diffusion, &
          & calc_intertile_diffusion_new, calc_intertile_diffusion_reduced
 
  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_process'

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of surface runoff
  !!
  PURE ELEMENTAL SUBROUTINE calc_runoff_transport(            &
    & dtime, source_pond_max, source_fract, target_fract, runoff_source,       &
    & wtr_latflow_source, wtr_latflow_target, weq_pond_source, weq_pond_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        source_pond_max,    & !< maximum surface water storage on source tile
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract          !< area fraction of water target tile
    REAL(wp), INTENT(inout) :: runoff_source,      & !< surface runoff on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        weq_pond_source,    & !< surface water storage on source tile
      &                        weq_pond_target       !< surface water storage on target tile

    REAL(wp)                :: overland_flx
    !> Update target ponds with surface runoff from source tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Surface runoff is scaled with the tile fractions and moved into the target tile's
    !!   pond storage
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    overland_flx    = MAX(0._wp, weq_pond_source + runoff_source * dtime / rhoh2o - source_pond_max)
    weq_pond_source = MIN(source_pond_max, weq_pond_source + runoff_source * dtime / rhoh2o)

    IF (target_fract > 1.0e-10_wp) THEN
      weq_pond_target    = weq_pond_target    + overland_flx * (source_fract / target_fract)
      wtr_latflow_target = wtr_latflow_target + overland_flx * rhoh2o / dtime * (source_fract / target_fract)
      wtr_latflow_source = wtr_latflow_source - overland_flx * rhoh2o / dtime
      runoff_source      = 0._wp
    END IF

  END SUBROUTINE calc_runoff_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile exchange of lateral drainage
  !!
  PURE ELEMENTAL SUBROUTINE calc_drainage_transport(            &
    & dtime, source_fract, target_fract, field_cap_target,      &
    & latdrain_source, drainage_source, wtr_latflow_source,     &
    & wtr_latflow_target, wtr_soil_target &
    )

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in)    :: dtime,              & !< timestep length
      &                        source_fract,       & !< area fraction of water source tile
      &                        target_fract,       & !< area fraction of water target tile
      &                        field_cap_target      !< soil field capacity of target tile
    REAL(wp), INTENT(inout) :: latdrain_source,    & !< lateral drainage on source tile
      &                        drainage_source,    & !< bottom drainage on source tile
      &                        wtr_latflow_source, & !< lateral water flow from source tile
      &                        wtr_latflow_target, & !< lateral water flow towards target tile
      &                        wtr_soil_target       !< soil water storage on target tile

    ! Local variables
    REAL(wp)                :: latflow_act           !< Actual lateral flow during time step
                                                     !! considering the target storage limits

    !> Update target soil with lateral drainage from source tile
    !!   *** DUMMY implementation which later on will be replaced with a physics based process ***
    !!   Lateral drainage is scaled with the tile fractions and moved into the target tile's
    !!   soil storage. It flows into the same layer as it originates from on the source tile.
    !!   Excess drainage is added to the total drainage of the source tile.
    !! @todo If a flow lag is required it needs to be implemented here, probably with an
    !!   temporary storage on the source tile, which then also needs to be part of the global
    !!   water balance

    IF (target_fract > 1.0e-10_wp) THEN
      ! Substract diagnosed lateral drainage from total drainage
      drainage_source    = drainage_source - latdrain_source
      ! Scale and add all lateral drainage to target tile soil layer
      latflow_act        = latdrain_source * (source_fract / target_fract)
      wtr_soil_target    = wtr_soil_target + latflow_act * dtime / rhoh2o
      ! Reduce lateral flow by upper soil water limit and correct for any excess soil moisture
      !pdv >
      !latflow_act        = latflow_act - MAX(wtr_soil_target - field_cap_target, 0._wp) * rhoh2o / dtime
      !wtr_soil_target    = MIN(wtr_soil_target, field_cap_target)
      !pdv <
      wtr_latflow_target = wtr_latflow_target + latflow_act
      ! Rescale the actually transported water back to the source tile
      latflow_act        = latflow_act * (target_fract / source_fract)
      wtr_latflow_source = wtr_latflow_source - latflow_act
      ! Add any residual lateral flow to the total drainage of the source tile
      !pdv >
      !drainage_source    = drainage_source + MAX(latdrain_source - latflow_act, 0._wp)
      !pdv <
      latdrain_source    = 0._wp
    END IF

  END SUBROUTINE calc_drainage_transport

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles
  !!
  PURE ELEMENTAL SUBROUTINE calc_intertile_diffusion(           &
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

  END SUBROUTINE calc_intertile_diffusion


  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles with different surface elevation (elevation/depression)
  !!
  SUBROUTINE calc_intertile_diffusion_new(nsoil, i_shm, dtime,   &
      &      lvt_fract, dpr_fract, dsoil_lvt,dsoil_dpr,          &
      &      field_cap_lvt, field_cap_dpr, wsat_lvt, wsat_dpr,   &
      &      pwp_lvt, pwp_dpr, mpot_sat_lvt, mpot_sat_dpr,       &
      &      k_sat_lvt, k_sat_dpr, bclapp_lvt, bclapp_dpr, ps_index_lvt, ps_index_dpr,             &
      &      ice_soil_lvt, ice_soil_dpr, wtr_scool_lvt, wtr_scool_dpr, weq_pond_lvt, weq_pond_dpr, &
      &      wtr_latflow_lvt, wtr_latflow_dpr, wtr_soil_lvt, wtr_soil_dpr         &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: get_conductivies_diffusivities_mix

    INTEGER, INTENT(in)     :: nsoil,               & !< number of soil layers
                               i_shm                  !< Choice of soil hydrological model

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
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


    REAL(wp), INTENT(inout) :: weq_pond_lvt,    & !< lateral water flow from elevation tile
      &                        weq_pond_dpr,    & !< lateral water flow towards depression tile
      &                        wtr_latflow_lvt, & !< lateral water flow from elevation tile
      &                        wtr_latflow_dpr    !< lateral water flow towards depression tile

    REAL(wp), INTENT(inout) :: wtr_soil_lvt(:),     & !< soil water storage on elevation tile
      &                        wtr_soil_dpr(:)        !< soil water storage on depression tile



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: varindex, nrgrd, i, i_lvt, i_dpr

    REAL(wp)                :: fr_lvt_rel,       &
      &                        fr_dpr_rel,       &
      &                        pond_flx_lvt,     &
      &                        pond_flx_dpr,     &
      &                        weq_pond_total,   &
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

    REAL(wp)                :: diffus_lvt(nsoil),        & !< Diffusion flux at elevation tile per time step
      &                        diffus_dpr(nsoil),        & !< Diffusion flux at depression tile per time step
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



!    REAL(wp), PARAMETER ::     Pi           = 3.1415927_wp,  &
!      &                        length_scale = 1.0_wp,        & ! assumed diameter of depressions (could be input field / namelist para)
!      &                        h_offset     = 0.1_wp,        & ! hight difference between ridge and depression (should later be input field)
!      &                        ck           = 8.0_wp           ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    REAL(wp) ::     Pi,              &
      &             length_scale,    & ! assumed diameter of depressions (could be input field / namelist para)
      &             h_offset,        & ! hight difference between ridge and depression (should later be input field)
      &             ck                 ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)


    Pi           = 3.1415927_wp  
    h_offset     = 0.1_wp           ! hight difference between ridge and depression (should later be input field)

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

    !Assumption of circular depression with A = 1m² with in a circular elevation
    !length_scale = ((1._wp/Pi)**0.5_wp+((1._wp/fr_dpr_rel)/Pi)**0.5_wp)/2._wp ! From center depression to half of outer rim
    length_scale = ((1._wp/fr_dpr_rel)/Pi)**0.5_wp                      ! From center to outer rim  

    ! Fist calculate the standing water at the surface of both tiles ... with the assumption that runoff is instant
    weq_pond_total = weq_pond_lvt * fr_lvt_rel + weq_pond_dpr * fr_dpr_rel

    IF (weq_pond_total / fr_dpr_rel <= h_offset ) THEN
      pond_flx_dpr =   weq_pond_total / fr_dpr_rel - weq_pond_dpr
      pond_flx_lvt = - weq_pond_lvt
      weq_pond_dpr =   weq_pond_total / fr_dpr_rel
      weq_pond_lvt =   0._wp
    ELSE
      pond_flx_dpr =   h_offset  + (weq_pond_total - h_offset * fr_dpr_rel) - weq_pond_dpr
      pond_flx_lvt =   weq_pond_total - h_offset * fr_dpr_rel - weq_pond_lvt
      weq_pond_dpr =   h_offset  + (weq_pond_total - h_offset * fr_dpr_rel)
      weq_pond_lvt =   weq_pond_total - h_offset * fr_dpr_rel
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


    varindex      = 1
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & ice_impede_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & ice_impede_dpr(:))


    varindex      = 2
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & k_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & k_dpr(:))


    varindex      = 4
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & mpot_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & mpot_dpr(:))


    
    ! Start setting up the interpolated grid
    i_lvt = 1
    z_lvt = dsoil_lvt(i_lvt)
    i_dpr = 1
    IF (h_offset > 0._wp) THEN
      i_dpr = 0
      IF (h_offset > weq_pond_dpr) THEN
        z_dpr = h_offset - weq_pond_dpr
      ELSE
        z_dpr = h_offset
      END IF
    ELSE
      z_dpr = dsoil_dpr(i_dpr)
    END IF

    z_rgrd_top = 0._wp
    z_rgrd_bot = MIN(z_lvt, z_dpr)

    len_inter = ((4._wp * Pi) ** 0.5_wp) * fr_dpr_rel ! base assumption is the circumference of circle with area 1m² within elevation area x
    
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

        IF (z_rgrd_bot <= h_offset - weq_pond_dpr) THEN ! Still above water level in depression
          ! Water flows out of the soil in the elevation tile into the pond reservoir of the depression tile
          wtr_flx_lat       = k_lvt_rgrd * a_inter
          wtr_avail         = MAX(0._wp,(wtr_soil_lvt_rgrd - pwp_lvt_rgrd) * fr_lvt_rel)
          wtr_lat           = MIN(wtr_flx_lat * dtime, wtr_avail)
          wtr_flx_lat       = wtr_lat / dtime
          weq_pond_dpr      = weq_pond_dpr + wtr_lat / fr_dpr_rel
          wtr_soil_lvt_rgrd = wtr_soil_lvt_rgrd - wtr_lat / fr_lvt_rel
          wtr_soil_dpr_rgrd = 0.0_wp

        ELSE IF (z_rgrd_bot <= h_offset ) THEN ! Within the ponded surface water in the depression
          ! Water flows from the pond reservoir of the depression tile into the soil of the elevation tile
!          wtr_flx_lat       = k_lvt_rgrd * a_inter
          wtr_flx_lat       = ((k_lvt_rgrd * k_sat_lvt_rgrd)**0.5_wp)*a_inter
!          free_space        = MAX(0._wp,(wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          free_space        = MAX(0._wp,(field_cap_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          wtr_avail         = MIN(h, weq_pond_dpr) * fr_dpr_rel
          wtr_lat           = MAX(-wtr_flx_lat * dtime, -free_space, -wtr_avail)
          wtr_flx_lat       = wtr_lat / dtime
          weq_pond_dpr      = weq_pond_dpr + wtr_lat / fr_dpr_rel
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


  END SUBROUTINE calc_intertile_diffusion_new

  ! ===============================================================================================================================
  !>
  !! @brief   Compute intertile lateral diffusion between tiles with different surface elevation (elevation/depression)
  !!
  SUBROUTINE calc_intertile_diffusion_reduced(nsoil, i_shm, dtime,   &
      &      lvt_fract, dpr_fract, dsoil_lvt,dsoil_dpr,          &
      &      field_cap_lvt, field_cap_dpr, wsat_lvt, wsat_dpr,   &
      &      pwp_lvt, pwp_dpr, mpot_sat_lvt, mpot_sat_dpr,       &
      &      k_sat_lvt, k_sat_dpr, bclapp_lvt, bclapp_dpr, ps_index_lvt, ps_index_dpr,             &
      &      ice_soil_lvt, ice_soil_dpr, wtr_scool_lvt, wtr_scool_dpr, weq_pond_lvt, weq_pond_dpr, &
      &      wtr_latflow_lvt, wtr_latflow_dpr, wtr_soil_lvt, wtr_soil_dpr         &
      &      )

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: get_conductivies_diffusivities_mix

    INTEGER, INTENT(in)     :: nsoil,               & !< number of soil layers
                               i_shm                  !< Choice of soil hydrological model

    REAL(wp), INTENT(in)    :: dtime,               & !< timestep length
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


    REAL(wp), INTENT(inout) :: weq_pond_lvt,    & !< lateral water flow from elevation tile
      &                        weq_pond_dpr,    & !< lateral water flow towards depression tile
      &                        wtr_latflow_lvt, & !< lateral water flow from elevation tile
      &                        wtr_latflow_dpr    !< lateral water flow towards depression tile

    REAL(wp), INTENT(inout) :: wtr_soil_lvt(:),     & !< soil water storage on elevation tile
      &                        wtr_soil_dpr(:)        !< soil water storage on depression tile



    ! Local variables ------------------------------------------------------------------------------------
    INTEGER                 :: varindex, nrgrd, i, i_lvt, i_dpr

    REAL(wp)                :: fr_lvt_rel,       &
      &                        fr_dpr_rel,       &
      &                        weq_pond_total,   &
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

    REAL(wp)                :: diffus_lvt(nsoil),        & !< Diffusion flux at elevation tile per time step
      &                        diffus_dpr(nsoil),        & !< Diffusion flux at depression tile per time step
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



!    REAL(wp), PARAMETER ::     Pi           = 3.1415927_wp,  &
!      &                        length_scale = 1.0_wp,        & ! assumed diameter of depressions (could be input field / namelist para)
!      &                        h_offset     = 0.1_wp,        & ! hight difference between ridge and depression (should later be input field)
!      &                        ck           = 8.0_wp           ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)

    REAL(wp) ::     Pi,              &
      &             length_scale,    & ! assumed diameter of depressions (could be input field / namelist para)
      &             h_offset,        & ! hight difference between ridge and depression (should later be input field)
      &             ck                 ! For matric potential formula from Zhang 2007 (https://doi.org/10.1175/JHM605.1)


    Pi           = 3.1415927_wp  
    h_offset     = 0.1_wp           ! hight difference between ridge and depression (should later be input field)

    nrgrd         = nsoil * 2 + 3

    ! Initilize fields for accumulating regridded output
    wtr_soil_lvt_updt(:) = 0._wp
    wtr_soil_dpr_updt(:) = 0._wp

    diffus_lvt(:)        = 0._wp
    diffus_dpr(:)        = 0._wp

    ! Determine relative fractions from absolute fractions -- don't forget to scale fluxes afterwards
    fr_lvt_rel = lvt_fract / (lvt_fract+dpr_fract)
    fr_dpr_rel = dpr_fract / (lvt_fract+dpr_fract)

    !Assumption of circular depression with A = 1m² with in a circular elevation
    length_scale = ((1._wp/fr_dpr_rel)/Pi)**0.5_wp                      ! From center to outer rim  

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


    varindex      = 1
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & ice_impede_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & ice_impede_dpr(:))


    varindex      = 2
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & k_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & k_dpr(:))


    varindex      = 4
    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_lvt(:), wtr_soil_lvt(:), ice_soil_lvt(:), wtr_scool_lvt(:),   &
                    & wsat_wi_lvt(:), field_cap_wi_lvt(:), pwp_wi_lvt(:), k_sat_lvt(:),   &
                    & mpot_sat_lvt(:), bclapp_lvt(:), ps_index_lvt(:),                    &
                    & mpot_lvt(:))

    call get_conductivies_diffusivities_mix(i_shm,varindex,                                   &
                    & dsoil_dpr(:), wtr_soil_dpr(:), ice_soil_dpr(:), wtr_scool_dpr(:),   &
                    & wsat_wi_dpr(:), field_cap_wi_dpr(:), pwp_wi_dpr(:), k_sat_dpr(:),   &
                    & mpot_sat_dpr(:), bclapp_dpr(:), ps_index_dpr(:),                    &
                    & mpot_dpr(:))

    len_inter = ((4._wp * Pi) ** 0.5_wp) * fr_dpr_rel ! base assumption is the circumference of circle with area 1m² within elevation area x
    
    DO i = 1, nsoil
      i_lvt = i
      i_dpr = i
      IF (i_lvt <= nsoil .AND. i_dpr <= nsoil .AND. dsoil_lvt(i_lvt) > 0._wp .AND. dsoil_dpr(MAX(i_dpr,1)) > 0._wp) THEN ! Above Bedrock in both tiles
        h                   = dsoil_lvt(i)
        a_inter             = len_inter * h
    
        wtr_soil_lvt_rgrd   = wtr_soil_lvt(i_lvt)      * h / dsoil_lvt(i_lvt)
        ice_soil_lvt_rgrd   = ice_soil_lvt(i_lvt)      * h / dsoil_lvt(i_lvt)
        wsat_lvt_rgrd       = wsat_wi_lvt(i_lvt)       * h / dsoil_lvt(i_lvt)
        field_cap_lvt_rgrd  = field_cap_wi_lvt(i_lvt)  * h / dsoil_lvt(i_lvt)
        pwp_lvt_rgrd        = pwp_wi_lvt(i_lvt)        * h / dsoil_lvt(i_lvt)
        mpot_lvt_rgrd       = mpot_lvt(i_lvt)
        k_sat_lvt_rgrd      = k_sat_lvt(i_lvt)
        k_lvt_rgrd          = k_lvt(i_lvt)

        wtr_soil_dpr_rgrd = wtr_soil_dpr(i_dpr)        * h / dsoil_dpr(i_dpr)
        ice_soil_dpr_rgrd = ice_soil_dpr(i_dpr)        * h / dsoil_dpr(i_dpr)
        wsat_dpr_rgrd     = wsat_wi_dpr(i_dpr)         * h / dsoil_dpr(i_dpr)
        field_cap_dpr_rgrd= field_cap_wi_dpr(i_dpr)    * h / dsoil_dpr(i_dpr)
        pwp_dpr_rgrd      = pwp_wi_dpr(i_dpr)          * h / dsoil_dpr(i_dpr)
        mpot_dpr_rgrd     = mpot_dpr(i_dpr)
        k_dpr_rgrd        = k_dpr(i_dpr)

        k_inter           =(k_dpr_rgrd * k_lvt_rgrd)**0.5_wp
        d_mpot            = mpot_dpr_rgrd - mpot_lvt_rgrd
        wtr_flx_lat       = - k_inter  * d_mpot / length_scale *  a_inter

        wtr_rel_qlbrm     = ((wtr_soil_lvt_rgrd+ice_soil_lvt_rgrd) * fr_lvt_rel       &
                          & +  (wtr_soil_dpr_rgrd+ice_soil_dpr_rgrd) * fr_dpr_rel)      &
                          & / (wsat_lvt_rgrd * fr_lvt_rel +wsat_dpr_rgrd * fr_dpr_rel )

        IF (wtr_flx_lat > 0._wp) THEN ! Flux is going from elevation to depression 

          wtr_2_qlbrm     = MAX(0._wp, &
                        &  (wtr_rel_qlbrm * wsat_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
!          free_space      = MAX(0._wp,(wsat_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
          free_space      = MAX(0._wp,(field_cap_dpr_rgrd - wtr_soil_dpr_rgrd - ice_soil_dpr_rgrd) * fr_dpr_rel)
          wtr_avail       = MAX(0._wp,(wtr_soil_lvt_rgrd - pwp_lvt_rgrd) * fr_lvt_rel)
          wtr_lat         = MIN(wtr_flx_lat * dtime, free_space, wtr_avail, wtr_2_qlbrm)
          wtr_lat         = MAX(0._wp, wtr_lat)
          wtr_flx_lat     = wtr_lat / dtime

        ELSE ! Flux is going from depression to elevation 

          wtr_2_qlbrm     = MAX(0._wp, &
                         & (wtr_rel_qlbrm  * wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
!          free_space      = MAX(0._wp, (wsat_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          free_space      = MAX(0._wp, (field_cap_lvt_rgrd - wtr_soil_lvt_rgrd - ice_soil_lvt_rgrd) * fr_lvt_rel)
          wtr_avail       = MAX(0._wp, (wtr_soil_dpr_rgrd - pwp_dpr_rgrd) * fr_dpr_rel)
          wtr_lat         = MAX(wtr_flx_lat * dtime, -free_space, -wtr_avail, -wtr_2_qlbrm)
          wtr_lat         = MIN(0._wp, wtr_lat)
          wtr_flx_lat     = wtr_lat / dtime

        END IF

        wtr_soil_lvt_rgrd = wtr_soil_lvt_rgrd - wtr_lat / fr_lvt_rel          
        wtr_soil_dpr_rgrd = wtr_soil_dpr_rgrd + wtr_lat / fr_dpr_rel

        wtr_soil_lvt_updt(i_lvt)  = wtr_soil_lvt_updt(i_lvt) + wtr_soil_lvt_rgrd
        diffus_lvt(i_lvt) = diffus_lvt(i_lvt) - wtr_lat / fr_lvt_rel

        wtr_soil_dpr_updt(i_dpr)  = wtr_soil_dpr_updt(i_dpr) + wtr_soil_dpr_rgrd
        diffus_dpr(i_dpr) = diffus_dpr(i_dpr) + wtr_lat / fr_dpr_rel

      ELSE

        wtr_soil_lvt_rgrd = 0.0_wp
        wtr_soil_dpr_rgrd = 0.0_wp

        wtr_soil_lvt_updt(i_lvt)  = 0.0_wp
        diffus_lvt(i_lvt) = 0.0_wp

        wtr_soil_dpr_updt(i_dpr)  = 0.0_wp
        diffus_dpr(i_dpr) = 0.0_wp

      END IF
    END DO
    wtr_soil_lvt(:)  = wtr_soil_lvt_updt(:)
    wtr_soil_dpr(:)  = wtr_soil_dpr_updt(:)

!    wtr_latflow_lvt = wtr_latflow_lvt + SUM(diffus_lvt(:)) * (lvt_fract+dpr_fract) * rhoh2o / dtime
!    wtr_latflow_dpr = wtr_latflow_dpr + SUM(diffus_dpr(:)) * (lvt_fract+dpr_fract) * rhoh2o / dtime

    wtr_latflow_lvt = wtr_latflow_lvt + SUM(diffus_lvt(:)) * rhoh2o / dtime
    wtr_latflow_dpr = wtr_latflow_dpr + SUM(diffus_dpr(:)) * rhoh2o / dtime


  END SUBROUTINE calc_intertile_diffusion_reduced



#endif
END MODULE mo_mix_process
