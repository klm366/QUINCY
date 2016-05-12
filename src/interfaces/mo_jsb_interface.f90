!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_util_string,         ONLY: int2string
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_jsbach
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_Task_queue
  USE mo_hsm_class,           ONLY: t_Hsm, t_Message, t_State
!  USE mo_jsb_controller,      ONLY: Run_process_task
  USE mo_atmland_interface,   ONLY: update_atm2land, update_land2atm
  USE mo_jsb_parallel_iface,  ONLY: get_my_global_mpi_id

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_interface, interface_full, jsbach_start_timestep, jsbach_finish_timestep

  INTERFACE jsbach_interface
    MODULE PROCEDURE interface_full
    MODULE PROCEDURE interface_inquire
  END INTERFACE jsbach_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_interface'

CONTAINS

  SUBROUTINE interface_full( &
    & model_id, iblk, ics, ice, dtime, steplen,                    &
    & t_air, q_air, rain, snow, wind_air, wind_10m,                &
    & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, &
    & press_srf, drag_srf,                                         &
    & t_acoef, t_bcoef, q_acoef, q_bcoef,                          &
    & pch, cos_zenith_angle,                                       &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                           &
    & fact_q_air, fact_qsat_srf,                                   &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,&
    & rough_h_srf, rough_m_srf,                                    &
    & tte_corr,                                                    &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,          &
    ! For lakes:
    & latent_hflx_wtr, sensible_hflx_wtr,                          & !< in
    & t_wtr,                                                       & !< out
    & latent_hflx_ice, sensible_hflx_ice,                          & !< in
    & t_ice, ice_fract_lake,                                       & !< out
    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr,           & !< out
    & alb_nir_dif_wtr, albedo_ice                                  & !< out
    & )

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_full'

    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice

    REAL(wp), INTENT(in) :: dtime, steplen
    REAL(wp), INTENT(in) ::    &
      & t_air             (:), &
      & q_air             (:), &
      & rain              (:), &
      & snow              (:), &
      & wind_air          (:), &
      & wind_10m          (:), &
      & lw_srf_down       (:), &
      & swvis_srf_down    (:), &
      & swnir_srf_down    (:), &
      & swpar_srf_down    (:), &
      & press_srf         (:), &
      & drag_srf          (:), &
      & t_acoef           (:), &
      & t_bcoef           (:), &
      & q_acoef           (:), &
      & q_bcoef           (:), &
      & pch               (:), &
      & cos_zenith_angle  (:), &
      ! For lakes:
      & latent_hflx_wtr   (:), &
      & sensible_hflx_wtr (:), &
      & latent_hflx_ice   (:), &
      & sensible_hflx_ice (:)
    REAL(wp), INTENT(out) ::   &
      & t_srf             (:), &
      & t_eff_srf         (:), &
      & qsat_srf          (:), &
      & s_srf             (:), &
      & fact_q_air        (:), &
      & fact_qsat_srf     (:), &
      & evapotrans        (:), &
      & latent_hflx       (:), &
      & sensible_hflx     (:), &
      & grnd_hflx         (:), &
      & grnd_hcap         (:), &
      & rough_h_srf       (:), &
      & rough_m_srf       (:), &
      & tte_corr          (:), &
      & alb_vis_dir       (:), &
      & alb_nir_dir       (:), &
      & alb_vis_dif       (:), &
      & alb_nir_dif       (:), &
      & t_wtr             (:), &
      & t_ice             (:), &
      & ice_fract_lake    (:), &
      & alb_vis_dir_wtr   (:), &
      & alb_vis_dif_wtr   (:), &
      & alb_nir_dir_wtr   (:), &
      & alb_nir_dif_wtr   (:), &
      & albedo_ice        (:)

    ! Local variables
    !
    CLASS(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    INTEGER :: nc, i

    CALL timer_start(timer_jsbach)

    nc = ice - ics + 1

    model => Get_model(model_id)

!!$    WRITE(message_text, '(A,I4,A,I6,A,I6)') 'Updating '//TRIM(model%name)// &
!!$      &                 ' for: block=', iblk, ', ics=', ics, ', ice=', ice
!!$    CALL message(TRIM(routine), message_text, all_print=.TRUE.)
!!$    p_pe = get_my_global_mpi_id()
!!$    CALL message( TRIM(routine), 'Updating '//TRIM(model%name)//' on PE '//int2string(p_pe)// &
!!$      &                          ', nblk='//int2string(iblk), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

    ! Update options
    CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, dtime=dtime, steplen=steplen)

    ! Put atmospheric forcing into the atm2land memory.
    tile => model%Get_top_tile()
    CALL update_atm2land(                           &
      & tile,                                       &
      & model%options,                              &
      & t_air              = t_air(:),              &
      & q_air              = q_air(:),              &
      & rain               = rain(:),               &
      & snow               = snow(:),               &
      & wind_air           = wind_air(:),           &
!!$      & wind_10m          = wind_10m(:),          &
      & wind_10m           = wind_air(:)*0.8_wp,    & ! Temporary
      & lw_srf_down        = lw_srf_down(:),        &
      & swvis_srf_down     = swvis_srf_down(:),     &
      & swnir_srf_down     = swnir_srf_down(:),     &
      & swpar_srf_down     = swpar_srf_down(:),     &
      & press_srf          = press_srf(:),          &
      & drag_srf           = drag_srf(:),           &
      & t_acoef            = t_acoef(:),            &
      & t_bcoef            = t_bcoef(:),            &
      & q_acoef            = q_acoef(:),            &
      & q_bcoef            = q_bcoef(:),            &
      & pch                = pch(:),                &
      & cos_zenith_angle   = cos_zenith_angle(:),   &
      ! For lakes:
      & latent_hflx_wtr    = latent_hflx_wtr(:),    &
      & sensible_hflx_wtr  = sensible_hflx_wtr(:),  &
      & latent_hflx_ice    = latent_hflx_ice(:),    &
      & sensible_hflx_ice  = sensible_hflx_ice(:)   &
      & )
    NULLIFY(tile)

    ! Run tasks in queue
    CALL model%Run_tasks()


!!$    ! Update wet skin and snow fractions
!!$    CALL  update_snow_and_skin_fraction(tile, iblk, ics, ice, nc, dtime, steplen)
!!$
!!$    ! Compute new surface temperatures and saturation surface specific humidity
!!$    CALL update_surface_energy(tile, iblk, ics, ice, nc, dtime, steplen)
!!$
!!$    ! Compute surface fluxes of heat and moisture based on the new surface temperature and qsat
!!$    CALL update_surface_fluxes(tile, iblk, ics, ice, nc, steplen)
!!$
!!$    ! First part of vegetation that is applied before the soil is updated,
!!$    ! in particular:
!!$    ! 1. Compute canopy conductance not limited by water stress
!!$    ! 2. Compute soil moisture in root zone
!!$    ! 3. Compute water stress factor
!!$    ! 4. Compute canopy conductance limited by water stress
!!$    ! 5. Compute skin reservoir and snow fractions on canopy
!!$    !    For the moment, this is done in update_soil_1 to be compatible with echam6
!!$    CALL update_vegetation_1(tile, iblk, ics, ice, nc)
!!$
!!$    ! Soil hydrology
!!$    CALL update_soil_1(tile, iblk, ics, ice, nc, dtime, steplen)
!!$
!!$#ifndef __NO_JSBACH_HD__
!!$    dsl4jsb_Get_config(HD)
!!$
!!$    IF (tile%Is_process_active(HD)) THEN
!!$      CALL update_hd(tile, iblk, ics, ice, nc, steplen)
!!$    END IF
!!$#endif
!!$
!!$    ! Re-calculate wet skin and snow fractions with updated moisture parameters
!!$    CALL  update_snow_and_skin_fraction(tile, iblk, ics, ice, nc, dtime, steplen)
!!$
!!$    ! Update effective evaporation weights for qsat and q_air
!!$    CALL update_humidity_scaling(tile, iblk, ics, ice, nc)
!!$
!!$    ! Second part of vegetation that is applied after the soil has been updated,
!!$    ! in particular:
!!$    ! 1. Update LAI and foliage projected cover
!!$    ! 2. Compute canopy albedo
!!$    CALL update_vegetation_2(tile, iblk, ics, ice, nc)
!!$
!!$    ! Second part of soil that is applied after the soil energy and water balance has been updated,
!!$    ! in particular:
!!$    ! 1. Compute ground albedo (at the moment this is done together with canopy in update_vegetation_2)
!!$    CALL update_soil_2(tile, iblk, ics, ice, nc)
!!$
!!$    ! Update roughness length
!!$    ! At the moment, roughness length for momentum is read in from ini file, only for heat is changed here
!!$    CALL update_roughness(tile, iblk, ics, ice, nc)
!!$
!!$    ! Apply Asselin time filter for surface temperature (t_srf)
!!$    CALL update_t_srf_asselin(tile, iblk, ics, ice, nc)
!!$
    ! Finally, give the required variables back to the atmosphere
    tile => model%Get_top_tile()
    CALL update_land2atm(                     &
      & tile,                                 &
      & model%options,                        &
      & t_srf           = t_srf(:),           &
      & t_eff_srf       = t_eff_srf(:),       &
      & qsat_srf        = qsat_srf(:),        &
      & s_srf           = s_srf(:),           &
      & fact_q_air      = fact_q_air(:),      &
      & fact_qsat_srf   = fact_qsat_srf(:),   &
      & evapotrans      = evapotrans(:),      &
      & latent_hflx     = latent_hflx(:),     &
      & sensible_hflx   = sensible_hflx(:),   &
      & grnd_hflx       = grnd_hflx(:),       &
      & grnd_hcap       = grnd_hcap,          &
      & rough_h_srf     = rough_h_srf(:),     &
      & rough_m_srf     = rough_m_srf(:),     &
      & tte_corr        = tte_corr(:),        &
      & alb_vis_dir     = alb_vis_dir(:),     &
      & alb_nir_dir     = alb_nir_dir(:),     &
      & alb_vis_dif     = alb_vis_dif(:),     &
      & alb_nir_dif     = alb_nir_dif(:),     &
      & t_wtr           = t_wtr(:),           &
      & t_ice           = t_ice(:),           &
      & ice_fract_lake  = ice_fract_lake(:),  &
      & alb_vis_dir_wtr = alb_vis_dir_wtr(:), &
      & alb_vis_dif_wtr = alb_vis_dif_wtr(:), &
      & alb_nir_dir_wtr = alb_nir_dir_wtr(:), &
      & alb_nir_dif_wtr = alb_nir_dif_wtr(:), &
      & albedo_ice      = albedo_ice(:)       &
      & )
    NULLIFY(tile)

!!$    IF (iblk == 1)                                                              &
!!$      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//          &
!!$      &                             ' on PE '//int2string(p_pe)//' completed.', &
!!$      &              all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

    CALL timer_stop(timer_jsbach)

  END SUBROUTINE interface_full

  SUBROUTINE interface_inquire( model_id, iblk, ics, ice,         &
    & t_srf, t_eff_srf, qsat_srf, zh_srf, zm_srf, alb_vis, alb_nir)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_inquire'

    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice

    REAL(wp), INTENT(inout) :: &
      & t_srf    (:), &
      & t_eff_srf(:), &
      & qsat_srf (:), &
      & zh_srf   (:), &
      & zm_srf   (:), &
      & alb_vis  (:), &
      & alb_nir  (:)

  END SUBROUTINE interface_inquire

  SUBROUTINE jsbach_start_timestep(model_id)
 
    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsbach_start_timestep'

    INTEGER,  INTENT(in) ::  model_id

    TYPE(t_jsb_model), POINTER :: model

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    model => Get_model(model_id)

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_start_timestep

  SUBROUTINE jsbach_finish_timestep(model_id, steplen)

#ifndef __NO_JSBACH_HD__
    USE mo_interface_hd_ocean, ONLY: interface_hd_ocean
    USE mo_hd_interface,       ONLY: hd_check_water_budget, hd_lateral_flow
#endif

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsbach_finish_timestep'

    INTEGER, INTENT(in)  :: model_id
    REAL(wp), INTENT(in) :: steplen

    TYPE(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    model => Get_model(model_id)
    tile => model%Get_top_tile()

#ifndef __NO_JSBACH_HD__
      CALL hd_lateral_flow(tile, model%options)
      CALL hd_check_water_budget(tile, model%options)

      CALL interface_hd_ocean(tile, model%options)
#endif

  END SUBROUTINE jsbach_finish_timestep

#endif
END MODULE mo_jsb_interface
