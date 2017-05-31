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

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_util_string,         ONLY: int2string
  USE mo_jsb_parallel,        ONLY: Get_omp_thread

  USE mo_jsb_control,         ONLY: debug_on, l_timer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_jsbach
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_Task_queue
  USE mo_hsm_class,           ONLY: t_Hsm, t_Message, t_State
!  USE mo_jsb_controller,      ONLY: Run_process_task
  USE mo_atmland_interface,   ONLY: update_atm2land, update_land2atm
  USE mo_jsb_parallel,        ONLY: get_my_global_mpi_id

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_interface, interface_full, jsbach_start_timestep, jsbach_finish_timestep

  INTERFACE jsbach_interface
    MODULE PROCEDURE interface_full
    ! MODULE PROCEDURE interface_standalone
    MODULE PROCEDURE interface_inquire
  END INTERFACE jsbach_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_interface'

CONTAINS

  SUBROUTINE interface_full( &
    & model_id,         &
    & iblk, ics, ice, dtime, steplen, &
    & t_air,            &
    & q_air,            &
    & rain,             &  
    & snow,             &
    & wind_air,         &
    & wind_10m,         &
    & lw_srf_down,      &
    & swvis_srf_down,   &
    & swnir_srf_down,   &
    & swpar_srf_down,   &
    & frac_par_diffuse, &
    & press_srf,        &
    & drag_srf,         &
    & t_acoef,          &
    & t_bcoef,          &
    & q_acoef,          &
    & q_bcoef,          &
    & pch,              &
    & cos_zenith_angle, &
    & CO2_air,          &
    & t_srf,            &
    & t_eff_srf,        &
    & qsat_srf,         &
    & s_srf,            &
    & fact_q_air,       &
    & fact_qsat_srf,    &
    & evapotrans,       &
    & latent_hflx,      &
    & sensible_hflx,    &
    & grnd_hflx,        &
    & grnd_hcap,        &
    & rough_h_srf,      &
    & rough_m_srf,      &
    & q_snocpymlt,      &
    & alb_vis_dir,      &
    & alb_nir_dir,      &
    & alb_vis_dif,      &
    & alb_nir_dif,      &
    ! For lakes (optional):
    & drag_wtr, drag_ice,                                    & !< in
    & t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,    & !< in
    & t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice,    & !< in
    & t_lwtr, evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr, & !< out
    & albedo_lwtr,                                           & !< out
    & t_lice, evapo_ice, latent_hflx_ice, sensible_hflx_ice, & !< out
    & albedo_lice,                                           & !< out
    & ice_fract_lake,                                        & !< out
    ! For standalone jsbach
    & evapopot                                               & !< out; optional
    & )

    USE mo_jsb_test,          ONLY: write_interface_vars, write_interface_variables

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_full'

    ! Arguments
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
      & frac_par_diffuse  (:), &
      & press_srf         (:), &
      & drag_srf          (:), &
      & t_acoef           (:), &
      & t_bcoef           (:), &
      & q_acoef           (:), &
      & q_bcoef           (:), &
      & pch               (:), &
      & cos_zenith_angle  (:), &
      & CO2_air           (:)
    ! For lakes:
    REAL(wp), INTENT(in), OPTIONAL ::    &
      & drag_wtr          (:), &
      & drag_ice          (:), &
      & t_acoef_wtr       (:), &
      & t_bcoef_wtr       (:), &
      & q_acoef_wtr       (:), &
      & q_bcoef_wtr       (:), &
      & t_acoef_ice       (:), &
      & t_bcoef_ice       (:), &
      & q_acoef_ice       (:), &
      & q_bcoef_ice       (:)
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
      & q_snocpymlt       (:), &
      & alb_vis_dir       (:), &
      & alb_nir_dir       (:), &
      & alb_vis_dif       (:), &
      & alb_nir_dif       (:)
    ! For lakes:
    REAL(wp), INTENT(out), OPTIONAL ::   &
      & t_lwtr            (:), &
      & evapo_wtr         (:), &
      & latent_hflx_wtr   (:), &
      & sensible_hflx_wtr (:), &
      & albedo_lwtr       (:), &
      & t_lice            (:), &
      & evapo_ice         (:), &
      & latent_hflx_ice   (:), &
      & sensible_hflx_ice (:), &
      & albedo_lice       (:), &
      & ice_fract_lake    (:)
    ! For standalone jsbach:
    REAL(wp), INTENT(out), OPTIONAL :: evapopot(:)

    ! Local variables
    CLASS(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    INTEGER :: nc
    INTEGER :: no_omp_thread

    ! ---------------------------
    ! Go

    no_omp_thread = Get_omp_thread()

    IF (write_interface_vars) THEN
      CALL write_interface_variables(ics, ice, iblk, t_air, q_air, rain, snow, wind_air, wind_10m, &
        & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf,    &
        & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle)
    END IF

    IF (l_timer) CALL timer_start(timer_jsbach)

    nc = ice - ics + 1

    model => Get_model(model_id)

!!$    WRITE(message_text, '(A,I4,A,I6,A,I6)') 'Updating '//TRIM(model%name)// &
!!$      &                 ' for: block=', iblk, ', ics=', ics, ', ice=', ice
!!$    CALL message(TRIM(routine), message_text, all_print=.TRUE.)
!!$    p_pe = get_my_global_mpi_id()
!!$    CALL message( TRIM(routine), 'Updating '//TRIM(model%name)//' on PE '//int2string(p_pe)// &
!!$      &                          ', nblk='//int2string(iblk), all_print=.TRUE.)
    IF (debug_on('basic') .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

    ! Update options
    CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, dtime=dtime, steplen=steplen)

    ! Put atmospheric forcing into the a2l memory.
    tile => model%Get_top_tile()
    CALL update_atm2land(                           &
      & tile,                                       &
      & model%options(no_omp_thread),               &
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
      & frac_par_diffuse   = frac_par_diffuse(:),   &
      & press_srf          = press_srf(:),          &
      & drag_srf           = drag_srf(:),           &
      & t_acoef            = t_acoef(:),            &
      & t_bcoef            = t_bcoef(:),            &
      & q_acoef            = q_acoef(:),            &
      & q_bcoef            = q_bcoef(:),            &
      & pch                = pch(:),                &
      & cos_zenith_angle   = cos_zenith_angle(:),   &
      & CO2_air            = CO2_air(:),            &
      ! For lakes:
      & drag_wtr           = drag_wtr,              &
      & drag_ice           = drag_ice,              &
      & t_acoef_wtr        = t_acoef_wtr,           &
      & t_bcoef_wtr        = t_bcoef_wtr,           &
      & q_acoef_wtr        = q_acoef_wtr,           &
      & q_bcoef_wtr        = q_bcoef_wtr,           &
      & t_acoef_ice        = t_acoef_ice,           &
      & t_bcoef_ice        = t_bcoef_ice,           &
      & q_acoef_ice        = q_acoef_ice,           &
      & q_bcoef_ice        = q_bcoef_ice            &
      & )
    NULLIFY(tile)

    ! Run tasks in queue
    CALL model%Run_tasks(debug_on('hsm') .AND. iblk == 1)

    ! Finally, give the required variables back to the atmosphere
    ! Note: evapopot is an optional parameter, so it has to be passed without the colon!
    tile => model%Get_top_tile()
    CALL update_land2atm(                        &
      & tile,                                    &
      & model%options(no_omp_thread),            &
      & t_srf              = t_srf(:),           &
      & t_eff_srf          = t_eff_srf(:),       &
      & qsat_srf           = qsat_srf(:),        &
      & s_srf              = s_srf(:),           &
      & fact_q_air         = fact_q_air(:),      &
      & fact_qsat_srf      = fact_qsat_srf(:),   &
      & evapopot           = evapopot,           &
      & evapotrans         = evapotrans(:),      &
      & latent_hflx        = latent_hflx(:),     &
      & sensible_hflx      = sensible_hflx(:),   &
      & grnd_hflx          = grnd_hflx(:),       &
      & grnd_hcap          = grnd_hcap,          &
      & rough_h_srf        = rough_h_srf(:),     &
      & rough_m_srf        = rough_m_srf(:),     &
      & q_snocpymlt        = q_snocpymlt(:),     &
      & alb_vis_dir        = alb_vis_dir(:),     &
      & alb_nir_dir        = alb_nir_dir(:),     &
      & alb_vis_dif        = alb_vis_dif(:),     &
      & alb_nir_dif        = alb_nir_dif(:),     &
      & t_lwtr             = t_lwtr,             &
      & evapo_wtr          = evapo_wtr,          &
      & latent_hflx_wtr    = latent_hflx_wtr,    &
      & sensible_hflx_wtr  = sensible_hflx_wtr,  &
      & t_lice             = t_lice,             &
      & evapo_ice          = evapo_ice,          &
      & latent_hflx_ice    = latent_hflx_ice,    &
      & sensible_hflx_ice  = sensible_hflx_ice,  &
      & ice_fract_lake     = ice_fract_lake,     &
      & albedo_lwtr        = albedo_lwtr,        &
      & albedo_lice        = albedo_lice         &
      & )
    NULLIFY(tile)

!!$    IF (iblk == 1)                                                              &
!!$      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//          &
!!$      &                             ' on PE '//int2string(p_pe)//' completed.', &
!!$      &              all_print=.TRUE.)
    IF (debug_on('basic') .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

    IF (l_timer) CALL timer_stop(timer_jsbach)

  END SUBROUTINE interface_full

!   SUBROUTINE interface_standalone( &
!     & model_id,         &
!     & iblk, ics, ice, dtime, steplen, &
!     & t_air,            &
!     & q_air,            &
!     & rain,             &  
!     & snow,             &
!     & wind_air,         &
!     & wind_10m,         &
!     & lw_srf_down,      &
!     & swvis_srf_net,    &
!     & swnir_srf_net,    &
!     & swpar_srf_down,   &
!     & frac_par_diffuse, &
!     & press_srf,        &
!     & drag_srf,         &
!     & t_acoef,          &
!     & t_bcoef,          &
!     & q_acoef,          &
!     & q_bcoef,          &
!     & pch,              &
!     & cos_zenith_angle, &
!     & CO2_air,          &
!     & t_srf,            &
!     & t_eff_srf,        &
!     & qsat_srf,         &
!     & s_srf,            &
!     & fact_q_air,       &
!     & fact_qsat_srf,    &
!     & evapotrans,       &
!     & latent_hflx,      &
!     & sensible_hflx,    &
!     & grnd_hflx,        &
!     & grnd_hcap,        &
!     & rough_h_srf,      &
!     & rough_m_srf,      &
!     & q_snocpymlt,      &
!     & alb_vis_dir,      &
!     & alb_nir_dir,      &
!     & alb_vis_dif,      &
!     & alb_nir_dif,      &
!     & evapopot                                                & !< out; optional
!     & )
! 
!     USE mo_jsb_test,          ONLY: write_interface_vars, write_interface_variables
! 
!     INTEGER :: p_pe
! 
!     ! Arguments
!     INTEGER, INTENT(in) :: &
!       & model_id,          &
!       & iblk,              &
!       & ics,               &
!       & ice
! 
!     REAL(wp), INTENT(in) :: dtime, steplen
!     REAL(wp), INTENT(in) ::    &
!       & t_air             (:), &
!       & q_air             (:), &
!       & rain              (:), &
!       & snow              (:), &
!       & wind_air          (:), &
!       & wind_10m          (:), &
!       & lw_srf_down       (:), &
!       & swvis_srf_net     (:), &
!       & swnir_srf_net     (:), &
!       & swpar_srf_down    (:), &
!       & frac_par_diffuse  (:), &
!       & press_srf         (:), &
!       & drag_srf          (:), &
!       & t_acoef           (:), &
!       & t_bcoef           (:), &
!       & q_acoef           (:), &
!       & q_bcoef           (:), &
!       & pch               (:), &
!       & cos_zenith_angle  (:), &
!       & CO2_air           (:)
!     REAL(wp), INTENT(out) ::   &
!       & t_srf             (:), &
!       & t_eff_srf         (:), &
!       & qsat_srf          (:), &
!       & s_srf             (:), &
!       & fact_q_air        (:), &
!       & fact_qsat_srf     (:), &
!       & evapotrans        (:), &
!       & latent_hflx       (:), &
!       & sensible_hflx     (:), &
!       & grnd_hflx         (:), &
!       & grnd_hcap         (:), &
!       & rough_h_srf       (:), &
!       & rough_m_srf       (:), &
!       & q_snocpymlt       (:), &
!       & alb_vis_dir       (:), &
!       & alb_nir_dir       (:), &
!       & alb_vis_dif       (:), &
!       & alb_nir_dif       (:)
! 
!     REAL(wp), INTENT(out), OPTIONAL :: evapopot(:)
! 
!     ! Local variables
!     CLASS(t_jsb_model),          POINTER :: model
!     CLASS(t_jsb_tile_abstract),  POINTER :: tile
! 
!     INTEGER :: nc
!     INTEGER :: no_omp_thread
! 
!     CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_standalone'
! 
!     ! ---------------------------
!     ! Go
! 
!     no_omp_thread = Get_omp_thread()
! 
!     IF (write_interface_vars) THEN
!       CALL write_interface_variables(ics, ice, iblk, t_air, q_air, rain, snow, wind_air, wind_10m, &
!         & lw_srf_down, swvis_srf_net, swnir_srf_net, swpar_srf_down, press_srf, drag_srf,    &
!         & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle)
!     END IF
! 
!     CALL timer_start(timer_jsbach)
! 
!     nc = ice - ics + 1
! 
!     model => Get_model(model_id)
! 
! !!$    WRITE(message_text, '(A,I4,A,I6,A,I6)') 'Updating '//TRIM(model%name)// &
! !!$      &                 ' for: block=', iblk, ', ics=', ics, ', ice=', ice
! !!$    CALL message(TRIM(routine), message_text, all_print=.TRUE.)
! !!$    p_pe = get_my_global_mpi_id()
! !!$    CALL message( TRIM(routine), 'Updating '//TRIM(model%name)//' on PE '//int2string(p_pe)// &
! !!$      &                          ', nblk='//int2string(iblk), all_print=.TRUE.)
!     IF (debug_on('basic') .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))
! 
!     ! Update options
!     CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, dtime=dtime, steplen=steplen)
! 
!     ! Put atmospheric forcing into the a2l memory.
!     tile => model%Get_top_tile()
!     IF (tile%contains_lake) CALL finish(TRIM(routine), 'Standalone JSBACH interface doesn"t work with lakes yet!')
!     CALL update_atm2land(                           &
!       & tile,                                       &
!       & model%options(no_omp_thread),               &
!       & t_air              = t_air(:),              &
!       & q_air              = q_air(:),              &
!       & rain               = rain(:),               &
!       & snow               = snow(:),               &
!       & wind_air           = wind_air(:),           &
! !!$      & wind_10m          = wind_10m(:),          &
!       & wind_10m           = wind_air(:)*0.8_wp,    & ! Temporary
!       & lw_srf_down        = lw_srf_down(:),        &
!       & swvis_srf_net      = swvis_srf_net(:),      &
!       & swnir_srf_net      = swnir_srf_net(:),      &
!       & swpar_srf_down     = swpar_srf_down(:),     &
!       & frac_par_diffuse   = frac_par_diffuse(:),   &
!       & press_srf          = press_srf(:),          &
!       & drag_srf           = drag_srf(:),           &
!       & t_acoef            = t_acoef(:),            &
!       & t_bcoef            = t_bcoef(:),            &
!       & q_acoef            = q_acoef(:),            &
!       & q_bcoef            = q_bcoef(:),            &
!       & pch                = pch(:),                &
!       & cos_zenith_angle   = cos_zenith_angle(:),   &
!       & CO2_air            = CO2_air(:)             &
!       & )
!     NULLIFY(tile)
! 
!     ! Run tasks in queue
!     CALL model%Run_tasks(debug_on('hsm') .AND. iblk == 1)
! 
!     ! Finally, give the required variables back to the atmosphere
!     ! Note: evapopot is an optional parameter, so it has to be passed without the colon!
!     tile => model%Get_top_tile()
!     CALL update_land2atm(                     &
!       & tile,                                 &
!       & model%options(no_omp_thread),         &
!       & t_srf           = t_srf(:),           &
!       & t_eff_srf       = t_eff_srf(:),       &
!       & qsat_srf        = qsat_srf(:),        &
!       & s_srf           = s_srf(:),           &
!       & fact_q_air      = fact_q_air(:),      &
!       & fact_qsat_srf   = fact_qsat_srf(:),   &
!       & evapopot        = evapopot,           &
!       & evapotrans      = evapotrans(:),      &
!       & latent_hflx     = latent_hflx(:),     &
!       & sensible_hflx   = sensible_hflx(:),   &
!       & grnd_hflx       = grnd_hflx(:),       &
!       & grnd_hcap       = grnd_hcap,          &
!       & rough_h_srf     = rough_h_srf(:),     &
!       & rough_m_srf     = rough_m_srf(:),     &
!       & q_snocpymlt     = q_snocpymlt(:),     &
!       & alb_vis_dir     = alb_vis_dir(:),     &
!       & alb_nir_dir     = alb_nir_dir(:),     &
!       & alb_vis_dif     = alb_vis_dif(:),     &
!       & alb_nir_dif     = alb_nir_dif(:)      &
!       & )
!     NULLIFY(tile)
! 
! !!$    IF (iblk == 1)                                                              &
! !!$      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//          &
! !!$      &                             ' on PE '//int2string(p_pe)//' completed.', &
! !!$      &              all_print=.TRUE.)
!     IF (debug_on('basic') .AND. iblk == 1) &
!       & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')
! 
!     CALL timer_stop(timer_jsbach)
! 
!   END SUBROUTINE interface_standalone

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

    INTEGER :: no_omp_thread

    no_omp_thread = Get_omp_thread()

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    model => Get_model(model_id)
    tile => model%Get_top_tile()

    ! Update options
    CALL model%Set_options(steplen=steplen)

#ifndef __NO_JSBACH_HD__
      CALL hd_lateral_flow(tile, model%options(no_omp_thread))
      CALL hd_check_water_budget(tile, model%options(no_omp_thread))

      CALL interface_hd_ocean(tile, model%options(no_omp_thread))
#endif

  END SUBROUTINE jsbach_finish_timestep

#endif
END MODULE mo_jsb_interface
