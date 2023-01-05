!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish !,message_text
  ! USE mo_util_string,         ONLY: int2string

  USE mo_jsb_control,         ONLY: debug_on, timer_on, l_timer_host, timer_jsbach
  USE mo_timer,               ONLY: timer_start, timer_stop
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
!  USE mo_jsb_controller,      ONLY: Run_process_task
  USE mo_atmland_interface,   ONLY: update_atm2land, update_land2atm
  USE mo_jsb_subset,          ONLY: ON_CHUNK
  USE mo_jsb_parallel,        ONLY: Get_omp_thread, Is_omp_inside_serial !, get_my_global_mpi_id
  USE mo_jsb_time,            ONLY: t_datetime, is_newyear

  USE mo_anthro2l_interface,  ONLY: update_anthro2l  ! <Rackhun Son>

  dsl4jsb_Use_processes CARBON_, ALCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_interface, jsbach_start_timestep, jsbach_finish_timestep
  PUBLIC :: jsbach_get_var

  INTERFACE jsbach_interface
    MODULE PROCEDURE interface_full
    ! MODULE PROCEDURE interface_inquire
  END INTERFACE jsbach_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_interface'

CONTAINS

  SUBROUTINE interface_full( &
    & model_id,         &
    & iblk, ics, ice,   &
    & current_datetime, &
    & dtime, steplen,   &
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
    & fract_par_diffuse, &
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
    & CO2_flux,         &
    & DEBUG_VAR,        &
    ! For lakes (optional):
    & drag_wtr, drag_ice,                                               & !< in
    & t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,               & !< in
    & t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice,               & !< in
    & t_lwtr, qsat_lwtr, s_lwtr,                                        & !< out
    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                    & !< out
    & albedo_lwtr,                                                      & !< out
    & t_lice, qsat_lice, s_lice,                                        &
    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                    & !< out
    & albedo_lice,                                                      & !< out
    & ice_fract_lake,                                                   & !< out
    ! For standalone jsbach
    & evapopot,                                                         & !< out; optional
    ! For anthropogenic forcing     < Rackhun Son >
    & a_temperature,                                                    &
    & a_rel_humidity,                                                   &
    & population,                                                       &
    & gdp,                                                              &
    & hdi,                                                              &
    & land_use_primf,                                                   &
    & land_use_c4per,                                                   &
    & land_use_urban,                                                   &
    & land_use_secdf,                                                   &
    & land_use_secma,                                                   &
    & land_use_primn,                                                   &
    & land_use_secdn,                                                   &
    & land_use_c4ann,                                                   &
    & land_use_pastr,                                                   &
    & land_use_range,                                                   &
    & land_use_c3ann,                                                   &
    & land_use_c3nfx,                                                   &
    & land_use_c3per,                                                   &
    & land_use_secmb)

#ifndef __ICON__
    USE mo_jsb_test,          ONLY: write_interface_vars, write_interface_variables
#endif

    ! INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_full'

    ! Arguments
    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice
    TYPE(t_datetime), POINTER, INTENT(in) :: &
      & current_datetime
    REAL(wp), INTENT(in) :: &
      & dtime,              &
      & steplen
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
      & fract_par_diffuse  (:), &
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
      & DEBUG_VAR         (:), &
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
      & alb_nir_dif       (:), &
      & CO2_flux          (:)
    ! For lakes:
    REAL(wp), INTENT(out), OPTIONAL ::   &
      & t_lwtr            (:), &
      & qsat_lwtr         (:), &
      & s_lwtr            (:), &
      & evapo_wtr         (:), &
      & latent_hflx_wtr   (:), &
      & sensible_hflx_wtr (:), &
      & albedo_lwtr       (:), &
      & t_lice            (:), &
      & qsat_lice         (:), &
      & s_lice            (:), &
      & evapo_ice         (:), &
      & latent_hflx_ice   (:), &
      & sensible_hflx_ice (:), &
      & albedo_lice       (:), &
      & ice_fract_lake    (:)
    ! For standalone jsbach:
    REAL(wp), INTENT(out), OPTIONAL :: evapopot(:)

    ! For anthropogenic forcing:     < Rackhun Son >
    REAL(wp), INTENT(in), OPTIONAL ::    &
      & a_temperature(:),                &
      & a_rel_humidity(:),               &
      & population(:),                   &
      & gdp(:),                          &
      & hdi(:),                          &
      & land_use_primf(:),               &
      & land_use_c4per(:),               &
      & land_use_urban(:),               &
      & land_use_secdf(:),               &
      & land_use_secma(:),               &
      & land_use_primn(:),               &
      & land_use_secdn(:),               &
      & land_use_c4ann(:),               &
      & land_use_pastr(:),               &
      & land_use_range(:),               &
      & land_use_c3ann(:),               &
      & land_use_c3nfx(:),               &
      & land_use_c3per(:),               &
      & land_use_secmb(:)
    

    ! Local variables
    CLASS(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    INTEGER :: nc
    INTEGER :: no_omp_thread

    ! ---------------------------
    ! Go

    IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_jsbach(model_id))

    no_omp_thread = Get_omp_thread()

    ! Note: Since write_interface_variables is not yet implemented for ICON, the call to it is suppressed here for ICON.
    !       This is an exception since we don't want to use conditional compiling, but otherwise there are lots
    !       of compiler warnings about unused dummy arguments.
#ifndef __ICON__
    IF (write_interface_vars) THEN
      CALL write_interface_variables(ics, ice, iblk, current_datetime, t_air, q_air, rain, snow, wind_air, wind_10m, &
        & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf,    &
        & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle)
    END IF
#endif

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
    CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, current_datetime=current_datetime, dtime=dtime, steplen=steplen)
    CALL model%Set_subset(type=ON_CHUNK, iblk=iblk, ics=ics, ice=ice)

    CALL model%Associate_var_pointers(ics, ice, iblk, iblk)

    ! Put atmospheric forcing into the a2l memory.

    CALL model%Get_top_tile(tile)

    CALL update_atm2land(                           &
      & tile,                                       &
      & model%options(no_omp_thread),               &
      & t_air              = t_air(:),              &
      & q_air              = q_air(:),              &
      & rain               = rain(:),               &
      & snow               = snow(:),               &
      & wind_air           = wind_air(:),           &
      & wind_10m           = wind_10m(:),           &
      & lw_srf_down        = lw_srf_down(:),        &
      & swvis_srf_down     = swvis_srf_down(:),     &
      & swnir_srf_down     = swnir_srf_down(:),     &
      & swpar_srf_down     = swpar_srf_down(:),     &
      & fract_par_diffuse   = fract_par_diffuse(:), &
      & press_srf          = press_srf(:),          &
      & drag_srf           = drag_srf(:),           &
      & t_acoef            = t_acoef(:),            &
      & t_bcoef            = t_bcoef(:),            &
      & q_acoef            = q_acoef(:),            &
      & q_bcoef            = q_bcoef(:),            &
      & pch                = pch(:),                &
      & cos_zenith_angle   = cos_zenith_angle(:),   &
      & CO2_air            = CO2_air(:),            &
      & DEBUG_VAR          = DEBUG_VAR,             &
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

    ! Anthropogegnic forcing     < Rackhun Son >
    CALL model%Get_top_tile(tile)
    CALL update_anthro2l(                           &
      & tile,                                       &
      & model%options(no_omp_thread),               &
      & a_temperature      = a_temperature(:),      &
      & a_rel_humidity     = a_rel_humidity(:),     &
      & population         = population(:),         &
      & gdp                = gdp(:),                &
      & hdi                = hdi(:),                &
      & land_use_primf     = land_use_primf(:),     &
      & land_use_c4per     = land_use_c4per(:),     &
      & land_use_urban     = land_use_urban(:),     &
      & land_use_secdf     = land_use_secdf(:),     &
      & land_use_secma     = land_use_secma(:),     &
      & land_use_primn     = land_use_primn(:),     &
      & land_use_secdn     = land_use_secdn(:),     &
      & land_use_c4ann     = land_use_c4ann(:),     &
      & land_use_pastr     = land_use_pastr(:),     &
      & land_use_range     = land_use_range(:),     &
      & land_use_c3ann     = land_use_c3ann(:),     &
      & land_use_c3nfx     = land_use_c3nfx(:),     &
      & land_use_c3per     = land_use_c3per(:),     &
      & land_use_secmb     = land_use_secmb(:),     &
      & DEBUG_VAR          = DEBUG_VAR              &
      & )
    NULLIFY(tile)
    !------------------------------------------------

    ! Run tasks in queue
    CALL model%Run_tasks(debug_on('hsm') .AND. iblk == 1)


    ! give the required variables back to the atmosphere
    !   Note: evapopot is an optional parameter, so it has to be passed without the colon!
    CALL model%Get_top_tile(tile)
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
      & CO2_flux           = CO2_flux(:),        &
      & t_lwtr             = t_lwtr,             &
      & qsat_lwtr          = qsat_lwtr,          &
      & s_lwtr             = s_lwtr,             &
      & evapo_wtr          = evapo_wtr,          &
      & latent_hflx_wtr    = latent_hflx_wtr,    &
      & sensible_hflx_wtr  = sensible_hflx_wtr,  &
      & t_lice             = t_lice,             &
      & qsat_lice          = qsat_lice,          &
      & s_lice             = s_lice,             &
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

    IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_jsbach(model_id))

  END SUBROUTINE interface_full

  !
  ! Not used currently
  ! SUBROUTINE interface_inquire( model_id, iblk, ics, ice,         &
  !   & t_srf, t_eff_srf, qsat_srf, zh_srf, zm_srf, alb_vis, alb_nir)

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':interface_inquire'

  !   INTEGER, INTENT(in) :: &
  !     & model_id,          &
  !     & iblk,              &
  !     & ics,               &
  !     & ice

  !   REAL(wp), INTENT(inout) :: &
  !     & t_srf    (:), &
  !     & t_eff_srf(:), &
  !     & qsat_srf (:), &
  !     & zh_srf   (:), &
  !     & zm_srf   (:), &
  !     & alb_vis  (:), &
  !     & alb_nir  (:)

  ! END SUBROUTINE interface_inquire

  SUBROUTINE jsbach_start_timestep(model_id, current_datetime, dtime)

    USE mo_alcc_init,           ONLY : read_land_use_data

    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(wp),                  INTENT(in) :: dtime

    TYPE(t_jsb_model),           POINTER :: model

    INTEGER :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_start_timestep'

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    no_omp_thread = Get_omp_thread()   ! Should always be 1
    IF (no_omp_thread /= 1) CALL finish(routine, 'should not be called inside OMP parallel region!')

    model => Get_model(model_id)

    CALL model%Set_options(current_datetime=current_datetime, dtime=dtime)

    ! Reset all tiles at the beginning of the timestep
    CALL model%Reset_tiles()

    !> For the ALCC process land use data has to be read once a year
    IF (is_newyear(current_datetime,dtime) .AND. model%processes(ALCC_)%p%config%active) THEN
      CALL read_land_use_data(model_id, current_datetime)
    ENDIF

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_start_timestep

  SUBROUTINE jsbach_finish_timestep(model_id, current_datetime, steplen)

    USE mo_util,        ONLY: report_memory_usage
    USE mo_jsb_control, ONLY: get_debug_memory_level
    USE mo_jsb_time,    ONLY: is_newday, is_newmonth, is_time_experiment_start
    USE mo_carbon_interface,  ONLY: yday_carbon_conservation_test

    USE mo_seb_interface,      ONLY: global_seb_diagnostics
    USE mo_hydro_interface,    ONLY: global_hydrology_diagnostics
    USE mo_pheno_interface,    ONLY: global_phenology_diagnostics
    USE mo_carbon_interface,   ONLY: global_carbon_diagnostics
    USE mo_pplcc_interface,    ONLY: global_pplcc_diagnostics
#ifndef __NO_JSBACH_HD__
    USE mo_interface_hd_ocean, ONLY: interface_hd_ocean
    USE mo_hd_interface,       ONLY: hd_check_water_budget, hd_lateral_flow
#endif

    dsl4jsb_Use_processes SEB_
    dsl4jsb_Use_processes HYDRO_
    dsl4jsb_Use_processes PHENO_
    dsl4jsb_Use_processes PPLCC_

    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(wp),                  INTENT(in) :: steplen

    TYPE(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_finish_timestep'

    REAL(wp) :: dtime
    LOGICAL :: l_trace_memory
    INTEGER :: no_omp_thread

    no_omp_thread = Get_omp_thread()

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    model => Get_model(model_id)

    ! Update options
    CALL model%Set_options(current_datetime=current_datetime, steplen=steplen)

    dtime = model%options(no_omp_thread)%dtime

    ! Carbon conservation test
    CALL model%Get_top_tile(tile)
    IF (model%processes(CARBON_)%p%config%active) THEN
      IF (is_newday(current_datetime, dtime) .AND. .NOT. is_time_experiment_start(current_datetime)) THEN
        CALL yday_carbon_conservation_test(tile)
      END IF
    END IF

#ifndef __NO_JSBACH_HD__
    CALL model%Get_top_tile(tile)
    CALL hd_lateral_flow(tile, model%options(no_omp_thread))
    CALL hd_check_water_budget(tile, model%options(no_omp_thread))

    CALL interface_hd_ocean(tile, model%options(no_omp_thread))
#endif

! #ifdef _OPENACC
!    CALL model%gpu_to_cpu()
! #endif

    ! Global diagnostics
    !TODO: go through all processes.
    IF (model%processes(SEB_)%p%config%active) THEN
      CALL global_seb_diagnostics(tile)
    END IF
    IF (model%processes(HYDRO_)%p%config%active) THEN
      CALL global_hydrology_diagnostics(tile)
    END IF
    IF (model%processes(PHENO_)%p%config%active) THEN
      CALL global_phenology_diagnostics(tile)
    END IF
    IF (model%processes(CARBON_)%p%config%active) THEN
      CALL global_carbon_diagnostics(tile)
    END IF
    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL global_pplcc_diagnostics(tile)
    END IF

    l_trace_memory = .FALSE.
    SELECT CASE (get_debug_memory_level())
    CASE (1)
      l_trace_memory = is_newmonth(current_datetime, steplen)
    CASE (2)
      l_trace_memory = is_newday(current_datetime, steplen)
    CASE (3)
      l_trace_memory = .TRUE.
    END SELECT

    IF (l_trace_memory) CALL report_memory_usage()

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_finish_timestep

  SUBROUTINE jsbach_get_var(var_name, model_id, ics, ice, iblk, ptr1d, ptr2d, ptr3d)

    USE mo_jsb_var_class,     ONLY: t_jsb_var, REAL2D, REAL3D
    USE mo_jsb_process_class, ONLY: Get_process_id

    CHARACTER(len=*),  INTENT(in) :: var_name
    INTEGER,           INTENT(in) :: model_id
    INTEGER, OPTIONAL, INTENT(in) :: &
      & iblk,              &
      & ics,               &
      & ice
    REAL(wp), OPTIONAL, POINTER   :: ptr1d(:), ptr2d(:,:), ptr3d(:,:,:)

    TYPE(t_jsb_model),           POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile
    CLASS(t_jsb_var),            POINTER :: var_ptr
    CHARACTER(len=:), ALLOCATABLE :: process_name, variable_name
    INTEGER :: process_id, pos
    LOGICAL :: l_on_chunk
    INTEGER :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_get_var'

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)

    pos = INDEX(var_name, '_')
    IF (pos < 1) CALL finish(routine, 'var_name must have format <process name>_<variable name>')
    process_name = var_name(1:pos-1)
    variable_name = var_name(pos+1:LEN(TRIM(var_name)))
    process_id = Get_process_id(process_name)
    var_ptr => tile%mem(process_id)%p%Get_var(variable_name)

    IF (PRESENT(ics) .OR. PRESENT(ice) .OR. PRESENT(iblk)) THEN
      IF (.NOT. (PRESENT(ics) .AND. PRESENT(ice) .AND. PRESENT(iblk))) THEN
        CALL finish(routine, 'Must provide either none or all of ics,ice,iblk arguments!')
      ELSE
        l_on_chunk = .TRUE.
        IF (var_ptr%type == REAL2d) THEN
          IF (.NOT. PRESENT(ptr1d)) CALL finish(routine, 'Must provide argument ptr1d')
        ELSE
          IF (.NOT. PRESENT(ptr2d)) CALL finish(routine, 'Must provide argument ptr2d')
        END IF
      END IF
    ELSE
      IF (.NOT. Is_omp_inside_serial()) THEN
        CALL finish(TRIM(routine), 'Should not be called inside parallel OMP region')
      END IF
      l_on_chunk = .FALSE.
      IF (var_ptr%type == REAL2d) THEN
        IF (.NOT. PRESENT(ptr2d)) CALL finish(routine, 'Must provide argument ptr2d')
      ELSE
        IF (.NOT. PRESENT(ptr3d)) CALL finish(routine, 'Must provide argument ptr3d')
      END IF
    END IF

    IF (l_on_chunk) THEN
      IF (var_ptr%type == REAL2d) THEN
        ptr1d => var_ptr%ptr2d(ics:ice,iblk)
      ELSE
        ptr2d => var_ptr%ptr3d(ics:ice,:,iblk)
      END IF
    ELSE
      IF (var_ptr%type == REAL2d) THEN
        ptr2d => var_ptr%ptr2d(:,:)
      ELSE
        ptr3d => var_ptr%ptr3d(:,:,:)
      END IF
    END IF

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE jsbach_get_var

#endif
END MODULE mo_jsb_interface
