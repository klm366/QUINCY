PROGRAM jsb4_driver

#ifndef __NO_JSBACH__

  ! USE mo_nonhydro_state, ONLY: p_nh_state, p_nh_state_lists, construct_nh_state, destruct_nh_state

#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  USE, INTRINSIC :: ieee_arithmetic
#endif

  ! Used everywhere
  USE mo_exception,                 ONLY: message, message_text, finish
  USE mo_kind,                      ONLY: wp
  USE mo_impl_constants,            ONLY: SUCCESS, MAX_CHAR_LENGTH

  ! General model configurations
  USE mo_run_config,                ONLY: configure_run,                                          &
    &                                     dtime, nshift, num_lev, nsteps,                         &
    &                                     iforcing,                                               &
    &                                     output_mode
    ! &                                     iqv, iqc, iqi, iqs, iqr, iqt, ico2, iqm_max, nqtendphy, &
    ! &                                     ltestcase, ldynamics, ltransport                        &
    ! &                                     lvert_nest, ntracer,                                    &
    ! &                                     ldass_lhn, msg_level,                                   &
    ! &                                     iqc, iqt,                                               &
    ! &                                     ico2, io3,                                              &
    ! &                                     number_of_grid_used,                                    &
    ! &                                     grid_generatingCenter,                                  & ! grid generating center     R: nur zum Debuggen
    ! &                                     grid_generatingSubcenter                                  ! grid generating subcenter  R: nur zum Debuggen
  USE mo_nonhydrostatic_config,     ONLY: configure_nonhydrostatic ! R: sonst fehler in time loop wg substepping
  !USE mo_atm_phy_nwp_config,        ONLY: atm_phy_nwp_config ! R: nur zum debuggen
  USE mo_initicon_config,           ONLY: configure_initicon
  USE mo_dynamics_config,           ONLY: iequations, nnow_rcf
  USE mo_upatmo_config,             ONLY: upatmo_config
  USE mo_impl_constants,            ONLY: inh_atmosphere, & ! kind of atmosphere. R: sonst fehler in time loop wg substepping
    &                                     max_dom

  ! Parallelisation and decomposition
  USE mo_mpi,                       ONLY: start_mpi , stop_mpi, my_process_is_global_root, my_process_is_io, &
    &                                     set_mpi_work_communicators,   &
    &                                     process_mpi_io_size, my_process_is_pref, process_mpi_pref_size
  USE mo_icon_comm_interface,       ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_parallel_config,           ONLY: p_test_run, num_test_pe, l_test_openmp, num_io_procs, num_prefetch_proc, pio_type, &
    &                                     nproma
  USE mo_jsb_parallel,              ONLY: Get_omp_thread
  USE mo_build_decomposition,       ONLY: build_decomposition
  USE mo_async_latbc,               ONLY: prefetch_main_proc
  USE mo_complete_subdivision,      ONLY: setup_phys_patches

  ! Horizontal grids/patches and interpolation
  USE mo_grid_config,               ONLY: n_dom, n_dom_start !, &
    !&                                     end_time, dynamics_grid_filename
  USE mo_gridref_config,            ONLY: configure_gridref
  USE mo_interpol_config,           ONLY: configure_interpolation
  USE mo_nonhydrostatic_config,     ONLY: ivctype  ! R: nur damit parse_variable_groups keinen Abbruch erzeugt.
                                                  !    Wird normalerweise über die nonhydrostatic_nml eingelesen und ist 1 oder 2. Default ist 2.
  !USE mo_extpar_config,             ONLY: itopo ! event. necessary for correct building of vaxis
  USE mo_model_domain,              ONLY: p_patch, p_patch_local_parent
  USE mo_alloc_patches,             ONLY: destruct_patches, destruct_comm_patterns
  USE mo_intp_lonlat_types,         ONLY: lonlat_grids
  USE mo_intp_lonlat,               ONLY: compute_lonlat_intp_coeffs
  USE mo_grf_intp_state,            ONLY: construct_2d_gridref_state
  USE mo_intp_data_strc,            ONLY: p_int_state_local_parent, &
    &                                     p_int_state              ! R: wg. ivctype, ! p_int_state wird bei unserem call unten nur von der Subroutine verwendet wenn ivctype=2 (SLEVE). 
                                                                !    Was der Fall ist, ansonsten könnte ich das nur Platzhalter ungefüllt verwenden!!!
  USE mo_grf_intp_data_strc,        ONLY: p_grf_state, p_grf_state_local_parent
  USE mo_grf_intp_state,            ONLY: transfer_grf_state, create_grf_index_lists
  USE mo_intp_state,                ONLY: construct_2d_interpol_state, &
    &                                     transfer_interpol_state
  !USE mo_limarea_config,            ONLY: latbc_config
  !USE mo_opt_diagnostics,           ONLY: compute_lonlat_area_weights
  !USE mo_async_latbc_types,         ONLY: t_latbc_data

  ! Vertical grids
  USE mo_vertical_coord_table,      ONLY: vct_a, vct_b, vct, allocate_vct_atmo !, &
    !&                                     apzero ! R: wg. ivctype
  USE mo_util_vgrid,                ONLY: construct_vertical_grid                      ! R: wg. ivctype
  USE mo_init_vgrid,                ONLY: nflatlev                                     ! R: wg. ivctype
  USE mo_zaxis_type,                ONLY: zaxisTypeList, t_zaxisTypeList
  USE mo_nh_pzlev_config,           ONLY: configure_nh_pzlev
                            
  ! Reading the patches
  USE mo_loopindices,               ONLY: get_indices_c
  USE mo_impl_constants,            ONLY: min_rlcell_int
  USE mo_impl_constants_grf,        ONLY: grf_bdywidth_c

  ! Namelists
  USE mo_master_nml,                ONLY: read_master_namelist
  USE mo_master_config,             ONLY: isRestart !, setRestart
  USE mo_nml_crosscheck,            ONLY: atm_crosscheck
  USE mo_nml_crosscheck,            ONLY: land_crosscheck

  ! Restart
  USE mo_restart,                   ONLY: detachRestartProcs, t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor
  USE mo_load_restart,              ONLY: read_restart_files, read_restart_header
  USE mo_restart_attributes,        ONLY: t_RestartAttributeList, getAttributesForRestarting
  USE mo_output_event_handler,      ONLY: get_current_jfile

  ! I/O
  USE mo_io_config,                 ONLY: restartWritingParameters, configure_io, inextra_2d, inextra_3d
  USE mo_name_list_output_config,   ONLY: use_async_name_list_io ! asynchronous namelist I/O scheme, default is .FALSE. 
  USE mo_impl_constants,            ONLY: pio_type_async, pio_type_cdipio
  USE mo_name_list_output,          ONLY: close_name_list_output
  USE mo_cdi,                       ONLY: FILETYPE_GRB, FILETYPE_GRB2, FILETYPE_NC2, FILETYPE_NC4 
  USE mo_name_list_output,          ONLY: write_name_list_output, istime4name_list_output, name_list_io_main_proc

  USE mo_name_list_output_init,     ONLY:  init_name_list_output,        &
    &                                      parse_variable_groups ,       &
    &                                      collect_requested_ipz_levels, &
    &                                      output_file, create_vertical_axes
  USE mo_derived_variable_handling, ONLY: init_statistics_streams, finish_statistics_streams
  USE mo_level_selection,           ONLY: create_mipz_level_selections
  !USE mo_var_list,                  ONLY: nvar_lists, var_lists, print_var_list ! only for printing varlists
  USE mo_derived_variable_handling, ONLY: update_statistics, reset_statistics ! R: this is necessary and not only for debugging
  USE mo_nonhydro_state,            ONLY: p_nh_state, p_nh_state_lists ! R: wird beides nicht benötigt nur zum Debuggen von deallocate fehler
                                                                      !    Doch wird für zeitlichen Durchschnitt wohl benötigt
  USE mo_initicon_utils,            ONLY: reinit_average_first_guess

  ! Timer for model performance
  !USE mo_timer,                     ONLY: init_timer, timer_start, timer_stop,                &
  !  &                                     timers_level, timer_model_init,                     &
  !  &                                     timer_domain_decomp, timer_compute_coeffs,          &
  !  &                                     timer_ext_data, print_timer, ltimer, timer_total

  ! Time and event management
  USE mo_time_config,               ONLY: time_config, dt_restart
  USE mo_initicon_config,           ONLY: is_avgFG_time ! , timeshift
  USE mtime,                        ONLY: datetimeToString, datetime, newDatetime,                                           &
    &                                     event, eventGroup, newEvent, addEventToEventGroup, isCurrentEventActive, NO_ERROR, &
    &                                     OPERATOR(+), OPERATOR(==), OPERATOR(>), OPERATOR(>=),                              &
    &                                     MAX_MTIME_ERROR_STR_LEN, MAX_DATETIME_STR_LEN, max_timedelta_str_len,              &
    &                                     mtime_strerror, timedeltaToString, deallocateDatetime, timedelta,                  &
    &                                     getPTStringFromSeconds, newTimedelta, deallocateTimedelta
    !&                                     getTotalSecondsTimedelta                                                    &
    !&                                     OPERATOR(-), OPERATOR(*),                                                          &
    !&                                     ASSIGNMENT(=), OPERATOR(>=), OPERATOR(/=),                                         &
    !&                                     MAX_TIMEDELTA_STR_LEN,                                                             &
    !&                                     , getTimedeltaFromDatetime
  USE mo_time_management,           ONLY: compute_timestep_settings, compute_restart_settings, compute_date_settings
  USE mo_util_mtime,                ONLY: t_datetime_ptr, getElapsedSimTimeInSeconds
  USE mo_event_manager,             ONLY: initEventManager
  USE mo_event_manager,             ONLY: addEventGroup, getEventGroup, printEventGroup
  USE mo_action,                    ONLY: reset_act
  USE mo_jsb4_forcing,              ONLY: delta_time, time_step_len
                                        !next_date_echamfmt    => next_date,    &
                                        !current_date_echamfmt => current_date, & !Do not confuse with ICON-A: USE mtime, ONLY: datetime  ;   TYPE(datetime) :: current_date
  USE mo_jsb_control,               ONLY: jsbach_is_restarted

  !USE mo_datetime,                 ONLY: t_datetime, print_datetime

  ! Time loop management
  USE mo_nonhydrostatic_config,     ONLY: ndyn_substeps, ndyn_substeps_var
  USE mo_nh_stepping,               ONLY: deallocate_nh_stepping, diag_for_output_dyn, & !, reset_to_initial_state
    &                                     set_ndyn_substeps
  USE mo_output_event_types,        ONLY: t_sim_step_info
  USE mo_nh_prepadv_types,          ONLY: jstep_adv ! , prep_adv, t_prepare_adv,
  USE mo_util_timer,                ONLY: util_walltime

  ! Radiation requirements
  USE mo_psrad_orbit,               ONLY: orbit_kepler, orbit_vsop87, get_orbit_times ! for observational forcing
  USE mo_psrad_solar_parameters,    ONLY: solar_parameters
  USE mo_echam_rad_config,          ONLY: echam_rad_config
  USE mo_math_constants,            ONLY: pi

  ! External data, physics
  USE mo_ext_data_state,            ONLY: ext_data, destruct_ext_data
  USE mo_ext_data_init,             ONLY: init_ext_data
  USE mo_echam_phy_cleanup,         ONLY: cleanup_echam_phy
  USE mo_echam_phy_config,          ONLY: eval_echam_phy_config, eval_echam_phy_tc, print_echam_phy_config

  ! JSBACH4
  USE mo_lnd_nwp_config,            ONLY: nlev_soil, nlev_snow !, sstice_mode
  USE mo_jsb_version,               ONLY: jsbach_label_run !, jsbach_init_version
  USE mo_jsb_model_class,           ONLY: t_jsb_model
  USE mo_jsb_class,                 ONLY: get_model
  USE mo_jsb_base,                  ONLY: jsbach_setup_models, jsbach_setup_tiles
  USE mo_jsb_model_init,            ONLY: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart
  USE mo_jsb_model_final,           ONLY: jsbach_finalize
  !USE mo_jsb_process_class,         ONLY: SEB_, TURB_, HYDRO_
  dsl4jsb_Use_processes SEB_, TURB_, HYDRO_
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(HYDRO_)
  USE mo_jsb_tile_class,            ONLY: t_jsb_tile_abstract
  USE mo_jsb_interface,             ONLY: jsbach_interface, jsbach_start_timestep, jsbach_finish_timestep

  USE mo_jsb_namelist_iface,        ONLY: POSITIONED, open_nml, position_nml
  USE mo_jsb_driver_iface,          ONLY: read_atmo_namelists_for_jsb_standalone, &
    &                                     read_atm_forcing_for_jsb_stdalone_from_modelout, forcing_in_modeltyp
  USE mo_jsb_control,               ONLY: forcing_synchron_factor, forcing_steps_per_file, l_forcing_observations, &
    &                                     jsbach_runs_standalone !, jsbach_is_restarted

  USE mo_jsb_parallel,              ONLY: my_process_is_stdio, p_bcast, my_process_is_mpi_parallel !, mpi_comm , p_io
  !USE mo_jsb_time_iface,            ONLY: get_time_dt !, get_time_current
  USE mo_jsb_grid,                  ONLY: Get_grid
  USE mo_jsb_grid_iface,            ONLY: get_lat, get_lon
  USE mo_jsb_grid_class,            ONLY: t_jsb_grid
  USE mo_jsb_domain_iface,          ONLY: t_patch, get_nproma, get_nblks
  USE mo_jsb_time,                  ONLY: is_newday, jstep

  USE mo_jsb4_forcing,              ONLY: read_atm_forcing_for_jsb_stdalone_from_observ => get_interface_variables_from_external_forcing, init_forcing
  !USE mo_jsb_echam_iface,           ONLY: mtime_forcing  ! always one model timestep after mtime_current
                                                                          
  !USE mo_jsb4_forcing,              ONLY: mtime_current => next_date

  ! R: used only for debugging:
  !USE mo_rttov_interface,           ONLY: rttov_finalize, rttov_initialize
  !USE mo_synsat_config,             ONLY: lsynsat
  !USE mo_netcdf,                    ONLY: nf_check, nf_close
  ! time control and orbit calc
  !USE mo_nwp_ww,                    ONLY: configure_ww


  IMPLICIT NONE

  ! mtime_current is for the first forcing file one model_time_step ahead of mtime_forcing!
  ! In mo_jsb_forcing the naming is shifted as compared to the model/js4-driver.
  TYPE(datetime), POINTER       :: mtime_forcing
  TYPE(datetime), POINTER       :: mtime_current


  ! Technical stuff
#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  TYPE(ieee_status_type)   :: saved_fpscr
  LOGICAL                  :: halting_mode,  current_flags(size(ieee_all))
  REAL(wp)                 :: r
#endif
#endif
#if defined (__xlC__)
  INTEGER                     :: core_dump_flag
  INTEGER                     :: signals(1)
  INTEGER                     :: iret
#endif
  
  ! Used everywhere
  CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
  INTEGER                                :: error_status
  INTEGER                                :: i

  ! Time and event management
  REAL(wp)                               :: zwtime
  REAL(wp)                               :: sim_time          !< elapsed simulation time
  REAL(wp), ALLOCATABLE                  :: elapsedTime(:)    ! time elapsed since last call of NWP physics routines. For restart purposes.
  TYPE(datetime), POINTER                :: mtime_forcing_ptr ! 
  CHARACTER(len=32)                      :: mtime_forcing_char
  TYPE(datetime), POINTER                :: mtime_old         ! copy of current datetime (mtime)
  TYPE(datetime), POINTER                :: mtime_old_jsb     ! copy of current datetime (mtime)
  TYPE(datetime), POINTER                :: eventStartDate    => NULL(), &
    &                                       eventEndDate      => NULL()
  TYPE(datetime), POINTER                :: checkpointRefDate => NULL(), &
    &                                       restartRefDate    => NULL()
  TYPE(timedelta), POINTER               :: eventInterval     => NULL()
  TYPE(timedelta), POINTER               :: model_time_step   => NULL()
  TYPE(t_datetime_ptr)                   :: datetime_current(max_dom), datetime_local(max_dom)
  TYPE(event), POINTER                   :: checkpointEvent   => NULL()
  TYPE(event), POINTER                   :: restartEvent      => NULL()
  INTEGER                                :: checkpointEvents
  TYPE(eventGroup), POINTER              :: checkpointEventGroup => NULL()
  TYPE(t_sim_step_info)                  :: sim_step_info
  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)    :: modelTimeStep
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dt_string, dstring                                                         
  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)   :: td_string
  LOGICAL                                :: lprint_timestep, l_nml_output
    !&                                      lwrite_checkpoint, 
  CLASS(t_RestartDescriptor), POINTER    :: restartDescriptor
  TYPE(t_RestartAttributeList), POINTER  :: restartAttributes
  LOGICAL                                :: lret, lwrite_checkpoint
  LOGICAL                                :: l_isStartdate, l_isExpStopdate, l_isCheckpoint, l_doWriteRestart !, jsbach_is_restarted
  !INTEGER                               :: num_steps,    &    !< number of time steps to be executed if we would have a nested domain
    !&                                      jsubstep
  INTEGER                                :: jstep0, jstep_shift, jgp, step_in_forcing_file !, loop_step ! counters for time loop

  LOGICAL                                :: l_first_forcing_file = .TRUE. , l_shift_one_timestep = .TRUE.
  ! Parallelisation and decomposition
  INTEGER                                :: iblk, nblks, model_id, ics, ice, nc !, no_omp_thread
  INTEGER                                :: nblks_c !, nblks_c_max
  INTEGER                                :: dedicatedRestartProcs
  INTEGER, PARAMETER                     :: atmo_process = 1 ! to set_mpi_work_communicators: root_buffer(1)

  ! Vertical grids
  TYPE(t_jsb_grid),  POINTER             :: hgrid    ! Horizontal grid
  ! infos from the patch used for jsbach
  TYPE(t_patch)                          :: patch
  INTEGER                                :: jg          !< grid index / domain
  INTEGER                                :: ncd         !< number of child domains of grid jg
  INTEGER                                :: rls, rle    !< row limits, start and end
  INTEGER                                :: jb          !< index of block loop
  INTEGER                                :: jbs, jbe    !< start and end indices of block   loop
  INTEGER                                :: jcs, jce    !< start and end indices of columns loop
  !INTEGER                               :: nlev        !< number of levels

  ! Observation forcing
  LOGICAL                                :: lyr_perp = .FALSE. ! This switch has no effect on the results but we need it as a dummmy.
  INTEGER                                :: yr_perp  = -99999  ! Year if (lyr_perp == .TRUE.). We need it as a dummy.
  REAL(wp)                               :: time_of_day, time_of_day_rt, orbit_date, orbit_date_rt, &
    &                                       rasc_sun, decl_sun, dist_sun, flx_ratio
  REAL(wp), ALLOCATABLE                  :: lon(:,:), lat(:,:), coslon(:,:), coslat(:,:), sinlon(:,:), sinlat(:,:)
  REAL(wp), ALLOCATABLE                  :: daylight_frc(:,:) ! daylight fraction (0 or 1) with diurnal cycle
  REAL(wp), POINTER                      :: cecc, cobld, clonp
  INTEGER , POINTER                      :: icosmu0
  LOGICAL , POINTER                      :: l_orbvsop87
  REAL(wp)                               :: dt_ext
  LOGICAL, POINTER                       :: ldiur, l_sph_symm_irr

  ! Structural rest
  INTEGER                               :: master_namelist_status
  !CHARACTER(len=MAX_CHAR_LENGTH)       :: restart_filename
  CHARACTER(len=MAX_CHAR_LENGTH)        :: attname   ! restart file attribute name
  CHARACTER(len=*), PARAMETER           :: routine = 'PROGRAM jsb4_icon-driver'
  TYPE(t_jsb_model), POINTER            :: model
  CLASS(t_jsb_tile_abstract),  POINTER  :: tile
  dsl4jsb_Def_memory(SEB_)
  dsl4jsb_Def_memory(TURB_)
  dsl4jsb_Def_memory(HYDRO_)
  INTEGER, ALLOCATABLE                  :: output_jfile(:) ! number of output file
  INTEGER                               :: remain, step_for_synchron
  REAL(wp)                              :: synchr_ratio 
  !TYPE(t_latbc_data)                   :: latbc !< data structure for async latbc prefetching

  ! JSB4 interface variables
  ! in vars
  REAL(wp), ALLOCATABLE :: DEBUG_VAR        (:,:), & ! only for debugging purpose
    &                      t_air            (:,:), &
    &                      q_air            (:,:), &
    &                      rain             (:,:), &
    &                      snow             (:,:), &
    &                      wind_air         (:,:), &
    &                      wind_10m         (:,:), &
    &                      lw_srf_down      (:,:), &
    &                      swvis_srf_down   (:,:), &
    &                      swnir_srf_down   (:,:), &
    &                      swpar_srf_down   (:,:), &
    &                      fract_par_diffuse(:,:), &
    &                      press_srf        (:,:), &
    &                      drag_srf         (:,:), &
    &                      t_acoef          (:,:), &
    &                      t_bcoef          (:,:), &
    &                      q_acoef          (:,:), &
    &                      q_bcoef          (:,:), &
    &                      pch              (:,:), &
    &                      cos_zenith_angle (:,:), & 
    &                      CO2_air          (:,:), & 
    &                      trig_rad_m1      (:,:), &
    &                      newday_Merker    (:,:)

  REAL(WP), ALLOCATABLE ::  l_day_reminder   (:,:)

  ! out vars, for standalone only dummy use but necessary!
  REAL(wp), ALLOCATABLE :: t_srf_fromJS4        (:,:), &
    &                      t_eff_srf         (:,:),   &
    &                      qsat_srf          (:,:),   &
    &                      s_srf             (:,:),   &
    &                      fact_q_air_fromJS4   (:,:), &
    &                      fact_qsat_srf_fromJS4(:,:), &
    &                      evapotrans        (:,:),   &
    &                      evapopot          (:,:),   & 
    &                      latent_hflx       (:,:),   &
    &                      sensible_hflx     (:,:),   &
    &                      grnd_hflx         (:,:),   &
    &                      grnd_hcap         (:,:),   &
    &                      rough_h_srf_fromJS4  (:,:), &
    &                      rough_m_srf_fromJS4  (:,:), &
    &                      q_snocpymlt       (:,:), &
    &                      alb_vis_dir       (:,:), &
    &                      alb_nir_dir       (:,:), &
    &                      alb_vis_dif       (:,:), &
    &                      alb_nir_dif       (:,:), &
    &                      CO2_flux          (:,:)

  ! in vars for lakes, optional
  REAL(wp), ALLOCATABLE :: drag_wtr         (:,:), &
    &                      drag_ice         (:,:), &
    &                      t_acoef_wtr      (:,:), &
    &                      t_bcoef_wtr      (:,:), &
    &                      q_acoef_wtr      (:,:), &
    &                      q_bcoef_wtr      (:,:), &
    &                      t_acoef_ice      (:,:), &
    &                      t_bcoef_ice      (:,:), &
    &                      q_acoef_ice      (:,:), &
    &                      q_bcoef_ice      (:,:)

  ! out vars for lakes, optional. Only dummy use for standalone but necessary!
  REAL(wp), ALLOCATABLE :: t_lwtr           (:,:), &
    &                      qsat_lwtr        (:,:), &
    &                      s_lwtr           (:,:), &
    &                      evapo_wtr        (:,:), &
    &                      latent_hflx_wtr  (:,:), &
    &                      sensible_hflx_wtr(:,:), &
    &                      albedo_lwtr      (:,:), &
    &                      t_lice           (:,:), &
    &                      qsat_lice        (:,:), &
    &                      s_lice           (:,:), &
    &                      evapo_ice        (:,:), &
    &                      latent_hflx_ice  (:,:), &
    &                      sensible_hflx_ice(:,:), &
    &                      albedo_lice      (:,:), &
    &                      ice_fract_lake   (:,:)

  ! perhaps we want to use later
  ! LOGICAL, ALLOCATABLE :: linit_dyn (:) ! Determines whether dynamics must be initialized on given patch.
                                          ! Only used when jg ==1 and l_limited_area is .true.


  !     Nicht in JULIAS VERSION:
  !         ! For lakes (optional):
  !         & drag_wtr, drag_ice,                                               & !< in
  !         & t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,               & !< in
  !         & t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice,               & !< in
  !         & t_lwtr, qsat_lwtr, s_lwtr,                                        &
  !         & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr, & !< out
  !         & albedo_lwtr,                         INTEGER                      :: counter                               & !< out
  !         & t_lice, qsat_lice, s_lice,                                        &
  !         & evapo_ice, latent_hflx_ice, sensible_hflx_ice, & !< out
  !         & albedo_lice,                                                      & !< out
  !         & ice_fract_lake,                                                   & !< out
  !         
  ! Reste:
  REAL(wp), ALLOCATABLE :: evapo_act2pot_fromJS4generated(:,:) !Observations: evapo_act2pot.
  !  REAL(wp), ALLOCATABLE :: CO2_concentration(:,:)  ! Observations: CO2_air in jsb_interface. Nicht mehr benötigt ist jetzt auch in mo_forcing CO2_air!
  !  REAL(wp), ALLOCATABLE :: t_air(:,:)            ! Observations: model-interface: t_air
  !  REAL(wp), ALLOCATABLE :: evapo_act_sum_fromJS4(:,:)
  !  REAL(wp), ALLOCATABLE :: evapo_pot_sum_fromJS4(:,:)


  ! R: used only for debugging:
  INTEGER                      :: counter

  ! -------------------------------------------------------------------
  ! Go

  ! -------------------------------------------------------------------------------------------------------------------------------
  ! 1. Technical stuff
  ! -------------------------------------------------------------------------------------------------------------------------------

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  CALL ieee_get_status(saved_fpscr)
  CALL ieee_set_halting_mode(ieee_all, .TRUE.)
#endif
#endif

#if defined (__INTEL_COMPILER) && ! defined (VARLIST_INITIZIALIZE_WITH_NAN)
  ! Important on Intel: disable underflow exceptions:
  CALL ieee_set_halting_mode(ieee_underflow, .FALSE.)
#endif
 

  ! -------------------------------------------------------------------------------------------------------------------------------
  ! 2. Setup the basic ICON-A framework:
  ! (Start and setup of MPI, read basic model infos from namelists, read restart header, setup time controls)
  ! -------------------------------------------------------------------------------------------------------------------------------

  CALL message(TRIM(routine), 'Starting the ICON JSB4 driver ...')

    ! Start MPI
  CALL start_mpi('JS4_driver')

  ! Set up signal trapping on IBM: export USE_SIGNAL_HANDLING=yes
#if defined (__xlC__)
  core_dump_flag = 0
  signals(1)     = 0
  iret = signal_trap(core_dump_flag, signals)
  IF (iret == -2) THEN
    CALL message('', 'Signal trapping disabled by environment')
  ELSE IF (iret == -1) THEN
    WRITE(message_text,'(a,i0)') 'Error: ', iret
    CALL message('', message_text)
  ELSE IF (iret == 0) THEN
    CALL message('', 'FPE trapping is not set')
  ELSE
    WRITE(message_text,'(a,i0)') 'FPE trapping mode =', iret
    CALL message('', message_text)
  END IF
#endif
#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  WRITE(message_text,'(a,l1)') ' IEEE standard supported: ', ieee_support_standard(r)
  CALL message('', message_text)
#endif
#endif
 
  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime()
!$OMP END MASTER
!$OMP END PARALLEL
 
  ! Read basic namelists for initial setup of the model
  ! Reads master namelists of ICON (for JS4 standalone case): 
  ! master_nml, master_time_control_nml, master_model_nml
  master_namelist_status = read_master_namelist('NAMELIST_JS4-standalone_master')

  ! Initialize global registry of lon-lat grids
  CALL lonlat_grids%init()

  ! If this is a restart read the header of the restart file
  restartAttributes => NULL()
  IF (isRestart()) THEN
    CALL message('','Read restart file meta data ...')
    ! Reads attributes and namelists for all available domains from restart file.
    ! These namelist values will overwrite the model default, and will later be overwritten if the user has
    ! specified something different for this integration.
    CALL read_restart_header("atm") ! atm is the modelType
    ENDIF

  ! Reads master namelists of ICON (for JS4 standalone case): time_nml
  ! Reads ICON-A atmo namelists (for JS4 standalone case): 
  ! parallel_nml, run_nml, io_nml, output_nml, grid_nml, gridref_nml and calls init_grid_configuration
  CALL read_atmo_namelists_for_jsb_standalone('NAMELIST_JS4-standalone_master')

  ! We want to avoid reading the whole dynamics_nml, but we need iequations in compute_restart_settings
  iequations = 3 ! 0: shallow water model
                ! 1: hydrostatic atmosphere, T
                ! 2: hydrostatic atm., θ·dp
                ! 3: non-hydrostatic atmosphere (this is default)
                !-1: hydrostatic ocean

  ! Set times and events
  CALL compute_timestep_settings()
  CALL compute_restart_settings()
  CALL compute_date_settings("atm", dt_restart, nsteps)
  CALL initEventManager(time_config%tc_exp_refdate)

  ! Check if compiled with JSBACH and if in echam_phy_config llake=true but JSBACH is false
  CALL land_crosscheck() ! wieder rein nehmen nur zum debuggen raus genommen

  ! Finish filling the run_config state.
  CALL configure_run( )

  ! This has to be done before configure_initicon
  ! We want to avoid configure_upatmo in the driver. As icon forces us we do it like this.
  ALLOCATE(upatmo_config(n_dom_start:n_dom), STAT=error_status)
  IF (error_status /= SUCCESS) CALL finish(TRIM(routine), "Allocation of upatmo_config failed.")


  ! Complete init icon config-state e.g. mtime
  CALL configure_initicon()

  ! Parallelization check
  ! CALL check_parallel_configuration()
  ! Reset num_prefetch_proc to zero if the model does not run in limited-area mode
  ! or if there are no lateral boundary data to be read
  ! For the JS4-driver this check is not necessary we just set:
  num_prefetch_proc = 0


  ! Initialize the MPI work groups
  CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs) ! out, get dedicatedRestartProcs.  
                                                                                ! opt_dedicatedProcCount is the number 
                                                                                ! of processes that are split off the
                                                                                ! MPI_Communicator to serve as 
                                                                                ! dedicated restart processes.
  CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, dedicatedRestartProcs, atmo_process, & ! inout
    &                             num_prefetch_proc, num_test_pe, pio_type)                                       ! in, optional

  ! Initialize various timers
  ! ltimer switches  on/off  the  output  of  some performance  related  information  (run time)
  ! IF (ltimer) CALL init_timer
  ! IF (timers_level > 3) CALL timer_start(timer_model_init)

  ! Initialize dynamic list of vertical axes
  zaxisTypeList = t_zaxisTypeList()

  ! -------------------------------------------------------------------------------------------------------------------------------
  ! 3. Setup basic JSBACH4:
  ! (Read namelists, configure JS4-models for each domain, add zaxes to zaxisTypeList)
  ! -------------------------------------------------------------------------------------------------------------------------------
  ! This has to be done after (!) the ICON zaxes have been created in the above line but before (!) the 
  ! restart PEs are detached a line below! The reason for that is that JSBACH4 adds its zaxes to zaxisTypeList.

  ! Reads master   namelists of JSBACH4: jsb_control_nml, jsb_model_nml.
  ! Reads inferior namelists of JSBACH4: NAMELIST_JS4-standalone_jsbach: jsb_model_nml and the process_nmls.
  CALL jsbach_setup_models('NAMELIST_JS4-standalone_master')


  ! -------------------------------------------------------------------------------------------------------------------------------
  ! 4. MPI: detach restart, prefetch and I/O procs, setup time controls for this simulation step, build_decomposition)
  ! -------------------------------------------------------------------------------------------------------------------------------

  ! MPI: detach the restart processes 
  CALL detachRestartProcs() ! Note, without MPI there are no dedicated restart processes and the model exits.
                            ! The dedicated restart processes separate here and start according to the nml a 
                            ! (1) sychron, (2) asynchron or (3) multifile restart-loop.
                            ! From here on they are waiting for a signal for writing the restart file or for 
                            ! the end signal of the model run.
    
  ! MPI: detatch the prefetching PEs if asynchronous input prefetching is enabled
  ! These prefetching PEs are used to load boundary conditions
  IF (process_mpi_pref_size > 0) THEN ! if we belong to the prefetching of data PEs 
    num_prefetch_proc = 1             ! in the moment not more than 1 possible!
    CALL message(routine, 'asynchronous input prefetching is enabled.')
    IF (my_process_is_pref()) CALL prefetch_main_proc ! Start prefetching loop (similar to the restart loop).
  ENDIF

  ! MPI: detach the I/O processes
  ! Note, process_mpi_io_size = 0 means this is no I/O PE
  IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_async) THEN ! If we want I/O PEs  AND  pio_type is asynchron
                                                                    ! (pio_type_async=1). This is the default case.
    ! Decide whether async vlist or name_list IO is to be used, only one of both may be enabled!
    IF (output_mode%l_nml) THEN ! if we want namelist output
      ! asynchronous I/O
      use_async_name_list_io = .TRUE.
      CALL message(routine, 'asynchronous namelist I/O scheme is enabled.')
      ! consistency check
      IF (my_process_is_io()) THEN  ! if we are on an I/O PE
        ! compute sim_start, sim_end
        CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start) ! e.g. "1980 1  1  0  0  0  0"   to 
                                                                                    !      "1980-01-01T00:00:00.000" (ISO_8601)
        CALL datetimeToString(time_config%tc_exp_stopdate,  sim_step_info%sim_end)
        CALL datetimeToString(time_config%tc_startdate,     sim_step_info%run_start)
        CALL datetimeToString(time_config%tc_stopdate,      sim_step_info%restart_time)
        sim_step_info%dtime = dtime
        jstep0 = 0

        restartAttributes => getAttributesForRestarting()
        IF (ASSOCIATED(restartAttributes)) THEN
          ! get start counter for time loop from restart file:
          jstep0 = restartAttributes%getInteger("jstep")
        END IF
        sim_step_info%jstep0    = jstep0
        CALL name_list_io_main_proc(sim_step_info)  ! here we detach the I/O processes!
      END IF
    ELSE IF (my_process_is_io()) THEN  ! if we are on an I/O PE but output_mode%l_nml = FALSE
      ! Shut down MPI
      CALL stop_mpi
      STOP
    ENDIF ! output_mode%l_nml
  ELSE IF (process_mpi_io_size > 0 .AND. pio_type == pio_type_cdipio) THEN  ! If we want I/O PEs  AND  pio_type 
                                                                            ! is cdi (pio_type_cdipio = 2).
    CALL finish(routine, 'CDI-PIO requested but not implemented in JS4 driver yet.')
  ELSE  ! if we are not on an I/O PE
    ! non-asynchronous I/O (performed by PE #0)
    IF (output_mode%l_nml) THEN
      CALL message(routine, 'synchronous namelist I/O scheme is enabled.')
      CALL message(routine, 'in combination with asynchronous input prefetching this will &
                          &  lead to a missing latbc_config%latbc_varnames_map_file....')
    ENDIF
  ENDIF

  ! Import patches, perform domain decomposition and define vertical coordinates
  CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .false.) ! all in


  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 5. Construct horizontal and vertical grid interpolation state, compute interpolation coefficients.
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! Prepare grid information
  CALL configure_interpolation( n_dom, p_patch(1:)%level, p_patch(1:)%geometry_info ) ! all in

  ! Necessary as soon as we will use more than one grid domain (n_dom>1)
  ! CALL configure_gridref(n_dom, p_patch(1:)%geometry_info%mean_characteristic_length) ! all in

  ! Allocate array for interpolation state
  ALLOCATE( p_int_state(n_dom_start:n_dom), &
    &       p_grf_state(n_dom_start:n_dom), &
    &       STAT=error_status )
  IF (error_status /= SUCCESS) THEN
    CALL finish(routine, 'allocation for ptr_int_state or p_grf_state failed')
  ENDIF

  ALLOCATE( p_int_state_local_parent(n_dom_start+1:n_dom), &
    &       p_grf_state_local_parent(n_dom_start+1:n_dom), &
    &       STAT=error_status )
  IF (error_status /= SUCCESS) THEN
    CALL finish(routine, 'allocation for p_int_state_local_parent or p_grf_state_local_parent failed')
  ENDIF

  ! Construct interpolation state
  ! Please note that for parallel runs the divided state is constructed here
  CALL construct_2d_interpol_state(p_patch, p_int_state)

  ! Transfer interpolation state to local parent
  DO jg = n_dom_start+1, n_dom
    jgp = p_patch(jg)%parent_id
    CALL transfer_interpol_state( p_patch(jgp), p_patch_local_parent(jg), &
      &  p_int_state(jgp), p_int_state_local_parent(jg) )
  ENDDO


  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 6. Construct grid refinment state, compute coefficients
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! For the NH model, the initialization routines called from
  ! construct_2d_gridref_state require the metric terms to be present

  ! Necessary as soon as we will use more than one grid domain (n_dom>1)
  ! IF (n_dom_start==0 .OR. n_dom > 1) THEN ! n_dom_start is 1 except no_of_radiation_grids > 0
  !   ! Construct gridref state
  !   ! For non-parallel runs (either 1 PE only or on the test PE) this is done as usual,
  !   ! for parallel runs, the main part of the gridref state is constructed on the
  !   ! local parent with the following call
  !   CALL construct_2d_gridref_state (p_patch, p_grf_state)
  !
  !   ! Transfer gridref state from local parent to p_grf_state
  !   DO jg = n_dom_start+1, n_dom
  !     jgp = p_patch(jg)%parent_id
  !     CALL transfer_grf_state(p_patch(jgp), p_patch_local_parent(jg),         &
  !       &                     p_grf_state(jgp), p_grf_state_local_parent(jg), &
  !       &                     p_patch(jg)%parent_child_index)
  !   ENDDO
  ! ENDIF


  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 7. MPI: Establish the communication
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! Initialize icon_comm_lib
  CALL construct_icon_communication(p_patch, n_dom) ! all inout (INTENT not defined)

  ! Setup the information for the physical patches
  CALL setup_phys_patches

  ! Constructing data for lon-lat interpolation
  CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:)) ! inout, inout

  ! Necessary as soon as we will use more than one grid domain (n_dom>1)
  ! IF (n_dom_start==0 .OR. n_dom > 1) THEN ! n_dom_start is 1 except no_of_radiation_grids > 0
  !   CALL create_grf_index_lists(p_patch, p_grf_state, p_int_state)
  ! ENDIF


  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 8. Construct vertical atmospheric grid
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! The following is needed for this driver only because SUBROUTINE construct_vertical_grid needs a topography
  ! Create and optionally read external data fields
  ALLOCATE (ext_data(n_dom), STAT=error_status)
  IF (error_status /= SUCCESS) THEN
    CALL finish(routine, 'allocation for ext_data failed')
  ENDIF
  ! allocate memory for atmospheric/oceanic external data and optionally read those data from netCDF file.
  CALL init_ext_data (p_patch(1:), p_int_state(1:), ext_data) ! inout, in, inout. Note, p_int_state is used only for iforcing == inwp. 
                                                              ! ext_data
  ! Import vertical grid/ define vertical coordinate
  CALL allocate_vct_atmo(p_patch(1)%nlevp1)
  DO jg = 1, n_dom
    ext_data(jg)%atm%topography_c(1:nproma,1:p_patch(jg)%nblks_c) = 0.0_wp
  ENDDO

  ! R: wäre eigentlich effizienter so: 
  ! initicon%const%topography_c => topography_c !< topographic height at cell centers [m], index1=1,nproma, index2=1,nblks_c
  ! ALLOCATE(dummy_target_topo(n_dom, nproma, nblks_c_max))
  ! dummy_target_topo(:, :,:) = 1
  ! DO jg = 1,n_dom
  !   nblks_c   = p_patch(jg)%nblks_c
  !   !ext_data(jg)%atm%topography_c => dummy_target_topo(jg, :,1:nblks_c) ! R: Dummy Wert setzten
  !   ext_data(jg)%atm%topography_c => dummy_target_topo(:,1:nblks_c) ! R: Dummy Wert setzten
        ! R: Problem: Pointer zeigt ins leere weil SUBROUTINE construct_atmo_nonhydrostatic nicht 
        !    aufgerufen wurde und dadurch SUBROUTINE init_icon nicht aufgerufen wurde.
        !    und dadurch wurde dies nicht ausgeführt: 
        !    DO jg = 1, n_dom
        !         initicon(jg)%const => initicon_const(jg)
        !         CALL construct_initicon(initicon(jg), p_patch(jg), ext_data(jg)%atm%topography_c, p_nh_state(jg)%metrics)
        !    END DO
        !    Die macht wiederum nur initicon%const%topography_c => topography_c
  ! END DO

  ! This is only necessary because the namelist output needs a atmospheric zaxis
  ! (would not be necessary for this driver and could only be changed in ICON)
  CALL construct_vertical_grid(p_patch(1:),       & ! inout
                              & p_int_state(1:),  & ! in
                              & ext_data,         & ! inout, is needed for CALL init_vert_coord to provide ext_data(jg)%atm%topography_c
                              & vct_a,            & ! inout, hight of atmo level in m.
                              & vct_b,            & ! inout
                              & vct,              & ! inout
                              & nflatlev)           ! inout, nflatlev = 0 even when num_lev=1 or =31


  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 9. Horizontal and vertical grid(s) are now defined. Assign values to derived variables in the configuration states
  !-------------------------------------------------------------------------------------------------------------------------------------

  ! CALL configure_diffusion( n_dom, dynamics_parent_grid_id,       &
  !   &                       p_patch(1)%nlev, vct_a, vct_b, apzero )
  ! CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_phys_dom)

  ! IF (iequations == inh_atmosphere) THEN
  DO jg =1,n_dom
    CALL configure_nonhydrostatic( jg, p_patch(jg)%nlev,     &
      &                            p_patch(jg)%nshift_total  )
    !    IF ( iforcing == inwp) THEN
    !      CALL configure_ww( time_config%tc_startdate, jg, p_patch(jg)%nlev, p_patch(jg)%nshift_total, 'ICON')
    !    END IF
  ENDDO
  ! ENDIF

  ! R: hier die ALLOCATES raus und nach unten geschoben
  ! Außerdem wurden die JSBACH4 Teile nach unten verschoben

    !   IF (.NOT. read_interface_vars) THEN ! R: forcing data not from ICON => das wird bereits vom JS4-driver für ECHAM6 erledigt!
    !    ...
    !   ENDIF

  !!! SUBROUTINE construct_atmo_nonhydrostatic:

  ! R: Anfang: wird beides nicht benötigt nur zum Debuggen von deallocate fehler:
    ALLOCATE (p_nh_state(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state failed')
    ENDIF

    ALLOCATE (p_nh_state_lists(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state_lists failed')
    ENDIF
  ! R: Ende: wird beides nicht benötigt nur zum Debuggen von deallocate fehler

  !R: Das hier ist wohl nötig um echam_phy_nml nicht nur zu lesen so wie dies in mo_jsb_driver_interface.f90: CALL process_echam_phy_nml  getan wird sondern auch bereit zu stellen:
  ! drivers/mo_atmo_nonhydrostatic.f90: SUBROUTINE construct_atmo_nonhydrostatic(latbc): Darin: CALL init_echam_phy_params       ( p_patch(1:) ) und darin:

    ! ECHAM physics time control
    ! --------------------------
    ! Evaluate the ECHAM physics configuration variables echam_phy_config(:)
    ! and the derived time control variables echam_phy_tc(:) on all grids
    ! and for all controled processes.
    !
    CALL eval_echam_phy_config
    CALL eval_echam_phy_tc
    CALL print_echam_phy_config



  !-------------------------------------------------------------------------------------------------------------------------------------
  ! 10. Setup further JSBACH4
  !-------------------------------------------------------------------------------------------------------------------------------------
  ! Note: we had to move this after eval_echam_phy_tc because we need echam_phy_tc

  ! Setup horizontal grids and HSM (tiles) for JSBACH
  DO jg=1,n_dom
      CALL jsbach_setup_grid( jg, p_patch(jg), type='icon')
      CALL jsbach_setup_tiles(jg)
  END DO

  ! Create queue tasks, create process memories and initialize them for each tile (read bc and ic files)
  DO jg=1,n_dom
    CALL jsbach_init(jg)
  END DO

  ! Print out informations
  CALL jsbach_label_run(jsbach_runs_standalone())

    ! Now we know the decomposition and we can allocate the interface variables
    ! R: diesen Block muß ich sobald ich die folgenden ALLOCATE los bin nach X1!
    jg = 1 ! R: wir haben nur ein Modell aber im weiteren Verlauf wird dennoch immer ein looü über jg gemacht, wir können die zeile hier also wieder löschen
    model_id = jg ! R: das mußte ich hinzufügen da dies normalerweise über CALL jsbach_interface ( jg, nblock, jcs, kproma, pdtime, pdtime,...) durchgeführt würde was hier im driver weg fällt



  ! Read restart file

  !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   IF Restart'

  IF (jsbach_is_restarted()) THEN ! isRestart() verwenden?
    l_shift_one_timestep = .FALSE.
    ! Brauch ich wohl nicht:
    ! CALL get_restart_attribute("sim_time_DOM01", sim_time)
    !IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
    !    ! get start counter for time loop from restart file:
    !    CALL get_restart_attribute("jstep", jstep0)
    ! END IF


    !     ! init LES
    !     DO jg = 1 , n_dom
    !       IF(atm_phy_nwp_config(jg)%is_les_phy) THEN
    !         CALL init_les_phy_interface(jg, p_patch(jg)       ,&  ! R: das macht nur: CALL init_vertical_grid_for_les(jg, p_patch, p_int_state, p_metrics)   weiß nicht ob ich das wirklich brauche.
    !            &                        p_int_state(jg)       ,&
    !            &                        p_nh_state(jg)%metrics)
    !       END IF
    !     END DO


    !DO jg = 1,n_dom
    !  IF (p_patch(jg)%ldom_active) THEN
    !    CALL read_restart_files( p_patch(jg), n_dom)
    !  END IF
    !END DO
    CALL read_restart_files( p_patch(jg), n_dom) ! n_dom = 1, from CALL init_grid_configuration() in read_atmo_namelists_for_jsb_standalone

    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   Restart'
    !DO jg = 1,n_dom
    !  IF (.NOT. p_patch(jg)%ldom_active) CYCLE
    !    IF (echam_phy_config(jg)%ljsb) THEN
    !      CALL jsbach_init_after_restart(jg)
    !  END IF
    !END DO
    CALL jsbach_init_after_restart(model_id)

    ! R: gibts nicht in ICON  CALL init_output
    ! R: gibts nicht in ICON  CALL open_output_streams

  ELSE
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   kein Restart'
    l_shift_one_timestep=.TRUE.
  ENDIF
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR  num_prefetch_proc:', num_prefetch_proc
  !IF ((num_prefetch_proc == 1) .AND. (latbc_config%itype_latbc > 0)) THEN ! itype_latbc=1: time-dependent lateral boundary conditions; setting for limited area mode
  !  CALL init_prefetch(latbc) ! R: ich glaube hier wird latbc ALLOCATET und mit Werten gefüllt! TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching
  !ENDIF

  !------------------------------------------------------------------
  ! Prepare output file
  !------------------------------------------------------------------

  CALL configure_io()   ! set n_chkpt and n_diag, which control
                        ! writing of restart files and tot_int diagnostics.

  ! Add a special metrics variable containing the area weights of
  ! the regular lon-lat grid.
  !CALL compute_lonlat_area_weights(lonlat_grids, p_nh_state_lists) ! R: wahrscheinlich optional

  ! Map the variable groups given in the output namelist onto the
  ! corresponding variable subsets:
  ! ATTENTION: all add_vars must be finished before calling this routine.

  ! In Namelist nonhydrostatic_nml: Type of vertical coordinate:
  ! 1: Gal-Chen hybrid
  ! 2: SLEVE (uses sleve_nml, this is default,). If used we have to add the sleve_nml !!!
  IF (output_mode%l_nml) THEN
    CALL parse_variable_groups()
  END IF
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver L'

  ! if output on z and/or p-levels is required do some config
  !
  ! Note that on the compute PEs we must call this *before* the
  ! initialization of the vertical interpolation,
  ! "pp_scheduler_init"
  IF (output_mode%l_nml) THEN
    CALL collect_requested_ipz_levels()
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver L2'
    DO jg = 1, n_dom
      WRITE(*,*) jg, nproma, p_patch(jg)%npromz_c,  p_patch(jg)%nblks_c
      CALL configure_nh_pzlev(jg, nproma, p_patch(jg)%npromz_c,  &
        &                     p_patch(jg)%nblks_c)
    ENDDO
  END IF

  IF (output_mode%l_nml) THEN
    ! compute sim_start, sim_end <--------- siehe oben, warum hier nochmal?
    CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start)
    CALL datetimeToString(time_config%tc_exp_stopdate, sim_step_info%sim_end)
    CALL datetimeToString(time_config%tc_startdate, sim_step_info%run_start)
    CALL datetimeToString(time_config%tc_stopdate, sim_step_info%restart_time)
    !WRITE(*,*) time_config%tc_exp_startdate, sim_step_info%sim_start
    sim_step_info%dtime  = dtime
    jstep0 = 0

    restartAttributes => getAttributesForRestarting()
    IF (ASSOCIATED(restartAttributes)) THEN
      ! get start counter for time loop from restart file:
      jstep0 = restartAttributes%getInteger("jstep")
    END IF

    sim_step_info%jstep0    = jstep0
    CALL init_statistics_streams ! R: sonst Abbruch
    !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver L5B jstep0', jstep0
    CALL init_name_list_output(sim_step_info)
    !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver L6'
    CALL create_mipz_level_selections(output_file)
    CALL create_vertical_axes(output_file)
  END IF

  !Anurag Dipankar, MPIM (2015-08-01): always call this routine
  !for LES simulation
  !DO jg = 1 , n_dom
  !  atm_phy_nwp_config(jg)%lcalc_moist_integral_avg  = atm_phy_nwp_config(jg)%lcalc_moist_integral_avg &
  !       .OR. atm_phy_nwp_config(jg)%is_les_phy
  !END DO

  !----------------------!
  !  Initialize actions  !
  !----------------------!

  ! Initialize reset-Action, i.e. assign variables to action object
  !CALL reset_act%initialize(ACTION_RESET)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aus perform_nh_stepping
  !LOGICAL :: is_les_phy          !>TRUE is turbulence is 3D 
  !                               !>FALSE otherwise
  ! --------------------
  ! AUS    CALL perform_nh_stepping( time_config%tc_current_date, latbc )
  ! R: in Namelist run_nml in runskript ist iforcing = 2. Das ist für JSBACH4 Vorraussetzung und die anderen Fälle brauche ich nicht berücksichtigen!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   AUS    CALL allocate_nh_stepping (mtime_current):
  ALLOCATE(jstep_adv(n_dom), STAT=error_status )
  IF (error_status /= SUCCESS) THEN
    CALL finish ('RR: perform_nh_stepping',           &
      &      'allocation for jstep_adv failed' )
  ENDIF

  !   ! allocate flow control variables for transport and slow physics calls
  !   ALLOCATE(linit_dyn(n_dom), STAT=error_status )
  !   IF (error_status /= SUCCESS) THEN
  !     CALL finish ('RR: perform_nh_stepping',           &
  !     &      'allocation for flow control variables failed' )
  !   ENDIF

  ! initialize
  restartAttributes => getAttributesForRestarting()
  IF (ASSOCIATED(restartAttributes)) THEN
    ! Get attributes from restart file
    DO jg = 1,n_dom
      WRITE(attname,'(a,i2.2)') 'ndyn_substeps_DOM',jg
      ndyn_substeps_var(jg) = restartAttributes%getInteger(TRIM(attname))
      WRITE(attname,'(a,i2.2)') 'jstep_adv_marchuk_order_DOM',jg
      jstep_adv(jg)%marchuk_order = restartAttributes%getInteger(TRIM(attname))
      ENDDO
    !linit_dyn(:)      = .FALSE.
  ELSE
    jstep_adv(:)%marchuk_order = 0
    !linit_dyn(:)      = .TRUE.
  ENDIF

  !IF (l_limited_area .AND. latbc_config%itype_latbc > 0 .AND. num_prefetch_proc == 0) THEN
  !  CALL prepare_latbc_data(p_patch(1), p_int_state(1), p_nh_state(1), ext_data(1))
  !ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! wieder aus    CALL perform_nh_stepping
  ! Compute diagnostics for output if necessary
  CALL diag_for_output_dyn () ! R: nur debuggen


  !  DO jg=1, n_dom
  !     CALL compute_airmass(p_patch(jg),                  &
  !       &                  p_nh_state(jg)%metrics,       &
  !       &                  p_nh_state(jg)%prog(nnow(jg)),&
  !       &                  p_nh_state(jg)%diag, itlev = 2)
  ! 
  !     ! initialize exner_pr if the model domain is active
  !     IF (p_patch(jg)%ldom_active .AND. .NOT. isRestart()) CALL init_exner_pr(jg, nnow(jg))
  !   ENDDO

  !------------------------------------------------------------------
  !  get and write out some of the initial values
  !------------------------------------------------------------------
  mtime_current => time_config%tc_current_date ! R: außschließlich in atm_phy_psrad/mo_ps_radiation_model.f90:130
  IF (.NOT.isRestart() .AND. (mtime_current >= time_config%tc_exp_startdate)) THEN

    !     !--------------------------------------------------------------------------
    !     ! loop over the list of internal post-processing tasks, e.g.
    !     ! interpolate selected fields to p- and/or z-levels
    !     simulation_status = new_simulation_status(l_first_step   = .TRUE.,                  &
    !       &                                       l_output_step  = .TRUE.,                  &
    !       &                                       l_dom_active   = p_patch(1:)%ldom_active, &
    !       &                                       i_timelevel_dyn= nnow, i_timelevel_phy= nnow_rcf)
    !     CALL pp_scheduler_process(simulation_status) 
    CALL update_statistics ! R: notwendig siehe Rechnertext
    ! R: only for debugging
    ! IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
    !   CALL update_opt_acc(p_nh_opt_diag(1)%acc,            &
    !     &                 p_nh_state(1)%prog(nnow_rcf(1)), &
    !     &                 p_nh_state(1)%prog(nnow(1))%rho, &
    !     &                 p_nh_state(1)%diag,              &
    !     &                 p_patch(1)%cells%owned,          &
    !     &                 p_patch(1)%nlev                  )
    ! END IF

    ! R: das gibt die initvalus vor der Zeitloop d.h. ohne Berechungen aus!!! Im Driver gibt es aber noch keine Initwerte daher sollte das wieder raus!!!
    IF (output_mode%l_nml) THEN
      CALL write_name_list_output(jstep=0)
    END IF

    !     IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
    !       CALL reset_opt_acc(p_nh_opt_diag(1)%acc)
    !     END IF
   
    CALL reset_statistics

    !     ! sample meteogram output
    !     DO jg = 1, n_dom
    !       IF (.NOT. output_mode%l_none .AND. &    ! meteogram output is not initialized for output=none
    !         & p_patch(jg)%ldom_active  .AND. &
    !         & meteogram_is_sample_step( meteogram_output_config(jg), 0 ) ) THEN
    !         CALL meteogram_sample_vars(jg, 0, time_config%tc_startdate)
    !       END IF
    !     END DO
    ! 
    !     !AD: Also output special diagnostics for LES on torus
    !     IF (atm_phy_nwp_config(1)%is_les_phy &
    !       .AND. sampl_freq_step>0)THEN
    !       CALL calculate_turbulent_diagnostics(                      &
    !                              & p_patch(1),                       & !in
    !                              & p_nh_state(1)%prog(nnow(1)),      &
    !                              & p_nh_state(1)%prog(nnow_rcf(1)),  & !in
    !                              & p_nh_state(1)%diag,                   & !in
    !                              & p_lnd_state(1)%prog_lnd(nnow_rcf(1)), &
    !                              & p_lnd_state(1)%diag_lnd,              &
    !                              & prm_nwp_tend(1),                      &
    !                              & prm_diag(1)                )     !inout
    ! 
    !       !write out time series
    !       CALL write_time_series(prm_diag(1)%turb_diag_0dvar, mtime_current)
    !       CALL write_vertical_profiles(prm_diag(1)%turb_diag_1dvar, mtime_current, 1)
    !       prm_diag(1)%turb_diag_1dvar = 0._wp
    !     END IF
  END IF ! not isRestart()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aus CALL perform_nh_timeloop (mtime_current, latbc)
  ! --------------------
  ! AUS    CALL perform_nh_timeloop (mtime_current, latbc) !Organizes nonhydrostatic time stepping
  ! IF (ltimer) CALL timer_start(timer_total)
  ! calculate elapsed simulation time in seconds
  sim_time = getElapsedSimTimeInSeconds(mtime_current) 


  ! allocate temporary variable for restarting purposes
  IF (output_mode%l_nml) THEN
    ALLOCATE(output_jfile(SIZE(output_file)), STAT=error_status)
    IF (error_status /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed!')
  ENDIF
 
  ! output_file ist allocatet hier!
  
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver P'

  jstep_shift = 0  ! 0 is default. R: only different from 0 if Incremental Analysis Update (IAU)
  !   IF (timeshift%dt_shift < 0._wp  .AND. .NOT. isRestart()) THEN
  !     jstep_shift = NINT(timeshift%dt_shift/dtime)
  !     WRITE(message_text,'(a,i6,a)') 'Model start shifted backwards by ', ABS(jstep_shift),' time steps'
  !     CALL message(TRIM(routine),message_text)
  !     !atm_phy_nwp_config(:)%lcalc_acc_avg = .FALSE.
  !   ELSE
  !     jstep_shift = 0
  !   ENDIF

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Q '
  mtime_old     => newDatetime(mtime_current) ! in
  mtime_old_jsb => newDatetime(mtime_current) ! in
  DO jg=1, n_dom
    datetime_current(jg)%ptr => newDatetime(mtime_current) ! R: für jedes Modell datetime_current setzen
  END DO

  restartDescriptor => createRestartDescriptor("atm")
  !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver R'

  ! R. get jstep0: Das wird oben schon 2 mal gemacht. ist das wirklich notwendig???
  jstep0 = 0
  restartAttributes => getAttributesForRestarting()
  IF (isRestart()) THEN
    ! get start counter for time loop from restart file:
    jstep0 = restartAttributes%getInteger("jstep")
  END IF

  ! R: gibt viel Output:
  !  DO i = 1, nvar_lists
  !    CALL print_var_list(var_lists(i), .FALSE.) ! .TRUE: for short varlist printing; .FALSE. for long varlist printing
  !  ENDDO
 
  ! init routine for mo_nh_supervise module (eg. opening of files)
  !CALL init_supervise_nh() ! Print-out of max winds to an ASCII file. if output_mode%l_maxwinds
 
  ! set events, group and the events
  CALL message('','')
  eventStartDate => time_config%tc_exp_startdate
  eventEndDate   => time_config%tc_exp_stopdate
  ! for debugging purposes the referenece (anchor) date for checkpoint
  ! and restart may be switched to be relative to current jobs start
  ! date instead of the experiments start date.
  IF (time_config%is_relative_time) THEN
    checkpointRefDate => time_config%tc_startdate
    restartRefDate    => time_config%tc_startdate
  ELSE
    checkpointRefDate => time_config%tc_exp_startdate
    restartRefDate    => time_config%tc_exp_startdate
  ENDIF
  ! --- create an event group for checkpointing and restart
  checkpointEvents =  addEventGroup('checkpointEventGroup')
  checkpointEventGroup => getEventGroup(checkpointEvents)
  ! --- --- create checkpointing event
  eventInterval  => time_config%tc_dt_checkpoint
  checkpointEvent => newEvent('checkpoint', checkpointRefDate, eventStartDate, eventEndDate, eventInterval, errno=error_status)
  IF (error_status /= NO_ERROR) THEN
    ! give an elaborate error message:
    CALL datetimeToString(checkpointRefDate, dt_string)
    WRITE (0,*) "event reference date: ",    dt_string
    CALL datetimeToString(eventStartDate,    dt_string)
    WRITE (0,*) "event start date    : ",    dt_string
    CALL datetimeToString(eventEndDate,      dt_string)
    WRITE (0,*) "event end date      : ",    dt_string
    CALL timedeltaToString(eventInterval,    td_string)
    WRITE (0,*) "event interval      : ",    td_string
    CALL mtime_strerror(error_status, errstring)
    CALL finish('perform_nh_timeloop', "event 'checkpoint': "//errstring)
  ENDIF
  lret = addEventToEventGroup(checkpointEvent, checkpointEventGroup)
  ! --- --- create restart event, ie. checkpoint + model stop
  eventInterval  => time_config%tc_dt_restart
  restartEvent => newEvent('restart', restartRefDate, eventStartDate, eventEndDate, eventInterval, errno=error_status)
  IF (error_status /= NO_ERROR) THEN
    ! give an elaborate error message:
    CALL datetimeToString(restartRefDate, dt_string)
    WRITE (0,*) "event reference date: ", dt_string
    CALL datetimeToString(eventStartDate, dt_string)
    WRITE (0,*) "event start date    : ", dt_string
    CALL datetimeToString(eventEndDate,   dt_string)
    WRITE (0,*) "event end date      : ", dt_string
    CALL timedeltaToString(eventInterval, td_string)
    WRITE (0,*) "event interval      : ", td_string
    CALL mtime_strerror(error_status, errstring)
    CALL finish('perform_nh_timeloop', "event 'restart': "//errstring)
  ENDIF
  lret = addEventToEventGroup(restartEvent, checkpointEventGroup)
  CALL printEventGroup(checkpointEvents)

  ! set time loop properties
  model_time_step => time_config%tc_dt_model
  
  CALL message('','')
  CALL datetimeToString(mtime_current, dstring)
  WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
  CALL message('',message_text)
  CALL datetimeToString(time_config%tc_stopdate, dstring)
  WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
  CALL message('',message_text)
  CALL message('','')

  jstep = jstep0+jstep_shift+1

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U jstep0', jstep0
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U jstep', jstep

  counter = 0 ! R: NUR ZUM DEBUGGEN

  ! -------------------------------------------------------------------------------------------------------------------------------
  ! Time loop
  ! -------------------------------------------------------------------------------------------------------------------------------
  !loop_step = 0 ! R: timestep counter within the loop
  TIME_LOOP: DO ! atm_dyn_icoham/mo_ha_stepping.f90: SUBROUTINE perform_nh_timeloop (mtime_current, latbc): Zeile 699 vor integrate_nh
    !loop_step = loop_step + 1
    !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U1 jstep in loop', jstep ! geht z.B. von 1 bis 72 sollte forcing_steps_per_file entsprechen????
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U1 mtime_current at start loop', mtime_current     ! z.B. für jstep 72:  1979           1           2           0           0    0           0
    ! WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U1 model_time_step in loop', model_time_step ! z.B. für jstep 72:     1 +  0           0           0           0          20        0           0
    ! optional memory loggin
    !CALL memory_log_add

    ! Check if a nested domain needs to be turned off
    !DO jg=2, n_dom
    !  IF (p_patch(jg)%ldom_active .AND. (sim_time >= end_time(jg))) THEN
    !    p_patch(jg)%ldom_active = .FALSE.
    !    WRITE(message_text,'(a,i2,a,f12.2)') 'domain ',jg,' stopped at time ',sim_time
    !    CALL message('in TIME_LOOP', TRIM(message_text))
    !  ENDIF
    !ENDDO

    ! update model date and time mtime based
    IF (l_first_forcing_file) THEN
      mtime_forcing => newDatetime(mtime_current) ! mtime_current is here the init time step before the model calculations.
                                                  ! The forcing input is at exp start one step before the calculation step
    ELSE
      mtime_forcing => newDatetime(mtime_current)
      mtime_forcing = mtime_forcing + model_time_step ! Changed this for debugging the synchron. now it is synchron!
    END IF
    mtime_old_jsb = mtime_current
    mtime_current = mtime_current + model_time_step
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U1 mtime_current after update at begin of loop', mtime_current
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U1 mtime_forcing after update at begin of loop', mtime_forcing
    !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver U3'
    ! store state of output files for restarting purposes
    IF (output_mode%l_nml .AND. jstep>=0 ) THEN
      DO i=1,SIZE(output_file)
        output_jfile(i) = get_current_jfile(output_file(i)%out_event)
      END DO
    ENDIF
 
    ! turn on calculation of averaged and accumulated quantities at the first regular time step
    !delete    IF (jstep-jstep0 == 1) atm_phy_nwp_config(:)%lcalc_acc_avg = .TRUE.

    ! R: Nur für limited area
    ! read boundary data if necessary
    !IF (l_limited_area .AND. latbc_config%itype_latbc > 0 .AND. num_prefetch_proc == 0) &
    !  CALL read_latbc_data_sync(p_patch(1), p_nh_state(1), ext_data(1), p_int_state(1), mtime_current)

    ! always print the first and the last time step
    lprint_timestep = .TRUE.
    lprint_timestep = lprint_timestep .OR. (jstep == jstep0+1) .OR. (jstep == jstep0+nsteps)
    IF (lprint_timestep) THEN
        WRITE(message_text,'(a,i8,a,i0,a,4(i2.2,a),i2.2)') &
              &             'Time step: ', jstep, ' model time ',                                & !  Time step: 1     model time 1979-01-01 00:20:00
              &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
              &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
              &             mtime_current%time%minute, ':', mtime_current%time%second
      CALL message('',message_text)
    ENDIF

    !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver V'
    ! set namelist output flag
    l_nml_output = output_mode%l_nml .AND. jstep >= 0 .AND. istime4name_list_output(jstep)
    ! Computation of diagnostic quantities may also be necessary for
    ! meteogram sampling:
    !DR Note that this may be incorrect for meteograms in case that
    !DR meteogram_output_config is not the same for all domains.
    !l_compute_diagnostic_quants = l_nml_output
    !DO jg = 1, n_dom
    !  l_compute_diagnostic_quants = l_compute_diagnostic_quants .OR. &
    !    &          meteogram_is_sample_step(meteogram_output_config(jg), jstep )
    !END DO
    !l_compute_diagnostic_quants = jstep >= 0 .AND. l_compute_diagnostic_quants .AND. &
    !  &                           .NOT. output_mode%l_none

    ! Calculations for enhanced sound-wave and gravity-wave damping during the spinup phase
    ! if mixed second-order/fourth-order divergence damping (divdamp_order=24) is chosen.
    ! Includes increased vertical wind off-centering during the first 2 hours of integration.
    !IF (divdamp_order==24 .AND. .NOT. isRestart()) THEN
    !  elapsed_time_global = (REAL(jstep,wp)-0.5_wp)*dtime
    !  IF (elapsed_time_global <= 7200._wp+0.5_wp*dtime .AND. .NOT. ltestcase) THEN
    !    CALL update_spinup_damping(elapsed_time_global)
    !  ENDIF
    !ELSE IF (divdamp_order==24) THEN
    !  divdamp_fac_o2 = 0._wp
    !ENDIF

    !--------------------------------------------------------------------------
    !
    ! dynamics stepping
    !
    ! R: ich sehe kaum etwas das evtl. aus dieser Prozedur gebraucht würde.
    ! R: In Runskript nichts angegeben für itime_scheme: ist 14 by default!
    !    In Runskript run_nml: ltimer = .FALSE. . Timer for monitoring the runtime of specific routines
    !CALL integrate_nh(datetime_current, 1, jstep-jstep_shift, iau_iter, dtime, model_time_step, 1, latbc) ! R: iau_iter=0 wenn iterate_iau = FALSE in initicon_nml, default ist false

    !!!!!!!!!!!!!!!!!!!!! integrate_nh START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate elapsed simulation time in seconds (local time for
    ! this domain!)
    datetime_local = datetime_current ! R: because we do not have nested domains now. , wird so in ICON übergeben
    sim_time = getElapsedSimTimeInSeconds(datetime_local(jg)%ptr) 

    ! R: Normalerweise wäre hier irgendwo in ICON-A eine model loop, da wir aber immer nur ein Modell verwenden fällt das weg

    !R: Wozu brauche ich jpg überhaut?
    ! Determine parent domain ID
    jg    = 1 ! ist im SRC Code für den 1. globalen aufruf so fixiert. Da sich integrat_nh aber rekursiv aufruft würde sich dies für jedes Childgitter ändern!
    IF ( jg > 1) THEN
      jgp = p_patch(jg)%parent_id
    ELSE IF (n_dom_start == 0) THEN
      jgp = 0
    ELSE
      jgp = 1
    ENDIF

    ! This executes one time step for the global domain and two steps for nested domains
    !num_steps = 2
    !DO jsubstep = 1, num_steps 

    ! simplified setting (may be removed lateron)
    jstep_adv(jg)%marchuk_order = jstep_adv(jg)%marchuk_order + 1
    !p_patch(jg)%n_childdom=1 ! R: ist 0. parallel_infrastructure/mo_setup_subdivision.f90:1565:      wrk_p_patch%n_childdom         = wrk_p_patch_pre%n_childdom

    ! Hier wird interface_iconam_echam gecallt!
    !!!!!! Ab hier: SUBROUTINE interface_iconam_echam
    !TYPE(datetime)          , POINTER    :: datetime_old      !< date and time at the beginning of this time step
    
    !datetime_local(jg)%ptr = datetime_local(jg)%ptr + mtime_dt_loc
    !CALL interface_iconam_echam( dt_loc ,datetime_local(jg)%ptr...) wird an datetime_new übergeben...
    !datetime_old      => newDatetime(datetime_new)
    !datetime_old      =  datetime_new + neg_dt_loc_mtime
    ! Berechungen des Orbits:
    !patch=p_patch(jg) ! R: wenn child domains dann ist das keine 1 hier mehr. Das muß rein wg. xdec, was ich bisher fixiert habe!
    !CALL echam_phy_bcs( patch        ,&! in
    !  &                 datetime_old ,&! in
    !  &                 dt_loc       ) ! in


    
    CALL jsbach_start_timestep(jg, mtime_old_jsb, dtime) ! mtime_current
    !delete WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver X'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! R: Das Folgende ist aus echam_phy_main
    !...
    ! R: Das Folgende ist aus echam_phy_main und dort aus CALL omp_loop_cell_prog(patch, interface_echam_vdf      ,&
    !            &                  is_in_sd_ed_interval, is_active ,&
    !            &                  datetime_old, pdtime            )

    ! START ALLOCATE BLOCK SOLL AUS TIMELOOP RAUS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Block loop 1 um nblks zu ermitteln (corresponding to SUBROUTINE omp_loop_cell_prog)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Get necessary informations from the patch of the actual thread for jsbach 
    patch=p_patch(1) ! R: Die Einträge in Patch können für jeden Prozessor unterschiedlich sein!
    jg  = patch%id ! see above R: kann wohl wieder raus

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1) ! R: jbs ist immer =1
    jbe = patch%cells%  end_blk(rle,ncd) ! R: jbe ist unterschiedlich. Z.B. 244744400, 251860448, 176309936, -2013259902 (und das bei 20480 Gitterpunkten!)

    !     !nlev= patch%nlev
    !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR  DRIVER jbs' , jbs !  1
    !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR  DRIVER jbe:', jbe !im std-gekoppelten-Lauf ist: jbe immer =39 , hier ist es aber 61
!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    nblks=0
    DO jb = jbs,jbe ! loop über den block index
      nblks = nblks + 1
      IF (jb == 1)  THEN        ! R: we want only jce max
        CALL get_indices_c(patch,  & ! in 
            &                   jb,  & ! in  index of block loop
            &                   jbs, & ! in  start indice of block loop
            &                   jbe, & ! in  end   indice of block loop
            &                   jcs, & ! out start indice of columns loop
            &                   jce, & ! out end   indice of columns loop
            &                   rls, & ! in  row limit start
            &                   rle)   ! in  row limit end
        nproma = jce
        !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   nproma', nproma, 'jce', jce ! beides immer 16!
      END IF
    END DO
!$OMP END PARALLEL DO
    !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   nproma', nproma, 'jce', jce ! 
    ! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR   ALLOCATEN muß nach der Ermittlung von nblks gemacht werden!!!
    !call model%Get_top_tile(tile)
    !hgrid => Get_grid(model%grid_id)
    !jsb_nproma = hgrid%nproma
    !nblks  = hgrid%nblks is the same as nblks  = get_nblks (hgrid)
    !nblks  = p_patch(1)%nblks_c ! R: ich nehme es direkt aus dem Hauptpatch ist wohl für jeden Prozessor unterschiedlich
    ! nblks wird nicht direkt hier oben bestimmt, muß ich woanderst her holen:
    !   nblks  = hgrid%nblks
    !nblks  = p_patch(1)%cell%nblks_c
    !nblks=20480
    ! R: 20480 Gitterpunkt bei R2B4. aproma = 16, es sollte also 1280 Blocks global geben. Auf 32 Prozessoren aufgeteilt wären das 40 Blocks pro Prozessor!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Block loop 1 um nblks zu ermitteln (corresponding to SUBROUTINE omp_loop_cell_prog)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (.NOT.(ALLOCATED(t_air))) THEN
      ! R: brauche ich wahrscheinlich nicht als lokale Variablen hier da sie von jsbach_init -> call tile%init angelegt werden. !!!!!
      ! Vars used in interface_full
      ! in vars (forcing) not necessary:
      ALLOCATE(DEBUG_VAR (nproma, nblks))
      ALLOCATE(t_air              (nproma, nblks)) ! R: e.g. (16, 50)
      ALLOCATE(q_air              (nproma, nblks))
      ALLOCATE(rain               (nproma, nblks))
      ALLOCATE(snow               (nproma, nblks))
      ALLOCATE(wind_air           (nproma, nblks))
      ALLOCATE(wind_10m           (nproma, nblks))
      ALLOCATE(lw_srf_down        (nproma, nblks))
      ALLOCATE(swvis_srf_down     (nproma, nblks))
      ALLOCATE(swnir_srf_down     (nproma, nblks))
      ALLOCATE(swpar_srf_down     (nproma, nblks))
      ALLOCATE(fract_par_diffuse  (nproma, nblks)) ! so far not used in jsb4 
      ALLOCATE(press_srf          (nproma, nblks))
      ALLOCATE(drag_srf           (nproma, nblks))
      ALLOCATE(t_acoef            (nproma, nblks))
      ALLOCATE(t_bcoef            (nproma, nblks))
      ALLOCATE(q_acoef            (nproma, nblks))
      ALLOCATE(q_bcoef            (nproma, nblks))
      ALLOCATE(pch                (nproma, nblks))
      ALLOCATE(cos_zenith_angle   (nproma, nblks))
      ALLOCATE(CO2_air            (nproma, nblks)) ! so far not used in jsb4
      ALLOCATE(trig_rad_m1        (nproma, nblks))
      ALLOCATE(newday_Merker        (nproma, nblks))
      ALLOCATE(l_day_reminder       (nproma, nblks))
      ! out vars, necessary for standalone:
      ALLOCATE(t_srf_fromJS4         (nproma, nblks))
      ALLOCATE(t_eff_srf    (nproma, nblks))
      ALLOCATE(qsat_srf     (nproma, nblks))
      ALLOCATE(s_srf        (nproma, nblks))
      ALLOCATE(fact_q_air_fromJS4    (nproma, nblks))
      ALLOCATE(fact_qsat_srf_fromJS4 (nproma, nblks))
      ALLOCATE(evapotrans   (nproma, nblks))
      ALLOCATE(latent_hflx  (nproma, nblks))
      ALLOCATE(sensible_hflx(nproma, nblks))
      ALLOCATE(grnd_hflx    (nproma, nblks))
      ALLOCATE(grnd_hcap    (nproma, nblks))
      ALLOCATE(rough_h_srf_fromJS4   (nproma, nblks))
      ALLOCATE(rough_m_srf_fromJS4   (nproma, nblks))
      ALLOCATE(q_snocpymlt  (nproma, nblks))
      ALLOCATE(alb_vis_dir  (nproma, nblks))
      ALLOCATE(alb_nir_dir  (nproma, nblks))
      ALLOCATE(alb_vis_dif  (nproma, nblks))
      ALLOCATE(alb_nir_dif  (nproma, nblks))
      ALLOCATE(CO2_flux     (nproma, nblks))
      ! in vars for lakes, optional. Not necessary:
      ALLOCATE( drag_wtr          (nproma, nblks))
      ALLOCATE( drag_ice          (nproma, nblks))
      ALLOCATE( t_acoef_wtr       (nproma, nblks))
      ALLOCATE( t_bcoef_wtr       (nproma, nblks))
      ALLOCATE( q_acoef_wtr       (nproma, nblks))
      ALLOCATE( q_bcoef_wtr       (nproma, nblks))
      ALLOCATE( t_acoef_ice       (nproma, nblks))
      ALLOCATE( t_bcoef_ice       (nproma, nblks))
      ALLOCATE( q_acoef_ice       (nproma, nblks))
      ALLOCATE( q_bcoef_ice       (nproma, nblks))
      ! out vars for lakes, necessary for standalone, optional:
      ALLOCATE( t_lwtr            (nproma, nblks))
      ALLOCATE( qsat_lwtr         (nproma, nblks))
      ALLOCATE( s_lwtr            (nproma, nblks))
      ALLOCATE( evapo_wtr         (nproma, nblks))
      ALLOCATE( latent_hflx_wtr   (nproma, nblks))
      ALLOCATE( sensible_hflx_wtr (nproma, nblks))
      ALLOCATE( albedo_lwtr       (nproma, nblks))
      ALLOCATE( t_lice            (nproma, nblks))
      ALLOCATE( qsat_lice         (nproma, nblks))
      ALLOCATE( s_lice            (nproma, nblks))
      ALLOCATE( evapo_ice         (nproma, nblks))
      ALLOCATE( latent_hflx_ice   (nproma, nblks))
      ALLOCATE( sensible_hflx_ice (nproma, nblks))
      ALLOCATE( albedo_lice       (nproma, nblks))
      ALLOCATE( ice_fract_lake    (nproma, nblks))
      ! For standalone echam-jsbach
      ALLOCATE(evapopot           (nproma, nblks)) ! R: warum nblks raus genommen? auch unten im call des interfaces macht das keinen Sinn...

      ! dummy output variables required for next timesteps forcing (if not reading interface vars)
      !  ALLOCATE(evapo_act_sum_fromJS4 (nproma, nblks))
      !  ALLOCATE(evapo_pot_sum_fromJS4 (nproma, nblks))

      !Set dummy values.   RRR Kann ich die zum Allocaten hoch schieben??? RRRRRR
      ! in vars: not necessary, set these dummy values only for testing purpose:
!       t_air             = 290        ! Kelvin
!       q_air             = 0.1        ! kg kg-1
!       rain              = 0.00004    ! kg/m2/s
!       snow              = 0.000001   ! kg/m2/s
!       wind_air          = 0.1        ! air_wind_speed, m/s
!       wind_10m          = 0.2        ! air_wind_speed, m/s
!       lw_srf_down       = 340        ! W/m2   
!       swvis_srf_down    = 120        ! W/m2
!       swnir_srf_down    = 100        ! W/m2
!       swpar_srf_down    = 60         ! fictional only for testing purpuses
!       fract_par_diffuse = 0.5        ! fictional only for testing purpuses
!       press_srf         = 101325.0   ! N/m2 =Pa
!       drag_srf          = 0          ! fictional only for testing purpuses
!       t_acoef           = 0.0000     ! fictional only for testing purpuses
!       t_bcoef           = 0.0000     ! fictional only for testing purpuses
!       q_acoef           = 0.0000     ! fictional only for testing purpuses
!       q_bcoef           = 0.0000     ! fictional only for testing purpuses
!       pch               = 0.0000     ! fictional only for testing purpuses
!       cos_zenith_angle  = 0.5        ! fictional only for testing purpuses
!       CO2_air           = 0.00052868 ! this is in mol: 0.00034800
!       trig_rad_m1       = 0          !
!       newday_Merker       = 4          !
!       l_day_reminder = 
      ! out vars, necessary but only dummy values for standalone:
      t_srf_fromJS4         = 0   ! fictional
      t_eff_srf          = 0   ! fictional
      qsat_srf           = 0   ! fictional
      s_srf              = 0   ! fictional
      fact_q_air_fromJS4    = 0   ! fictional
      fact_qsat_srf_fromJS4 = 0   ! fictional
      evapotrans         = 0   ! fictional! R: land2atmo Umrechung
      latent_hflx        = 0   ! fictional! R: land2atmo Umrechung
      sensible_hflx      = 0   ! fictional! R: land2atmo Umrechung
      grnd_hflx          = 0   ! fictional
      grnd_hcap          = 0   ! fictional
      rough_h_srf_fromJS4   = 0   ! fictional
      rough_m_srf_fromJS4   = 0   ! fictional
      q_snocpymlt        = 0   ! fictional! R: land2atmo Umrechung
      alb_vis_dir        = 0   ! fictional
      alb_nir_dir        = 0   ! fictional
      alb_vis_dif        = 0   ! fictional
      alb_nir_dif        = 0   ! fictional
      CO2_flux           = 0   ! fictional! R: land2atmo Umrechung
      ! in vars for lakes, optional: Not necessary only for testing purpose:
!       drag_wtr          =  0.0000  ! fictional only for testing purpuses
!       drag_ice          =  0.0000  ! fictional only for testing purpuses
!       t_acoef_wtr       =  0.0000  ! fictional only for testing purpuses
!       t_bcoef_wtr       =  0.0000  ! fictional only for testing purpuses
!       q_acoef_wtr       =  0.0000  ! fictional only for testing purpuses
!       q_bcoef_wtr       =  0.0000  ! fictional only for testing purpuses
!       t_acoef_ice       =  0.0000  ! fictional only for testing purpuses
!       t_bcoef_ice       =  0.0000  ! fictional only for testing purpuses
!       q_acoef_ice       =  0.0000  ! fictional only for testing purpuses
!       q_bcoef_ice       =  0.0000  ! fictional only for testing purpuses
      ! out vars for lakes, necessary for standalone, but only dummy values, optional:
      t_lwtr            =  0.0000  ! fictional
      qsat_lwtr         =  0.0000  ! fictional
      s_lwtr            =  0.0000  ! fictional
      evapo_wtr         =  0.0000  ! fictional! R: land2atmo Umrechung
      latent_hflx_wtr   =  0.0000  ! fictional! R: land2atmo Umrechung
      sensible_hflx_wtr =  0.0000  ! fictional! R: land2atmo Umrechung
      albedo_lwtr       =  0.0000  ! fictional
      t_lice            =  0.0000  ! fictional
      qsat_lice         =  0.0000  ! fictional
      s_lice            =  0.0000  ! fictional
      evapo_ice         =  0.0000  ! fictional
      latent_hflx_ice   =  0.0000  ! fictional! R: land2atmo Umrechung
      sensible_hflx_ice =  0.0000  ! fictional ! R: land2atmo Umrechung
      albedo_lice       =  0.0000  ! fictional
      ice_fract_lake    =  0.0000  ! fictional
      evapopot          =  0.0000  ! fictional ! R: land2atmo Umrechung

      ! WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver t_air:'       , forcing_in_modeltyp%t_air%data_read       (:, iblk, step_in_forcing_file), 'step_in_forcing_file', step_in_forcing_file, 'iblk', iblk
      ! WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver q_acoef_ice:' ,  forcing_in_modeltyp%q_acoef_ice%data_read(:, iblk, step_in_forcing_file)
      ! WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver q_bcoef_ice:' ,  forcing_in_modeltyp%q_bcoef_ice%data_read(:, iblk, step_in_forcing_file)

      IF ( l_forcing_observations ) THEN
        !IF ( model_time_step ) THEN
        !  CALL finish('The time step of the model must not be larger than one day if l_forcing_observations is true.')
        !END IF

        CALL init_forcing(model_id)
        write(*,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        ALLOCATE(lon(nproma,nblks))
        ALLOCATE(lat(nproma,nblks))

        ALLOCATE(coslon(nproma,nblks))
        ALLOCATE(coslat(nproma,nblks))
        ALLOCATE(sinlon(nproma,nblks))
        ALLOCATE(sinlat(nproma,nblks))
        ALLOCATE(daylight_frc(nproma, nblks))

        !ALLOCATE(CO2_concentration (nproma, nblks)) ! Observations: CO2_air in jsb_interface. Nicht mehr benötigt ist jetzt auch in mo_forcing CO2_air!
        ALLOCATE(evapo_act2pot_fromJS4generated (nproma, nblks)) ! für Obsevations
        !ALLOCATE(t_air (nproma, nblks))

        !init for first timestep/day (according to jsb3: mo_jsbalone: init_driving)
        t_srf_fromJS4 = 273._wp
        fact_q_air_fromJS4 = 0.5_wp
        fact_qsat_srf_fromJS4 = 0.5_wp
        rough_h_srf_fromJS4 = 1._wp
        rough_m_srf_fromJS4 = 1._wp
        evapo_act2pot_fromJS4generated = 1._wp

        ! evapo_act_sum_fromJS4 = 0._wp
        ! evapo_pot_sum_fromJS4 = 0._wp
      ENDIF

    END IF
    
    ! Allocate the forcing input vars
    IF (.NOT.(ALLOCATED(forcing_in_modeltyp%t_air%data_read))) THEN
      ALLOCATE(forcing_in_modeltyp%t_air%data_read             (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_air%data_read             (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%rain%data_read              (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%snow%data_read              (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%wind_air%data_read          (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%wind_10m%data_read          (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%lw_srf_down%data_read       (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%swvis_srf_down%data_read    (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%swnir_srf_down%data_read    (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%swpar_srf_down%data_read    (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%fract_par_diffuse%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%press_srf%data_read         (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%drag_srf%data_read          (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_acoef%data_read           (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_bcoef%data_read           (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_acoef%data_read           (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_bcoef%data_read           (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%pch%data_read               (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%cos_zenith_angle%data_read  (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%co2_air%data_read           (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%trig_rad_m1%data_read       (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%newday_Merker%data_read     (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%l_day_reminder%data_read    (nproma, nblks, forcing_steps_per_file))
      !for lakes
      ALLOCATE(forcing_in_modeltyp%drag_wtr%data_read    (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%drag_ice%data_read    (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_acoef_wtr%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_bcoef_wtr%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_acoef_wtr%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_bcoef_wtr%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_acoef_ice%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%t_bcoef_ice%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_acoef_ice%data_read (nproma, nblks, forcing_steps_per_file))
      ALLOCATE(forcing_in_modeltyp%q_bcoef_ice%data_read (nproma, nblks, forcing_steps_per_file))

!      ! R: debugging: set values before they are set from observational or model forcing. Delete later
!      forcing_in_modeltyp%t_air%data_read(:,:,:)             = SPREAD ( t_air            , DIM=3, NCOPIES=forcing_steps_per_file)
!      forcing_in_modeltyp%q_air%data_read(:,:,:)             = SPREAD ( q_air            , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%rain%data_read(:,:,:)              = SPREAD ( rain             , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%snow%data_read(:,:,:)              = SPREAD ( snow             , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%wind_air%data_read(:,:,:)          = SPREAD ( wind_air         , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%wind_10m%data_read(:,:,:)          = SPREAD ( wind_10m         , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%lw_srf_down%data_read(:,:,:)       = SPREAD ( lw_srf_down      , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%swvis_srf_down%data_read(:,:,:)    = SPREAD ( swvis_srf_down   , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%swnir_srf_down%data_read(:,:,:)    = SPREAD ( swnir_srf_down   , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%swpar_srf_down%data_read(:,:,:)    = SPREAD ( swpar_srf_down   , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%fract_par_diffuse%data_read(:,:,:) = SPREAD ( fract_par_diffuse, DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%press_srf%data_read(:,:,:)         = SPREAD ( press_srf        , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%drag_srf%data_read(:,:,:)          = SPREAD ( drag_srf         , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_acoef%data_read(:,:,:)           = SPREAD ( t_acoef          , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_bcoef%data_read(:,:,:)           = SPREAD ( t_bcoef          , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_acoef%data_read(:,:,:)           = SPREAD ( q_acoef          , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_bcoef%data_read(:,:,:)           = SPREAD ( q_bcoef          , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%pch%data_read(:,:,:)               = SPREAD ( pch              , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%cos_zenith_angle%data_read(:,:,:)  = SPREAD ( cos_zenith_angle , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%CO2_air%data_read(:,:,:)           = SPREAD ( CO2_air          , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%trig_rad_m1%data_read(:,:,:)       = SPREAD ( trig_rad_m1      , DIM=3, NCOPIES=forcing_steps_per_file) 
!         !
!      forcing_in_modeltyp%drag_wtr%data_read(:,:,:)    = SPREAD ( drag_wtr    , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%drag_ice%data_read(:,:,:)    = SPREAD ( drag_ice    , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_acoef_wtr%data_read(:,:,:) = SPREAD ( t_acoef_wtr , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_bcoef_wtr%data_read(:,:,:) = SPREAD ( t_bcoef_wtr , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_acoef_wtr%data_read(:,:,:) = SPREAD ( q_acoef_wtr , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_bcoef_wtr%data_read(:,:,:) = SPREAD ( q_bcoef_wtr , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_acoef_ice%data_read(:,:,:) = SPREAD ( t_acoef_ice , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%t_bcoef_ice%data_read(:,:,:) = SPREAD ( t_bcoef_ice , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_acoef_ice%data_read(:,:,:) = SPREAD ( q_acoef_ice , DIM=3, NCOPIES=forcing_steps_per_file) 
!      forcing_in_modeltyp%q_bcoef_ice%data_read(:,:,:) = SPREAD ( q_bcoef_ice , DIM=3, NCOPIES=forcing_steps_per_file) 
    END IF

! END ALLOCATE BLOCK SOLL AUS TIMELOOP RAUS

    !IF (my_process_is_stdio()) THEN ! führt dazu daß der output ewig geschrieben wird. Eigentlich sollte io prozess hier nie ankommen sondern oben schon abgetrennt werden.
    model => get_model(model_id)

    ! Read vars that are input for observations and are taken dynamically from the JS4 output
    ! Therefore we have to read them after the restart and after the var allocation.
    ! And we have to get them from the JS4 memory!
    IF ( l_forcing_observations ) THEN
      IF(jsbach_is_restarted()) THEN
        CALL model%Get_top_tile(tile)

        WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR Einlesen fuer Observations1...'
        dsl4jsb_Get_memory(SEB_) ! R: in SEB ist nur dann etwas drin wenn es aus Restart File gelesen wurde ansonsten sind die init Werte drin.
        dsl4jsb_Get_memory(TURB_)
        dsl4jsb_Get_memory(HYDRO_)
        !R: Debugg: Die Felder hier passen nicht zusammen, muß korrigiert werden!
        CALL model%Set_options(iblk=jb, ics=jcs, ice=jce)
        t_srf_fromJS4(:,jb)           = dsl4jsb_var2D_onChunk(SEB_, t)
        fact_q_air_fromJS4(:,jb)      = dsl4jsb_var2D_onChunk(TURB_, fact_q_air)
        fact_qsat_srf_fromJS4(:,jb)   = dsl4jsb_var2D_onChunk(TURB_, fact_qsat_srf)
        rough_h_srf_fromJS4(:,jb)     = dsl4jsb_var2D_onChunk(TURB_, rough_h)
        rough_m_srf_fromJS4(:,jb)     = dsl4jsb_var2D_onChunk(TURB_, rough_m)
        evapotrans(:,jb)              = dsl4jsb_var2D_onChunk(HYDRO_, evapotrans)
        evapopot(:,jb)                = dsl4jsb_var2D_onChunk(HYDRO_, evapopot)
        ! TODO: if really required calculate daily sums of evapotrans and evapopot in jsbach and get them here (see comment below)
        evapo_act2pot_fromJS4generated(:,jb) =  MIN(evapotrans(:,jb),-EPSILON(1._wp)) & ! wird nur für IF (.NOT. read_interface_vars) THEN bei Julia gebraucht!
          &                                     / MIN(evapopot(:,jb),  -EPSILON(1._wp))
        WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR Einlesen fuer Observations 2...'
      ENDIF
    ENDIF
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FORCING lesen mit dem oben ermittelten nblks
    remain = mod((jstep),( forcing_steps_per_file * forcing_synchron_factor ) )
    IF (  jstep == 1 .OR.  remain == 0  ) THEN  ! Check if jstep is the first time step or a mulible of the time steps * forcing_synchron_factor to read the forcing
      step_in_forcing_file = 1

      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver jstep:' , jstep, 'remain', remain, 'forcing_synchron_factor', forcing_synchron_factor
      WRITE (mtime_forcing_char,'(i4.4,2(i2.2),a,3(i2.2),a)')                                                 &
      &                     mtime_forcing%date%year,mtime_forcing%date%month,mtime_forcing%date%day, 'T',   &
      &                     mtime_forcing%time%hour,mtime_forcing%time%minute,mtime_forcing%time%second, 'Z'
      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Drivermtime_forcing_char:' , mtime_forcing_char ! 19790101T002000Z


      IF ( l_forcing_observations ) THEN
        write(*,*) 'RRRRRRRRRRRRRRRRRRR l_forcing_observations', l_forcing_observations
        !----------------Julias treiber
        !Wenn modellforcing benutzt wird:  liest den current timestep für jeden timestep in der loop ein!!!
        !Wenn observations  benutzt wird:  liest aus dem forcing file das maximal wohl 1 Tag enthalten darf (routine only implemented for constant, daily or timestep frequency.)
        !                                  alle Zeitschritte bis zum aktuellen Zeitschritt in der time loop ein
        !                                     get_interface_variables_from_external_forcing ruft (für jeden Timestep in loop) die SUBROUTINE get_data_and_update_field auf.
        !                                     Darin: 
        !                                           Schleife die currentTimeStep hochzählt: currentTimeStep = currentTimeStep + 1. Diese Schleife vergleicht in jeder Iteration, ob 
        !                                           currentTimeStep noch kleiner gleich der Anzahl der Zeitschritte im forcing file entspricht.
        !                                           Wenn das ydate (integer aus year+month+day, also der Tagesangabe) von mo_time_control: next_date entspricht,
        !                                             dann wird thisDayEnd = currentTimeStep gesetzt.
        !                                             ansonsten wird die schleife verlassen.
        !                                           Die Schleife wird also so lange hoch gezählt bis entweder der Zeitschritt im File zuende ist oder der next_date Tag im File erreicht wird.
        !                                           Außerdem wir in der Schleife thisDayStart immer auf 1 d.h. den ersten currentTimeStep festgelegt.
        !                                           Erst nach dieser Schleife wird das File eingelesen mit: 
        !                                              return_pointer => input_file%Read_2d_time(forcing_input_variable%variable_name, &
        !                                                & start_time_step=thisDayStart, end_time_step=thisDayEnd, fill_array=forcing_input_variable%data_read)
        !
        !                                           => Es wird also immer genau 1 Tag (mit beliebig vielen Zeitschritten darin) eingelesen.
        !                                              Würde die Routine öfters gecallt würde das dennoch funktionieren, dann wird der Tag eben nochmal eingelesen.
        !                                              Forcing Files sind bei Climate_dCRUNCEP_T63_1975.nc jährliche Files mit täglichen Zeitschritten.
        !                                              Mein Driver kann auch Zeitschritte kleiner als einen Tag haben.
        !
        ! Meine Organisation ist freier: es wird im Runskipt festgelegt: 
        ! forcing_filename_base="forcing365_interface_vars_for_JS4-standalone_ml"                         => Climate_dCRUNCEP_R2B4 
        ! forcing_start_date="1979-01-01T00:00:00Z" # ISO 8601, in the format "1979-01-01T00:00:00Z"      => 1979-01-01T12:00:00Z
        !                                           # It is assumed that the first forcing file name starts with this date and uses a format specification like "19790101T000000Z"
        ! forcing_end_date="1980-01-01T00:00:00Z"   # The last forcing file name not necessarily uses this date, but the file has to include this date.    => 1980-01-01T00:00:00Z
        ! forcing_modelTimeStep="P1T" # e.g. PT20M only necessary to calclulate the date for the name of the forcing input files => einen Tag verwenden
        ! forcing_steps_per_file="72" # z.B. 72 um bei einem Modellzeitschritt von 20 Minuten 24 Stunden in einem File zusammen zu packen     => 365 um Tagesschritte für ein Jahr zu erhalten
        ! forcing_synchron_factor="1"     #  => erst mal 1 beibehalten und forcing_modelTimeStep dafür auf einen Tag setzen
        !

        ! Find out time_of_day and orbit_date of the current time step
        ! Note: In contrast to ECHAM5 orbit_date is not used anymore in ICON-A!
        !       Instead orbit_date_rt is used for both subsequent calls to find 
        !       out decl_sun and dist_sun. However, we need time_of_day because
        !       this is differentiated from time_of_day_rt when we call solar_parameters.

        ! Read in forcing file once an day
        ! We could read it more but not less frequent

        exit ! R: debug!!
        IF ( is_newday(mtime_old_jsb, dtime) ) THEN ! mtime_current
          mtime_forcing_ptr => mtime_forcing

          l_orbvsop87   => echam_rad_config(jg)% l_orbvsop87 ! Default: .TRUE.
          cecc          => echam_rad_config(jg)% cecc
          cobld         => echam_rad_config(jg)% cobld
          clonp         => echam_rad_config(jg)% clonp

          CALL get_orbit_times(l_orbvsop87, mtime_forcing_ptr, lyr_perp, yr_perp, & ! in
            &                  time_of_day, orbit_date)            ! out. 

          ! Find out time_of_day_rt and orbit_date_rt of the radiation transfer of the current time step
          CALL get_orbit_times(l_orbvsop87, mtime_forcing_ptr, lyr_perp, yr_perp, & ! in
            &                  time_of_day_rt, orbit_date_rt)      ! out
          
          ! Find out decl_sun and dist_sun
          IF (l_orbvsop87) THEN
            CALL orbit_vsop87(orbit_date_rt, rasc_sun, decl_sun, dist_sun) ! in, out, out, out
          ELSE
            CALL orbit_kepler(cecc, cobld, clonp, orbit_date_rt, & ! in 
              &               rasc_sun, decl_sun, dist_sun)        ! out
          END IF

          ! Get the solar angels
          !TODO: get from jsb grid? but how to name such a function get_lon_from_grid?
          !      and if, shouldn't then also get_nproma should be called get_nproma_from_grid?
          ! Note, this can be done here because we get lon lat for the whole horizontal grid!
          lon = get_lon(patch)
          lat = get_lat(patch)
          coslon = COS(pi * lon / 180._wp)
          sinlon = SIN(pi * lon / 180._wp)
          coslat = COS(pi * lat / 180._wp)
          sinlat = SIN(pi * lat / 180._wp)

          ! With the informations above we can determine the cos_zenith_angle of the current time step
          icosmu0        => echam_rad_config(jg)% icosmu0        ! defines if and how cosmu0 is defined for extended sunlit areas
          ldiur          => echam_rad_config(jg)% ldiur          ! .TRUE.
          l_sph_symm_irr => echam_rad_config(jg)% l_sph_symm_irr ! .FALSE.
          dt_ext         = 0.0_wp                                ! time interval over which the insolated area is extended, for us =0

          CALL solar_parameters(decl_sun,         time_of_day,     & ! in
            &                   icosmu0,          dt_ext,          & ! in
            &                   ldiur,            l_sph_symm_irr,  & ! in
            &                   patch,                             & ! in 
            &                   cos_zenith_angle, daylight_frc)      ! out

          ! We first need to set values for the variables USEd in mo_jsb_forcing that do not exist in ICON-A except for mo_jsb_forcing.
          ! We have to transfer these values from ICON-A/this driver.

          delta_time    = dtime
          time_step_len = dtime

          ! Now, as we have the cos_zenith_angle we will be abel to calculate swvis_srf_down and we can compile the forcing
          ! ToDo diese Variablen entweder gleich mit       forcing_in_modeltyp%t_air%data_read(jcs:nproma, jb, step_for_synchron)   usw. ersetzten oder
          ! nach dem call explizit übergeben!
          ! Dabei beachten: t_air_K entspricht t_air beim Model-forcing. Da im mo_forcing aber beide namen auftauchen kann man es dort nicht einfach ersetzen! Daher müßte es beim Einlesen der observations beibehalten werden.
          !                 Es müßte dann vor dem Aufruf von jsb_interface an t_air übergeben werden! Das möchte ich aber so nicht machen. 
          !                 => habe es in mo_forcing so geändert, daß es statt t_air ein t_air_C gibt und dort das neue t_air =t_air_K ist.
          ! Dabei beachten: CO2_air im jsb_interface wird von model-forcing bei mir so bereit gestellt, bei observations entspricht es CO2_concentration, was von mo_forcing als out geliefert wird! 
          !                 Dieses wird in mo_forcing aber zunächst aus foring files als mol eingelesen und wird erst dort in die vom jsb_interface und dem model-forcing verwendeten Einheit umgerechnet.
          !                 Man benötigt es in mo_forcing aber nicht mit diesem Namen weil es nicht mit diesen Variablennamen einlesen wird! 
          !                 Wahrscheinlich hat sie die Benennung von mo_forcing übernommen weil sie diese Variable beim model-forcing nicht einliest! Sie hat sich auch nicht mit Werten gefüllt weil sie damals wohl nicht
          !                 von JS4 verwendet wurde. => habe sie in mo_forcing umbenannt!
          !                 
          CALL read_atm_forcing_for_jsb_stdalone_from_observ(1, 1, nproma, nblks, mtime_forcing, mtime_current, sinlat, coslat, lon, evapo_act2pot_fromJS4generated, &
          & t_srf_fromJS4, fact_q_air_fromJS4, fact_qsat_srf_fromJS4, rough_h_srf_fromJS4, rough_m_srf_fromJS4, cos_zenith_angle,       &
          & CO2_air, t_air, q_air, rain, snow, wind_air, wind_10m,                                           &
          & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, fract_par_diffuse,                              &
          & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

          ! R: Put results in forcing data structure
          forcing_in_modeltyp%t_air%data_read             (jcs:nproma, jb, step_for_synchron) =   t_air            (jcs:nproma, jb) ! R: jb mit nblks ersetzen
          forcing_in_modeltyp%q_air%data_read             (jcs:nproma, jb, step_for_synchron) =   q_air            (jcs:nproma, jb) 
          forcing_in_modeltyp%rain%data_read              (jcs:nproma, jb, step_for_synchron) =   rain             (jcs:nproma, jb) 
          forcing_in_modeltyp%snow%data_read              (jcs:nproma, jb, step_for_synchron) =   snow             (jcs:nproma, jb) 
          forcing_in_modeltyp%wind_air%data_read          (jcs:nproma, jb, step_for_synchron) =   wind_air         (jcs:nproma, jb) 
          forcing_in_modeltyp%wind_10m%data_read          (jcs:nproma, jb, step_for_synchron) =   wind_10m         (jcs:nproma, jb) 
          forcing_in_modeltyp%lw_srf_down%data_read       (jcs:nproma, jb, step_for_synchron) =   lw_srf_down      (jcs:nproma, jb) 
          forcing_in_modeltyp%swvis_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swvis_srf_down   (jcs:nproma, jb) 
          forcing_in_modeltyp%swnir_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swnir_srf_down   (jcs:nproma, jb) 
          forcing_in_modeltyp%swpar_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swpar_srf_down   (jcs:nproma, jb) 
          forcing_in_modeltyp%fract_par_diffuse%data_read (jcs:nproma, jb, step_for_synchron) =   fract_par_diffuse(jcs:nproma, jb) 
          forcing_in_modeltyp%press_srf%data_read         (jcs:nproma, jb, step_for_synchron) =   press_srf        (jcs:nproma, jb) 
          forcing_in_modeltyp%drag_srf%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_srf         (jcs:nproma, jb) 
          forcing_in_modeltyp%t_acoef%data_read           (jcs:nproma, jb, step_for_synchron) =   t_acoef          (jcs:nproma, jb) 
          forcing_in_modeltyp%t_bcoef%data_read           (jcs:nproma, jb, step_for_synchron) =   t_bcoef          (jcs:nproma, jb) 
          forcing_in_modeltyp%q_acoef%data_read           (jcs:nproma, jb, step_for_synchron) =   q_acoef          (jcs:nproma, jb) 
          forcing_in_modeltyp%q_bcoef%data_read           (jcs:nproma, jb, step_for_synchron) =   q_bcoef          (jcs:nproma, jb) 
          forcing_in_modeltyp%pch%data_read               (jcs:nproma, jb, step_for_synchron) =   pch              (jcs:nproma, jb) 
          forcing_in_modeltyp%cos_zenith_angle%data_read  (jcs:nproma, jb, step_for_synchron) =   cos_zenith_angle (jcs:nproma, jb) 
          forcing_in_modeltyp%CO2_air%data_read           (jcs:nproma, jb, step_for_synchron) =   CO2_air          (jcs:nproma, jb) 
! Not present in the observational forcing
!             IF (model%config%use_lakes) THEN
!               forcing_in_modeltyp%drag_wtr%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_wtr         (jcs:nproma, jb) 
!               forcing_in_modeltyp%drag_ice%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_ice         (jcs:nproma, jb) 
!               forcing_in_modeltyp%t_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   t_acoef_wtr      (jcs:nproma, jb) 
!               forcing_in_modeltyp%t_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   t_bcoef_wtr      (jcs:nproma, jb) 
!               forcing_in_modeltyp%q_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   q_acoef_wtr      (jcs:nproma, jb) 
!               forcing_in_modeltyp%q_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   q_bcoef_wtr      (jcs:nproma, jb) 
!               forcing_in_modeltyp%t_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   t_acoef_ice      (jcs:nproma, jb) 
!               forcing_in_modeltyp%t_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   t_bcoef_ice      (jcs:nproma, jb) 
!               forcing_in_modeltyp%q_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   q_acoef_ice      (jcs:nproma, jb) 
!               forcing_in_modeltyp%q_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   q_bcoef_ice      (jcs:nproma, jb) 
!             END IF
        END IF ! is_newday
      ELSE
        counter = counter + 1
        ! IF (counter < 2) THEN ! R: NUR ZUM DEBUGGEN
        CALL read_atm_forcing_for_jsb_stdalone_from_modelout(model%config%use_lakes, forcing_synchron_factor, forcing_steps_per_file, mtime_forcing_char, nblks, nproma, model%grid_id) ! ICON-A forcing einlesen
        ! END IF 
      END IF ! l_forcing_observations
    ! ELSE
      !'forcing_steps_per_file', forcing_steps_per_file, 'forcing_synchron_factor', forcing_synchron_factor
      IF (.NOT. l_first_forcing_file) THEN
        l_shift_one_timestep = .FALSE.
      END IF
      l_first_forcing_file= .FALSE.
      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver: jstep:' , jstep, 'remain:', remain, 'counter:', counter, 'l_shift_one_timestep:', l_shift_one_timestep
    END IF  ! remain
    !END IF

    !Block loop 2

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs,jbe ! loop über den block index
      CALL get_indices_c(patch,  & ! in 
        &                   jb,  & ! in  index of block loop
        &                   jbs, & ! in  start indice of block loop
        &                   jbe, & ! in  end   indice of block loop
        &                   jcs, & ! out start indice of columns loop
        &                   jce, & ! out end   indice of columns loop
        &                   rls, & ! in  row limit start
        &                   rle)   ! in  row limit end
      ! WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR in loop   jbs', jbs, 'jbe', jbe
      IF (jcs>jce) CYCLE
      !write(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR jg:', jg, 'jb:', jb, 'jcs:', jcs, 'jce:', jce
      !counter=jb + 1
      !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR  DRIVER jb:' ,   jb, 'counter' , counter ! im std-gekoppelten-Lauf ist: jb max 39 , counter max 40, d.h. von 2 bis 40 =39
                                                                              !                    hier ist: jb start 1 , counter start 2, d.h. von 2 bis 2 =1
      ! CALL routine(jg,jb,jcs,jce        ,& ! R: das würde nun SUBROUTINE interface_echam_vdf callen
      !   &       nproma,nlev,ntracer  ,&
      !   &       is_in_sd_ed_interval ,&
      !   &       is_active            ,&
      !   &       datetime_old         ,&
      !   &       pdtime               )

      ! Anstatt die ganzen Prozeduren aufzurufen hier die direkte Übergabe zu den entsprechendne interface_full variablen:
      !iblk= jb
      !ics = jcs
      !ice = jce
      !nc =  jce - jcs + 1
      !CALL model%Set_options(iblk=jb, ics=jcs, ice=jce, nc=nc, dtime=dtime, steplen=dtime) ! wird in interface_full gemacht.

      !IF ( ics /= 1)      write(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRRRR ics:', ics ! should be 1
      !IF ( ice /= nproma) write(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRRRR ice:', ice ! should be e.g. 16

      !no_omp_thread = Get_omp_thread()
      !iblk = model%options(no_omp_thread)%iblk
      !ics = model%options(no_omp_thread)%ics
      !ice = model%options(no_omp_thread)%ice
      !dtime = model%options(no_omp_thread)%dtime

      !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver model%options(no_omp_thread)%iblk:', model%options(no_omp_thread)%iblk , 'model%options(no_omp_thread)%ics:'  , model%options(no_omp_thread)%ics , 'model%options(no_omp_thread)%ice' , model%options(no_omp_thread)%ice , 'model%options(no_omp_thread)%dtime' , model%options(no_omp_thread)%dtime

      synchr_ratio = step_in_forcing_file / forcing_synchron_factor
      step_for_synchron =nint(synchr_ratio)
      IF (l_shift_one_timestep) THEN
        step_for_synchron = step_for_synchron + 1
      END IF
      ! R: used only for debugging
      !  t_air            (jcs:nproma, jb) = forcing_in_modeltyp%t_air%data_read             (jcs:nproma, jb, step_for_synchron)
      !  q_air            (jcs:nproma, jb) = forcing_in_modeltyp%q_air%data_read             (jcs:nproma, jb, step_for_synchron)
      !  rain             (jcs:nproma, jb) = forcing_in_modeltyp%rain%data_read              (jcs:nproma, jb, step_for_synchron)
      !  snow             (jcs:nproma, jb) = forcing_in_modeltyp%snow%data_read              (jcs:nproma, jb, step_for_synchron)
      !  wind_air         (jcs:nproma, jb) = forcing_in_modeltyp%wind_air%data_read          (jcs:nproma, jb, step_for_synchron)
      !  wind_10m         (jcs:nproma, jb) = forcing_in_modeltyp%wind_10m%data_read          (jcs:nproma, jb, step_for_synchron)
      !  lw_srf_down      (jcs:nproma, jb) = forcing_in_modeltyp%lw_srf_down%data_read       (jcs:nproma, jb, step_for_synchron)
      !  swvis_srf_down   (jcs:nproma, jb) = forcing_in_modeltyp%swvis_srf_down%data_read    (jcs:nproma, jb, step_for_synchron)
      !  swnir_srf_down   (jcs:nproma, jb) = forcing_in_modeltyp%swnir_srf_down%data_read    (jcs:nproma, jb, step_for_synchron)
      !  swpar_srf_down   (jcs:nproma, jb) = forcing_in_modeltyp%swpar_srf_down%data_read    (jcs:nproma, jb, step_for_synchron)
      !  fract_par_diffuse(jcs:nproma, jb) = forcing_in_modeltyp%fract_par_diffuse%data_read (jcs:nproma, jb, step_for_synchron)
      !  press_srf        (jcs:nproma, jb) = forcing_in_modeltyp%press_srf%data_read         (jcs:nproma, jb, step_for_synchron)
      !  drag_srf         (jcs:nproma, jb) = forcing_in_modeltyp%drag_srf%data_read          (jcs:nproma, jb, step_for_synchron)
      !  t_acoef          (jcs:nproma, jb) = forcing_in_modeltyp%t_acoef%data_read           (jcs:nproma, jb, step_for_synchron)
      !  t_bcoef          (jcs:nproma, jb) = forcing_in_modeltyp%t_bcoef%data_read           (jcs:nproma, jb, step_for_synchron)
      !  q_acoef          (jcs:nproma, jb) = forcing_in_modeltyp%q_acoef%data_read           (jcs:nproma, jb, step_for_synchron)
      !  q_bcoef          (jcs:nproma, jb) = forcing_in_modeltyp%q_bcoef%data_read           (jcs:nproma, jb, step_for_synchron)
      !  pch              (jcs:nproma, jb) = forcing_in_modeltyp%pch%data_read               (jcs:nproma, jb, step_for_synchron)
      !  cos_zenith_angle (jcs:nproma, jb) = forcing_in_modeltyp%cos_zenith_angle%data_read  (jcs:nproma, jb, step_for_synchron)
      !  CO2_air          (jcs:nproma, jb) = forcing_in_modeltyp%CO2_air%data_read           (jcs:nproma, jb, step_for_synchron)
      !  trig_rad_m1      (jcs:nproma, jb) = forcing_in_modeltyp%trig_rad_m1%data_read       (jcs:nproma, jb, step_for_synchron)
      ! optional
      !  drag_wtr         (jcs:nproma, jb) = forcing_in_modeltyp%drag_wtr%data_read          (jcs:nproma, jb, step_for_synchron) 
      !  drag_ice         (jcs:nproma, jb) = forcing_in_modeltyp%drag_ice%data_read          (jcs:nproma, jb, step_for_synchron)
      !  t_acoef_wtr      (jcs:nproma, jb) = forcing_in_modeltyp%t_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron)
      !  t_bcoef_wtr      (jcs:nproma, jb) = forcing_in_modeltyp%t_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron)
      !  q_acoef_wtr      (jcs:nproma, jb) = forcing_in_modeltyp%q_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron)
      !  q_bcoef_wtr      (jcs:nproma, jb) = forcing_in_modeltyp%q_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron)
      !  t_acoef_ice      (jcs:nproma, jb) = forcing_in_modeltyp%t_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron)
      !  t_bcoef_ice      (jcs:nproma, jb) = forcing_in_modeltyp%t_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron)
      !  q_acoef_ice      (jcs:nproma, jb) = forcing_in_modeltyp%q_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron)
      !  q_bcoef_ice      (jcs:nproma, jb) = forcing_in_modeltyp%q_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron)

      IF (jb == 1) THEN
        WRITE(*,*) 'jcs, nproma, jb:' , jcs, nproma, jb, 'step_for_synchron:', step_for_synchron,  'jstep:', jstep
      END IF
      ! R: to check if data from forcing that arrives here is the same as in the forcing file
      ! DEBUG_VAR(jcs:nproma,jb) = forcing_in_modeltyp%snow%data_read(jcs:nproma, jb, step_for_synchron) ! - lw_srf_down(jcs:nproma,jb)

      !t_air     (jcs:nproma,jb) = forcing_in_modeltyp%t_air%data_read            (jcs:nproma, jb, step_for_synchron)
      !lw_srf_down(jcs:nproma,jb) = forcing_in_modeltyp%lw_srf_down%data_read      (jcs:nproma, jb, step_for_synchron)
      !lw_srf_down       = 103.265228271484        ! W/m2
      !forcing_in_modeltyp%lw_srf_down%data_read(jcs:nproma, jb, step_for_synchron) = lw_srf_down (jcs:nproma,jb)
      !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver local t_air MIN:',   MINVAL(t_air(jcs:nproma, jb)), 'MAX:' , MAXVAL(t_air(jcs:nproma, jb))
      !WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver forcing_in_modeltyp%t_air%data_read(jcs:nproma, jb, step_for_synchron) MIN:',  &
      !  & MINVAL(forcing_in_modeltyp%t_air%data_read(jcs:nproma, jb, step_for_synchron)), &
      !  & 'MAX:' , MAXVAL(forcing_in_modeltyp%t_air%data_read(jcs:nproma, jb, step_for_synchron))

      ! R: work around forcing data problem
  ! R: debugging: set values AFTER they are set from observational or model forcing to keep the model running to get DEBUG_VAR. Delete later
  !      !if (step_for_synchron < 5) then
  !         forcing_in_modeltyp%t_air%data_read             (jcs:nproma, jb, step_for_synchron) =   t_air            (jcs:nproma, jb)
  !         forcing_in_modeltyp%q_air%data_read             (jcs:nproma, jb, step_for_synchron) =   q_air            (jcs:nproma, jb) 
  !         forcing_in_modeltyp%rain%data_read              (jcs:nproma, jb, step_for_synchron) =   rain             (jcs:nproma, jb) 
  !         forcing_in_modeltyp%snow%data_read              (jcs:nproma, jb, step_for_synchron) =   snow             (jcs:nproma, jb) 
  !         forcing_in_modeltyp%wind_air%data_read          (jcs:nproma, jb, step_for_synchron) =   wind_air         (jcs:nproma, jb) 
  !         forcing_in_modeltyp%wind_10m%data_read          (jcs:nproma, jb, step_for_synchron) =   wind_10m         (jcs:nproma, jb) 
  !         forcing_in_modeltyp%lw_srf_down%data_read       (jcs:nproma, jb, step_for_synchron) =   lw_srf_down      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%swvis_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swvis_srf_down   (jcs:nproma, jb) 
  !         forcing_in_modeltyp%swnir_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swnir_srf_down   (jcs:nproma, jb) 
  !         forcing_in_modeltyp%swpar_srf_down%data_read    (jcs:nproma, jb, step_for_synchron) =   swpar_srf_down   (jcs:nproma, jb) 
  !         forcing_in_modeltyp%fract_par_diffuse%data_read (jcs:nproma, jb, step_for_synchron) =   fract_par_diffuse(jcs:nproma, jb) 
  !         forcing_in_modeltyp%press_srf%data_read         (jcs:nproma, jb, step_for_synchron) =   press_srf        (jcs:nproma, jb) 
  !         forcing_in_modeltyp%drag_srf%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_srf         (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_acoef%data_read           (jcs:nproma, jb, step_for_synchron) =   t_acoef          (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_bcoef%data_read           (jcs:nproma, jb, step_for_synchron) =   t_bcoef          (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_acoef%data_read           (jcs:nproma, jb, step_for_synchron) =   q_acoef          (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_bcoef%data_read           (jcs:nproma, jb, step_for_synchron) =   q_bcoef          (jcs:nproma, jb) 
  !         forcing_in_modeltyp%pch%data_read               (jcs:nproma, jb, step_for_synchron) =   pch              (jcs:nproma, jb) 
  !         forcing_in_modeltyp%cos_zenith_angle%data_read  (jcs:nproma, jb, step_for_synchron) =   cos_zenith_angle (jcs:nproma, jb) 
  !         forcing_in_modeltyp%CO2_air%data_read           (jcs:nproma, jb, step_for_synchron) =   CO2_air          (jcs:nproma, jb) 
  !         forcing_in_modeltyp%trig_rad_m1%data_read       (jcs:nproma, jb, step_for_synchron) =   trig_rad_m1          (jcs:nproma, jb) 
  ! !         !
  !      !END IF
  !         forcing_in_modeltyp%drag_wtr%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_wtr         (jcs:nproma, jb) 
  !         forcing_in_modeltyp%drag_ice%data_read          (jcs:nproma, jb, step_for_synchron) =   drag_ice         (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   t_acoef_wtr      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   t_bcoef_wtr      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_acoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   q_acoef_wtr      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_bcoef_wtr%data_read       (jcs:nproma, jb, step_for_synchron) =   q_bcoef_wtr      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   t_acoef_ice      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%t_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   t_bcoef_ice      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_acoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   q_acoef_ice      (jcs:nproma, jb) 
  !         forcing_in_modeltyp%q_bcoef_ice%data_read       (jcs:nproma, jb, step_for_synchron) =   q_bcoef_ice      (jcs:nproma, jb) 


      !WRITE(*,*) 'RRRRRRRRRRRR: Rein geht forcing_in_modeltyp%newday_Merker%data_read (jcs, jb, step_for_synchron)' , forcing_in_modeltyp%newday_Merker%data_read (jcs, jb, step_for_synchron)
      !  WRITE(*,*) 'RRRRRRRRRRRR: Interface wird aufgerufen'

      CALL jsbach_interface (            &
        ! kontext:
        & 1,                             & ! model_id
        & jb, jcs, nproma, &
        & mtime_old_jsb, & ! mtime_old
        & dtime, dtime, & ! R: Note, the first dtime corresponds to delta_time and the second dtime corresponds to time_step_len (time_step_len = 2.0_dp*delta_time) when js4 runs with echam.
                                            !    Note, normally nproma should be ice in this call! Because ice can be smaller than nproma on a block/jb/iblks. However, the nml Output can not be written with ice.
        ! in vars:
        & forcing_in_modeltyp%t_air%data_read            (jcs:nproma, jb, step_for_synchron), &  ! kg kg-1
        & forcing_in_modeltyp%q_air%data_read            (jcs:nproma, jb, step_for_synchron), &  ! kg kg-1
        & forcing_in_modeltyp%rain%data_read             (jcs:nproma, jb, step_for_synchron), &  ! kg/m2/s
        & forcing_in_modeltyp%snow%data_read             (jcs:nproma, jb, step_for_synchron), &  ! kg/m2/s
        & forcing_in_modeltyp%wind_air%data_read         (jcs:nproma, jb, step_for_synchron), &  ! air_wind_speed, m/s
        & forcing_in_modeltyp%wind_10m%data_read         (jcs:nproma, jb, step_for_synchron), &  ! air_wind_speed, m/s
        & forcing_in_modeltyp%lw_srf_down%data_read      (jcs:nproma, jb, step_for_synchron), &  ! W/m2   
        & forcing_in_modeltyp%swvis_srf_down%data_read   (jcs:nproma, jb, step_for_synchron), &  ! W/m2
        & forcing_in_modeltyp%swnir_srf_down%data_read   (jcs:nproma, jb, step_for_synchron), &  ! W/m2
        & forcing_in_modeltyp%swpar_srf_down%data_read   (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%fract_par_diffuse%data_read(jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%press_srf%data_read        (jcs:nproma, jb, step_for_synchron), &  ! N/m2 =Pa
        & forcing_in_modeltyp%drag_srf%data_read         (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%t_acoef%data_read          (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%t_bcoef%data_read          (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%q_acoef%data_read          (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%q_bcoef%data_read          (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%pch%data_read              (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%cos_zenith_angle%data_read (jcs:nproma, jb, step_for_synchron), &  ! 
        & forcing_in_modeltyp%CO2_air%data_read          (jcs:nproma, jb, step_for_synchron), &  ! mol: 0.00034800
        ! out vars, only dummy use for standalone:
        & t_srf_fromJS4        (jcs:nproma, jb),   & !
        & t_eff_srf    (jcs:nproma, jb),             & !
        & qsat_srf     (jcs:nproma, jb),             & !
        & s_srf        (jcs:nproma, jb),             & !
        & fact_q_air_fromJS4   (jcs:nproma, jb),   & !
        & fact_qsat_srf_fromJS4(jcs:nproma, jb),   & !
        & evapotrans   (jcs:nproma, jb),             & !
        & latent_hflx  (jcs:nproma, jb),             & !
        & sensible_hflx(jcs:nproma, jb),             & !
        & grnd_hflx    (jcs:nproma, jb),             & !
        & grnd_hcap    (jcs:nproma, jb),             & !
        & rough_h_srf_fromJS4  (jcs:nproma, jb),   & !
        & rough_m_srf_fromJS4  (jcs:nproma, jb),   & !
        & q_snocpymlt  (jcs:nproma, jb),             & !
        & alb_vis_dir  (jcs:nproma, jb),             & !
        & alb_nir_dir  (jcs:nproma, jb),             & !
        & alb_vis_dif  (jcs:nproma, jb),             & !
        & alb_nir_dif  (jcs:nproma, jb),             & !
        & CO2_flux     (jcs:nproma, jb),             & !
        ! in vars for lakes, optional:
        & DEBUG_VAR   = DEBUG_VAR (jcs:nproma, jb),  &  
        & drag_wtr          =  forcing_in_modeltyp%drag_wtr%data_read    (jcs:nproma, jb, step_for_synchron),   &
        & drag_ice          =  forcing_in_modeltyp%drag_ice%data_read    (jcs:nproma, jb, step_for_synchron),   &
        & t_acoef_wtr       =  forcing_in_modeltyp%t_acoef_wtr%data_read (jcs:nproma, jb, step_for_synchron),   &
        & t_bcoef_wtr       =  forcing_in_modeltyp%t_bcoef_wtr%data_read (jcs:nproma, jb, step_for_synchron),   &
        & q_acoef_wtr       =  forcing_in_modeltyp%q_acoef_wtr%data_read (jcs:nproma, jb, step_for_synchron),   &
        & q_bcoef_wtr       =  forcing_in_modeltyp%q_bcoef_wtr%data_read (jcs:nproma, jb, step_for_synchron),   &
        & t_acoef_ice       =  forcing_in_modeltyp%t_acoef_ice%data_read (jcs:nproma, jb, step_for_synchron),   &
        & t_bcoef_ice       =  forcing_in_modeltyp%t_bcoef_ice%data_read (jcs:nproma, jb, step_for_synchron),   &
        & q_acoef_ice       =  forcing_in_modeltyp%q_acoef_ice%data_read (jcs:nproma, jb, step_for_synchron),   &
        & q_bcoef_ice       =  forcing_in_modeltyp%q_bcoef_ice%data_read (jcs:nproma, jb, step_for_synchron),   &
        ! out vars for lakes, only dummy use for standalone, optional:
        & t_lwtr            =  t_lwtr           (jcs:nproma, jb),  &  !
        & qsat_lwtr         =  qsat_lwtr        (jcs:nproma, jb),  &  !
        & s_lwtr            =  s_lwtr           (jcs:nproma, jb),  &  !
        & evapo_wtr         =  evapo_wtr        (jcs:nproma, jb),  &  !
        & latent_hflx_wtr   =  latent_hflx_wtr  (jcs:nproma, jb),  &  !
        & sensible_hflx_wtr =  sensible_hflx_wtr(jcs:nproma, jb),  &  !
        & albedo_lwtr       =  albedo_lwtr      (jcs:nproma, jb),  &  !
        & t_lice            =  t_lice           (jcs:nproma, jb),  &  !
        & qsat_lice         =  qsat_lice        (jcs:nproma, jb),  &  !
        & s_lice            =  s_lice           (jcs:nproma, jb),  &  !
        & evapo_ice         =  evapo_ice        (jcs:nproma, jb),  &  !
        & latent_hflx_ice   =  latent_hflx_ice  (jcs:nproma, jb),  &  !
        & sensible_hflx_ice =  sensible_hflx_ice(jcs:nproma, jb),  &  !
        & albedo_lice       =  albedo_lice      (jcs:nproma, jb),  &  !
        & ice_fract_lake    =  ice_fract_lake   (jcs:nproma, jb),  &  !
        ! out vars for standalone ECHAM-jsbach, optional:
        & evapopot          = evapopot(jcs:nproma, jb),            &  !
        ! in var for standalone ICON-jsbach, optional:
        & trig_rad_m1       = forcing_in_modeltyp%trig_rad_m1%data_read (jcs:nproma, jb, step_for_synchron),    &  !
        & newday_Merker     = forcing_in_modeltyp%newday_Merker%data_read (jcs:nproma, jb, step_for_synchron),  &  !
        & l_day_reminder    = forcing_in_modeltyp%l_day_reminder%data_read (jcs:nproma, jb, step_for_synchron)  &  !
        & )
    END DO

    ! Only for observations as input to calculate vapour_pressure_from_evapor
    IF ( l_forcing_observations ) THEN
      ! TODO: if really required (only ERA-interim forcing data):
      ! To ensure restart identity, the daily evapo_act_sum and evapo_pot_sum would need to be calculated in jsbach!
      ! => for now the instantaneous values are used instead of the daily sums, to ensure restart identity
      evapo_act2pot_fromJS4generated(:, jb) = MIN(evapotrans(:, jb),-EPSILON(1._wp)) / MIN(evapopot(:, jb),  -EPSILON(1._wp))
      ! evapo_act_sum_fromJS4(:,jb) = evapo_act_sum(:,jb) + MIN(evapotrans(:),-EPSILON(1._wp))
      ! evapo_pot_sum_fromJS4(:,jb) = evapo_pot_sum(:,jb) + MIN(evapopot(:),  -EPSILON(1._wp))
    ENDIF

!$OMP END PARALLEL DO

! R:echam_phy_main ENDE

    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Y'
    CALL jsbach_finish_timestep(jg, mtime_old_jsb, dtime) ! all in  mtime_current
    !!!!!! ENDE: SUBROUTINE interface_iconam_echam

    !!!!!!!!!!!!!!!!!!!!! integrate_nh ENDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! END DO ! für jsubstep = 1, num_steps, hab ich hier her gesetzt muß evtl. woanderst hin...  -------------------------------------------------------------------------------

    ! Compute diagnostics for output if necessary
    !delete    !IF (l_compute_diagnostic_quants .OR. iforcing==iecham .OR. iforcing==inoforcing) THEN
    CALL diag_for_output_dyn ()
    !delete    !ENDIF

    !--------------------------------------------------------------------------
    ! loop over the list of internal post-processing tasks, e.g.
    ! interpolate selected fields to p- and/or z-levels
    !
    ! Mean sea level pressure needs to be computed also at
    ! no-output-steps for accumulation purposes; set by l_accumulation_step
    !simulation_status = new_simulation_status(l_output_step  = l_nml_output,             &
    !  &                                       l_last_step    = (jstep==(nsteps+jstep0)), &
    !  &                                       l_accumulation_step = (iforcing == iecham),&
    !  &                                       l_dom_active   = p_patch(1:)%ldom_active,  &
    !  &                                       i_timelevel_dyn= nnow, i_timelevel_phy= nnow_rcf)
    !CALL pp_scheduler_process(simulation_status)

    ! update accumlated values
      CALL update_statistics ! 
    ! R: only for debugging
    ! IF (p_nh_opt_diag(1)%acc%l_any_m) THEN
    !   CALL update_opt_acc(p_nh_opt_diag(1)%acc,            &
    !     &                 p_nh_state(1)%prog(nnow_rcf(1)), &
    !     &                 p_nh_state(1)%prog(nnow(1))%rho, &
    !     &                 p_nh_state(1)%diag,              &
    !     &                 p_patch(1)%cells%owned,          &
    !     &                 p_patch(1)%nlev)
    ! IF (l_nml_output) CALL calc_mean_opt_acc(p_nh_opt_diag(1)%acc)
    ! END IF

    ! output of results
    ! note: nnew has been replaced by nnow here because the update

    IF (l_nml_output) THEN
      CALL write_name_list_output(jstep)
    ENDIF

    CALL reset_statistics  ! R: nur zum deallocate debuggen
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z'

    ! Diagnostics: computation of total integrals
    !
    ! Diagnostics computation is not yet properly MPI-parallelized
    !
    !IF (output_mode%l_totint .AND. is_totint_time(current_step =jstep,   &
    !  &                                           restart_step = jstep0, &
    !  &                                           n_diag       = n_diag, &
    !  &                                           n_steps      = nsteps) ) THEN

    !  kstep = jstep-jstep0

! #ifdef NOMPI
!       IF (my_process_is_mpi_all_seq()) &
! #endif
!         CALL supervise_total_integrals_nh( kstep, p_patch(1:), p_nh_state, p_int_state(1:), &
!         &                                  nnow(1:n_dom), nnow_rcf(1:n_dom), jstep == (nsteps+jstep0))
!     ENDIF


    ! re-initialize MAX/MIN fields with 'resetval'
    ! must be done AFTER output

    CALL reset_act%execute(slack=dtime, mtime_date=mtime_current) !R: zeigt auf mo_action.f90: SUBROUTINE reset_kernel

    ! re-initialization for FG-averaging. Ensures that average is centered in time.
    IF (is_avgFG_time(mtime_old_jsb)) THEN ! mtime_current
      IF (p_nh_state(1)%diag%nsteps_avg(1) == 0) THEN
        CALL reinit_average_first_guess(p_patch(1), p_nh_state(1)%diag, p_nh_state(1)%prog(nnow_rcf(1)))
      END IF
    ENDIF

    !--------------------------------------------------------------------------
    ! Write restart file
    !--------------------------------------------------------------------------
    ! check whether time has come for writing restart file
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z1'
    CALL message('','')
    !
    ! default is to assume we do not write a checkpoint/restart file
    lwrite_checkpoint = .FALSE.
    ! if thwe model is not supposed to write output, do not write checkpoints
    IF (.NOT. output_mode%l_none ) THEN
      ! to clarify the decision tree we use shorter and more expressive names:

      l_isStartdate    = (time_config%tc_startdate == mtime_current)
      l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
      ! jsbach_is_restarted      = isCurrentEventActive(restartEvent, mtime_current) ! now from mo_jsb_control
      l_isCheckpoint   = isCurrentEventActive(checkpointEvent, mtime_current)
      l_doWriteRestart = time_config%tc_write_restart

      IF ( &
            !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
            &         (jsbach_is_restarted() .OR. l_isCheckpoint)                     &
            &  .AND.                                                        &
            !  and the current date differs from the start date
            &        .NOT. l_isStartdate                                    &
            &  .AND.                                                        &
            !  and end of run has not been reached or restart writing has been disabled
            &        (.NOT. l_isExpStopdate .OR. l_doWriteRestart)          &
            & ) THEN
        lwrite_checkpoint = .TRUE.
      END IF
    END IF
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z2'

    IF (lwrite_checkpoint) THEN
      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z3'
      ! apply nest boundary filling; this has no impact on the correctness (in the sense of reproducibility)
      ! of the restart but facilitates debugging
      CALL diag_for_output_dyn () ! R: Textoutput des Outputs als Tabelle

      DO jg = 1, n_dom

        ! get elapsed time since last call of NWP physics processes and write it into 
        ! the restart file

        CALL restartDescriptor%updatePatch(p_patch(jg), &
          & opt_t_elapsed_phy          = elapsedTime,                &
          & opt_ndyn_substeps          = ndyn_substeps_var(jg),      & ! R: ist 5
          & opt_jstep_adv_marchuk_order= jstep_adv(jg)%marchuk_order,&
          & opt_depth_lnd              = nlev_soil,                  &
          & opt_nlev_snow              = nlev_snow,                  &
          & opt_ndom                   = n_dom)

      ENDDO
      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z4'
      ! finally write restart file
      CALL restartDescriptor%writeRestart(mtime_current, jstep, opt_output_jfile = output_jfile)

      WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z5'
      IF (ALLOCATED(elapsedTime)) THEN
        DEALLOCATE(elapsedTime, STAT=error_status)
        IF (error_status /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed!')
      ENDIF
    END IF  ! lwrite_checkpoint

    IF (mtime_current >= time_config%tc_stopdate) THEN
      ! leave time loop
      EXIT TIME_LOOP
    END IF

    jstep = jstep + 1
    step_in_forcing_file = step_in_forcing_file + 1
    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z6'
    sim_time = getElapsedSimTimeInSeconds(mtime_current) 

  END DO TIME_LOOP

  ! -------------------------------------------------------------------------------------------------------------------------------
  ! Clean up
  ! -------------------------------------------------------------------------------------------------------------------------------
  ! deallocate in- and output arrays
  !   DEALLOCATE(t_air, q_air, rain, snow, wind_air, wind_10m)
  !   DEALLOCATE(lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, fract_par_diffuse)
  !   DEALLOCATE(press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle, CO2_air, trig_rad_m1, newday_Merker, l_day_reminder)
  ! 
  !   DEALLOCATE(t_eff_srf, qsat_srf, s_srf)
  !   DEALLOCATE(evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap)
  !   DEALLOCATE(q_snocpymlt, alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
  !   DEALLOCATE(CO2_flux)

  ! clean-up routine for mo_nh_supervise module (eg. closing of files)
  !CALL finalize_supervise_nh()
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver Z7'
  CALL deleteRestartDescriptor(restartDescriptor)

  !IF (ltimer) CALL timer_stop(timer_total)

  ! clean up
  IF (output_mode%l_nml) THEN
    DEALLOCATE(output_jfile, STAT=error_status)
    IF (error_status /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed!')
  ENDIF

  CALL deallocateDatetime(mtime_forcing)
  CALL deallocateDatetime(mtime_old)
  CALL deallocateDatetime(mtime_old_jsb)
  DO jg=1,n_dom
    IF (ASSOCIATED(datetime_current(jg)%ptr)) &
      &  CALL deallocateDatetime(datetime_current(jg)%ptr)
  END DO

  ! --------------------
  ! WIEDER ZURÜCK: AUS    CALL perform_nh_timeloop
  !                              also in perform_nh_stepping

  ! Aus CALL deallocate_nh_stepping (latbc):

  !-----------------------------------------------------------------------
  !
  ! deallocate auxiliary fields for tracer transport and rcf
  !
  !   DO jg = 1, n_dom
  !     DEALLOCATE( prep_adv(jg)%mass_flx_me, prep_adv(jg)%mass_flx_ic,     &
  !       &         prep_adv(jg)%vn_traj, prep_adv(jg)%topflx_tra, STAT=error_status )
  !     IF (error_status /= SUCCESS) THEN
  !       CALL finish ( modname//': perform_nh_stepping',            &
  !         &    'deallocation for mass_flx_me, mass_flx_ic, vn_traj,' // &
  !         &    'topflx_tra failed' )
  !     ENDIF
  !   ENDDO

  !   DEALLOCATE( prep_adv, STAT=error_status )
  !   IF (error_status /= SUCCESS) THEN
  !     CALL finish ( modname//': perform_nh_stepping',              &
  !       &    'deallocation for prep_adv failed' )
  !   ENDIF

  DEALLOCATE( jstep_adv, STAT=error_status )
  IF (error_status /= SUCCESS) THEN
    CALL finish ( 'RR: perform_nh_stepping',              &
      &    'deallocation for jstep_adv failed' )
  ENDIF
  
  !
  ! deallocate flow control variables
  !
  !   DEALLOCATE( linit_dyn, STAT=error_status )
  !   IF (error_status /= SUCCESS) THEN
  !     CALL finish ( 'RR: perform_nh_stepping',    &
  !       &    'deallocation for linit_dyn failed' )
  !   ENDIF

  ! --------------------
  ! WIEDER ZURÜCK: AUS    CALL atmo_nonhydrostatic

  !---------------------------------------------------------------------
  ! 6. Integration finished. Clean up.
  !---------------------------------------------------------------------
  !Aus CALL destruct_atmo_nonhydrostatic(latbc) -> MODULE mo_atmo_nonhydrostatic: SUBROUTINE destruct_atmo_nonhydrostatic


  ! Delete optional diagnostics
  !CALL destruct_opt_diag()

  ! Delete state variables

  !     CALL destruct_nh_state( p_nh_state, p_nh_state_lists )  ! R: wird beides nicht benötigt nur zum Debuggen von deallocate fehler
  !      DEALLOCATE (p_nh_state, STAT=error_status)
  !      IF (error_status /= SUCCESS) THEN
  !        CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
  !      ENDIF
  !      DEALLOCATE (p_nh_state_lists, STAT=error_status)
  !      IF (error_status /= SUCCESS) THEN
  !        CALL finish(TRIM(routine),'deallocation for p_nh_state_lists failed')
  !      ENDIF

  CALL cleanup_echam_phy
  ! CALL SLEEP(240) ! R: 4 Minuten warten um output klarer zeitlich zu trennen!
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 0000000000000000000000000000000000000000000000000000000000000000000000'
  ! Delete output variable lists
  IF (output_mode%l_nml) THEN
    CALL close_name_list_output  ! R: In ICON-A this would be done independent of the io type synch/asynch. 
                                  !    However, in the asynch case this is done also in name_list_io_main_proc (in the background).
                                  !    Then it is not necessary and produces a wondrous "allocatable array or pointer is not allocated" 
                                  !    error when deallocating the last (!) nml output_file entry in the command DEALLOCATE(output_file)
                                  !    in SUBROUTINE close_name_list_output. 
    CALL finish_statistics_streams
  END IF

!     ! finalize meteogram output
!     IF (output_mode%l_nml) THEN
!       DO jg = 1, n_dom
!         IF (meteogram_output_config(jg)%lenabled) THEN
!           CALL meteogram_finalize(jg)
!         END IF
!       END DO
!       DO jg = 1, max_dom
!         DEALLOCATE(meteogram_output_config(jg)%station_list)
!       END DO
!     END IF

   !Close LES diag files
   
!     IF (ldass_lhn) THEN 
!       ! deallocate ext_data array
!       DEALLOCATE(radar_data, STAT=ist)
!       IF (ist /= SUCCESS) THEN
!         CALL finish(TRIM(routine), 'deallocation of radar_data')
!       ENDIF
!       CALL destruct_lhn (lhn_fields)
!     ENDIF

  CALL message(TRIM(routine),'clean-up destruct_atmo_nonhydrostatic finished')

 ! --------------------
 ! WIEDER ZURÜCK: AUS    CALL atmo_model
  ! print performance timers:
 !    IF (ltimer) CALL print_timer

  ! 13. Integration finished. Carry out the shared clean-up processes
  !---------------------------------------------------------------------
  !Aus mo_atmo_model: CALL destruct_atmo_model ():
  CALL destruct_ext_data

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 11111111111111111111111111'
  CALL lonlat_grids%finalize() !USE mo_intp_lonlat_types,       ONLY: lonlat_grids    CALL lonlat_grids%init()

  !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 22222222222222222222222222'
  CALL destruct_comm_patterns( p_patch, p_patch_local_parent ) ! USE mo_model_domain,        ONLY: t_patch, t_grid_cells, p_patch_local_parent, p_patch
  !     ! Deallocate grid patches
  !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 33333333333333333333333333'
  CALL destruct_patches( p_patch )

  !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 44444444444444444444444444'
  CALL destruct_patches( p_patch_local_parent )

  !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 55555555555555555555555555'
  ! !    CALL message('','destruct_patches is done')
  CALL destruct_icon_communication()

  !  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 66666666666666666666666666'
  CALL jsbach_finalize()

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 77777777777777777777777777'

  ! --------------------
  ! WIEDER ZURÜCK: AUS   PROGRAM icon

  IF (ASSOCIATED(time_config%tc_exp_stopdate) .AND. ASSOCIATED(time_config%tc_stopdate)) THEN
    ! write the control.status file

    WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 88888888888888888888888888'

    IF (my_process_is_global_root()) THEN
      OPEN (500, FILE="finish.status")
      IF ((time_config%tc_exp_stopdate > time_config%tc_stopdate) .AND. time_config%tc_write_restart) THEN
        WRITE(500,*) "RESTART"
      ELSE
        WRITE(500,*) "OK"
      ENDIF
      CLOSE(500)
    END IF
  END IF

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 999999999999999999999999999999'

  ! Shut down MPI
  CALL stop_mpi
  !STOP
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 1010110101010101010101010101010'

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver 11 11 11 11 11 11' ! R: hier kommt er nicht rein
  CALL ieee_set_status(saved_fpscr)
#endif
#endif

  WRITE(*,*) 'RRRRRRRRRRRRRRRRRRRRRRRRRR   JS4-Driver ffffffffffffffffffffffffffffffff'

#endif
END PROGRAM jsb4_driver
