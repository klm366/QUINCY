!> ICON driver for the JSBACH standalone model
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
!>#### Main routine for running ICON-Land standalone
!>
!> `jsbach_model()` is called from the ICON program `src/drivers/icon:icon`
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsbach_model

#ifndef __NO_JSBACH__

  USE mo_exception,   ONLY: message, message_text, finish
  USE mo_kind,        ONLY: lp, wp
  USE mo_util_string, ONLY: int2string, real2string

  USE mtime,        ONLY: datetime, newDatetime, datetimeToString, timedelta, newTimedelta, getPTStringFromSeconds, &
    &                     mtime_strerror, timedeltaToString, deallocateDatetime, deallocateTimedelta, &
    &                     OPERATOR(+), OPERATOR(==), OPERATOR(<), OPERATOR(>), OPERATOR(>=), &
    &                     event, eventGroup, newEvent, addEventToEventGroup, &
    &                     NO_ERROR, MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, MAX_MTIME_ERROR_STR_LEN
  USE mo_timer,     ONLY: init_timer, timer_start, timer_stop, timers_level,    &
    &                     timer_model_init, print_timer,                        &
    &                     timer_domain_decomp, timer_compute_coeffs,            &
    &                     timer_read_restart, timer_total
  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: get_model

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: jsbach_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsbach_model'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> main routine for running ICON-Land standalone using ICON infrastructure (without atmosphere)
  !>
  !>  The standalone model is selected with `model_type=5` (`jsbach_process`) in `master_model_nml` of the
  !>  ICON master namelist. It runs either on a global (or regional) grid or for one or more single grid cells
  !>  (site-level).
  !>
  !>  'jsbach_model()' is called from the ICON program `src/drivers/icon:icon`
  !>
  !>  Structure:
  !>
  !>  1. configuration of the ICON infrastructure used by ICON-Land
  !>  2. configuration and initialization of ICON-Land
  !>  3. preparations for the time loop
  !>  4. time loop
  !>  5. cleanup and finalization
  !>
  SUBROUTINE jsbach_model(jsb_namelist_filename, shr_namelist_filename)

    USE mo_master_control,          ONLY: jsbach_process, get_my_process_name
    USE mo_master_config,           ONLY: isRestart

    USE mo_load_restart,            ONLY: read_restart_header, read_restart_files
    USE mo_restart_nml_and_att,     ONLY: getAttributesForRestarting
    USE mo_key_value_store,         ONLY: t_key_value_store

    USE mo_read_namelists_iface,    ONLY: read_infrastructure_namelists_for_icon_land, print_infrastructure_namelists_for_icon_land
    USE mo_namelist,                ONLY: open_nml_output, close_nml_output
    USE mo_nml_annotate,            ONLY: log_nml_settings
    USE mo_dynamics_config,         ONLY: nold, nnow, nnew
    USE mo_name_list_output_init,   ONLY: read_name_list_output_namelists, init_name_list_output, parse_variable_groups, &
      &                                   output_file, create_vertical_axes
    USE mo_level_selection,         ONLY: create_mipz_level_selections

    USE mo_time_management,          ONLY: compute_timestep_settings,                        &
      &                                    compute_restart_settings,                         &
      &                                    compute_date_settings
    USE mo_time_config,              ONLY: time_config, dt_restart
    USE mo_run_config,               ONLY: configure_run, dtime, nsteps, ltimer, output_mode, &
      &                                    nshift, num_lev, msg_level, grid_generatingCenter, grid_generatingSubcenter
    USE mo_event_manager,            ONLY: initEventManager, addEventGroup, getEventGroup, printEventGroup

    USE mo_io_config,               ONLY: restartWritingParameters, configure_io

    USE mo_memory_log,              ONLY: memory_log_add

    USE mo_mpi,                     ONLY: stop_mpi, my_process_is_io, my_process_is_stdio, &
      &                                   set_mpi_work_communicators, process_mpi_io_size
    USE mo_parallel_config,         ONLY: p_test_run, num_test_pe, l_test_openmp,                              &
      &                                   num_io_procs, proc0_shift,                                           &
      &                                   num_prefetch_proc, nproma, proc0_offloading, check_parallel_configuration

    USE mo_zaxis_type,              ONLY: zaxisTypeList, t_zaxisTypeList
    USE mo_jsb_base,                ONLY: jsbach_setup_models, jsbach_setup_tiles
    USE mo_jsb_time,                ONLY: timesteps_per_day
    USE mo_jsb_model_init,          ONLY: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart
    USE mo_jsb_model_final,         ONLY: jsbach_finalize
    USE mo_jsb_control,             ONLY: get_no_of_models, force_from_observations
    USE mo_jsb4_forcing,            ONLY: setup_forcing, init_forcing, finalize_external_forcing
    USE mo_jsb_parallel,            ONLY: Get_omp_no_of_threads

    USE mo_build_decomposition,     ONLY: build_decomposition

    ! I/O
    USE mo_restart,                 ONLY: t_RestartDescriptor, createRestartDescriptor, deleteRestartDescriptor, detachRestartProcs
    USE mo_name_list_output,        ONLY: name_list_io_main_proc, write_name_list_output, close_name_list_output, &
      &                                   istime4name_list_output
    USE mo_output_event_handler,    ONLY: get_current_jfile
    USE mo_icon_output_tools,       ONLY: init_io_processes
    USE mo_name_list_output_config, ONLY: use_async_name_list_io
    USE mo_output_event_types,      ONLY: t_sim_step_info
    USE mo_impl_constants,          ONLY: SUCCESS
    USE mo_memory_log,              ONLY: memory_log_terminate
    USE mo_util_mtime,              ONLY: getElapsedSimTimeInSeconds, is_event_active

    USE mo_grid_config,             ONLY: n_dom, n_dom_start
    USE mo_model_domain,            ONLY: p_patch
    USE mo_model_domain,            ONLY: p_patch_local_parent
    USE mo_icon_comm_interface,     ONLY: construct_icon_communication, destruct_icon_communication
    USE mo_complete_subdivision,    ONLY: setup_phys_patches
    USE mo_gribout_config,          ONLY: configure_gribout

    USE mo_vertical_coord_table,    ONLY: vct_a, vct_b, vct, allocate_vct_atmo
    ! Horizontal interpolation
    USE mo_intp_lonlat_types,       ONLY: lonlat_grids
    USE mo_intp_data_strc,          ONLY: p_int_state, p_int_state_local_parent
    USE mo_interpol_config,         ONLY: configure_interpolation
    USE mo_intp_state,              ONLY: construct_2d_interpol_state,                          &
      &                                   destruct_2d_interpol_state !, transfer_interpol_state
    USE mo_intp_lonlat,             ONLY: compute_lonlat_intp_coeffs
    USE mo_gridref_config,          ONLY: configure_gridref

    USE mo_alloc_patches,           ONLY: destruct_patches, destruct_comm_patterns
    USE mo_action,                  ONLY: ACTION_RESET, reset_act
    USE mo_var_list_register_utils, ONLY: vlr_print_groups, vlr_print_vls

    USE mo_derived_variable_handling, ONLY: update_statistics
#ifndef __NO_RAGNAROK__
    USE mo_ragnarok_setup,          ONLY: init_ragnarok_model, finalize_ragnarok_model
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(LEN=*), INTENT(in) :: jsb_namelist_filename  !< ICON-Land master namelist file
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename  !< shared namelist file
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER   :: model

    TYPE(t_key_value_store),      POINTER :: restartAttributes
    CLASS(t_RestartDescriptor),   POINTER :: restartDescriptor
    INTEGER :: jstep0, jstep, dedicatedRestartProcs
    INTEGER :: no_models, jg, no_omp_threads
    INTEGER :: i
    CHARACTER(len=MAX_MTIME_ERROR_STR_LEN) :: errstring
    INTEGER                                :: error_status

    ! Time and event management
    TYPE(datetime), POINTER       :: mtime_current
    TYPE(datetime), POINTER       :: mtime_old
    REAL(lp) :: sim_time          ! elapsed simulation time
    TYPE(datetime),  POINTER :: eventStartDate    => NULL(), &
      &                         eventEndDate      => NULL()
    TYPE(datetime),  POINTER :: checkpointRefDate => NULL(), &
      &                         restartRefDate    => NULL()
    TYPE(event), POINTER                   :: checkpointEvent   => NULL()
    TYPE(event), POINTER                   :: restartEvent      => NULL()
    INTEGER                                :: checkpointEvents
    TYPE(eventGroup), POINTER              :: checkpointEventGroup => NULL()
    TYPE(timedelta), POINTER               :: eventInterval     => NULL()
    TYPE(timedelta), POINTER               :: model_time_step   => NULL()
    TYPE(t_sim_step_info) :: sim_step_info
    LOGICAL                                :: lret, lwrite_checkpoint
    LOGICAL                                :: l_isStartdate, l_isRestart, l_isExpStopdate, l_isCheckpoint, l_doWriteRestart
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: dt_string, dstring
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)   :: td_string
    LOGICAL                                :: lprint_timestep, l_nml_output
    LOGICAL                                :: is_sitelevel

    INTEGER, ALLOCATABLE :: output_jfile(:) ! number of output files

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_model'
    ! ----------------------------------------------------------------------------------------------------- !

    CALL message('======== ', 'Starting the ICON-Land standalone model ...')
    CALL message('', '')
    CALL message('Namelist files', TRIM(jsb_namelist_filename) // ', ' // TRIM(shr_namelist_filename))
    CALL message('', '')

    ! n_dom: number of domains (ICON uses several domains for nested grids); ICON-Land supports only one
    ! (site-level runs mimic one domain per site, see build_sitelevel_decomposition)
    IF (n_dom > 1) CALL finish(routine, 'Only one atmospheric domain allowed')

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 1.0 Configuration of the ICON infrastructure used with ICON-Land
    !>
    !>   only 'CALL jsbach_setup_models()' (1.10) is specific to ICON-Land
    !>
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.1 Initialize global registry of lon-lat grids
    !>
    !>   ICON-Land does not use regular lon-lat grids, but the registry is part of the ICON infrastructure
    !>   (e.g. for output on lon-lat grids)
    !>
    CALL lonlat_grids%init()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.2 If this is a resumed or warm-start run, read the restart file meta data
    !>
    restartAttributes => NULL()
    IF (isRestart()) THEN
      CALL message('','Read restart file meta data ...')
      ! Reads attributes and namelists for all available domains from restart file.
      ! These namelist values will overwrite the model default, and will later be overwritten if the user has
      ! specified something different for this integration.
      CALL read_restart_header(TRIM(get_my_process_name()))
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.3 read namelists required for the ICON infrastructure and fill the
    !>     corresponding sections of the ICON configuration states
    !>
    ! Create a new file ('<namelist>_output') in which all the namelist variables and their
    ! actual values used in the model run will be stored.
    IF (my_process_is_stdio()) CALL open_nml_output(TRIM(jsb_namelist_filename)//'_output')

    ! Namelists of the ICON infrastructure (parallelization, run control, time, I/O, ...) needed by ICON-Land
    CALL read_infrastructure_namelists_for_icon_land(TRIM(jsb_namelist_filename), TRIM(shr_namelist_filename))

    ! Output namelists. The following can't be called from the previous subroutine since it's causing a
    ! circular dependency.
    CALL read_name_list_output_namelists (TRIM(jsb_namelist_filename))

    ! Close the file in which all the namelist variables and their
    ! actual values were stored and write an annotated table of all namelist settings to a text file
    IF (my_process_is_stdio()) THEN
      CALL close_nml_output()
      CALL log_nml_settings(TRIM(jsb_namelist_filename)//".log")
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.4 Cross-check the namelist setup of the parallelization (`parallel_nml`)
    !>
    CALL check_parallel_configuration()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.5 Set times and events, e.g. time step, start and stop dates, restart interval
    !>
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("jsbach", dt_restart, nsteps)
    CALL initEventManager(time_config%tc_exp_refdate)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.6 Call ICON configure_run to finish filling the run_config state.
    !>
    !>    This needs to be done very early (but anyway after atm_crosscheck)
    !>    because some components of the state, e.g., num_lev, may be
    !>    modified in this subroutine which affects the following CALLs.
    !>
    CALL configure_run( )

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.7 Initialize the MPI work groups, incl. processes dedicated to I/O and to writing restart files
    !>
    num_prefetch_proc = 0
    CALL restartWritingParameters(opt_dedicatedProcCount = dedicatedRestartProcs)
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, &
         &                          num_io_procs, dedicatedRestartProcs, &
         &                          jsbach_process, num_prefetch_proc, num_test_pe,      &
         &                          num_dio_procs=proc0_shift)

#ifndef __NO_RAGNAROK__
    CALL init_ragnarok_model()
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.8 Initialize the timers
    !>
    !>   Timers measure the CPU time of specific parts of the model (e.g. subroutines), which is written
    !>   to the log file at the end of the run (see module `mo_timer`, `ltimer` and `timers_level` in
    !>   `run_nml`). Some timers are started below in this module (`CALL timer_start()`).
    !>
    IF (ltimer) CALL init_timer
    ! start timer for the model initialization
    IF (timers_level > 1) CALL timer_start(timer_model_init)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.9 initialize dynamic list of vertical axes (CDI zaxes for output, also used for the ICON-Land
    !>     vertical grids)
    !>
    zaxisTypeList = t_zaxisTypeList()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.10 Setup the model in ICON-Land: read namelists, configure models for each domain
    !>
    !>   @NOTE This has to be after (!) the ICON zaxes have been created in the above line but
    !>         before (!) the restart PEs are detached a few lines below since ICON-Land
    !>         adds its zaxes to zaxisTypeList
    !>
    CALL jsbach_setup_models(shr_namelist_filename)
    no_models = get_no_of_models()

    ! Site-level runs: one model per site, all models have to be site-level; only one OpenMP thread
    is_sitelevel = .FALSE.
    DO jg=1,no_models
      model => Get_model(jg)
      IF (model%config%is_single_site) THEN
        is_sitelevel = .TRUE.
        EXIT
      END IF
    END DO
    no_omp_threads = Get_omp_no_of_threads()
    DO jg=1,no_models
      model => Get_model(jg)
      IF (model%config%is_single_site .AND. no_omp_threads /= 1) &
          & CALL finish(routine, 'Only one OMP thread allowed for site-level models!')
      IF (is_sitelevel .AND. .NOT. model%config%is_single_site) &
          & CALL finish(routine, 'Must not mix site-level and non-site-level models!')
    END DO
    IF (is_sitelevel) THEN
      CALL message('', 'This is a site-level run with '//int2string(no_models)//' sites')
      CALL message('', '')
      DO jg=1,no_models
        model => Get_model(jg)
        CALL message('---', 'Site '//int2string(jg))
        CALL message('Name', model%config%single_site_name)
        CALL message('Short name', model%config%single_site_shortname)
        CALL message('Longitude', real2string(model%config%single_site_lonlat(1)))
        CALL message('Latitude', real2string(model%config%single_site_lonlat(2)))
      END DO
    END IF
    CALL print_infrastructure_namelists_for_icon_land(no_models)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.11 I/O initialization
    !>
    ! The dedicated restart processes are detached here: they don't RETURN, but start their main loop
    ! in which they only write restart files when the restart event is due.
    CALL detachRestartProcs(ltimer)
    CALL init_io_processes()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.12 Import patches, perform domain decomposition
    !>
    !>   define the horizontal and vertical grids since they are aready
    !>   needed for some derived control parameters. This includes
    !>   - patch import
    !>   - domain decompistion
    !>   - vertical coordinates
    !>
    IF (timers_level > 4) CALL timer_start(timer_domain_decomp)
    IF (is_sitelevel) THEN
      CALL build_sitelevel_decomposition(num_lev, nshift, no_models)
    ELSE
      CALL build_decomposition(num_lev, nshift, is_ocean_decomposition = .false.)
    END IF
    IF (timers_level > 4) CALL timer_stop(timer_domain_decomp)

    p_patch(:)%ldom_active = .TRUE.

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.13 Construct interpolation state, compute interpolation coefficients.
    !>
    IF (timers_level > 4) CALL timer_start(timer_compute_coeffs)
    CALL configure_interpolation( n_dom, p_patch(1:)%level, &
                                  p_patch(1:)%geometry_info )

    ! Not actually used in ICON-Land, but needed by ICON (gridref refers to nested setups)
    CALL configure_gridref(n_dom, p_patch(1:)%geometry_info%mean_characteristic_length)

    ! Allocate array for interpolation state
    ALLOCATE( p_int_state(n_dom_start:n_dom), &
            & STAT=error_status)
    IF (error_status /= SUCCESS) THEN
      CALL finish(routine, 'allocation for ptr_int_state failed')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.14 Initialize icon_comm_lib
    !>
    CALL construct_icon_communication(p_patch, n_dom)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.15 Setup the information for the physical patches
    !>
    CALL setup_phys_patches

    ! Copy needed parts of the patches to GPU
    DO jg=1,no_models
      !$ACC ENTER DATA &
      !$ACC   COPYIN(p_patch(jg)%cells%decomp_info%owner_mask)
    END DO

    ! Constructing data for lon-lat interpolation
    ! CALL compute_lonlat_intp_coeffs(p_patch(1:), p_int_state(1:)) ! inout, inout

    CALL allocate_vct_atmo(p_patch(1)%nlevp1)

    IF (timers_level > 4) CALL timer_stop(timer_compute_coeffs)

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 2.0 Configuration and initialization of ICON-Land
    !>
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.1 Main setup of ICON-Land, incl. setup and initialization of the external forcing
    !>     (e.g. observation based data or output of another model simulation)
    !>
    ! setup forcing (only set up some variables)
    IF (force_from_observations()) CALL setup_forcing(no_models)  ! @TODO this might always be TRUE with icon-land standalone
    ! loop over ICON-Land models (sites in site-level runs)
    DO jg=1,no_models
      model => Get_model(jg)
      ! collect information from the ICON patch (per grid cell) into the ICON-Land grid structure
      CALL jsbach_setup_grid(jg, p_patch(jg), type='icon')
      ! read usecase and lctlib, and set up the tile structure of the usecase
      CALL jsbach_setup_tiles(jg)
      ! main ICON-Land model init, incl. reading ic/bc files and setting up the task queue
      CALL jsbach_init(jg)
      ! init forcing and read the forcing data
      IF (force_from_observations()) CALL init_forcing(jg, model%namelist_filename, REAL(dtime,KIND=lp))
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.2 restart a simulation: read model state from restart file(s)
    !>
    IF (isRestart()) THEN
      IF (timers_level > 4) CALL timer_start(timer_read_restart)
      ! n_dom: one domain per site in site-level runs
      DO jg = 1,n_dom
        IF (p_patch(jg)%ldom_active) THEN
          CALL read_restart_files( p_patch(jg), n_dom)
        END IF
      END DO
      CALL message(TRIM(routine),'normal exit from read_restart_files')

      IF (timers_level > 4) CALL timer_stop(timer_read_restart)

      ! Initializations that depend on the state read from the restart file (e.g. for land cover change)
      DO jg=1,no_models
        CALL jsbach_init_after_restart(jg)
      END DO

    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.3 Prepare output
    !>
    !>   ICON-Land does not use GRIB output, but `configure_gribout` is needed by the ICON infrastructure
    !>
    CALL configure_gribout(grid_generatingCenter, grid_generatingSubcenter, n_dom)
    CALL configure_io()

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.4 Map the variable groups given in the output namelist onto the corresponding variable subsets:
    !>
    !>   @NOTE all add_vars must be finished before calling this routine (done in jsbach_init above).
    !>   Each ICON-Land process has its own groups, see `mo_jsb_memory_class`.
    !>
    IF (output_mode%l_nml) THEN
      CALL parse_variable_groups()
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.5 Initialize the output (`output_nml`)
    !>
    !>   If async IO is in effect, init_name_list_output is a collective call
    !>   with the IO procs and effectively starts async IO
    !>
    IF (output_mode%l_nml) THEN
      ! compute sim_start, sim_end
      sim_step_info%sim_start    = time_config%tc_exp_startdate
      sim_step_info%sim_end      = time_config%tc_exp_stopdate
      sim_step_info%run_start    = time_config%tc_startdate
      sim_step_info%restart_time = time_config%tc_stopdate

      sim_step_info%dtime  = dtime
      sim_step_info%jstep0 = 0

      CALL getAttributesForRestarting(restartAttributes)
      ! get start counter for time loop from restart file:
      IF (restartAttributes%is_init) &
        & CALL restartAttributes%get("jstep", sim_step_info%jstep0)
      ! Set time levels to 1 in order to avoid confusing finish message about <var_name>.TL0
      ! for non-existent variable <var_name> in output_nml
      nold(:) = 1; nnow(:) = 1; nnew(:) = 1
      CALL init_name_list_output(sim_step_info)
      CALL create_mipz_level_selections(output_file)
      ! let CDI know about the vertical axes
      CALL create_vertical_axes(output_file)

    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.6 Initialize reset-Action, i.e. assign variables to action object
    !>
    !>   Actions are an ICON concept for something done at certain times; the reset-Action re-initializes
    !>   variables with their reset value (e.g. accumulated or MAX/MIN fields) after output
    !>
    CALL reset_act%initialize(ACTION_RESET)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.7 Optional: for `msg_level >= 15` print the variables of each output group to the log file
    !>
    IF (my_process_is_stdio() .AND. (msg_level >= 15)) THEN
      DO jg=1,n_dom
        CALL vlr_print_groups(idom=jg)
      END DO
    END IF

    IF (timers_level > 1) CALL timer_stop(timer_model_init)

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> End of model initialization
    !>
    !> 3.0 Preparations for the time loop
    !>
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.1 get and write out some of the initial values
    !>
    mtime_current => time_config%tc_current_date
    ! at the start of the experiment (no restart)
    IF (.NOT. isRestart() .AND. (mtime_current >= time_config%tc_exp_startdate)) THEN
      CALL update_statistics(lacc=.TRUE.)
      ! write initial values to output
      IF (output_mode%l_nml) THEN
          CALL write_name_list_output(jstep=0, lacc=.TRUE.)
      END IF

    END IF ! not isRestart()

    ! start timer for the model run without initialization
    IF (ltimer) CALL timer_start(timer_total)

    ! calculate elapsed simulation time in seconds
    ! TODO: sim_time and mtime_old (below) are not used anymore and could be removed
    sim_time = getElapsedSimTimeInSeconds(mtime_current)

    ! allocate temporary variable for restarting purposes
    ! (the restart files store the state of the output files, see time loop)
    IF (output_mode%l_nml) THEN
      ALLOCATE(output_jfile(SIZE(output_file)), STAT=error_status)
      IF (error_status /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed!')
    ENDIF

    mtime_old     => newDatetime(mtime_current) ! in

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.2 init restart descriptor for ICON-Land ("jsbach")
    !>
    restartDescriptor => createRestartDescriptor("jsbach")
    jstep0 = 0
    CALL getAttributesForRestarting(restartAttributes)
    IF (isRestart()) THEN
      ! get start counter for time loop from restart file:
      CALL restartAttributes%get("jstep", jstep0)
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.3 for debug purposes print var lists: for msg_level >= 13 short and for >= 20 long format
    !>
    IF  (msg_level >= 13) CALL vlr_print_vls(lshort=(msg_level < 20))

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.4 set up the events for writing checkpoint and restart files
    !>
    !>   events define when actions are done, e.g. writing output or restart files
    !>
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
    ! print the events to the log file
    CALL printEventGroup(checkpointEvents)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.5 set time loop properties
    !>
    ! length of the time step
    model_time_step => time_config%tc_dt_model

    CALL message('','')
    CALL datetimeToString(mtime_current, dstring)
    WRITE(message_text,'(a,a)') 'Start date of this run: ', dstring
    CALL message('',message_text)
    CALL datetimeToString(time_config%tc_stopdate, dstring)
    WRITE(message_text,'(a,a)') 'Stop date of this run:  ', dstring
    CALL message('',message_text)
    CALL message('','')

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 4.0 Time loop
    !>
    ! ----------------------------------------------------------------------------------------------------- !

    ! time step counter (jstep0 /= 0 for restarted runs)
    jstep = jstep0+1

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.1 time loop for running ICON-Land
    !>
    !>   Note: mtime_current is the date and time at the end of the time step
    !>
    TIME_LOOP: DO
      ! optional memory logging
      CALL memory_log_add

      mtime_current = mtime_current + model_time_step

      ! store state of output files for restarting purposes
      IF (output_mode%l_nml .AND. jstep>=0 ) THEN
        DO i=1,SIZE(output_file)
          output_jfile(i) = get_current_jfile(output_file(i)%out_event)
        END DO
      ENDIF

      ! print the time step to the log file (each time step for msg_level > 2, else once a day)
      IF (msg_level > 2) THEN
        lprint_timestep = .TRUE.
      ELSE
        lprint_timestep = MOD(jstep,timesteps_per_day(REAL(dtime,KIND=lp),1)) == 0
      ENDIF
      ! always print the first and the last time step
      lprint_timestep = lprint_timestep .OR. (jstep == jstep0+1) .OR. (jstep == jstep0+nsteps)

      IF (lprint_timestep) THEN
        WRITE(message_text,'(a,i8,a,i0,a,5(i2.2,a),i3.3)') &
              &             'Time step: ', jstep, ' model time ',                                &
              &             mtime_current%date%year,   '-', mtime_current%date%month,    '-',    &
              &             mtime_current%date%day,    ' ', mtime_current%time%hour,     ':',    &
              &             mtime_current%time%minute, ':', mtime_current%time%second,   '.',    &
              &             mtime_current%time%ms
        CALL message('',message_text)
      ENDIF

      ! is output to be written in this time step?
      l_nml_output = output_mode%l_nml .AND. jstep >= 0 .AND. istime4name_list_output(jstep)

      ! --------------------------------------------------------------------------------------------------- !
      !> 4.1.1 run ICON-Land for one timestep
      !>
      DO jg = 1,no_models
        CALL run_one_timestep( &
          & jg, &                   ! model
          & p_patch(jg)%nblks_c, &  ! number of blocks (last dimension of ICON-Land 2d variables)
          & REAL(dtime,KIND=lp), &  ! length of timestep in seconds
          & mtime_current)          ! date and time at the end of the time step
      END DO

      ! update accumulated values (e.g. for mean values in output)
      CALL update_statistics(lacc=.TRUE.)

      IF (l_nml_output) THEN
        CALL write_name_list_output(jstep, lacc=.TRUE.)
      ENDIF

      ! re-initialize MAX/MIN fields with 'resetval'
      !   @NOTE must be done AFTER output
      CALL reset_act%execute(slack=dtime, mtime_date=mtime_current)

      ! --------------------------------------------------------------------------------------------------- !
      !> 4.1.2 Write restart file
      !>
      !>  check whether time has come for writing restart file
      !>
      ! default is to assume we do not write a checkpoint/restart file
      lwrite_checkpoint = .FALSE.
      ! if the model is not supposed to write output, do not write checkpoints
      IF (.NOT. output_mode%l_none ) THEN
        ! to clarify the decision tree we use shorter and more expressive names:
        l_isStartdate    = (time_config%tc_startdate == mtime_current)
        l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
        l_isRestart      = is_event_active(restartEvent, mtime_current, proc0_offloading)
        l_isCheckpoint   = is_event_active(checkpointEvent, mtime_current, proc0_offloading)
        l_doWriteRestart = time_config%tc_write_restart
        ! ...
        IF ( &
              !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
              &         (l_isRestart .OR. l_isCheckpoint)                     &
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
      IF (lwrite_checkpoint) THEN
        ! update patch information in the restart descriptor (ICON infrastructure)
        DO jg=1,n_dom
          CALL restartDescriptor%updatePatch(p_patch(jg), opt_ndom=n_dom)
        END DO
        ! finally write restart file
        CALL restartDescriptor%writeRestart(mtime_current, jstep, opt_output_jfile = output_jfile)
      END IF  ! lwrite_checkpoint
      ! stop date of this run (job) reached?
      IF (mtime_current >= time_config%tc_stopdate) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

      ! increase timestep by one
      jstep = jstep + 1
      ! step_in_forcing_file = step_in_forcing_file + 1
      sim_time = getElapsedSimTimeInSeconds(mtime_current)
    END DO TIME_LOOP

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 5.0 Cleanup and finalization (mostly ICON infrastructure)
    !>
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.1 cleanup restart descriptor and stop timer for the model run
    !>
    CALL deleteRestartDescriptor(restartDescriptor)

    IF (ltimer) CALL timer_stop(timer_total)

    ! temporary variable for the state of the output files
    IF (output_mode%l_nml) THEN
      DEALLOCATE(output_jfile, STAT=error_status)
      IF (error_status /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed!')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.2 deallocation of structure and variables
    !>
    ! @TODO CALL deallocateDatetime(mtime_forcing)
    CALL deallocateDatetime(mtime_old)

    ! Deallocate global registry for lon-lat grids
    CALL lonlat_grids%finalize()

    ! Destruction fails for site-level setup
    IF (.NOT. is_sitelevel) THEN
      ! Destruct communication patterns
      CALL destruct_comm_patterns( p_patch, p_patch_local_parent )
    END IF

    ! Deallocate grid patches
    DO jg=1,no_models
      !$ACC EXIT DATA &
      !$ACC   DELETE(p_patch(jg)%cells%decomp_info%owner_mask)
    END DO
    IF (is_sitelevel) THEN
      CALL destruct_patches( p_patch(1:1) )
    ELSE
      CALL destruct_patches( p_patch )
    END IF

    DEALLOCATE( p_patch, STAT=error_status )
    IF (error_status/=SUCCESS) THEN
      CALL finish(routine, 'deallocate for patch array failed')
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.3 close memory logging files and finalize model
    !>
    CALL memory_log_terminate
    CALL destruct_icon_communication()
    ! ICON-Land cleanup and finalization
    CALL finalize_external_forcing(no_models)
    CALL jsbach_finalize()
#ifndef __NO_RAGNAROK__
    CALL finalize_ragnarok_model()
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.4 Delete variable lists
    !>
    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.5 finish message and print performance timers
    !>
    CALL message(routine, 'clean-up finished')
    IF (ltimer) CALL print_timer

    ! print performance timers:
    IF (ltimer) CALL print_timer

  END SUBROUTINE jsbach_model

  ! ======================================================================================================= !
  !>
  !> run ICON-Land for one timestep
  !>
  !> `run_one_timestep` does not call the full interface (`mo_jsb_interface:interface_full`) used in coupled
  !> runs, but most of the routines called from there. The main difference is where the A2L variables are
  !> set: in coupled runs they are passed from the atmosphere and set by calling `update_atm2land`; here they
  !> are directly set by the routine that reads the external forcing and computes the forcing for the time step.
  !>
  SUBROUTINE run_one_timestep(model_id, nblks, dtime, datetime_new)

    USE mo_jsb_control,               ONLY: jsbach_is_restarted, force_from_observations, &
      &                                     debug_on, timer_on, l_timer_host, timer_jsbach, timer_forcing
    USE mo_jsb_time,                  ONLY: is_time_experiment_start, is_time_restart, get_secs_of_day, get_year, is_newyear
    USE mo_jsb_model_class,           ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
    USE mo_jsb_class,                 ONLY: get_model
    USE mo_jsb_grid,                  ONLY: Get_grid
    USE mo_jsb_grid_class,            ONLY: t_jsb_grid
    USE mo_jsb_tile_class,            ONLY: t_jsb_tile_abstract
    USE mo_jsb_interface,             ONLY: jsbach_start_timestep, jsbach_finish_timestep
    USE mo_jsb4_forcing,              ONLY: forcing_options, get_standalone_driver
    USE mo_orbit_solar,               ONLY: compute_cos_zenith_angle
    USE mo_jsb_physical_constants,    ONLY: molarMassDryAir, molarMassCO2
    USE mo_jsb_parallel,              ONLY: Get_omp_thread
    USE mo_jsb_subset,                ONLY: ON_CHUNK

    USE mo_iq_atm2land_process,       ONLY: update_local_time_and_daytime_counter
    USE mo_quincy_model_config,       ONLY: QCANOPY
    USE mo_isotope_util,              ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C
    dsl4jsb_Use_processes A2L_, HYDRO_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                INTENT(in)            :: model_id       !< model ID (number of the site in site-level runs)
    INTEGER,                INTENT(in)            :: nblks          !< number of blocks
    REAL(lp)              , INTENT(in)            :: dtime          !< time step
    TYPE(datetime)        , POINTER               :: datetime_new   !< date and time at the end of this time step
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=max_timedelta_str_len) :: dtime_string     !< time delta as string
    TYPE(timedelta)         , POINTER    :: dtime_mtime      !< time delta as mtime variable
    TYPE(datetime)          , POINTER    :: datetime_old     !< date and time at the beginning of this time step
    TYPE(datetime)          , POINTER    :: datetime_old_old !< date and time one time step before datetime_old
    TYPE(datetime)          , POINTER    :: datetime_forcing !< TODO: not used anymore, can be removed

    INTEGER  :: nproma                         !< number of grid cells per block
    INTEGER  :: iblk                           !< index of block loop
    INTEGER  :: startblk, endblk               !< start and end indices of block   loop
    INTEGER  :: startidx(nblks), endidx(nblks) !< start and end indices of columns loop (per block)

    TYPE(t_jsb_model), POINTER           :: model
    TYPE(t_jsb_grid),  POINTER           :: grid       ! Horizontal grid
    CLASS(t_jsb_tile_abstract), POINTER  :: tile

    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real2D_onDomain :: evapotrans
    dsl4jsb_Real2D_onDomain :: evapopot
    dsl4jsb_Real2D_onDomain :: co2_air
    dsl4jsb_Real2D_onDomain :: co2_air_mol
    ! quincy
    dsl4jsb_Real2D_onDomain :: co2_mixing_ratio
    dsl4jsb_Real2D_onDomain :: CO2_mixing_ratio_c13
    dsl4jsb_Real2D_onDomain :: CO2_mixing_ratio_c14
    dsl4jsb_Real2D_onDomain :: swpar_srf_down
    dsl4jsb_Real2D_onDomain :: daytime_counter
    dsl4jsb_Real2D_onDomain :: daylength_prev_day
    dsl4jsb_Real2D_onDomain :: local_time_day_seconds

    REAL(lp), POINTER     :: lon(:,:), lat(:,:), coslat(:,:), sinlat(:,:)
    REAL(lp), ALLOCATABLE :: evapo_act2pot(:,:), co2_air_dc13(:,:), co2_air_Dc14(:,:)

    REAL(lp) :: hlp1   ! Temporary scalar for mixing ratio calculations
    LOGICAL :: use_external_forcing               ! use forcing from external forcing (e.g., observation based data)
    LOGICAL :: l_first_forcing_file    = .TRUE.   ! first forcing file in the experiment
    LOGICAL :: l_shift_one_timestep    = .TRUE.   ! TODO: not used anymore, can be removed
    LOGICAL :: l_first_step_in_restart = .FALSE.  ! first time step after a restart

    INTEGER :: nc, ics, ice, ic
    INTEGER :: global_seconds_day, current_year   ! global_seconds_day: seconds since midnight (UTC)

    ! IQJ spin-up accelerator model configurations
    INTEGER :: accelerator_max_executions !< bookkeeping configurations for the iq spin-up accelerator: max number of executions
    INTEGER :: accelerator_frequency      !< bookkeeping configurations for the iq spin-up accelerator: frequency of executions
    INTEGER :: accelerator_start_year     !< bookkeeping configurations for the iq spin-up accelerator: start year of executions
    ! IQJ: further configuration
    LOGICAL :: flag_snow = .TRUE. !< quincy boolean indicating if simulations should include snow (else read snow is added to rain)

    CHARACTER(len=*), PARAMETER :: routine = modname//':run_one_timestep'
    ! ----------------------------------------------------------------------------------------------------- !

    ! use forcing from external forcing, such as, observation based data (e.g., CRUJRA)
    !                                             or output from another model simulation
    ! see namelist 'jsb_control_nml' and 'mo_jsb_control:l_force_from_obs'
    use_external_forcing = force_from_observations()

    IF (jsbach_is_restarted()) THEN
      l_first_forcing_file = .FALSE.
      l_shift_one_timestep = .FALSE.
    ELSE
      l_shift_one_timestep = .TRUE.
    END IF

    ! Note: datetime_new is the date and time at the end of the time step
    CALL getPTStringFromSeconds(-dtime, dtime_string) ! ATTN: dtime with negative sign!
    dtime_mtime  => newTimedelta(dtime_string)
    datetime_old => newDatetime(datetime_new)
    datetime_old =  datetime_new + dtime_mtime ! Substract dtime from datetime_new to get datetime_old at beginning of the time step
    ! TODO old_old is not used anymore and should be deleted
    datetime_old_old => newDatetime(datetime_old)
    datetime_old_old =  datetime_old + dtime_mtime
    CALL deallocateTimedelta(dtime_mtime)

    IF (is_time_restart(datetime_old)) THEN
      l_first_step_in_restart = .TRUE.
    END IF

    ! TODO: datetime_forcing is not used anymore, the following IF block can be removed
    datetime_forcing => newDatetime(datetime_old)
    IF (l_first_forcing_file) THEN
      IF (forcing_options(model_id)%forcing_steps_per_day == 1  .AND. .NOT. use_external_forcing ) THEN
        ! Go to next time step of next file
        CALL getPTStringFromSeconds(forcing_options(model_id)%forcing_synchron_factor * dtime, dtime_string)
        dtime_mtime  => newTimedelta(dtime_string)
        datetime_forcing  =  datetime_forcing + dtime_mtime
        CALL deallocateTimedelta(dtime_mtime)
      END IF
    ELSE
      IF (.NOT. l_first_step_in_restart  .AND. .NOT. use_external_forcing) THEN
        datetime_forcing = datetime_new
      END IF
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !> start timestep
    !>
    CALL jsbach_start_timestep(model_id, datetime_old, REAL(dtime,KIND=wp)) ! datetime_old: start of the time step

    model => get_model(model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()

    CALL model%Get_top_tile(tile)

    ! In coupled runs, iblk, ics and ice are provided by the atmosphere; here they are taken from the grid
    startblk = grid%Get_blk_start()
    endblk   = grid%Get_blk_end()
    DO iblk = startblk, endblk
      startidx(iblk) = grid%Get_col_start(iblk)
      endidx  (iblk) = grid%Get_col_end(iblk)
    END DO
    !$ACC DATA COPYIN(startidx, endidx) ASYNC(1)

    IF (use_external_forcing) THEN

      dsl4jsb_Get_memory(A2L_)

      ! lon and lat are required for the weather generator
      lon => grid%lon
      lat => grid%lat
      coslat => grid%coslat
      sinlat => grid%sinlat

      dsl4jsb_Get_var2d_onDomain(A2L_, co2_air_mol)
      dsl4jsb_Get_var2d_onDomain(A2L_, co2_air)

      CALL compute_cos_zenith_angle(datetime_old, grid%patch, dsl4jsb_var2D_onDomain(A2L_, cos_zenith_angle))
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(A2L_, cos_zenith_angle)) ASYNC(1)
      ! Note: inquire_declination should only be called each time step after compute_cos_zenith_angle

      ! TODO: if really required calculate daily sums of evapotrans and evapopot in jsbach and get them here (see comment below)
      ! -> would only be required if no Qair forcing is given and qair needs to be calculated from temperature!
      ALLOCATE(evapo_act2pot(nproma, nblks))
      !$ACC ENTER DATA CREATE(evapo_act2pot)

      dsl4jsb_Get_memory(HYDRO_)
      dsl4jsb_Get_var2d_onDomain(HYDRO_, evapotrans)
      dsl4jsb_Get_var2d_onDomain(HYDRO_, evapopot)

      IF (is_time_experiment_start(datetime_old)) THEN

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            evapo_act2pot(ic, iblk) = 1._lp
          END DO
        END DO
        !$ACC END PARALLEL LOOP

      ELSE

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            evapo_act2pot(ic, iblk) =   MIN(evapotrans(ic, iblk), -EPSILON(1._lp)) &
              &                       / MIN(evapopot  (ic, iblk), -EPSILON(1._lp))
          END DO
        END DO
        !$ACC END PARALLEL LOOP

      END IF

      IF (model%config%model_scheme == MODEL_QUINCY) THEN
       ! Allocate and enter data for C13 and C14 based on configuration flags
        ALLOCATE(co2_air_dc13(nproma, nblks))
        !$ACC ENTER DATA CREATE(co2_air_dc13)
        ALLOCATE(co2_air_Dc14(nproma, nblks))
        !$ACC ENTER DATA CREATE(co2_air_Dc14)
      END IF

      ! timer for reading and computing the forcing
      IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_forcing(model%id))

      ! --------------------------------------------------------------------------------------------------- !
      !> get the forcing variables and directly store them in the A2L_ memory
      !>   In coupled runs, these variables are provided by the atmosphere via the model interface
      !>   (see 'mo_jsb_interface:interface_full') and set by update_atm2land, which therefore does not need
      !>   to be called in standalone runs.
      !>

      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_JSBACH)
        CALL get_standalone_driver( &
          ! in
          & model_id,                                                                  &
          & model%config%model_scheme,                                                 &
          & flag_snow,                                                                 &
          & startblk, endblk, startidx, endidx, datetime_old, datetime_new,            &
          & dsl4jsb_var2D_onDomain(HYDRO_, elevation),                                 &
          & sinlat, coslat, lon, evapo_act2pot,                                        &
          & dsl4jsb_var2D_onDomain(A2L_,  cos_zenith_angle),                           &
          ! out
          & dsl4jsb_var2D_onDomain(A2L_, CO2_air),                                     &
          & dsl4jsb_var2D_onDomain(A2L_, t_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, q_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, rain),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, snow),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, wind_air),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, wind_10m),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, lw_srf_down),                                 &
          & dsl4jsb_var2D_onDomain(A2L_, swvis_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swnir_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swpar_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, fract_par_diffuse),                           &
          & dsl4jsb_var2D_onDomain(A2L_, press_srf)                                    &
          & )
      CASE (MODEL_QUINCY)
        CALL get_standalone_driver( &
          ! in
          & model_id,                                                                  &
          & model%config%model_scheme,                                                 &
          & flag_snow,                                                                 &
          & startblk, endblk, startidx, endidx, datetime_old, datetime_new,            &
          & dsl4jsb_var2D_onDomain(HYDRO_, elevation),                                 &
          & sinlat, coslat, lon, evapo_act2pot,                                        &
          & dsl4jsb_var2D_onDomain(A2L_,  cos_zenith_angle),                           &
          ! out
          & dsl4jsb_var2D_onDomain(A2L_, CO2_air),                                     &
          & dsl4jsb_var2D_onDomain(A2L_, t_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, q_air),                                       &
          & dsl4jsb_var2D_onDomain(A2L_, rain),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, snow),                                        &
          & dsl4jsb_var2D_onDomain(A2L_, wind_air),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, wind_10m),                                    &
          & dsl4jsb_var2D_onDomain(A2L_, lw_srf_down),                                 &
          & dsl4jsb_var2D_onDomain(A2L_, swvis_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swnir_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, swpar_srf_down),                              &
          & dsl4jsb_var2D_onDomain(A2L_, fract_par_diffuse),                           &
          & dsl4jsb_var2D_onDomain(A2L_, press_srf),                                   &
          & co2_air_dc13,                                                              &
          & co2_air_Dc14                                                               &
          & )
      END SELECT

      IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_forcing(model%id))

      ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk = 1, nblks
        DO ic = 1, nproma
          co2_air_mol(ic, iblk) = co2_air(ic, iblk) * molarMassDryAir / molarMassCO2
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP

      ! --------------------------------------------------------------------------------------------------- !
      !> code specific to the QUINCY model
      !>
      SELECT CASE (model%config%model_scheme)
      CASE (MODEL_QUINCY)

#ifdef _OPENACC
        IF (model%config%qmodel_id /= QCANOPY) &
          & CALL finish(routine, 'QUINCY is currently only ported to GPU for the CANOPY mode. '//model%config%quincy_model)
#endif

        ! convert CO2 from "molar ratio (volume)" to "co2 mixing ratio ppmv"
        dsl4jsb_Get_var2d_onDomain(A2L_, co2_mixing_ratio)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            co2_mixing_ratio(ic, iblk) = co2_air_mol(ic, iblk) * 1000000._lp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
        dsl4jsb_Get_var2d_onDomain(A2L_, CO2_mixing_ratio_c13)
        dsl4jsb_Get_var2d_onDomain(A2L_, CO2_mixing_ratio_c14)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) PRIVATE(hlp1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            hlp1 = calc_mixing_ratio_C13C12(co2_air_dc13(ic, iblk))
            CO2_mixing_ratio_c13(ic, iblk)  = co2_mixing_ratio(ic, iblk) / (1._lp + 1._lp / hlp1)
            hlp1 = calc_mixing_ratio_C14C(co2_air_dc13(ic, iblk), co2_air_Dc14(ic, iblk))
            CO2_mixing_ratio_c14(ic, iblk) = co2_mixing_ratio(ic, iblk) * hlp1
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC WAIT(1)
        !$ACC EXIT DATA DELETE(co2_air_dc13, co2_air_Dc14)
        DEALLOCATE(co2_air_dc13, co2_air_Dc14)

        ! Update local time and daytime counter (daytime is determined from radiation); these variables are
        ! used in several QUINCY processes to determine local midnight or noon and to track variables
        ! while there is shortwave radiation
        dsl4jsb_Get_var2d_onDomain(A2L_, daytime_counter)
        dsl4jsb_Get_var2d_onDomain(A2L_, daylength_prev_day)
        dsl4jsb_Get_var2d_onDomain(A2L_, local_time_day_seconds)
        dsl4jsb_Get_var2d_onDomain(A2L_, swpar_srf_down)
        global_seconds_day = NINT(get_secs_of_day(datetime_new)) ! TODO time steps w/ fractional seconds not supported by quincy, yet!

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            CALL update_local_time_and_daytime_counter( &
              &     global_seconds_day, dtime, lon(ic, iblk), swpar_srf_down(ic, iblk), &
              &     daytime_counter(ic, iblk), daylength_prev_day(ic, iblk), local_time_day_seconds(ic, iblk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP

      END SELECT

      !$ACC WAIT(1)

      !$ACC EXIT DATA DELETE(evapo_act2pot)
      DEALLOCATE(evapo_act2pot)

    END IF ! IF (use_external_forcing) THEN

    ! The following is also done in the interface in coupled runs
!$OMP PARALLEL DO PRIVATE(iblk,ics,ice,nc) IF(nblks > 1)
    ! ----------------------------------------------------------------------------------------------------- !
    !> loop over the chunks of the block
    !>
    DO iblk = startblk, endblk
      ! index of gridcell at start of the chunk
      ics = startidx(iblk)
      ! index of gridcell at end of the chunk
      ice = endidx  (iblk)

      IF (ics > ice) CYCLE

      IF (timer_on() .OR. l_timer_host) CALL timer_start(timer_jsbach(model%id))

      ! calc size of chunk, i.e., number of gridcells
      nc = ice - ics + 1

      IF (debug_on('basic') .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

      ! Update options
      CALL model%Set_options(iblk=iblk, ics=ics, ice=ice, nc=nc, current_datetime=datetime_old, &
        &                    dtime=dtime, alpha=1._lp)
      CALL model%Set_subset(type=ON_CHUNK, iblk=iblk, ics=ics, ice=ice)

      ! Run tasks in queue
      CALL model%Run_tasks(debug_on('tasks') .AND. iblk == 1)

      ! @TODO cleanup / remove ?
      !!$    IF (iblk == 1)                                                              &
      !!$      & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//          &
      !!$      &                             ' on PE '//int2string(p_pe)//' completed.', &
      !!$      &              all_print=.TRUE.)

      IF (debug_on('basic') .AND. iblk == 1) &
        & CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

      IF (timer_on() .OR. l_timer_host) CALL timer_stop(timer_jsbach(model%id))

    END DO
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> finish timestep
    !>
    CALL jsbach_finish_timestep(model_id, datetime_old, REAL(dtime,KIND=wp))

    CALL deallocateDatetime(datetime_old)
    CALL deallocateDatetime(datetime_old_old)
    CALL deallocateDatetime(datetime_forcing)

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE run_one_timestep

  ! ====================================================================================================== !
  !
  !> Build decomposition for site-level setup
  !
  ! This is a stripped-down version of `build_decomposition` from ICON for site-level setup.
  !
  ! Only one grid (patch) with one cell is read in (n_dom=1), but this patch is then replicated
  ! to n_dom=nsites domains
  !
  SUBROUTINE build_sitelevel_decomposition(num_lev,nshift, nsites)

    USE mo_model_domain,         ONLY: t_patch
    USE mo_pre_model_domain,     ONLY: t_pre_patch
    USE mo_model_domimp_patches, ONLY: import_pre_patches, complete_patches
    USE mo_setup_subdivision,    ONLY: decompose_domain
    USE mo_complete_subdivision, ONLY: finalize_decomposition, complete_parallel_setup
    USE mo_communication,        ONLY: setup_comm_gather_pattern
    USE mo_sync,                 ONLY: decomposition_statistics
    USE mo_grid_config,          ONLY: n_dom, n_dom_start
    USE mo_model_domain,         ONLY: p_patch
    USE mo_impl_constants,       ONLY: SUCCESS, max_dom
    USE mo_mpi,                  ONLY: my_process_is_mpi_parallel, proc_split, p_comm_work, p_pe_work


    INTEGER, INTENT(in) :: num_lev(max_dom),       &
      &                    nshift(max_dom),        &
      &                    nsites

    ! local variables
    TYPE(t_pre_patch), ALLOCATABLE :: p_patch_pre(:)
    TYPE(t_patch),     ALLOCATABLE :: p_patch_temp(:)
    INTEGER                        :: error_status, jg
    !> If .true., read fields related to grid refinement from separate  grid files:
    LOGICAL                        :: lsep_grfinfo
    LOGICAL                        :: is_ocean_decomposition

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':build_sitelevel_decomposition'

    is_ocean_decomposition = .TRUE. ! Just needed for `complete_patches` and `finalize_decomposition`

    ! Check p_patch allocation status
    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(routine), 'p_patch already allocated')
    END IF

    ALLOCATE(p_patch(1:nsites), stat=error_status)
    IF (error_status/=success) CALL finish(routine, 'allocation of p_patch failed')

    n_dom_start = 1

    ! -------------------------------------------------------------------------------------------------------
    !
    ! Read in one patch
    !
    n_dom = 1
    ALLOCATE(p_patch_temp(1:1))
    ALLOCATE(p_patch_pre(1:1))
    CALL import_pre_patches(p_patch_pre,num_lev,nshift,lsep_grfinfo)
    CALL decompose_domain(p_patch_temp, p_patch_pre)
    DEALLOCATE(p_patch_pre)

    ! Complete information which is not yet read or calculated
    CALL complete_patches( p_patch_temp, is_ocean_decomposition, lsep_grfinfo)

    CALL finalize_decomposition(p_patch_temp, is_ocean_decomposition)

    ! -------------------------------------------------------------------------------------------------------
    !
    ! Now mimic n_dom=nsites domains
    !
    n_dom = nsites

    DO jg=1,n_dom
      p_patch(jg) = p_patch_temp(1)
      p_patch(jg)%id = jg
      p_patch(jg)%cells%phys_id = jg
      p_patch(jg)%edges%phys_id = jg
      p_patch(jg)%verts%phys_id = jg
      p_patch(jg)%comm_pat_gather_c = p_patch_temp(1)%comm_pat_gather_c
      p_patch(jg)%comm_pat_c => p_patch_temp(1)%comm_pat_c
    END DO

    ! -------------------------------------------------------------------------------------------------------
    !
    ! This is a stripped-down version of `complete_parallel_setup` and the routines `set_patch_communicators`
    ! and `set_comm_pat_gather` called from there
    !
    DO jg=1,n_dom
      ! Set communication patterns for gathering on proc 0
      CALL setup_comm_gather_pattern( &
        p_patch(jg)%cells%decomp_info%glb2loc_index%global_size, &
        p_patch(jg)%cells%decomp_info%owner_local, &
        p_patch(jg)%cells%decomp_info%glb_index, p_patch(jg)%comm_pat_gather_c)
      CALL setup_comm_gather_pattern( &
        p_patch(jg)%verts%decomp_info%glb2loc_index%global_size, &
        p_patch(jg)%verts%decomp_info%owner_local, &
        p_patch(jg)%verts%decomp_info%glb_index, p_patch(jg)%comm_pat_gather_v)
      CALL setup_comm_gather_pattern( &
        p_patch(jg)%edges%decomp_info%glb2loc_index%global_size, &
        p_patch(jg)%edges%decomp_info%owner_local, &
        p_patch(jg)%edges%decomp_info%glb_index, p_patch(jg)%comm_pat_gather_e)
    END DO

    proc_split = .FALSE.
#ifndef NOMPI
    p_patch(:)%comm   = p_comm_work
    p_patch(:)%rank   = p_pe_work
#endif
    ! -------------------------------------------------------------------------------------------------------

    IF (my_process_is_mpi_parallel()) THEN ! the call below hangs in test mode
      ! Print diagnostic information about domain decomposition
      DO jg = 1, n_dom
        CALL decomposition_statistics(p_patch(jg))
      END DO
    END IF

    CALL message("build_sitelevel_decomposition", "Done.")

  END SUBROUTINE build_sitelevel_decomposition

#endif
END MODULE mo_jsbach_model
