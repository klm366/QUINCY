!> Driver for the QUINCY standalone model
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### main routine for running QUINCY standalone
!>

PROGRAM quincy_qs_model_driver

  ! echam
  USE mo_kind
  USE mo_util,                    ONLY: int2string
  USE mo_qs_time_util,            ONLY: calc_and_report_estimated_remaining_simulation_time, &
    &                                   open_quincy_time_output_files, close_quincy_time_output_files
  USE mo_exception,               ONLY: message_text, finish, message
  USE mo_io_units,                ONLY: nlog, nerr                ! LUN of log & error file
  ! land model
  USE mo_qs_model_interface,      ONLY: land_model_interface 
  ! base
  USE mo_jsb_math_constants,      ONLY: one_year, eps4, seconds_per_day
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_init,          ONLY: jsbach_setup_grid, &      ! interface to e.g. setup_grid_standalone
                                        quincy_init, &            ! JSBACH4: jsbach_init
                                        time_day, time_year, ntimestep, itimestep, timestep_prop_day, &
                                        flux_trans_timestep_out_state
  USE mo_jsb_class,               ONLY: get_model
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_process_class,       ONLY: A2L_, Q_SYL_
  USE mo_quincy_model_config,     ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
  USE mo_jsb_base,                ONLY: setup_model_vgrid_procconf, quincy_proc_parameter_init, jsbach_setup_tiles
  USE mo_jsb_time_iface,          ONLY: get_is_time_experiment_start, set_is_time_experiment_start
  USE mo_lookup_tables_constants, ONLY: init_satpres_coeffs
  ! forcing
  USE mo_qs_forcing,              ONLY: forcing_options, calc_output_time_interval_pool
  ! atmland
  USE mo_qs_atmland_forcing,      ONLY: update_forcing, open_forcing_input_file, close_forcing_input_file ! valid for land plant canopy and soil
  ! output
  USE mo_qs_output_txt,           ONLY: open_pool_output_files, &
    &                                   open_flux_output_files, &
    &                                   close_output_files, &
    &                                   open_output_files_soil_forcing, &
    &                                   close_output_files_soil_forcing 
  USE mo_qs_output_txt_util,      ONLY: init_output_txt_config, print_output_txt_config
  !> NOTE: writing of soil forcing still required if 'l_output_as_ascii_text = .FALSE.', i.e. if 'normal ascii' output is written to netcdf
  USE mo_qs_output_netcdf4,       ONLY: assert_possible_output_interval, init_quincy_nc4_output, define_nc4_output_files, &
    &                                   prepare_nc4_output_writing, close_nc4_output_files, deallocate_quincy_output

  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_config(Q_SYL_)

  IMPLICIT NONE

  !! 0.2 -----------------------------------------------------------------------------------------------
  !! 
  !! Local variables
  !! 
  LOGICAL       :: l_forcing                          ! IF .TRUE. THEN run forcing routines (e.g. update_forcing())
  INTEGER       :: first_timestep                     ! number of the first timestep simulated in the below DO loop
  LOGICAL       :: is_first_timestep                  ! is TRUE during the first timestep, i.e., experiment start
  INTEGER       :: transient_spinup_itimestep_period  ! count the number of timesteps covering one repeated spinup period
  INTEGER       :: nblks, iblk, nproma, ics, ice, nc  !< *slmtmp may be an input from JSBACH4 later in SLM dev
  INTEGER       :: npromz                             ! number of points in the last chunk (in QUINCY we only have one chunk)
  INTEGER       :: n_output_years                     !< number of years to be written into a to be prepared netcdf file
  INTEGER       :: nproca                             ! number of processors for jsbach
  INTEGER       :: nprocb                             ! has to be one
  INTEGER       :: nprocio                            ! number of processors for I/O server
  INTEGER       :: npedim                             ! Working dimension for blocks in each domain.
                                                      ! Default (-1): each domain is processed in one call
  INTEGER       :: model_id
  ! output control
  LOGICAL       :: l_output, &                        ! output is enabled (flux accumulation is done and l_output_interval_match is checked for being TRUE)
                   l_output_as_ascii_text, &          ! if output to txt or (if false) to nc4
                   l_output_interval_match            ! the specific output interval is a match with current time, hence output is written
  REAL(wp)      :: output_interval_pool_act, &        ! the actual interval number used at the current timestep
                   output_interval_pool, &            ! output interval of pools as a number (the same is used for flux output), converted from character
                   output_interval_pool_spinup        ! output interval of pools as a number valid for spinup period
  REAL(wp)      :: transient_spinup_timesteps         ! timesteps of spinup (transient forcing) as positive number, needed for output control
  REAL(wp)      :: elapsed_seconds                    ! seconds elapsed since start of this simulation (separate for spin-up, transient run, timestep output)

  CHARACTER(len=20) :: start_date_char
    !< start date of the simulation - used for the time unit in the netcdf4 files -- yyyy-mm-dd hh:mm
  INTEGER           :: nr_forcing_years
    !< number of forcing years available from the forcing file - required for calculating the start date char for static forcings

  INTEGER :: buffer_nc_index !< position in the netcdf output buffer to write to next (if writing to netcdf)
  INTEGER :: buffer_nc_size  !< size of the buffer used for netcdf output (if writing to netcdf)

  dsl4jsb_Def_config(Q_SYL_)

  ! A2L_
  dsl4jsb_Def_memory(A2L_)
  dsl4jsb_Real2D_onChunk :: t_air
  dsl4jsb_Real2D_onChunk :: q_air
  dsl4jsb_Real2D_onChunk :: rain
  dsl4jsb_Real2D_onChunk :: snow
  dsl4jsb_Real2D_onChunk :: wind_air
  dsl4jsb_Real2D_onChunk :: wind_10m
  dsl4jsb_Real2D_onChunk :: lw_srf_down
  dsl4jsb_Real2D_onChunk :: swvis_srf_down
  dsl4jsb_Real2D_onChunk :: swnir_srf_down
  dsl4jsb_Real2D_onChunk :: swpar_srf_down
  dsl4jsb_Real2D_onChunk :: fract_par_diffuse
  dsl4jsb_Real2D_onChunk :: press_srf
  dsl4jsb_Real2D_onChunk :: cos_zenith_angle
  dsl4jsb_Real2D_onChunk :: co2_mixing_ratio
  dsl4jsb_Real2D_onChunk :: co2_mixing_ratio_C13
  dsl4jsb_Real2D_onChunk :: co2_mixing_ratio_C14
  dsl4jsb_Real2D_onChunk :: nhx_deposition
  dsl4jsb_Real2D_onChunk :: noy_deposition
  dsl4jsb_Real2D_onChunk :: nhx_n15_deposition
  dsl4jsb_Real2D_onChunk :: noy_n15_deposition
  dsl4jsb_Real2D_onChunk :: p_deposition

  REAL(wp), ALLOCATABLE, DIMENSION(:) :: time_day_seconds
  CHARACTER(len=*), PARAMETER :: modname = 'quincy_qs_model_driver'

  !! 1.0 -----------------------------------------------------------------------------------------------
  !! 
  !! setup model, configure and init processes and grids
  !! 
  TYPE(t_jsb_model)         , POINTER :: model  !< local pointer to model
  CLASS(t_jsb_tile_abstract), POINTER :: tile   !< local data structure for dynamic vegetation and soil variables
                                                !! in jsbach4 tile is created in mo_jsb_model_usecases

  !! open log & error file
  OPEN (UNIT = nlog, FILE = 'quincy_standalone.log', ACTION="write", STATUS="unknown")
  OPEN (UNIT = nerr, FILE = 'quincy_standalone.err', ACTION="write", STATUS="unknown")

  !! output can currently go to text files or to netcdf4 files: decided upon compiling via cpp-key
  l_output_as_ascii_text = .TRUE.
#ifdef __QUINCY_WITH_NETCDF__
  l_output_as_ascii_text = .FALSE.
#endif

  PRINT*," "
  PRINT*,"#########################"
  PRINT*," Start QUINCY model Init "

  !! 1.1  define default values
  nproca   = 1      ! ?
  nprocb   = 1      ! ?
  nprocio  = 0      ! ?
  npedim   = -1     ! ?
  model_id = 1      ! must be consistent with model_idx set in setup_model_vgrid_procconf()
  nblks    = 1                ! must be always 1 in QUINCY (number of chunk in the "2D block of chunks" created from the 1D domain - information provided by ECHAM or ICON to JSBACH/QUINCY)
  iblk     = 1                ! must be always 1 in QUINCY (number of chunk in the "2D block of chunks" created from the 1D domain - information provided by ECHAM or ICON to JSBACH/QUINCY)
  nproma   = 1                ! number of points in the chunk
  npromz   = nproma           ! number of points in the last chunk (in QUINCY we only have one chunk anyway)
  ics      = 1                ! must be 1 in QUINCY (and is also always 1 in "jsb4 standalone with ECHAM") - mo_jsb_memory_class:Add_var calculates ice using nproma & based on the assumption of ics=1 
  ice      = ics + nproma - 1 ! ice depends on ics and nproma; ics is always 1, hence, ice = nproma always
  nc       = nproma           ! in QUINCY npromz = nproma always

  buffer_nc_index = 1

  !! 1.2  set up the model with vgrid and process configs (jsbach%models(model_id)%m => model  with 'model' of TYPE(t_jsb_model))
  CALL setup_model_vgrid_procconf("qs.namelist")
  model => get_model(model_id)

  !! 1.3 init of horizontal grid
  CALL jsbach_setup_grid(model%id, "hgrid0", nproma, npromz, nblks)  ! model%id, horizontal grid name (arbitrary name), number of points

  !! 1.4 init parameter of NWP atmosphere physics used for calculation of 'qsat' (saturation vapor pressure)
  CALL init_satpres_coeffs(2)   ! apply "more accurate IFS coefficients"

  !! 1.5 init process parameters (and modify parameter values for param-sensi)
  CALL quincy_proc_parameter_init(model%id)

  !! 1.6 set up tiles, read lctlib (and modify parameter values for param-sensi) and init usecase
  CALL jsbach_setup_tiles(model%id)

  !! 1.7 quincy-model init
  !!     and init of: forcing, SB_ & VEG_ pools, tiles, time, fluxnet diagnostics
  CALL quincy_init(model%id, iblk, nproma, nblks, ics, ice, nc, &           ! in
                   l_forcing)                                               ! out

  !! 1.8 init ascii txt and netcdf4 simulation output files
  IF (l_output_as_ascii_text) THEN
    ! init config
    CALL init_output_txt_config(TRIM(model%config%output_txt_file_set))
    ! print config to logfile
    CALL print_output_txt_config()
  ELSE
    ! assert possible output intervals
    ! spinup
    output_interval_pool = calc_output_time_interval_pool(model%config%output_interval_pool_spinup, model%config%dtime)
    CALL assert_possible_output_interval(output_interval_pool, model%config%output_interval_pool_spinup, modname, 'spin-up')
    ! transient simulation
    output_interval_pool = calc_output_time_interval_pool(model%config%output_interval_pool, model%config%dtime)
    CALL assert_possible_output_interval(output_interval_pool, model%config%output_interval_pool, modname)
    ! init netcdf4 simulation output
    CALL init_quincy_nc4_output(model%id)
  END IF

  !! 1.9 get the tile
  !!     note: in this QUINCY model driver just one tile is used
  CALL model%Get_top_tile(tile)

  !! 1.10 define the netcdf4 output files and their output variables
  !!      depending on the `model%config%qmodel_id`
  IF (.NOT. l_output_as_ascii_text) THEN
    CALL define_nc4_output_files(model%id, tile)

    IF ( (TRIM(model%config%output_interval_pool) .NE. TRIM(model%config%output_interval_flux)) .OR. &
        & (TRIM(model%config%output_interval_pool_spinup) .NE. TRIM(model%config%output_interval_flux_spinup)) ) THEN
      WRITE(message_text,*) 'Quincy nc4 output module expects output intervals to be identical for fluxes and pools, please check!'
      CALL finish(TRIM(modname),message_text)
    END IF
  END IF


  !! 2.0 -----------------------------------------------------------------------------------------------
  dsl4jsb_Get_config(Q_SYL_)
  dsl4jsb_Get_memory(A2L_)
  dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
  dsl4jsb_Get_var2D_onChunk(A2L_, q_air)
  dsl4jsb_Get_var2D_onChunk(A2L_, rain)
  dsl4jsb_Get_var2D_onChunk(A2L_, snow)
  dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)
  dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)
  dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)
  dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)
  dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)
  dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
  dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)
  dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
  dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)
  dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)
  dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C13)
  dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C14)
  dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
  dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
  dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
  dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
  dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
  !! 
  !! allocate forcing variable
  ALLOCATE(time_day_seconds(nc))

  !! 2.1 -----------------------------------------------------------------------------------------------
  !! 
  !! open forcing files
  !! 
  ! open forcing input files (and read the header into dummy variables); see mo_qs_atmland_forcing
  ! l_forcing is set in init_slm() depending on the specific SLM model
  ! for SLM soil different input is used than for land plant canopy
  IF (l_forcing) THEN
    CALL open_forcing_input_file(model%config%qmodel_id,tile)
  END IF


  !! 3.0 -----------------------------------------------------------------------------------------------
  !! 
  !! settings for the static/transient forcing_mode
  !! 
  ! get output interval of pools as a number of days
  output_interval_pool_act    = 0.0_wp
  output_interval_pool        = calc_output_time_interval_pool(model%config%output_interval_pool, &
    &                           model%config%dtime)
  output_interval_pool_spinup = calc_output_time_interval_pool(model%config%output_interval_pool_spinup, &
    &                           model%config%dtime)
  ! init with default values
  model%config%l_transient_forcing_mode          = .FALSE.
  model%config%l_transient_spinup                = .FALSE.
  model%config%l_transient_spinup_rewind_forcing = .FALSE.
  model%config%l_static_forcing_mode             = .FALSE.
  itimestep                         = 0             ! default start of time loop
  first_timestep                    = itimestep + 1 ! the first timestep in the below time loop is "itimestep + 1"
  transient_spinup_itimestep_period = 0  
  transient_spinup_timesteps        = 0             ! is used only during l_transient_spinup
  ! set simulation values
  IF ( TRIM(forcing_options%forcing_mode) == "static" ) THEN
    model%config%l_static_forcing_mode     = .TRUE.
    output_interval_pool_act               = output_interval_pool   ! output interval in days
  ELSEIF ( TRIM(forcing_options%forcing_mode) == "transient" ) THEN
    ! set logicals to TRUE
    model%config%l_transient_forcing_mode  = .TRUE.
    model%config%l_transient_spinup        = .TRUE.  ! transient forcing_mode starts with spinup; is set to FALSE during runtime when spinup is done
    ! output interval for spinup
    output_interval_pool_act          = output_interval_pool_spinup
    ! number of spinup timesteps
    itimestep                         = itimestep - forcing_options%transient_spinup_timesteps_all
    first_timestep                    = itimestep + 1             ! the first timestep in the below time loop is "itimestep + 1"
    transient_spinup_timesteps        = itimestep * (-1.0_wp)
    ! init spinup period counter
    transient_spinup_itimestep_period = 0
  ENDIF
  is_first_timestep = .TRUE.        ! TRUE for the first timestep, i.e., experiment start
  nr_forcing_years  = model%config%forcing_file_last_yr - model%config%forcing_file_start_yr + 1


  !! 4.0 -----------------------------------------------------------------------------------------------
  !! 
  !! open output files
  !! need to be done after 'output_interval_pool_act' is defined above
  !! test_canopy & test_radiation: special case, runs for just one timestep and pool/flux output routines are not called
  !! 
  !> 4.1 ascii pool & flux output files
  !>
  ! Initialise output and flux accumulation, and open output files; see mo_qs_output_txt
  ! the output routine for SLM canopy is located in mo_qs_model_interface -> canopy_model_write_data()
  SELECT CASE(model%config%qmodel_id)
    CASE(QLAND, QPLANT, QSOIL, QCANOPY)
      ! open output files for static mode or for transient spinup if output may be created (i.e. output interval is not 'none')
      IF (model%config%l_static_forcing_mode) THEN
        IF (l_output_as_ascii_text) THEN
          CALL open_pool_output_files(TRIM(model%config%output_interval_pool), '')
          CALL open_flux_output_files(tile, model%options, TRIM(model%config%output_interval_flux), '')
        ELSE
          ! In static forcing mode we cannot trust the simulation year offset as for longer
          ! simluations the model cycles the forcing years
          PRINT *, "Warning: Static forcing might produce wrong output years in netcdf output."
          PRINT *, "If you want to ensure to have the correct years associated with the output select transient forcing."

          ! Assert that the simulation length of static simulations is given in years 
          ! because this as assumed in the calculation of the startdate for the netcdf output
          IF(TRIM(forcing_options%simulation_length_unit) /= 'y') THEN
            WRITE(message_text,*) "Determination of startdate in netcdf output of static simulations presumes that "  &
              &                 // "the simulation length is specified in years ('y'). However, specified unit was '" &
              &                 // TRIM(forcing_options%simulation_length_unit)//"'. Please check!"
            CALL finish(TRIM(modname),message_text)
          ENDIF

          IF (forcing_options%simulation_length_number <= nr_forcing_years) THEN
            ! In this case the start date of the netcdf output file is set to the first year of forcing data
            start_date_char = TRIM(int2string(model%config%forcing_file_start_yr))//'-01-01 00:00'
          ELSE
            ! else the start date of the netcdf output file needs to be shifted
            WRITE(message_text,*) 'The start date of the netcdf output is dated back as compared to the start date ', &
              & 'of the input forcing, in order to have the last years matching between forcing data and simulation output!'
            CALL message(TRIM(modname), message_text)
            WRITE(message_text,*) 'Input forcing data starts at:  ', &
              &                   TRIM(int2string(model%config%forcing_file_start_yr)), '-01-01 00:00'
            CALL message(TRIM(modname), message_text)
            ! calc netcdf start date
            start_date_char = TRIM(int2string(model%config%forcing_file_start_yr &
              &                 - (forcing_options%simulation_length_number - nr_forcing_years)))//'-01-01 00:00'
            ! docu and info
            CALL message(TRIM(modname), 'Simulation output is written relative to the start date: '//TRIM(start_date_char))
          ENDIF
          ! prepare file
          n_output_years = model%config%output_end_last_day_year - model%config%output_start_first_day_year + 1
          CALL prepare_nc4_output_writing(tile, iblk, n_output_years, model%config%dtime, output_interval_pool_act, &
            &                             TRIM(model%config%output_interval_flux), 'static_', start_date_char)
        ENDIF
      ELSEIF (model%config%l_transient_forcing_mode .AND. output_interval_pool_act > 0.0_wp) THEN
        IF (l_output_as_ascii_text) THEN
          CALL open_pool_output_files(TRIM(model%config%output_interval_pool_spinup), 'spinup_')
          CALL open_flux_output_files(tile, model%options, TRIM(model%config%output_interval_flux_spinup), 'spinup_')
        ELSE
          ! docu and info
          WRITE(message_text,*) 'The start date of the netcdf output is dated back such that ', &
            &   'the spinup output ends in the year before the transient simulation output starts!'
          CALL message(TRIM(modname), message_text)
          ! calc netcdf start date and prepare netcdf output files
          start_date_char = TRIM(int2string(forcing_options%transient_simulation_start_year &
            &                                - forcing_options%transient_spinup_years))//'-01-01 00:00'
          n_output_years = forcing_options%transient_spinup_years
          CALL prepare_nc4_output_writing(tile, iblk, n_output_years, model%config%dtime, output_interval_pool_act, &
            &                             TRIM(model%config%output_interval_flux_spinup), 'spinup_', start_date_char)
          ! docu and info
          CALL message(TRIM(modname), 'Simulation output is written relative to the start date: '//TRIM(start_date_char))
        ENDIF
      ENDIF
      ! create soil forcing (for quincy-jsm standalone): output-control is in mo_shared:output:write_soil_forcing
      !  static:    output is created from a static simulation, for one year after static-spinup
      !  transient: output is created from a transient simulation with fluxnet/3 sites, for the last full forcing-period at the end of transient simulation
      IF (model%config%flag_generate_soil_forcing_stat .OR. model%config%flag_generate_soil_forcing_trans) THEN
        CALL open_output_files_soil_forcing(tile, iblk, ics, ice, nc)
      ENDIF
    CASE(Q_TEST_CANOPY, Q_TEST_RADIATION)
      ntimestep = 1     ! simulation starts with one, and loop asks for 'IF (itimestep .GT. ntimestep) EXIT'
  END SELECT
  !> 4.2 time related output files
  !>
  ! output of simulation runtime estimation
  CALL open_quincy_time_output_files()

  !! 5.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Time loop calling the land model interface
  !! 
  PRINT*," "
  PRINT*," Start QUINCY model Simulation"
  PRINT*," timeloop from timestep", itimestep + 1, "to", ntimestep
  elapsed_seconds = 0._wp
  DO
    ! exits the DO loop after completing the timestep with "itimestep == ntimestep"
    IF (itimestep == ntimestep) EXIT

    !! first timestep, experiment start
    !! 
    !! is set to FALSE at the end of the first timestep
    !! 
    IF (is_first_timestep) THEN
      CALL set_is_time_experiment_start(.TRUE.)
    END IF

    !! 5.1 advance timestep of this loop
    !! init itimestep with 0, the timestep 1 means running from 0 to 1, equally to time_day
    itimestep       = itimestep + 1             ! increment DO loop counter

    !! 5.2 update time information 
    !! advance simulation-timer and check timer for new year (year counter)
    !! 
    !! for timing data of the forcing input see forcing_options\%forcing_*
    !!
    !! Note: \n 
    !!  for each timestep the model simulates what happens during the period:
    !!  time_day (end of prevous timestep)  -->  time_day  + timestep_prop_day
    !!  time_day is initialized with 0.0
    !! 
    ! seconds passed in the simulation since simulation start
    elapsed_seconds = elapsed_seconds + model%config%dtime  ! add "timestep length in seconds"
    ! day and year of the simulation
    time_day  = time_day  + timestep_prop_day               ! time in days
    IF (time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_jsb_math_constants
      time_year = time_year + 1   

      IF (forcing_options%transient_spinup_years <= 20 .AND. forcing_options%simulation_length_number <= 20) THEN
        PRINT *,"New sim year ", time_year, " using forcing from ", forcing_options%forcing_year
      END IF
    END IF
    ! calculate the seconds passed on this day (subtract whole days and convert day to seconds)
    time_day_seconds(:) = (time_day - REAL(INT(time_day),wp)) * seconds_per_day

    ! update model options
    model%options%current_datetime%time_day = time_day
    model%options%current_datetime%time_year = time_year

    ! calculate and report the real time that has passed over this simulation and 
    !   estimate remaining time
    CALL calc_and_report_estimated_remaining_simulation_time(first_timestep, itimestep, ntimestep)

    !! 5.3 control for transient forcing_mode
    !! 
    !! 
    IF (model%config%l_transient_forcing_mode) THEN
      model%config%l_transient_spinup_rewind_forcing = .FALSE.
      ! control spinup
      IF (model%config%l_transient_spinup) THEN
        transient_spinup_itimestep_period = transient_spinup_itimestep_period + 1
        ! if one spinup period (i.e. repetitive part of the spinup) has finished
        IF (transient_spinup_itimestep_period .GT. forcing_options%transient_spinup_timesteps_repeat) THEN
          transient_spinup_itimestep_period              = 1
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
        END IF
        ! if spinup is completed and the simulation should start
        IF (itimestep .EQ. 1) THEN
          model%config%l_transient_spinup_rewind_forcing = .TRUE.
          model%config%l_transient_spinup                = .FALSE.           ! do not enter the "control spinup" IF statement anymore
          time_day                                       = timestep_prop_day ! equals 1.0 at the last simulated timestep of the day
          time_year                                      = 1                 ! the 1st year of a simulation == 1 (not zero)
          ! update model options
          model%options%current_datetime%time_day = time_day
          model%options%current_datetime%time_year = time_year
          ! close spinup output files if output was created (i.e. spinup output interval is not 'none')
          IF (output_interval_pool_act > 0.0_wp) THEN
            IF (l_output_as_ascii_text) THEN
              CALL close_output_files
            ELSE
              CALL close_nc4_output_files(model%id)
            ENDIF
          ENDIF
          ! change output interval
          output_interval_pool_act = output_interval_pool
          ! open new output files
          IF (l_output_as_ascii_text) THEN
            CALL open_pool_output_files(TRIM(model%config%output_interval_pool), '')
            CALL open_flux_output_files(tile, model%options, TRIM(model%config%output_interval_flux), '')
          ELSE
            start_date_char = TRIM(int2string(forcing_options%transient_simulation_start_year))//'-01-01 00:00'

            IF(model%config%fluxnet_type_transient_timestep_output .AND. model%config%output_start_first_day_year /= 1) THEN
              PRINT *, "Warning: with fluxnet type output following namelist settings are ignored: " &
                & // "output_start_first_day_year and output_end_last_day_year."
            ENDIF

            IF (model%config%fluxnet_type_transient_timestep_output) THEN
              n_output_years = forcing_options%simulation_length_number &
                &              - (model%config%fluxnet_static_forc_last_yr - model%config%fluxnet_static_forc_start_yr + 1)
            ELSE
              n_output_years = model%config%output_end_last_day_year - model%config%output_start_first_day_year + 1
              IF (forcing_options%simulation_length_number < model%config%output_start_first_day_year) THEN
                WRITE(message_text,*) 'With a simuation length of ('//TRIM(int2string(forcing_options%simulation_length_number)) &
                  &                 // ') no output would be created because the first output year would not be reached: ('      &
                  &                 // TRIM(int2string(model%config%output_start_first_day_year))//'). Please check!'
                CALL finish(TRIM(modname), message_text)
              ELSE IF (forcing_options%simulation_length_number < n_output_years) THEN
                n_output_years = forcing_options%simulation_length_number
              END IF
            END IF
            CALL prepare_nc4_output_writing(tile, iblk, n_output_years, model%config%dtime, output_interval_pool_act, &
              &                             TRIM(model%config%output_interval_flux), 'transient_', start_date_char)
          ENDIF
        END IF
      END IF
    END IF

    !! 5.4 set orbital parameters and update time-variant forcing
    !!
    !! this function also reads the "forcing year & doy & hour" into the forcing_options object \n
    !! l_forcing is set in init_slm() depending on the specific SLM model
    !! 
    IF (l_forcing) THEN
      CALL update_forcing( &
        tile      , & ! necessary for soil: veg%litterfall vegsoil%moisture
        ics       , &
        ice       , &
        nc        , &
        iblk      , &
        itimestep , &
        time_day  , & 
        timestep_prop_day, &
        time_day_seconds, &
        model%config%dtime, &
        !! forcing (read from file and passed to a2l memory via the land_model_interface)
        !! in update_forcing() defined as: REAL(wp), OPTIONAL, INTENT(out) :: VAR(:)  ! assumed shape
        t_air                   = t_air(:)              , &
        q_air                   = q_air(:)              , &
        rain                    = rain(:)               , &
        snow                    = snow(:)               , &
        wind_air                = wind_air(:)           , &
        wind_10m                = wind_10m(:)           , &
        lw_srf_down             = lw_srf_down(:)        , &
        swvis_srf_down          = swvis_srf_down(:)     , &
        swnir_srf_down          = swnir_srf_down(:)     , &
        swpar_srf_down          = swpar_srf_down(:)     , &
        fract_par_diffuse       = fract_par_diffuse(:)  , &
        press_srf               = press_srf(:)          , &
        cos_zenith_angle        = cos_zenith_angle(:)   , &
        co2_mixing_ratio        = co2_mixing_ratio(:)     , &
        co2_mixing_ratio_C13    = co2_mixing_ratio_C13(:) , &
        co2_mixing_ratio_C14    = co2_mixing_ratio_C14(:) , &
        nhx_deposition          = nhx_deposition(:)     , &
        noy_deposition          = noy_deposition(:)     , &
        nhx_n15_deposition      = nhx_n15_deposition(:) , &
        noy_n15_deposition      = noy_n15_deposition(:) , &
        p_deposition            = p_deposition(:)         &
        &  )
    END IF

    !! 5.5 output forcing & timestep info (not for QSOIL)
    !!   first timestep of the simulation
    !!   first timestep of transient simulation
    ! the QSOIL model mode has no time information included in the forcing files
    SELECT CASE(model%config%qmodel_id)
    !! land, plant, canopy, test_canopy, test_radiation
    CASE(QLAND, QPLANT, QCANOPY)
      IF (itimestep == first_timestep) THEN
        PRINT *," The simulation starts with the timestep", itimestep, " and forcing from [YEAR DOY HOUR]:"
        PRINT *,forcing_options%forcing_year, forcing_options%forcing_doy, forcing_options%forcing_hour

        ! Assert: the first year in the forcing file should equal the forcing file start year as specified in the namelist
        IF(forcing_options%forcing_year /= model%config%forcing_file_start_yr) THEN
          WRITE(message_text,*) 'First year read from the forcing file ('//TRIM(int2string(forcing_options%forcing_year))      &
            &                   // ') not equal to forcing start year as specified in the namelist or else set as DEFAULT ('   &
            &                   // TRIM(int2string(model%config%forcing_file_start_yr))//'). Please check!'
          CALL finish(TRIM(modname),message_text)
        ENDIF

        IF (model%config%l_transient_forcing_mode) THEN
          ! Assert: the first year in the forcing file should equal the transient spin-up start year
          IF(forcing_options%forcing_year /= forcing_options%transient_spinup_start_year) THEN
            WRITE(message_text,*) 'First year read from the forcing file ('//TRIM(int2string(forcing_options%forcing_year))    &
              &                   // ') not equal to start year of transient spinup simulation as specified in the namelist (' &
              &                   // TRIM(int2string(forcing_options%transient_spinup_start_year))//'). Please check!'
            CALL finish(TRIM(modname),message_text)
          ENDIF
        ENDIF

        IF (dsl4jsb_Config(Q_SYL_)%flag_stand_harvest) THEN
          IF (.NOT. dsl4jsb_Config(Q_SYL_)%harvest_active_in_qs) THEN
            PRINT *, "Warning: flag_stand_harvest is TRUE but harvest_active_in_qs is FALSE, thus there will not be any harvest"
          END IF
          ! When simulating with stand replacing harvest: check if the climate forcing is cycled and if so warn that therefore multiple
          ! harvest events are possible
          IF (forcing_options%simulation_length_number > nr_forcing_years) THEN
            PRINT *, "Warning: Beware, since the climate forcing is cycled, it is possible that the stand replacing harvest " &
              & // "is executed several times."
          END IF
        END IF
      END IF

      IF (model%config%l_transient_forcing_mode) THEN
        IF (itimestep == 1) THEN
          PRINT *," Spinup completed. Starting transient simulation with timestep 1 and forcing [YEAR DOY HOUR]:"
          PRINT *,forcing_options%forcing_year, forcing_options%forcing_doy, forcing_options%forcing_hour
          elapsed_seconds = 0._wp

          IF(forcing_options%forcing_year /= forcing_options%transient_simulation_start_year) THEN
            WRITE(message_text,*) 'First year read from the forcing file ('//TRIM(int2string(forcing_options%forcing_year)) &
              &                   // ') not equal to start year of transient simulation as specified in the namelist ('     &
              &                   // TRIM(int2string(forcing_options%transient_simulation_start_year))//'). Please check!'
            CALL finish(TRIM(modname),message_text)
          ENDIF
        ENDIF
      ENDIF
    END SELECT  ! select QUINCY model mode

    !! 5.6 control of timestep output (last full forcing period) for a transient simulation with a fluxnet-type site-set
    !!  timestep-output for the "last full forcing period" after spinup has finished, 
    !!  which is identical to the period static forcing is available 
    !! 
    !! this does work only with met-forcing, because the soil-forcing has no year/doy/hour information 
    !!  (hence, this functionality is disabled in the prepare_sim_quincy.sh see section that evaluates the option: fluxnetTypeTransientTimestepOutput)
    !! 
    ! if timestep output for the last full forcing period (of transient simulation) is enabled for a fluxnet-type site-set
    IF (model%config%fluxnet_type_transient_timestep_output .AND. .NOT. model%config%l_transient_spinup) THEN
      ! if first timestep of "start year of static forcing" of this site
      ! split into the two test "year" & "doy & hour" to check the latter only for one year (reduce CPU time)
      IF (forcing_options%forcing_year == model%config%fluxnet_static_forc_start_yr) THEN
        IF (forcing_options%forcing_doy == 1 .AND. forcing_options%forcing_hour < eps4) THEN
          ! update state: now writing timestep-output data
          flux_trans_timestep_out_state = 2 
          ! close pool & flux output files (opened above)
          IF (l_output_as_ascii_text) THEN
            CALL close_output_files 
          ELSE
            CALL close_nc4_output_files(model%id)
          ENDIF
          elapsed_seconds                     = 0._wp
          model%config%output_interval_pool   = 'timestep' 
          model%config%output_interval_flux   = 'timestep' 
          output_interval_pool                = calc_output_time_interval_pool(model%config%output_interval_pool, &
            &                                   model%config%dtime) ! character to int
          output_interval_pool_act            = output_interval_pool
          IF (l_output_as_ascii_text) THEN
            CALL open_pool_output_files(TRIM(model%config%output_interval_pool), '')
            CALL open_flux_output_files(tile, model%options, TRIM(model%config%output_interval_flux), '')
          ELSE
            start_date_char = int2string(model%config%fluxnet_static_forc_start_yr)//'-01-01 00:00'
            n_output_years = model%config%fluxnet_static_forc_last_yr - model%config%fluxnet_static_forc_start_yr + 1

            IF (TRIM(model%config%sim_siteset_specific) == "eucface_trans") THEN
              CALL prepare_nc4_output_writing(tile, iblk, n_output_years, model%config%dtime, output_interval_pool_act, &
                &                             TRIM(model%config%output_interval_flux), 'eucface_', start_date_char)
            ELSE
              CALL prepare_nc4_output_writing(tile, iblk, n_output_years, model%config%dtime, output_interval_pool_act, &
                &                             TRIM(model%config%output_interval_flux), 'fluxnetdata_', start_date_char)
            ENDIF
          ENDIF

          ! print info to default output
          IF (TRIM(model%config%sim_siteset_specific) == "eucface_trans") THEN
            PRINT*," EucFACE transient simulation: timestep-output starts with forcing [YEAR DOY HOUR]:"
          ELSE
            PRINT*," Fluxnet-type siteset transient simulation: timestep-output starts with forcing [YEAR DOY HOUR]:"
          ENDIF
          PRINT*,forcing_options%forcing_year, forcing_options%forcing_doy, forcing_options%forcing_hour
        ENDIF
      ! stop timestep-output at the 1st timestep after the last year of static forcing, i.e., the "last full forcing period" is over
      ! actually, the QUINCY model should stop before this happend, as the simulation end is reached with the prior timestep
      ELSE IF (forcing_options%forcing_year == model%config%fluxnet_static_forc_last_yr + 1) THEN
        IF (forcing_options%forcing_doy == 1 .AND. forcing_options%forcing_hour < eps4) THEN
          flux_trans_timestep_out_state = 3                       ! status update: writing output data finished
          ! print info to default output
          IF (TRIM(model%config%sim_siteset_specific) == "eucface_trans") THEN
            PRINT*," EucFACE transient simulation: timestep-output stopped before simulating the timestep:", itimestep
          ELSE
            PRINT*," Fluxnet diagnostics simulation: timestep-output stopped before simulating the timestep:", itimestep
          ENDIF
          PRINT*," Associated to this timestep is the forcing [YEAR DOY HOUR]:"
          PRINT*,forcing_options%forcing_year, forcing_options%forcing_doy, forcing_options%forcing_hour
        ENDIF
      ENDIF
    ENDIF

    !! 5.7 output control
    !!  does consider:  static & transient
    !!                  spinup & actual simulation
    !!  
    l_output = .FALSE.      ! is generally set to FALSE tested below if TRUE at this timestep
    !! static
    IF (model%config%l_static_forcing_mode) THEN
      ! output period
      IF (time_year >= model%config%output_start_first_day_year .AND. time_year <= model%config%output_end_last_day_year) THEN
        l_output = .TRUE.
      ENDIF        
    !! transient
    ELSEIF (model%config%l_transient_forcing_mode) THEN
      ! spinup
      IF (model%config%l_transient_spinup) THEN
        ! if output may be generated during spinup (i.e. output interval is not 'none')
        IF (output_interval_pool_act > 0.0_wp) THEN
          l_output = .TRUE.
        ENDIF
      ! actual simulation after spinup
      ELSE
        ! differ between a simple transient simulation and one with timestep output for the last full forcing period of the transient simulation
        IF (model%config%fluxnet_type_transient_timestep_output) THEN
          ! if active and not finished yet (finished == after the last full forcing period no further output is created)
          IF (flux_trans_timestep_out_state < 3.0_wp) THEN
            l_output = .TRUE.
          ENDIF
        ! a simple transient simulation, create output only for a given start/end output period
        ELSE
          ! output period
          IF (time_year >= model%config%output_start_first_day_year .AND. time_year <= model%config%output_end_last_day_year) THEN
            l_output = .TRUE.
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    !!  match output interval
    !! 
    l_output_interval_match = .FALSE.      ! is generally set to FALSE tested below if TRUE at this timestep
    IF (l_output) THEN
      ! if interval is smaller than 'daily' then set to TRUE always, because then it can only be timestep
      IF (output_interval_pool_act < 1.0_wp) THEN
        l_output_interval_match = .TRUE.
      ! transient spinup
      ! during spinup the timesteps are a negative number, which is compensated by adding the spinup length to time_day
      ELSEIF (model%config%l_transient_spinup) THEN
        IF (MODULO(time_day + transient_spinup_timesteps, output_interval_pool_act) < eps4 .OR. &
          ABS(output_interval_pool_act - MODULO(time_day + transient_spinup_timesteps, output_interval_pool_act)) < eps4) THEN
          l_output_interval_match = .TRUE.
        ENDIF
      ! static simulation or transient simulation after spinup
      ELSE
        IF (MODULO(time_day, output_interval_pool_act) < eps4 .OR. &
          ABS(output_interval_pool_act - MODULO(time_day, output_interval_pool_act)) < eps4) THEN
          l_output_interval_match = .TRUE.
        ENDIF
      ENDIF
    ENDIF

    ! Determine the size of the nc buffer for the current netcdf output files
    IF (output_interval_pool_act > 0.0_wp) THEN
      ! the buffer size depends on the number of entries written into one file
      buffer_nc_size = NINT(one_year/output_interval_pool_act)
    ELSE
      ! 1 if there is no output
      buffer_nc_size = 1
    ENDIF

    !! 5.8 land model interface
    !! 
    CALL land_model_interface(tile, &
                              buffer_nc_size, &
                              buffer_nc_index, &
                              l_forcing, &
                              l_output, &
                              l_output_as_ascii_text, &
                              l_output_interval_match, &
                              flux_trans_timestep_out_state, &
                              output_interval_pool_act, &      ! output interval in days, used for flux output: output_interval_flux
                              elapsed_seconds, &
                              time_day, &
                              time_year &
                              &  ) ! land_model_interface()

    !! end first timestep, experiment start
    !! 
    !! set to FALSE
    !! 
    IF (is_first_timestep) THEN
      is_first_timestep = .FALSE.
      CALL set_is_time_experiment_start(.FALSE.)
    END IF

  END DO  ! end time loop 

  !! 6.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Close output, forcing, log files, and time related output files
  !! 
  PRINT*," Simulation stopped after the timestep", itimestep, "has completed."
  PRINT*," Forcing data used for this last timestep [YEAR DOY HOUR]:"
  PRINT*,forcing_options%forcing_year, forcing_options%forcing_doy, forcing_options%forcing_hour
  IF (l_output_as_ascii_text) THEN
    CALL close_output_files
  ELSE
    CALL close_nc4_output_files(model%id)
  ENDIF
  IF (l_forcing) CALL close_forcing_input_file(model%config%qmodel_id, tile) ! only if forcing is enabled at all (land plant canopy soil)
  IF (model%config%flag_generate_soil_forcing_stat .OR. model%config%flag_generate_soil_forcing_trans) THEN
    CALL close_output_files_soil_forcing
  END IF
  CALL close_quincy_time_output_files()
  PRINT*," All output files were closed."
  ! close log & error file
  CLOSE(nlog)
  CLOSE(nerr)

  !! 7.0 -----------------------------------------------------------------------------------------------
  !! 
  !! DEALLOCATE variables
  DEALLOCATE(time_day_seconds)
  IF (.NOT. l_output_as_ascii_text) THEN
    CALL deallocate_quincy_output(model%id)
  ENDIF
  PRINT*," End QUINCY model"
  PRINT*,"#########################"

END PROGRAM quincy_qs_model_driver
