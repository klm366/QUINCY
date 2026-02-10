!> Methods for JSBACH initializations
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
!>#### Initialize JSBACH instance with grid, task queue, tile structure, initial/boundary conditions etc.
!>
MODULE mo_jsb_model_init

#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_grid, quincy_init

  PUBLIC :: time_day, time_year, ntimestep, itimestep, timestep_prop_day
  PUBLIC :: flux_trans_timestep_out_state

  ! variables used in the land model driver
  REAL(wp) :: time_day                      !< displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                     !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep                     !< number of timesteps the model runs
  INTEGER  :: itimestep                     !< current timestep i, looping over timesteps ntimestep
  REAL(wp) :: timestep_prop_day             !< proportion of one day that is simulated with one timestep (dtime/one_day)
  INTEGER  :: flux_trans_timestep_out_state !< state of the fluxnet diagnostics simulation [0,3], see the below init for details

  INTERFACE jsbach_setup_grid
    MODULE PROCEDURE setup_grid_standalone
    ! MODULE PROCEDURE setup_grid_patch
    ! MODULE PROCEDURE setup_grid_decomposition
  END INTERFACE jsbach_setup_grid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> set up the horizontal grid for the standalone model
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE setup_grid_standalone(model_id, grid_name, nproma, npromz, nblks)

    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,              ONLY: Register_grid

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                    INTENT(in) :: model_id
    CHARACTER(len=*),           INTENT(in) :: grid_name
    INTEGER,                    INTENT(in) :: nproma
    INTEGER,                    INTENT(in) :: npromz
    INTEGER,                    INTENT(in) :: nblks
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),  POINTER   :: model
    TYPE(t_jsb_grid),   POINTER   :: grid
    CHARACTER(len=*), PARAMETER   :: routine = modname//':setup_grid_standalone'
    ! ---------------------------
    ! 0.5 Get Memory
    model => Get_model(model_id)


    !< 0.9
    !! 
    CALL message(TRIM(routine), 'Standalone model currently only works with ECHAM')


    !< 1.0 create horizontal grid
    !! 
    grid        => new_grid(model%config%lon, model%config%lat) ! calls via an interface: FUNCTION new_grid_from_lonlat(lon, lat) RESULT(grid)
    grid%id     = 1                     ! *quincytmp -> must be identical with the position of this grid in mo_jsb_grid:jsbach_grids(:)
    grid%name   = grid_name             ! *quincytmp
    grid%nproma = nproma                ! *quincytmp
    grid%npromz = npromz                ! *quincytmp
    grid%nblks  = nblks                 ! *quincytmp
    CALL register_grid(grid)            ! add to mo_jsb_grid:jsbach_grids(:), actuall to jsbach_grids(1)
    model%grid_id = grid%id             ! tell the model the id of the hgrid, so that one can get the correct hgrid per model

    ALLOCATE(grid%lsm(grid%nproma, grid%nblks))
    ALLOCATE(grid%lsf(grid%nproma, grid%nblks))
    ALLOCATE(grid%area(grid%nproma, grid%nblks))
    grid%lsm(:,:)  = .TRUE.
    grid%lsf(:,:)  = 1.0_wp
    grid%area(:,:) = 1.0_wp

    NULLIFY(model)    ! NULLIFY pointer to model

  END SUBROUTINE setup_grid_standalone


  !-----------------------------------------------------------------------------------------------------
  !> init quincy model
  !! 
  !! and init of: forcing, SB_ & VEG_ pools, tiles, time, fluxnet diagnostics simulation state
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE quincy_init(model_id, iblk, nproma, nblks, ics, ice, nc, l_forcing)
    ! base
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_quincy_model_config,   ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_lct_class,         ONLY: Contains_lct
    USE mo_jsb_model_usecases,    ONLY: init_usecase
    USE mo_qs_forcing,            ONLY: Init_forcing_qs, forcing_options
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_jsb_process_class,     ONLY: VEG_, SSE_, Q_SYL_, Q_AGR_
    USE mo_jsb_math_constants,    ONLY: one_year, one_day
    ! icon_adapter
    USE mo_time_config,           ONLY: time_config
    ! shared 
    USE mo_util,                  ONLY: one_of, int2string
    ! quincy_standalone
    USE mo_jsb4_forcing,          ONLY: setup_forcing, init_forcing   ! setup icon-land forcing_options(:)
    ! process init
    USE mo_veg_init,              ONLY: veg_init
    USE mo_spq_init,              ONLY: spq_init
    USE mo_hydro_init,            ONLY: hydro_init
    USE mo_sse_init,              ONLY: sse_init
    USE mo_seb_init,              ONLY: seb_init
    USE mo_sb_init,               ONLY: sb_init
    USE mo_rad_init,              ONLY: rad_init
    USE mo_q_assimi_init,         ONLY: q_assimi_init
    USE mo_phyd_init,             ONLY: phyd_init
    USE mo_turb_init,             ONLY: turb_init
    USE mo_q_pheno_init,          ONLY: q_pheno_init
    USE mo_dfire_init,            ONLY: dfire_init
    USE mo_q_agr_init,            ONLY: q_agr_init
    USE mo_q_syl_init,            ONLY: q_syl_init
    ! process configurations (t_PROC_config)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SSE_)
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER ,INTENT(in)                       :: model_id
    INTEGER ,INTENT(in)                       :: iblk, nproma, nblks  !< block number, number of points per chunk, number of blocks   *slmtmp
    INTEGER ,INTENT(in)                       :: ics, ice, nc         !< first & last & number of points in the chunk the model currently works on  *slmtmp
    LOGICAL ,INTENT(out)                      :: l_forcing            !< IF .TRUE. THEN update_forcing()
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),          POINTER       :: model
    CLASS(t_jsb_tile_abstract), POINTER       :: tile
    TYPE(t_jsb_vgrid),          POINTER       :: soil_w           ! vertical grid HYDRO_
    TYPE(t_jsb_vgrid),          POINTER       :: soil_e           ! vertical grid SSE_
    TYPE(t_jsb_vgrid),          POINTER       :: snow_e           ! vertical grid SSE_
    INTEGER                                   :: no_omp_thread
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':quincy_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SSE_)
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    ! 
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SSE_)

    !>0.9 init local variables
    !>
    CALL message(TRIM(routine), 'starting model init')
    no_omp_thread = 1     ! no threading with OpenMPI in QUINCY


    !>1.0 QS forcing config
    !>
    CALL Init_forcing_qs(model_id)


    !>2.0 model options
    !>
    model%options%dtime        = model%config%dtime !< length of one time step
    model%options%alpha        = 1.0_wp
      !< in ICON-Land, 'alpha' is derived as alpha = MERGE(1._wp, cvdifts, jsbach_runs_standalone()))
      !< with cvdifts comming from ICON where it is set to 1.5; i.e. standalone: 1.0 // coupled: 1.5
    model%options%iblk         = iblk               !< 
    model%options%nproma       = nproma             !< 
    model%options%ics          = ics                !< 
    model%options%ice          = ice                !< 
    model%options%nc           = nc                 !< 

    ALLOCATE(model%options%current_datetime)
    ALLOCATE(time_config%tc_startdate)
    ALLOCATE(time_config%tc_stopdate)
    ALLOCATE(time_config%tc_current_date)   ! current not used, but should replace 'model%options%current_datetime' soon

    !>3.0 Init tiles
    !>
    ! At this point, the model configuration and tile structure is set up.
    ! We continue with initializing the memory structure and different model processes.
    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    !>  3.1 loop over tiles
    !>
    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      !>  3.2 tile init
      !>
      ! @TODO model_init -> init lcts of the tile via PROCEDURE Construct_jsb_lct  *slmtemp  ?
      tile%owner_model_id = model%id                      ! *slmtemp
      !! Set PFT ID, matching the lctlib (in ICON-Land this is done with 'init_usecase()')
      tile%lcts(1)%lib_id    = dsl4jsb_Config(VEG_)%pft_id   ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)

      !>  3.3 Initialize tile memory and read tile fractions
      !> mo_jsb_tile:Init_tile() -> init process memories
      !>
      ! Initialize tile process memories and read tile fractions
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name),                &
        & grid_id = model%grid_id,                                              &
        ! & in_var_groups=one_of(TRIM(tile%name), model%config%output_tiles) > 0)
        & in_var_groups=.TRUE.)

      ! setup bgc material store (if running with processes that have bgcms)
      CALL tile%Count_and_classify_bgc_materials()
      CALL tile%Collect_bgc_materials()

      !! QUINCY specific:
      !!   set lct of the tile (all is by default .FALSE.; see t_jsb_tile_abstract)
      !!   part of these values are also set with 'Add_lct()'
      tile%is_land              = .TRUE.
      tile%is_vegetation        = .TRUE.
      tile%contains_vegetation  = .TRUE.
      tile%contains_soil        = .TRUE.
      CALL message(TRIM(routine), 'tile%is_vegetation        = .TRUE.')
      CALL message(TRIM(routine), 'tile%contains_vegetation  = .TRUE.')
      CALL message(TRIM(routine), 'tile%contains_soil        = .TRUE.')
      ALLOCATE(tile%fract(nproma, nblks))
      ALLOCATE(tile%fract_old(nproma, nblks))
      tile%fract(:,:)           = 1.0_wp
      tile%fract_old(:,:)       = 1.0_wp

      !! this option is outdated; forcing is always used, but not for test_canopy, as the below CASE indicates :)
      !!
      !! enable/ disable forcing
      !! note: the soil model gets a different forcing than land plant canopy, see mo_qs_atmland_forcing
      SELECT CASE(model%config%qmodel_id)
        CASE(QLAND, QPLANT, QCANOPY, QSOIL)
          l_forcing = .TRUE.
        CASE DEFAULT
          l_forcing = .FALSE.
      END SELECT

      !>  3.4 run process init
      !>
      !> one_of() returns -1 if NONE of the QUINCY model names are found, but the 1st found SLM model name get the number 1 (no enum just loop)
      !> 
      !> a certain sequence of proc_init (i.e. the below sequence) is needed to correctly setup the processes, because of some 
      !>  dependencies across two or three processes
      !>
      ! store all bgcms into their matrices to work on the matrices within the init routines
      IF (ASSOCIATED(tile%bgcm_store)) THEN
        CALL tile%bgcm_store%Store_bgcms_in_matrices(tile%name, ics, ice, iblk)
      ENDIF
      ! dfire_init
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT ]) > 0) THEN
        CALL dfire_init(tile)
      END IF
      ! pheno_init 
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QCANOPY ]) > 0) THEN
        CALL q_pheno_init(tile)
      END IF
      ! pyhdro_init 
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT ]) > 0) THEN
        CALL phyd_init(tile)
      END IF
      ! turb_init 
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QCANOPY ]) > 0) THEN
        CALL turb_init(tile)
      END IF
      ! q_syl_init - sylviculture
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT]) > 0) THEN
        IF (model%processes(Q_SYL_)%p%config%active) THEN
          CALL q_syl_init(tile)
        END IF
      END IF
      ! q_agr - agriculture
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QCANOPY ]) > 0) THEN
        IF (model%processes(Q_AGR_)%p%config%active) THEN
          CALL q_agr_init(tile)
        END IF
      END IF
      ! assimi_init and rad_init (initialised for the same models)
      IF (one_of(model%config%qmodel_id, &
          & [ QLAND, QPLANT, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION ]) > 0) THEN
        CALL q_assimi_init(tile)
        CALL rad_init(tile)
      END IF
      ! sse_init
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QSOIL, QCANOPY ]) > 0) THEN
        CALL sse_init(tile)
      END IF
      ! seb_init
      ! initialize SEB_ before HYDRO_
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QSOIL, QCANOPY ]) > 0) THEN
        CALL seb_init(tile)
      END IF
      ! hydro_init
      ! initialize HYDRO_ before VEG_ and SPQ_
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QSOIL, QCANOPY ]) > 0) THEN
        CALL hydro_init(tile)
      END IF
      ! spq_init
      IF (.NOT. model%config%use_soil_phys_jsbach) THEN
        IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QSOIL, QCANOPY ]) > 0) THEN
          CALL spq_init(tile)
        END IF
      END IF
      ! veg_init
      IF (one_of(model%config%qmodel_id, &
          & [ QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION ]) > 0) THEN
        CALL veg_init(tile)
      END IF
      ! sb_init (needs spq_init to be finished)
      IF (one_of(model%config%qmodel_id, [ QLAND, QPLANT, QSOIL ]) > 0) THEN
        CALL sb_init(tile)
      END IF
      !!
      ! write back the matrices to the bgcms
      IF (ASSOCIATED(tile%bgcm_store)) THEN
        CALL tile%bgcm_store%Write_stored_matrices_to_bgcms(tile%name, ics, ice, iblk)
      ENDIF
      ! next tile
      tile => model%Goto_next_tile()
    END DO

    !>  3.5 Iterate over tiles again and set lct fractions (after ALL tile fractions have been set above)
    !>
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    !>4.0 Make number of vertical layers available from the options type
    !>
    ! soil water
    soil_w => Get_vgrid('soil_depth_water')
    model%options%nsoil_w = soil_w%n_levels
    ! soil energy
    soil_e => Get_vgrid('soil_depth_energy')
    model%options%nsoil_e = soil_e%n_levels
    ! snow (SSE_ vgrid)
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e => Get_vgrid('snow_depth_energy')
      model%options%nsnow_e = snow_e%n_levels
    END IF
    ! crosscheck process configs
    IF ((model%options%nsoil_e /= model%options%nsoil_w)) THEN
      CALL finish(routine, &
        & 'The number of soil layers used for water and energy currently need to be identical:'  &
        & //' nsoil_w = '//TRIM(int2string(model%options%nsoil_w))//', nsoil_e = '//TRIM(int2string(model%options%nsoil_e)))
    END IF

    !>4.0 Time init for e.g. the time loop in the quincy_qs_model_driver
    !>
    time_day          = 0.0_wp                                        !! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                                             !! the 1st year of a simulation == 1 (not zero)
    ntimestep         = forcing_options%Timesteps(model%config%dtime) !! number of time steps (considers the simulation_mode static & transient)
    timestep_prop_day = model%config%dtime / one_day                  !! proportion of one day that is simulated with one timestep

    model%options%current_datetime%time_day  = time_day
    model%options%current_datetime%time_year = time_year

    time_config%tc_startdate%time_day        = time_day                           ! just zero, i.e., start of the year
    time_config%tc_startdate%time_year       = model%config%forcing_file_start_yr ! first year in the forcing file

    !>5.0 Setup and init icon-land "quincy_standalone/mo_jsb4_forcing:forcing_options"
    !>
    ! setup 'forcing_options()' (in quincy standalone model%id is set to 1)
    CALL setup_forcing(model%id)
    ! init 'forcing_options()' with default values
    CALL init_forcing( &
      &               model%id, &
      &               "qs.namelist", &
      &               model%config%dtime)

    !>6.0 Init a simulation setup specific for a site-set
    !>

    !>  6.1 EucFACE transient simulations
    !>  using the fluxnet diagnostics-simulations infrastructure
    !>
    IF (TRIM(model%config%sim_siteset_specific) == "eucface_trans") THEN
      model%config%fluxnet_type_transient_timestep_output = .TRUE.
      PRINT*, TRIM(routine), ' Init EucFACE transient simulations.'
      PRINT*, TRIM(routine), ' Set fluxnet_type_transient_timestep_output=T; use fluxnet diagnostics infrastructure.'
    ENDIF

    !>  6.2 Fluxnet diagnostics simulation
    !>    set the state to "active"
    !>    states: \n
    !>     0: not a fluxnet diagnostics simulation
    !>     1: fluxnet diagnostics simulation; spinup or output prior to the "last full forcing period" (static output available)
    !>     2: fluxnet diagnostics simulation: timestep-output over the "last full forcing period"
    !>     3: fluxnet diagnostics simulation finished, first timestep after the "last full forcing period", no output at all
    !>
    IF (model%config%fluxnet_type_transient_timestep_output) THEN
      !! activate fluxnet diagnostics output
      flux_trans_timestep_out_state = 1
    ELSE
      flux_trans_timestep_out_state = 0
    ENDIF

  END SUBROUTINE quincy_init

#endif
END MODULE mo_jsb_model_init
