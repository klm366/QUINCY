!>
!! @brief Contains methods for JSBACH initializations.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-22); starting with custom SLM init, moving to default JSBACH4 init later
!! Second version by Jan Engel (20161213): moving veg process init routines to mo_veg_init
!! Third version by Steffen Richter & JE (June 2018): split model init into init_model_standalone & init_model_common according to JSBACH4 code
!!
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_base,         ONLY: jsbach
  USE mo_jsb_model_class,  ONLY: t_jsb_model

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_model

  PUBLIC :: time_day, time_year, ntimestep, itimestep, timestep_prop_day
  PUBLIC :: fluxnet_output_start, fluxnet_output_stop, fluxnet_output_status

  ! variables used in the land model driver
  REAL(wp) :: time_day                     !< displays the 'real time' in days; time_day == 1.0 after one day has passed
  INTEGER  :: time_year                    !< displays the 'real time' in years; the 1st year of a simulation == 1 (not zero)
  INTEGER  :: ntimestep                    !< number of timesteps the model runs
  INTEGER  :: itimestep                    !< current timestep i, looping over timesteps ntimestep
  REAL(wp) :: timestep_prop_day            !< proportion of one day that is simulated with one timestep (dtime/one_day)
  ! special case: fluxnet output
  INTEGER  :: fluxnet_output_length        !< number of timesteps data may be written out to file, only for fluxnet sites
  INTEGER  :: fluxnet_output_start         !< write output at each time step starting with this timestep, only for fluxnet sites
  INTEGER  :: fluxnet_output_stop          !< write output at each time step stopping with this timestep, only for fluxnet sites
  INTEGER  :: fluxnet_output_status        !< a simple status var to inform other processes about it
                                           !! 0: not activiated at all [default], 1: active, 2: currently writing, 3: finished

  INTERFACE init_model
    MODULE PROCEDURE init_model_standalone
  END INTERFACE init_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> standalone model init
  !! 
  !! creates model and grid \n 
  !! runs init_model_common() for all other init and config tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_model_standalone(nproca, nprocb, npedim, model_id, l_forcing, nsoil_sb, iblk, ics, ice, nc, model_name)

    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                    INTENT(in) :: nproca, nprocb, &
                                              npedim, model_id
    LOGICAL,INTENT(out)                    :: l_forcing            !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
    INTEGER,INTENT(out)                    :: nsoil_sb             !< number of soil layers used in SB_ - equal to nsoil_water
                                                                   !!  necessary for open_pool_output_files() & init_output_writing_pool()
    INTEGER,INTENT(in)                     :: iblk, ics, ice, nc   !< block number and chunk size
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),  POINTER            :: model
    TYPE(t_jsb_grid),   POINTER            :: grid
    INTEGER                                :: nlon, nlat
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_standalone'


    !> 1.0 model
    !! 
    ! create one model
    model => jsbach%models(model_id)%m
    ! check if this model got the correct namelist (if model name is provided)
    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    !> 2.0 grid
    !! 
    ! create grid, currently for only one point
    grid => new_grid(1.0_wp,1.0_wp)
    ! register grid 
    CALL register_grid(grid)
    model%grid_id = grid%id

    !> 3.0 run model init and config
    !! 
    CALL init_model_common(model, l_forcing, nsoil_sb, iblk, ics, ice, nc, model_id)
    ! delete local pointer 
    NULLIFY(model)
  END SUBROUTINE init_model_standalone


  !-----------------------------------------------------------------------------------------------------
  !> init_model_common
  !! 
  !! subtask of init_model_standalone() \n 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_model_common(model, l_forcing, nsoil_sb, iblk, ics, ice, nc, model_id)

    ! base
    USE mo_jsb_process_factory,   ONLY: max_no_of_processes, Create_process
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_lct_class,         ONLY: Contains_lct
    USE mo_jsb_model_usecases,    ONLY: init_usecase
    USE mo_jsb4_forcing,          ONLY: Init_forcing, forcing_options
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid
    USE mo_jsb_grid,              ONLY: Get_Grid
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: Read_lctlib
    USE mo_jsb_config_class,      ONLY: new_model_config
    USE mo_base_config_class,     ONLY: base_config
    ! shared 
    USE mo_util,                  ONLY: one_of
    USE mo_shared_set_parameters, ONLY: parameters_from_file    !< TYPE t_parameters_from_file
    ! process memory
    USE mo_sb_memory_class,       ONLY: Init_sb_pools_elements_a
    USE mo_veg_memory_class,      ONLY: Init_veg_pools_elements_a
    ! process init
    USE mo_veg_init,              ONLY: veg_init
    USE mo_sse_init,              ONLY: sse_init
    USE mo_sb_init,               ONLY: sb_init
    USE mo_rad_init,              ONLY: rad_init
    USE mo_assimi_init,           ONLY: assimi_init
    USE mo_pheno_init,            ONLY: pheno_init
    ! process constants
    USE mo_base_constants         ! e.g. dtime
    USE mo_shared_constants,      ONLY: init_shared_constants
    USE mo_atmland_constants,     ONLY: init_atmland_constants
    USE mo_rad_constants,         ONLY: init_rad_constants
    USE mo_assimi_constants,      ONLY: init_assimi_constants
    USE mo_pheno_constants,       ONLY: init_pheno_constants
    USE mo_veg_constants,         ONLY: init_veg_constants
    USE mo_sse_constants,         ONLY: init_sse_constants
    USE mo_sb_constants,          ONLY: init_sb_constants
    ! process config
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(VEG_)        
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
    LOGICAL,INTENT(out)                       :: l_forcing            !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
    INTEGER,INTENT(out)                       :: nsoil_sb             !< number of soil layers used in SB_ - equal to nsoil_water
                                                                      !!  necessary for open_pool_output_files() & init_output_writing_pool()
    INTEGER,INTENT(in)                        :: iblk, ics, ice, nc   !< *slmtmp iBlock, iChunkStart, iChunkEnd, number of chunks
    INTEGER,INTENT(in)                        :: model_id
    ! ---------------------------
    ! 0.2 Local
    CLASS(t_jsb_tile_abstract), POINTER       :: tile, child          !< data structure for dynamic vegetation and soil variables
    INTEGER                                   :: iproc, ilct, jlct, &
                                                 i, no_children
    REAL(wp),                   ALLOCATABLE   :: fract_tmp(:,:)
    INTEGER                                   :: no_omp_thread
    TYPE(t_jsb_grid),           POINTER       :: grid
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_common'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SSE_)
    ! the according Get_config is called below


    !> 0.9 init local variables
    !! 
    no_omp_thread = 1


    !> 1.0 configuration
    !!
    !> 1.1 model config
    !!
    model%config => new_model_config(model%namelist_filename)
    !> 1.2 process config
    !!
    ALLOCATE(model%processes(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
    END DO
    base_config%namelist_filename = model%namelist_filename
    CALL base_config%Init                            ! includes base_config%slm_model_name => land plant soil canopy test_canopy
    CALL model%Configure_processes()
    
    
    !> 2.0 forcing config
    !! 
    CALL Init_forcing(model_id)


    !> 3.0 model init 
    !! 
    ! get process configs needed for model_init
    dsl4jsb_Get_config(SSE_)
    ! Set model options
    model%options%dtime        = base_config%dtime  ! length of one time step
    model%options%steplen      = base_config%dtime  ! ..
    model%options%iblk         = iblk               ! Number of current chunk that consists of multiple grid boxes
    model%options%ics          = ics                ! Index of chunk start
    model%options%ice          = ice                ! Index of chunk end
    model%options%nc           = nc                 ! Length of chunk  (end-number - start-number + 1)
    model%options%nsoil_energy = dsl4jsb_Config(SSE_)%nsoil_energy  ! number of soil layers for soil energy calculations
    model%options%nsoil_water  = dsl4jsb_Config(SSE_)%nsoil_water   ! number of soil layers for soil moisture calculations
    model%options%nsoil_sb     = dsl4jsb_Config(SSE_)%nsoil_water   ! number of soil layers used in SB_ - equals nsoil_water
    nsoil_sb                   = dsl4jsb_Config(SSE_)%nsoil_water   ! SUBROUTINE output for land_model_driver -> open_pool_output_files()


    !> 4.0 Grid init
    !!
    grid => Get_grid(model%grid_id)
    grid%lat = base_config%lat
    grid%lon = base_config%lon


    !> 5.0 Constants init
    !!
    CALL init_base_constants          ! may be 1st of the 'init constants'
    CALL init_shared_constants        ! may be 2nd of the 'init constants'
    CALL init_rad_constants           ! may be 3rd of the 'init constants'
    CALL init_atmland_constants       ! uses hourday                  defined in mo_base_constants
    CALL init_assimi_constants        ! uses molar_mass_N, Dwv, Dco2  defined in mo_shared_constants
    CALL init_veg_constants           ! uses Tzero                    defined in mo_shared_constants
    CALL init_sse_constants
    CALL init_sb_constants
    CALL init_pheno_constants

    !> 5.1 Constants modification for parameter sensitivity analysis
    !!
    IF (base_config%set_parameter_values_from_file) CALL parameters_from_file%set_param


    !> 6.0 Memory init -> pool and element structures in veg_ and sb_ memory with ALLOCATABLE variables with operator overload
    !!
    CALL Init_sb_pools_elements_a(ics, ice, nc, iblk, nsoil_sb)
    CALL Init_veg_pools_elements_a(ics, ice, nc, iblk)


    !> 7.0 lctlib init and read
    !!
    !! Important! Read_lctlib() needs config and constants init to be finished !
    !!
    model%lctlib => Read_lctlib('lctlib_slm_nlct8.def')
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_slm_nlct8.def')


    !> 8.0 Init tiles
    !!
    ! get process configs needed for tile init
    dsl4jsb_Get_config(VEG_)
    ! Shortcut model configuration through usecase
    CALL init_usecase(model)
    ! At this point, the model configuration and tile structure is set up.
    ! We continue with initializing the memory structure and different model processes.
    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    !> 8.1 loop over tiles
    !! 
    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      !> 8.2 Initialize tile memory and read tile fractions
      !! @todo: once we have land use change or dynveg, l_fixed_fractions needs to be adapted <- JSBACH4
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name), l_fixed_fractions=.TRUE.)

      !> 8.3 tile init
      !! @todo model_init -> init lcts of the tile via PROCEDURE Construct_jsb_lct  *slmtemp  ?
      tile%owner_model_id = model_id                      ! *slmtemp
      ! lcts
      tile%lcts%id        = 2                             ! VEG_TYPE     2 == Vegetated
      tile%lcts%name      = "veg"                         ! Short name for lct
      tile%lcts%long_name = "slm vegetation"              ! Long name for lct
      tile%lcts%lib_id    = dsl4jsb_Config(VEG_)%pft_id   ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)

      !! this option is outdated; forcing is always used, but not for test_canopy, as the below CASE indicates :)
      !!
      !! enable/ disable forcing
      !! note: the soil model gets a different forcing than land plant canopy, see mo_atmland_forcing
      SELECT CASE(TRIM(base_config%slm_model_name))
        CASE('land','plant','canopy','soil')
          l_forcing = .TRUE.
        CASE DEFAULT
          l_forcing = .FALSE.
      END SELECT

      !> 8.4 run process init
      !! 
      !! one_of() returns -1 if NONE of the SLM model names are found, but the 1st found SLM model name get the number 1 (no enum just loop)
      !! 
      !! a certain sequence of proc_init (i.e. the below sequence) is needed to correctly setup the processes, because of some 
      !!  dependencies across two or three processes
      ! test model name
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','soil','canopy','test_canopy' ]) == -1) THEN
        CALL finish(routine,'Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
      END IF
      ! pheno_init
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','canopy' ]) > 0) THEN
        CALL pheno_init(tile, model%options)
      END IF
      ! assimi_init
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','canopy','test_canopy' ]) > 0) THEN
        CALL assimi_init(tile, model%options)
      END IF
      ! rad_init
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','canopy','test_canopy' ]) > 0) THEN
        CALL rad_init(tile, model%options)
      END IF
      ! veg_init
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','soil','canopy','test_canopy' ]) > 0) THEN
        CALL veg_init(tile, model%options)
      END IF
      ! sse_init
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','soil','canopy' ]) > 0) THEN
        CALL sse_init(tile, model%options)
      END IF
      ! sb_init (needs sse_init to be finished)
      IF (one_of(base_config%slm_model_name, [character(len=20) :: 'land','plant','soil' ]) > 0) THEN
        CALL sb_init(tile, model%options)
      END IF
      ! next tile
      tile => model%Goto_next_tile()
    END DO

    !> 8.5 Iterate over tiles again and set lct fractions (after ALL tile fractions have been set above)
    !! 
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')


    !> 9.0 Time init for e.g. the time loop in the land_model_driver
    !!
    time_day          = 0.0_wp                      !! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                           !! the 1st year of a simulation == 1 (not zero)
    dtime             = base_config%dtime           !! time step length in seconds, default: 1800.0 == 30 minutes, mo_base_constants
    ntimestep         = forcing_options%Timesteps() !! number of time steps (considers the simulation_mode static & transient)
    timestep_prop_day = dtime/one_day               !! proportion of one day that is simulated with one timestep


    !> 10.0 Initialise: Fluxnet specific simulation output
    !!
    !! the last full "forcing period" is written into extra output files with half-hourly write-out of values
    !!
    fluxnet_output_status = 0
    IF (base_config%fluxnet_output) THEN
      fluxnet_output_status = 1
      ! check whether simulation length covers at least one complete fluxnet forcing period (i.e. fluxnet_output_length)
      fluxnet_output_length = NINT(base_config%fluxnet_output_length * one_year * one_day / base_config%dtime)
      IF (fluxnet_output_length > ntimestep) THEN
        CALL finish('land_model_driver-init_routine','Simulation is too short for fluxnet output of one full forcing period')
      END IF

      ! calc timestep when fluxnet output starts (get the before last full output period that starts within ntimestep)
      fluxnet_output_start  = fluxnet_output_length * FLOOR( (REAL(ntimestep) / REAL(fluxnet_output_length)) -1.0_wp ) 
      ! calc timestep when fluxnet output stops
      fluxnet_output_stop   = fluxnet_output_start + fluxnet_output_length

      ! check whether general output start and stop covers the 'Fluxnet specific simulation output'
      IF ( (base_config%output_end_last_day_year * one_year * one_day / base_config%dtime) < fluxnet_output_stop) THEN
        CALL finish('land_model_driver-init_routine','Simulation output does not cover the last full forcing period')
      END IF
      IF ( (((base_config%output_start_first_day_year - 1) * one_year * one_day + 1) / base_config%dtime) > fluxnet_output_start) THEN
        CALL finish('land_model_driver-init_routine','Simulation output does not cover the last full forcing period')
      END IF
    END IF
  END SUBROUTINE init_model_common

#endif
END MODULE mo_jsb_model_init
