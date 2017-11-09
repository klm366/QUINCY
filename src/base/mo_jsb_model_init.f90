!>
!! @brief Contains routines for model init
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-22); starting with custom SLM init, moving to default JSBACH4 init later
!! Second version by Jan Engel (20161213): moving veg process init routines to mo_veg_init
!!
MODULE mo_jsb_model_init

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_exception,                       ONLY: finish, message, message_text
  ! driver
  USE mo_jsb4_forcing,                    ONLY: Init_forcing, forcing_options
  ! base
  USE mo_base_constants                     ! e.g. dtime, init_base_constants
  USE mo_jsb_class,                       ONLY: Get_model, jsbach
  USE mo_base_memory,                     ONLY: t_slm_tile, grid
  USE mo_jsb_model_class,                 ONLY: t_jsb_model, new_model
  USE mo_base_config_class,               ONLY: base_config
  USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element, Read_lctlib
  USE mo_jsb_memory_class,                ONLY: t_jsb_memory
  USE mo_jsb_process_class,               ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, &
                                                PHENO_, SB_, Get_process_name
  USE mo_jsb_process_factory,             ONLY: Create_process, Create_process_memory, max_no_of_processes
  ! shared
  USE mo_shared_constants,                ONLY: init_shared_constants
  USE mo_shared_set_parameters,           ONLY: parameters_from_file    !< TYPE t_parameters_from_file
  ! atmland
  USE mo_atmland_constants,               ONLY: init_atmland_constants
  ! radiation
  USE mo_rad_constants,                   ONLY: init_rad_constants
  USE mo_rad_init,                        ONLY: rad_init
  ! assimilation (photosynthesis)
  USE mo_assimi_constants,                ONLY: init_assimi_constants
  USE mo_assimi_init,                     ONLY: assimi_init
  ! phenology
  USE mo_pheno_constants,                 ONLY: init_pheno_constants
  USE mo_pheno_init,                      ONLY: pheno_init
  ! vegetation
  USE mo_veg_constants,                   ONLY: init_veg_constants
  USE mo_veg_config_class,                ONLY: t_veg_config
  USE mo_veg_init,                        ONLY: veg_init
  USE mo_veg_memory_class,                ONLY: Init_veg_pools_elements_a
  ! soil
  USE mo_soil_constants,                  ONLY: init_soil_constants
  USE mo_soil_config_class,               ONLY: t_soil_config
  USE mo_soil_init,                       ONLY: soil_init
  ! soil_biogeochemistry
  USE mo_sb_constants,                    ONLY: init_sb_constants
  USE mo_sb_init,                         ONLY: sb_init
  USE mo_sb_memory_class,                 ONLY: Init_sb_pools_elements_a

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
    MODULE PROCEDURE init_slm
    !MODULE PROCEDURE init_model_standalone
  END INTERFACE init_model


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> init all processes and memory of SLM
  !!
  !! @todo mo_jsb_model_init: processes may init themselves
  !! 
  !! @todo mo_jsb_model_init: move to default JSBACH4 init style
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm(tile, l_forcing, nsoil_sb, iblk, ics, ice, nc)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)              :: tile           !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(out)                         :: l_forcing      !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
    INTEGER,INTENT(out)                         :: nsoil_sb       !< number of soil layers used in SB_ - equal to nsoil_water
                                                                  !!  necessary for open_pool_output_files() & init_output_writing_pool()
    INTEGER,INTENT(in)                          :: iblk, &        !< iBlock, iChunkStart, iChunkEnd, number of chunks
                                                   ics, &
                                                   ice, &
                                                   nc 
    ! ---------------------------
    ! 0.2 Local
    CLASS(t_jsb_memory), POINTER                :: tmp_ptr        !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: process_name   !! for Memory init and alloc
    INTEGER                                     :: iproc          !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: prefix         !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: suffix         !! for Memory init and alloc
    INTEGER                                     :: no_of_models,& !! SLM: we use just one model
                                                   model_id       !! SLM: per default model_id should be 1
    TYPE(t_jsb_model), POINTER                  :: current => NULL() ! model allocation to jsbach struct - in jsbach4 this is looped over x models
    TYPE(t_jsb_model), POINTER                  :: model
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SOIL_)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_slm'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science / model init
    ! ------------------------------------------------------------------------------------------------------------


    !! 1 Model and Process init 
    !!
    !! 
    !! model
    no_of_models = 1            ! *slmtemp 
    model_id     = 1            ! *slmtemp 
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models
    current => new_model(model_id, "slm", "slm", base_config%slm_model_name, "namelist.slm") ! model_id, name, shortname, description, namelist_filename
    ! "current" is used for looping if no_of_models > 1 ; but looping would need to include process configs per model (see section 2)
    jsbach%models(model_id)%m => current
    ! point model to jsbach%models(1)%m
    model => Get_model(model_id)  !     model => Get_model(model_id)   *slmtemp hardcoded model_id !!!
    !! processes
    ALLOCATE(model%processes(max_no_of_processes))
    DO iproc=1,max_no_of_processes      ! *slmtmp currently only creating: SOIL_ VEG_ ASSIMI_ PHENO_ SB_  (mo_jsb_process_factory:Create_process)
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
    END DO


    !! 2 Config init, run PROC_config%Init routines and Forcing init
    !!
    !!  
    !! process config
    base_config%namelist_filename       = model%namelist_filename
    CALL base_config%Init               ! includes base_config%slm_model_name => land plant soil canopy test_canopy
    CALL model%Configure_processes()    ! mo_jsb_model_class:Configure_model_processes CALL mo_jsb_process_class:Configure_process CALL config%Init()
    !! forcing config
    CALL Init_forcing(model_id)         ! model_id = 1, see above new_model()
    !! get process configs needed for model_init
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SOIL_)


    !! 3 Set model options
    !!
    !! 
    model%options%dtime         = base_config%dtime  ! length of one time step
    model%options%steplen       = base_config%dtime  ! .. 
    model%options%iblk          = iblk               ! Number of current chunk that consists of multiple grid boxes
    model%options%ics           = ics                ! Index of chunk start
    model%options%ice           = ice                ! Index of chunk end
    model%options%nc            = nc                 ! Length of chunk  (end-number - start-number + 1)
    model%options%nsoil_energy  = dsl4jsb_Config(SOIL_)%nsoil_energy  ! number of soil layers for soil energy calculations
    model%options%nsoil_water   = dsl4jsb_Config(SOIL_)%nsoil_water   ! number of soil layers for soil moisture calculations
    model%options%nsoil_sb      = dsl4jsb_Config(SOIL_)%nsoil_water   ! number of soil layers used in SB_ - equals nsoil_water
    nsoil_sb                    = dsl4jsb_Config(SOIL_)%nsoil_water   ! SUBROUTINE output for land_model_driver -> open_pool_output_files()
    !model%lctlib => Read_lctlib()             ! not possible here - need init constants init !!


    !! 4 Tile init
    !! 
    !! 
    !! @todo model_init -> init lcts of the tile via PROCEDURE Construct_jsb_lct  *slmtemp
    tile%owner_model_id = current%id                    ! *slmtemp
    ! lcts
    tile%lcts%id        = 2                             ! VEG_TYPE     2 == Vegetated
    tile%lcts%name      = "veg"                         ! Short name for lct
    tile%lcts%long_name = "slm vegetation"              ! Long name for lct
    tile%lcts%lib_id    = dsl4jsb_Config(VEG_)%pft_id   ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)


    !! 5 Grid init
    !! 
    !! 
    grid%lat            = base_config%lat
    grid%lon            = base_config%lon


    !! 6 Constants
    !! 

    !! 6.1 Constants init
    !! 
    !! 
    CALL init_base_constants          ! may be 1st of the 'init constants'
    CALL init_shared_constants        ! may be 2nd of the 'init constants'
    CALL init_rad_constants           ! may be 3rd of the 'init constants'
    CALL init_atmland_constants       ! uses hourday                  defined in mo_base_constants
    CALL init_assimi_constants        ! uses molar_mass_N, Dwv, Dco2  defined in mo_shared_constants
    CALL init_veg_constants           ! uses Tzero                    defined in mo_shared_constants
    CALL init_soil_constants
    CALL init_sb_constants
    CALL init_pheno_constants

    !! 6.2 Constants set (using param.list file)
    !! 
    !! 
    IF (base_config%set_parameter_values_from_file)  CALL parameters_from_file%set_param


    !! 7 Memory init and alloc
    !! 
    !! >> this is work in progress ! and, thus, not very nice code, but it works for now (Jan 20160819) <<
    !! 
    !! process memory allocation - borrowed from JSBACH4 git 20160818 mo_jsb_tile 

    ! Loop over all processes and initialize memory for this tile
    ALLOCATE(tile%mem(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      tile%mem(iproc)%p => NULL()
      SELECT CASE(iproc)
        ! SLM: extend this CASE to all SLM processes that work with JSBACH4 memory structure; see mo_jsb_process_class for enumeration
        CASE (RAD_, ASSIMI_, VEG_, PHENO_, SOIL_, SB_, A2L_, L2A_)    
          prefix = 'no-need'
          suffix = 'no-need'
          process_name = Get_process_name(iproc)
          tmp_ptr => Create_process_memory(iproc) ! see mo_jsb_process_factory ==> uncomment lines here for new process memory !!
          SELECT TYPE (tmp_ptr)
          CLASS IS (t_jsb_memory)
            tile%mem(iproc)%p => tmp_ptr
            tile%mem(iproc)%p%id = iproc
            !tile%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
            !tile%mem(iproc)%p%grid_id = tile%grid_id
            tile%mem(iproc)%p%owner_model_id = tile%owner_model_id
            tile%mem(iproc)%p%owner_proc_id = iproc  ! process ID: A2L_=1, L2A_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10 , PHENO_=11
            !ALLOCATE(tile%mem(iproc)%p%owner_tile_path(SIZE(tile%path)))
            !tile%mem(iproc)%p%owner_tile_path(:) = tile%path(:)
            !tile%mem(iproc)%p%owner_tile_name = TRIM(tile%name)
!!            IF (prefix /= '') THEN
!!              CALL tile%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(suffix), tile%lcts(:)%id, tile%owner_model_id)
!!            ELSE
!!              CALL tile%mem(iproc)%p%Init(TRIM(process_name), TRIM(suffix), tile%lcts(:)%id, tile%owner_model_id)
!!            END IF
            ! call init routine of memory
            CALL tile%mem(iproc)%p%Init( TRIM(process_name), TRIM(suffix), tile%lcts%id, tile%owner_model_id)
        END SELECT  ! processes with memory (A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_)
        NULLIFY(tmp_ptr)    
      END SELECT    ! iproc
    END DO          ! iproc


    !! 7.1 Memory init -> pool and element structures in veg_ and sb_ memory with ALLOCATABLE variables with operator overload
    !! 
    CALL Init_sb_pools_elements_a(ics, ice, nc, iblk, nsoil_sb)
    CALL Init_veg_pools_elements_a(ics, ice, nc, iblk)


    !! 8 lctlib init and read
    !!
    !! Important! Read_lctlib() needs config and constants init to be finished !
    !! 
    model%lctlib => Read_lctlib("lctlib_slm_nlct8.def")


    !! 9 Time init
    !! 
    !! 
    time_day          = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                           ! the 1st year of a simulation == 1 (not zero)
    dtime             = base_config%dtime           ! time step length in seconds, default: 1800.0 == 30 minutes, mo_base_constants
    ntimestep         = forcing_options%Timesteps() ! number of time steps (considers the simulation_mode static & transient)
    timestep_prop_day = dtime/one_day               ! proportion of one day that is simulated with one timestep


    !! 10 Initialise: Fluxnet specific simulation output
    !! 
    !! the last full "forcing period" is written into extra output files with half-hourly write-out of values
    fluxnet_output_status = 0
    IF (base_config%fluxnet_output) THEN
      fluxnet_output_status = 1
      ! check whether simulation length covers at least one complete fluxnet forcing period (i.e. fluxnet_output_length)
      fluxnet_output_length = NINT(base_config%fluxnet_output_length * one_year * one_day / base_config%dtime)
      IF (fluxnet_output_length > ntimestep) THEN
        CALL finish('land_model_driver-init_routine','Simulation is too short for fluxnet output of one full forcing period')
      END IF

      ! calc timestep when fluxnet output starts
      fluxnet_output_start  = fluxnet_output_length * FLOOR((ntimestep / fluxnet_output_length ) - 1.0_wp)
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


    !! 11 run model-specific process init routine  (land plant soil canopy test_canopy)
    !! 
    !! 
    ! enable/ disable forcing
    ! note: the soil model gets a different forcing than land plant canopy, see mo_atmland_forcing
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy','soil')
        l_forcing = .TRUE.
      CASE DEFAULT
        l_forcing = .FALSE.
    END SELECT
    ! run model init
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        CALL pheno_init   (tile, model%options)
        CALL assimi_init  (tile, model%options)
        CALL rad_init     (tile, model%options)
        CALL veg_init     (tile, model%options)
        CALL soil_init    (tile, model%options)
        CALL sb_init      (tile, model%options)
      CASE('soil')
        CALL veg_init     (tile, model%options)       ! strictly speaking not necessary, but to be safe
        CALL soil_init    (tile, model%options)
        CALL sb_init      (tile, model%options)
      CASE('canopy')
        CALL veg_init     (tile, model%options)
        CALL assimi_init  (tile, model%options)
        CALL rad_init     (tile, model%options)
        CALL soil_init    (tile, model%options)
      CASE('test_canopy')
        CALL veg_init     (tile, model%options)
        CALL assimi_init  (tile, model%options)
        CALL rad_init     (tile, model%options)
        ! init is done within the run_test_canopy_model(), which is reasonable because this model runs for only one timestep
      CASE DEFAULT
        CALL finish('land_model_driver-init_routine','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
     END SELECT

  END SUBROUTINE init_slm

END MODULE mo_jsb_model_init
