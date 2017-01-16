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
    ! base
    USE mo_base_constants                     ! e.g. dtime, init_base_constants
    USE mo_jsb_class,                       ONLY: Get_model, jsbach
    USE mo_base_memory,                     ONLY: t_slm_tile, t_model, grid, new_model
    USE mo_base_config_class,               ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element, init_lctlib
    USE mo_jsb_memory_class,                ONLY: t_jsb_memory
    USE mo_jsb_process_class,               ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, &
                                                  PHENO_, SB_, Get_process_name
    USE mo_jsb_process_factory,             ONLY: Create_process_memory, max_no_of_processes
    ! shared
    USE mo_shared_constants,                ONLY: init_shared_constants
    USE mo_shared_set_parameters,           ONLY: parameters_from_file    !< TYPE t_parameters_from_file
    ! atmland
    USE mo_atmland_constants,               ONLY: init_atmland_constants
    ! radiation
    USE mo_rad_constants,                   ONLY: init_rad_constants
    ! assimilation (photosynthesis)
    USE mo_assimi_constants,                ONLY: init_assimi_constants
    USE mo_assimi_config_class,             ONLY: assimi_config
    ! phenology
    USE mo_pheno_constants,                 ONLY: init_pheno_constants
    USE mo_pheno_config_class,              ONLY: pheno_config
    USE mo_pheno_init,                      ONLY: pheno_init
    ! vegetation
    USE mo_veg_constants,                   ONLY: init_veg_constants
    USE mo_veg_config_class,                ONLY: veg_config
    USE mo_veg_init,                        ONLY: veg_init
    ! soil
    USE mo_soil_constants,                  ONLY: init_soil_constants
    USE mo_soil_config_class,               ONLY: soil_config
    USE mo_soil_init,                       ONLY: soil_init
    ! soil_biogeochemistry
    USE mo_sb_constants,                    ONLY: init_sb_constants
    USE mo_sb_config_class,                 ONLY: sb_config
    USE mo_sb_init,                         ONLY: sb_init

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_model
  PUBLIC :: time_day, time_year, ntimestep, itimestep, timestep_prop_day
  PUBLIC :: fluxnet_output_start, fluxnet_output_stop, fluxnet_output_status


    ! variables for the driver
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
  SUBROUTINE init_slm(tile, l_forcing) !(model_id, model_name)


    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)              :: tile       !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(out)                         :: l_forcing  !< IF .TRUE. THEN update_forcing(), update_atm2land(), update_land2atm()
    ! 0.2 Local
    CLASS(t_jsb_memory), POINTER                :: tmp_ptr        !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: process_name   !! for Memory init and alloc
    INTEGER                                     :: iproc          !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: prefix         !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: suffix         !! for Memory init and alloc
    INTEGER                                     :: no_of_models
    TYPE(t_model), POINTER                      :: current => NULL() ! model allocation to jsbach struct - in jsbach4 this is looped over x models
    TYPE(t_model), POINTER                      :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_slm'

    ! check number of tiles, exit if more than output routines can handle
    IF(ntile > 1) CALL finish('land_model_driver-init_routine','Number of tiles is too large, i.e. larger than 1')
    ! 



    !! 1 Read config paramters from namelist.plant section
    !!
    !!  
    base_config%namelist_filename       = 'namelist.slm'
    assimi_config%namelist_filename     = 'namelist.slm'
    veg_config%namelist_filename        = 'namelist.slm'
    soil_config%namelist_filename       = 'namelist.slm'
    sb_config%namelist_filename         = 'namelist.slm'
    pheno_config%namelist_filename      = 'namelist.slm'
    CALL base_config%Init              ! includes base_config%slm_model_name => land plant soil canopy test_canopy
    CALL assimi_config%Init
    CALL veg_config%Init
    CALL soil_config%Init
    CALL sb_config%Init
    CALL pheno_config%Init


    !! 2 Model init including lctlib init
    !!
    !! see jsbach4 mo_jsb_base
    !! 
    no_of_models = 1
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models
    current => new_model(1, "slm", "slm", "slm_veg_model", "not_available.txt") ! id, name, shortname, description, namelist_filename
    current%options%dtime   = base_config%dtime  ! length of one time step
    current%options%steplen = base_config%dtime  ! .. 
    current%options%iblk    = 1               ! Number of current block (chunk)
    current%options%ics     = 1               ! Index of chunk start
    current%options%ice     = 1               ! Index of chunk end
    current%options%nc      = current%options%ice - current%options%ics + 1   ! Length of chunk
    !current%lctlib => init_lctlib()           ! not possible here - need init constants init !!
    jsbach%models(1)%m => current
    ! point model to jsbach%models(1)%m
    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!


    !! 3 Tile init
    !! 
    !! 
    !! @todo model_init -> init lcts of the tile via PROCEDURE Construct_jsb_lct  *slmtemp
    tile%owner_model_id = current%id         ! *slmtemp
    ! lcts
    tile%lcts%id        = 2                  ! VEG_TYPE     2 == Vegetated
    tile%lcts%name      = "veg"              ! Short name for lct
    tile%lcts%long_name = "slm vegetation"   ! Long name for lct
    tile%lcts%lib_id    = veg_config%pft_id  ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)


    !! 4 Grid init
    !! 
    !! 
    grid%lat            = base_config%lat
    grid%lon            = base_config%lon


    !! 5 Constants
    !! 

    !! 5.1 Constants init
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

    !! 5.2 Constants set (using param.list file)
    !! 
    !! 
    IF (base_config%set_parameter_values_from_file)  CALL parameters_from_file%set_param


    !! 6 Memory init and alloc
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
        CASE (RAD_, ASSIMI_, VEG_, PHENO_, SOIL_, SB_, ATM2LAND_, LAND2ATM_)    
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
            tile%mem(iproc)%p%owner_proc_id = iproc  ! process ID: ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10 , PHENO_=11
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
        END SELECT  ! processes with memory (ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_)
        NULLIFY(tmp_ptr)    
      END SELECT    ! iproc
    END DO          ! iproc


    !! 7 lctlib init
    !!
    !! Important! init_lctlib() needs config and constants init to be finished !
    !! 
    model%lctlib => init_lctlib()         ! initialize the lctlib with parameters values across PFT's


    !! 8 Time init
    !! 
    !! 
    time_day          = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                           ! the 1st year of a simulation == 1 (not zero)
    dtime             = base_config%dtime              ! time step length in seconds, default: 1800.0 == 30 minutes, mo_base_constants
    ntimestep         = base_config%Timesteps()     ! number of time steps
    timestep_prop_day = dtime/one_day                  ! proportion of one day that is simulated with one timestep


    !! 9 Initialise: Fluxnet specific simulation output
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


    !! 10 run model-specific init routine  (land plant soil canopy test_canopy)
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
    ! for land the same init routine as plant is run
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        CALL pheno_init(tile)     ! pheno 
        CALL veg_init(tile)       ! veg 
        CALL sb_init(tile)        ! sb 
        CALL soil_init(tile)      ! soil
      CASE('soil')
        CALL veg_init(tile)       ! veg, strictly speaking not necessary, but to be safe
        CALL sb_init(tile)        ! sb 
        CALL soil_init(tile)      ! soil, strictly speaking not necessary, but to be safe
      CASE('canopy')
        CALL veg_init(tile)       ! veg 
        CALL soil_init(tile)      ! soil
      CASE('test_canopy')
        ! init is done within the run_test_canopy_model(), which is reasonable because this model runs for only one timestep
      CASE DEFAULT
        CALL finish('land_model_driver-init_routine','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
     END SELECT

  END SUBROUTINE init_slm

END MODULE mo_jsb_model_init
