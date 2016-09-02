!>
!! @brief Contains routines for model init
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-22); starting with custom SLM init, moving to default JSBACH4 init later
!!
MODULE mo_jsb_model_init

    USE mo_kind
    USE mo_exception,                       ONLY: finish, message, message_text
    ! base
    USE mo_base_constants
    USE mo_base_memory,                     ONLY: t_slm_tile, t_model, grid 
    USE mo_base_config_class,               ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element, init_lctlib
    USE mo_jsb_memory_class,                ONLY: t_jsb_memory
    USE mo_jsb_process_class,               ONLY: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, HD, RAD, ASSIMI, PHENO, &
                                                  Get_process_name
    USE mo_jsb_process_factory,             ONLY: Create_process_memory, max_no_of_processes
    ! shared
    USE mo_shared_isotope,                  ONLY: calc_fractionation 
    USE mo_shared_constants,                ONLY: init_shared_constants, molar_mass_N, Dwv, Dco2, Tzero
    ! atmland
    USE mo_atmland_constants
    ! radiation
    USE mo_rad_constants,                   ONLY: init_rad_constants
    ! assimilation (photosynthesis)
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps 
    USE mo_assimi_constants
    USE mo_assimi_config_class,             ONLY: assimi_config
    ! vegetation
    USE mo_veg_constants
    USE mo_veg_config_class,                ONLY: veg_config


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
    REAL(wp) :: timestep_prop_day            !< proportion of one day that is simulated with one timestep (dt/one_day)

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
  SUBROUTINE init_slm(model, tile, l_forcing) !(model_id, model_name)


    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_model),INTENT(inout)                 :: model      !< data structure for dynamic vegetation and soil variables
    TYPE(t_slm_tile),INTENT(inout)              :: tile       !< data structure for dynamic vegetation and soil variables
    LOGICAL,INTENT(out)                         :: l_forcing  !< IF .TRUE. THEN update_forcing(), update_atm2land(), update_land2atm()
    ! 0.2 Local
    CLASS(t_jsb_memory), POINTER                :: tmp_ptr        !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: process_name   !! for Memory init and alloc
    INTEGER                                     :: iproc          !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: prefix         !! for Memory init and alloc
    CHARACTER(LEN=15)                           :: suffix         !! for Memory init and alloc


    ! check number of tiles, exit if more than output routines can handle
    IF(ntile > 1) CALL finish('land_model_driver-init_routine','Number of tiles is too large, i.e. larger than 1')
    ! 



    !! 1 Read config paramters from namelist.plant section
    !!
    !!  
    base_config%namelist_filename       = 'namelist.slm'
    assimi_config%namelist_filename     = 'namelist.slm'
    veg_config%namelist_filename        = 'namelist.slm'
    CALL base_config%Init              ! includes base_config%slm_model_name => land plant soil canopy test_canopy
    CALL assimi_config%Init
    CALL veg_config%Init


    !! 2 Model init
    !!
    !!  
    !! 
    model%id = 0                                    ! this is just an arbitrary value consistent with tile%owner_model_id
    model%shortname = "slm"                         ! arbitrary name
    model%options%dtime   = base_config%dt  ! length of one time step
    model%options%steplen = base_config%dt  ! .. 
    model%options%iblk    = 1               ! Number of current block (chunk)
    model%options%ics     = 1               ! Index of chunk start
    model%options%ice     = 1               ! Index of chunk end
    model%options%nc      = model%options%ice - model%options%ics + 1   ! Length of chunk


    !! 3 Tile init
    !! 
    !! 
    !! @todo land model driver: update Intel compiler to enable: init lcts of the tile via PROCEDURE Construct_jsb_lct
    ! lcts
    tile%lcts%id        = 2                  ! VEG_TYPE     2 == Vegetated
    tile%lcts%name      = "veg"              ! Short name for lct
    tile%lcts%long_name = "slm vegetation"   ! Long name for lct
    tile%lcts%lib_id    = veg_config%pft_id  ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)
    tile%ics            = model%options%ics  ! simple (temporary) solution to make it available across Tasks
    tile%ice            = model%options%ice  ! simple (temporary) solution to make it available across Tasks
    tile%iblk           = model%options%iblk ! simple (temporary) solution to make it available across Tasks

    !! 4 Grid init
    !! 
    !! 
    grid%lat            = base_config%lat
    grid%lon            = base_config%lon


    !! 5 Memory init and alloc
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
        CASE (RAD, ASSIMI, VEG, ATM2LAND, LAND2ATM)    
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
            tile%mem(iproc)%p%owner_proc_id = iproc
            !ALLOCATE(tile%mem(iproc)%p%owner_tile_path(SIZE(tile%path)))
            !tile%mem(iproc)%p%owner_tile_path(:) = tile%path(:)
            !tile%mem(iproc)%p%owner_tile_name = TRIM(tile%name)
!!            IF (prefix /= '') THEN
!!              CALL tile%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(suffix), tile%lcts(:)%id, tile%owner_model_id)
!!            ELSE
!!              CALL tile%mem(iproc)%p%Init(TRIM(process_name), TRIM(suffix), tile%lcts(:)%id, tile%owner_model_id)
!!            END IF
            ! call init routine of memory here & added model%options -> changed for SLM compatibility
            CALL tile%mem(iproc)%p%Init( TRIM(process_name), TRIM(suffix), tile%lcts%id, tile%owner_model_id, &
                                        model%options%ics, model%options%ice, model%options%iblk) ! model%options -> for SLM compatibility
        END SELECT  ! processes with memory (ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, HD, RAD, ASSIMI, PHENO)
        NULLIFY(tmp_ptr)    
      END SELECT    ! iproc
    END DO          ! iproc


    !! 6 Constants init
    !! 
    !! 
    CALL init_base_constants          ! may be 1st of the 'init constants'
    CALL init_shared_constants        ! may be 2nd of the 'init constants'
    CALL init_rad_constants           ! may be 3rd of the 'init constants'
    CALL init_atmland_constants       ! uses hourday                  defined in mo_base_constants
    CALL init_assimi_constants        ! uses molar_mass_N, Dwv, Dco2  defined in mo_shared_constants
    CALL init_veg_constants           ! uses Tzero                    defined in mo_shared_constants
      ! CALL init_soil_constants


    !! 7 Plant init
    !! 
    !! 
    CALL tile%plant%init                  ! mo_veg_memory  ==> may move to veg init (or model init ?)
    CALL tile%plant%flux%set_zero         ! mo_veg_memory  ==> may move to veg init (or model init ?)
    model%lctlib => init_lctlib()         ! initialize the lctlib with parameters values across PFT's


    !! 8 Time init
    !! 
    !! 
    time_day          = 0.0_wp                      ! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                           ! the 1st year of a simulation == 1 (not zero)
    dt                = base_config%dt              ! time step length in seconds, default: 1800.0 == 30 minutes
    ntimestep         = base_config%Timesteps()     ! number of time steps
    timestep_prop_day = dt/one_day                  ! proportion of one day that is simulated with one timestep


    !! 9 Initialise: Fluxnet specific simulation output
    !! 
    !! the last full "forcing period" is written into extra output files with half-hourly write-out of values
    fluxnet_output_status = 0
    IF (base_config%fluxnet_output) THEN
      fluxnet_output_status = 1
      ! check whether simulation length covers at least one complete fluxnet forcing period (i.e. fluxnet_output_length)
      fluxnet_output_length = NINT(base_config%fluxnet_output_length * one_year * one_day / base_config%dt)
      IF (fluxnet_output_length > ntimestep) THEN
        CALL finish('land_model_driver-init_routine','Simulation is too short for fluxnet output of one full forcing period')
      END IF

      ! calc timestep when fluxnet output starts
      fluxnet_output_start  = fluxnet_output_length * FLOOR((ntimestep / fluxnet_output_length ) - 1.0_wp)
      ! calc timestep when fluxnet output stops
      fluxnet_output_stop   = fluxnet_output_start + fluxnet_output_length

      ! check whether genral output start and stop covers the 'Fluxnet specific simulation output'
      IF ( (veg_config%output_end_last_day_year * one_year * one_day / base_config%dt) < fluxnet_output_stop) THEN
        CALL finish('land_model_driver-init_routine','Simulation output does not cover the last full forcing period')
      END IF
      IF ( (((veg_config%output_start_first_day_year - 1) * one_year * one_day + 1) / base_config%dt) > fluxnet_output_start) THEN
        CALL finish('land_model_driver-init_routine','Simulation output does not cover the last full forcing period')
      END IF
    END IF


    !! 10 run model-specific init routine  (land plant soil canopy test_canopy)
    !! 
    !! 
    ! enable/ disable forcing
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy') ! ,'soil'
        l_forcing = .TRUE.
      CASE DEFAULT
        l_forcing = .FALSE.
    END SELECT
    ! run model init
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land')
        ! currently same as plant -> calling simply 'init_slm_plant'!
        CALL init_slm_plant(model, tile)
      CASE('plant')
        CALL init_slm_plant(model, tile)
      CASE('soil')
        CALL finish('land_model_driver-init_routine','Nice try, but the soil model does not exist yet :) ')
      CASE('canopy')
        CALL init_slm_canopy(model, tile)
      CASE('test_canopy')
!         init is done within the run_test_canopy_model(), which is reasonable because this model runs for only one timestep
      CASE DEFAULT
        CALL finish('land_model_driver-init_routine','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
     END SELECT

  END SUBROUTINE init_slm


  !-----------------------------------------------------------------------------------------------------
  ! init routine specific to the land model
  ! 
  ! tbc ...
  !-----------------------------------------------------------------------------------------------------




  !-----------------------------------------------------------------------------------------------------
  !> init routine specific to the plant model
  !!
  !! called from init_slm
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm_plant(model, tile) !(model_id, model_name)

    USE mo_shared_memory, ONLY: t_element

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_model),INTENT(inout),TARGET           :: model      !< data structure for dynamic vegetation and soil variables
    TYPE(t_slm_tile),INTENT(inout)               :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_lctlib_element),POINTER               :: lctlib     !< land-cover-type library - parameter across pft's
    TYPE(t_element)                              :: init_pool  !< land-cover-type library - parameter across pft's
    REAL(wp)                                     :: hlp1

    !  bi-directional pointer
    lctlib => model%lctlib(tile%lcts%lib_id)

    ! accumulation
    CALL tile%plant%daily_average%init      !< initially sets all variables to zero or .FALSE.
    CALL tile%plant%moving_average%init     !< init variables with reasonable default values and sets accumulation periods


      tile%soil%moisture                 = veg_config%soil_awc_prescribe
      tile%soil%theta                    = veg_config%soil_theta_prescribe 
      tile%soil%ammonium                 = 0.1_wp
      tile%soil%phosphate                = 0.01_wp
      tile%soil%pool%litter%nitrogen     = 5.0_wp
      tile%plant%pftid                   = veg_config%pft_id
      tile%plant%target_cn_leaf          = lctlib%cn_leaf
      tile%plant%target_np_leaf          = lctlib%np_leaf

      IF(veg_config%flag_start_bareground)THEN
         init_pool%carbon        = lctlib%seed_size * 10.0_wp
         tile%plant%dens_ind     = 0.0_wp 
      ELSE
         init_pool%carbon        = 50.0_wp 
         tile%plant%dens_ind     = 1.0_wp 
      ENDIF

      init_pool%nitrogen  = init_pool%carbon / &
                            (1._wp + fresp_growth) / lctlib%cn_leaf
      init_pool%phosphorus = init_pool%nitrogen / &
                            lctlib%np_leaf

      ! init_pool%carbon13
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ELSE
         hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ENDIF
      init_pool%carbon13 = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13,hlp1) * init_pool%carbon 


      ! init_pool%carbon14
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ELSE
         hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ENDIF
      init_pool%carbon14 = calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C14,hlp1) * init_pool%carbon

      ! start bareground ?
      IF(veg_config%flag_start_bareground)THEN
         tile%plant%pool%seed_bed = init_pool 
      ELSE
         tile%plant%pool%reserve  = init_pool 
      ENDIF
 
  END SUBROUTINE init_slm_plant


  !-----------------------------------------------------------------------------------------------------
  !> init routine specific to the canopy model
  !!
  !! called from init_slm
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm_canopy(model, tile) !(model_id, model_name)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_model),INTENT(inout),TARGET           :: model      !< data structure for dynamic vegetation and soil variables
    TYPE(t_slm_tile),INTENT(inout)               :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_lctlib_element),POINTER               :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                                     :: ref_lai    !< LAI
    REAL(wp)                                     :: ref_nleaf  !< N leaf per unit leaf area in gN/m2

    !  bi-directional pointer
    lctlib => model%lctlib(tile%lcts%lib_id)

    ref_lai   = 6._wp
    ref_nleaf = 2.0_wp

    ! accumulation
    CALL tile%plant%daily_average%init      !< initially sets all variables to zero or .FALSE.
    CALL tile%plant%moving_average%init     !< init variables with reasonable default values and sets accumulation periods

    !! 1.3 Initialise tile
    !! 
    !! Initialisation of soil and vegetation, so far relevant
    tile%soil%moisture   = veg_config%soil_awc_prescribe
    tile%soil%theta      = veg_config%soil_theta_prescribe
    tile%plant%lai       = ref_lai 
    IF(veg_config%pft_id==0) THEN
      tile%plant%pool%leaf%nitrogen   = ref_nleaf / molar_mass_N * tile%plant%lai
    ELSE
      tile%plant%pool%leaf%nitrogen   = tile%plant%lai / & 
                                                lctlib%sla / lctlib%cn_leaf
      tile%plant%pool%leaf%phosphorus = tile%plant%pool%leaf%nitrogen / & 
                                                lctlib%np_leaf
      tile%plant%pftid                = veg_config%pft_id
    ENDIF

  END SUBROUTINE init_slm_canopy


END MODULE mo_jsb_model_init
