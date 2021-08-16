!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief
!!  Contains methods for model initializations.
!!
!! @author
!!  Reiner Schnur, MPI-Met
!!
!! @par Revision History
!!  the JSBACH4 model-init module was tailored to fit QUINCY by Jan Engel MPI-BGC
!! 
!!
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

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
  SUBROUTINE setup_grid_standalone(model_id, grid_name, nproma)

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
    CALL register_grid(grid)            ! add to mo_jsb_grid:jsbach_grids(:), actuall to jsbach_grids(1)
    model%grid_id = grid%id             ! tell the model the id of the hgrid, so that one can get the correct hgrid per model

    NULLIFY(model)    ! NULLIFY pointer to model

  END SUBROUTINE setup_grid_standalone


  !-----------------------------------------------------------------------------------------------------
  !> init quincy model
  !! 
  !! and init of: forcing, SB_ & VEG_ pools, tiles, time, fluxnet diagnostics simulation state
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE quincy_init(model_id, iblk, nproma, ics, ice, nc, l_forcing)

    ! base
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_lct_class,         ONLY: Contains_lct
    USE mo_jsb_model_usecases,    ONLY: init_usecase
    USE mo_jsb4_forcing,          ONLY: Init_forcing, forcing_options
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_math_constants,    ONLY: dtime, one_year, one_day
    ! shared 
    USE mo_util,                  ONLY: one_of
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
    ! process configurations (t_PROC_config)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SSE_)
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER ,INTENT(in)                       :: model_id
    INTEGER ,INTENT(in)                       :: iblk, nproma         !< block number, number of points per chunk   *slmtmp
    INTEGER ,INTENT(in)                       :: ics, ice, nc         !< first & last & number of points in the chunk the model currently works on  *slmtmp
    LOGICAL ,INTENT(out)                      :: l_forcing            !< IF .TRUE. THEN update_forcing(), update_a2l(), update_l2a()
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER                :: model
    CLASS(t_jsb_tile_abstract), POINTER       :: tile
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


    !> 0.9 init local variables
    !! 
    CALL message(TRIM(routine), 'starting model init')
    no_omp_thread = 1     ! no threading with OpenMPI in QUINCY


    !> 1.0 forcing config
    !! 
    CALL Init_forcing(model_id)


    !> 2.0 model options
    !! 
    model%options%dtime        = model%config%dtime ! length of one time step
    model%options%steplen      = model%config%dtime ! 
    model%options%iblk         = iblk               ! 
    model%options%nproma       = nproma             ! 
    model%options%ics          = ics                ! 
    model%options%ice          = ice                ! 
    model%options%nc           = nc                 ! 


    !> 3.0 SB_ & VEG_ pools (pool and element structures in veg_ and sb_ memory with ALLOCATABLE variables with operator overload)
    !!
    CALL Init_sb_pools_elements_a(ics, ice, nc, iblk, dsl4jsb_Config(SSE_)%nsoil_water)
    CALL Init_veg_pools_elements_a(ics, ice, nc, iblk)


    !> 4.0 Init tiles
    !!
    ! At this point, the model configuration and tile structure is set up.
    ! We continue with initializing the memory structure and different model processes.
    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    !> 4.1 loop over tiles
    !! 
    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      !> 4.2 Initialize tile memory and read tile fractions
      !! mo_jsb_tile:Init_tile() -> init process memories
      !! @todo: once we have land use change or dynveg, l_fixed_fractions needs to be adapted <- JSBACH4
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name), l_fixed_fractions=.TRUE.)

      !> 4.3 tile init
      !! @todo model_init -> init lcts of the tile via PROCEDURE Construct_jsb_lct  *slmtemp  ?
      tile%owner_model_id = model%id                      ! *slmtemp
      ! lcts
      tile%lcts%id        = 2                             ! VEG_TYPE     2 == Vegetated
      tile%lcts%name      = "veg"                         ! Short name for lct
      tile%lcts%long_name = "quincy_vegetation"           ! Long name for lct
      tile%lcts%lib_id    = dsl4jsb_Config(VEG_)%pft_id   ! pftID == column number in model%lctlib == ID of lct in lctlib (JSBACH4)

      !! this option is outdated; forcing is always used, but not for test_canopy, as the below CASE indicates :)
      !!
      !! enable/ disable forcing
      !! note: the soil model gets a different forcing than land plant canopy, see mo_atmland_forcing
      SELECT CASE(TRIM(model%config%slm_model_name))
        CASE('land','plant','canopy','soil')
          l_forcing = .TRUE.
        CASE DEFAULT
          l_forcing = .FALSE.
      END SELECT

      !> 4.4 run process init
      !! 
      !! one_of() returns -1 if NONE of the SLM model names are found, but the 1st found SLM model name get the number 1 (no enum just loop)
      !! 
      !! a certain sequence of proc_init (i.e. the below sequence) is needed to correctly setup the processes, because of some 
      !!  dependencies across two or three processes
      ! test model name
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','soil','canopy','test_canopy','test_radiation' ]) == -1) THEN
        CALL finish(routine,'Invalid SLM model name. Try either land, plant, soil, canopy, test_canopy or test_radiation.')
      END IF
      ! pheno_init
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','canopy' ]) > 0) THEN
        CALL pheno_init(tile)
      END IF
      ! assimi_init
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','canopy','test_canopy','test_radiation' ]) > 0) THEN
        CALL assimi_init(tile)
      END IF
      ! rad_init
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','canopy','test_canopy','test_radiation' ]) > 0) THEN
        CALL rad_init(tile)
      END IF
      ! veg_init
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','soil','canopy','test_canopy','test_radiation' ]) > 0) THEN
        CALL veg_init(tile, model%options)
      END IF
      ! sse_init
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','soil','canopy' ]) > 0) THEN
        CALL sse_init(tile, model%options)
      END IF
      ! sb_init (needs sse_init to be finished)
      IF (one_of(model%config%slm_model_name, [character(len=20) :: &
          'land','plant','soil' ]) > 0) THEN
        CALL sb_init(tile, model%options)
      END IF
      ! next tile
      tile => model%Goto_next_tile()
    END DO

    !> 4.5 Iterate over tiles again and set lct fractions (after ALL tile fractions have been set above)
    !! 
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')


    !> 5.0 Time init for e.g. the time loop in the land_model_driver
    !!
    time_day          = 0.0_wp                                        !! equals 1.0 after one day has passed (i.e. start with zero)
    time_year         = 1                                             !! the 1st year of a simulation == 1 (not zero)
    dtime             = model%config%dtime                            !! time step length in seconds, default: 1800.0 == 30 minutes, mo_jsb_math_constants
    ntimestep         = forcing_options%Timesteps(model%config%dtime) !! number of time steps (considers the simulation_mode static & transient)
    timestep_prop_day = dtime/one_day                                 !! proportion of one day that is simulated with one timestep


    !> 6.0 Init a simulation setup specific for a site-set
    !! 

    !> 6.1 EucFACE transient simulations
    !! using the fluxnet diagnostics-simulations infrastructure
    IF (TRIM(model%config%sim_siteset_specific) == "eucface_trans") THEN
      model%config%fluxnet_type_transient_timestep_output = .TRUE.
      PRINT*, TRIM(routine), ' Init EucFACE transient simulations.'
      PRINT*, TRIM(routine), ' Set fluxnet_type_transient_timestep_output=T; use fluxnet diagnostics infrastructure.'
    ENDIF

    !> 6.2 Fluxnet diagnostics simulation
    !! set the state to "active"
    !! states: \n
    !!  0: not a fluxnet diagnostics simulation
    !!  1: fluxnet diagnostics simulation; spinup or output prior to the "last full forcing period" (static output available)
    !!  2: fluxnet diagnostics simulation: timestep-output over the "last full forcing period"
    !!  3: fluxnet diagnostics simulation finished, first timestep after the "last full forcing period", no output at all
    IF (model%config%fluxnet_type_transient_timestep_output) THEN
      !! activate fluxnet diagnostics output
      flux_trans_timestep_out_state = 1
    ELSE
      flux_trans_timestep_out_state = 0
    ENDIF

  END SUBROUTINE quincy_init

#endif
END MODULE mo_jsb_model_init
