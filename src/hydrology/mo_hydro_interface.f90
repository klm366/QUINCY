!>
!! @brief Contains the interfaces to the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB_, SSE_, TURB_, PHENO_, A2L_, HYDRO_, ASSIMI_, MIX_

  ! Use of process configurations
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(MIX_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks, global_hydrology_diagnostics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

  !> Type definition for surface_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_hydrology
  END TYPE tsk_surface_hydrology

  !> Constructor interface for surface_hydrology
  INTERFACE tsk_surface_hydrology
    PROCEDURE Create_task_surface_hydrology
  END INTERFACE tsk_surface_hydrology

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_properties
  END TYPE tsk_soil_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_properties
    PROCEDURE Create_task_soil_properties
  END INTERFACE tsk_soil_properties

  !> Type definition for soil_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_hydrology
  END TYPE tsk_soil_hydrology

  !> Constructor interface for soil_hydrology
  INTERFACE tsk_soil_hydrology
    PROCEDURE Create_task_soil_hydrology
  END INTERFACE tsk_soil_hydrology

  !> Type definition for evaporation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_evaporation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_evaporation
    PROCEDURE, NOPASS :: Aggregate => aggregate_evaporation
  END TYPE tsk_evaporation

  !> Constructor interface for evaporation task
  INTERFACE tsk_evaporation
    PROCEDURE Create_task_evaporation
  END INTERFACE tsk_evaporation

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for snow_hydrology task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_ice_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_ice_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_ice_hydrology
  END TYPE tsk_snow_and_ice_hydrology

  !> Constructor interface for snow_hydrology task
  INTERFACE tsk_snow_and_ice_hydrology
    PROCEDURE Create_task_snow_and_ice_hydrology
  END INTERFACE tsk_snow_and_ice_hydrology

  !> Type definition for snow_and_skin_fraction task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_skin_fraction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_skin_fraction    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_skin_fraction !< Aggregates computed task variables
  END TYPE tsk_snow_and_skin_fraction

  !> Constructor interface for snow_and_skin_fraction task
  INTERFACE tsk_snow_and_skin_fraction
    PROCEDURE Create_task_snow_and_skin_fraction                       !< Constructor function for task
  END INTERFACE tsk_snow_and_skin_fraction

  !> Type definition for water_balance task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_balance
  END TYPE tsk_water_balance

  !> Constructor interface for water balance task
  INTERFACE tsk_water_balance
    PROCEDURE Create_task_water_balance
  END INTERFACE tsk_water_balance

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_hydrology"
  !!
  FUNCTION Create_task_surface_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_hydrology::return_ptr)
    CALL return_ptr%Construct(name='surface_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_properties"
  !!
  FUNCTION Create_task_soil_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_properties', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_properties

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_hydrology"
  !!
  FUNCTION Create_task_soil_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_hydrology::return_ptr)
    CALL return_ptr%Construct(name='soil_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_unstressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_unstressed"
  !!
  FUNCTION Create_task_canopy_cond_unstressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_unstressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_stress"
  !!
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_stressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_stressed"
  !!
  FUNCTION Create_task_canopy_cond_stressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_stressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for evaporation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "evaporation"
  !!
  FUNCTION Create_task_evaporation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_evaporation::return_ptr)
    CALL return_ptr%Construct(name='evaporation', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_evaporation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snow_hydrology"
  !!
  FUNCTION Create_task_snow_and_ice_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_ice_hydrology::return_ptr)
    CALL return_ptr%Construct(name='snow_ice_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_ice_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow_and_skin_fraction task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_snow_and_skin_fraction"
  !!
  FUNCTION Create_task_snow_and_skin_fraction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_skin_fraction::return_ptr)
    CALL return_ptr%Construct(name='snow_and_skin_fraction', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_balance task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_balance"
  !!
  FUNCTION Create_task_water_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_balance::return_ptr)
    CALL return_ptr%Construct(name='water_balance', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_balance

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for hydrology process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of HYDRO_ process class
  !!
  SUBROUTINE Register_hydro_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task( tsk_surface_hydrology      (model_id))
    CALL process%Register_task( tsk_soil_properties        (model_id))
    CALL process%Register_task( tsk_soil_hydrology         (model_id))
    CALL process%Register_task( tsk_water_stress           (model_id))
    CALL process%Register_task( tsk_canopy_cond_unstressed (model_id))
    CALL process%Register_task( tsk_canopy_cond_stressed   (model_id))
    CALL process%Register_task( tsk_evaporation            (model_id))
    !CALL process%Register_task( tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task( tsk_snow_and_skin_fraction (model_id))
    CALL process%Register_task( tsk_water_balance          (model_id))

  END SUBROUTINE Register_hydro_tasks
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: calc_surface_hydrology_land, calc_surface_hydrology_glacier
    USE mo_jsb_physical_constants, ONLY: dens_snow

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: q_snocpymlt
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: wind_10m
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: wtr_skin
    dsl4jsb_Real2D_onChunk :: fract_wet
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real2D_onChunk :: fract_pond
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_max
    dsl4jsb_Real2D_onChunk :: wtr_pond
    dsl4jsb_Real2D_onChunk :: ice_pond
    dsl4jsb_Real2D_onChunk :: pond_freeze
    dsl4jsb_Real2D_onChunk :: pond_melt
    dsl4jsb_Real2D_onChunk :: weq_snow
    dsl4jsb_Real2D_onChunk :: snow_accum
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: weq_snow_can
    dsl4jsb_Real2D_onChunk :: evapotrans_soil
    dsl4jsb_Real2D_onChunk :: evapo_skin
    dsl4jsb_Real2D_onChunk :: evapo_deficit
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: transpiration
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapo_snow
    dsl4jsb_Real2D_onChunk :: evapo_pond
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real2D_onChunk :: infilt
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: runoff_glac
    dsl4jsb_Real2D_onChunk :: water_to_soil
    dsl4jsb_Real2D_onChunk :: weq_glac


    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & skinres_max,        &
      & skinres_canopy_max, &
      & trans_tmp

    LOGICAL  :: lstart
    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime, w_skin_max_config

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    IF (tile%is_lake) THEN
      ! compute runoff (P-E) for lakes here and exit
      ! Note: for the purpose of the water balance it is assumed here that all rain and snow fall is going into runoff
      ! @todo: compute actual water balance (snow/ice melt) in lake model and use for runoff. The corresponding energy
      ! fluxes and snow/ice budgets are already computed in the lake model.
      dsl4jsb_Get_var2D_onChunk(A2L_,   rain)               ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,   snow)               ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)           ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)         ! out

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_0 tile_name=tile%name
      runoff  (:) = rain(:) + snow(:) + evapotrans(:)
      infilt(:)   = 0._wp
      drainage(:) = 0._wp
      weq_snow(:) = 0._wp
      evapo_snow(:) = evapopot(:)
      !$claw end expand

      RETURN
    END IF

    !$claw acc data pcreate(skinres_max, skinres_canopy_max, trans_tmp)

    dsl4jsb_Get_memory(SSE_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_10m)       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   rain)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)           ! in

    dsl4jsb_Get_var2D_onChunk(SSE_,   hcap_grnd)      ! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens) ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_pond)     ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond_max)   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pond_freeze)    ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pond_melt)      ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_pond)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, q_snocpymlt)    ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)       ! out

    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)       ! in

    IF (tile%is_glacier) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_glac)        ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff_glac)     ! inout
    ELSE
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_to_soil)   ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_wet)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)        ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    snow_accum)      ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans_soil) ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_skin)      ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_deficit)   ! out
      IF (tile%is_vegetation .OR. tile%is_land) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_snow_can)   ! inout
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   transpiration)  ! in
        dsl4jsb_Get_var2D_onChunk(PHENO_,   lai)            ! in
        dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max)  ! in
      END IF
    END IF

    IF (is_time_experiment_start(options%current_datetime)) THEN          ! Start of experiment
      lstart = .TRUE.
      weq_snow(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    w_skin_max_config = dsl4jsb_Config(HYDRO_)%w_skin_max
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_1 tile_name=tile%name
      skinres_canopy_max(:) = w_skin_max_config * lai(:) * fract_fpc_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_2 tile_name=tile%name
      skinres_canopy_max(:) = 0._wp
    END IF
    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_3 tile_name=tile%name
      skinres_max(:) =  w_skin_max_config + skinres_canopy_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_4 tile_name=tile%name
      skinres_max(:) = 0._wp
    END IF

    ! weq_snow_soil_old(:) = weq_snow_soil(:) ! Remember old snow on ground for use in snow model

    IF (tile%is_glacier) THEN
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_surface_hydrology_glacier_1 tile_name=tile%name
      CALL calc_surface_hydrology_glacier( &
        & lstart,                          & !< in
        & dtime,                           & !< in
        & t_unfilt(:),                     & !< in
        & fract_snow(:),                   & !< in
        & hcap_grnd(:),                    & !< in
        & evapotrans(:),                   & !< in
        & evapopot(:),                     & !< in
        & rain(:),                         & !< in
        & snow(:),                         & !< in
        & weq_glac(:),                     & !< inout
        & q_snocpymlt(:),                  & !< out
        & snowmelt(:),                     & !< out
        & runoff_glac(:),                  & !< out
        & pme_glacier = runoff(:)          & !< out, P-E ... used as runoff for HD model
        & )
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_5 tile_name=tile%name
      infilt(:)   = 0._wp
      drainage(:) = 0._wp
      weq_snow_soil(:) = 0._wp
      evapo_snow(:) = evapopot(:)
      !$claw end expand
      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_6 tile_name=tile%name
        snow_soil_dens(:) = 330._wp
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_7 tile_name=tile%name
        snow_soil_dens(:) = dens_snow
      END IF
    ELSE
      IF (tile%is_vegetation) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_8 tile_name=tile%name
        trans_tmp(:) = transpiration(:)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_9 tile_name=tile%name
        trans_tmp(:) = 0._wp
      END IF
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_surface_hydrology_land_1 tile_name=tile%name
      CALL calc_surface_hydrology_land(   &
        & lstart,                         & !< in
        & dtime,                          & !< in
        & dsl4jsb_Config(SSE_)%l_dynsnow, & !< in
        & dsl4jsb_Config(HYDRO_)%snow_depth_max, & !< in
        & t_unfilt(:),                    & !< in
        & wind_10m(:),                    & !< in
        & t_air(:),                       & !< in
        & skinres_canopy_max(:),          & !< in
        & skinres_max(:),                 & !< in
        & weq_pond_max(:),                & !< in
        & fract_snow(:),                  & !< in
        & fract_wet(:),                   & !< in
        & hcap_grnd(:),                   & !< in
        & evapotrans(:),                  & !< in
        & evapopot(:),                    & !< in
        & trans_tmp(:),                   & !< in
        & rain(:),                        & !< in
        & snow(:),                        & !< in
        & wtr_skin(:),                    & !< inout
        & weq_snow_soil(:),               & !< inout
        & weq_snow_can(:),                & !< inout
        & snow_soil_dens(:),              & !< inout
        & fract_pond(:),                  & !< inout
        & weq_pond(:),                    & !< inout
        & wtr_pond(:),                    & !< inout
        & ice_pond(:),                    & !< inout
        & q_snocpymlt(:),                 & !< out
        & snow_accum(:),                  & !< out
        & snowmelt(:),                    & !< out
        & pond_freeze(:),                 & !< out
        & pond_melt(:),                   & !< out
        & evapotrans_soil(:),             & !< out
        & evapo_skin(:),                  & !< out
        & evapo_snow(:),                  & !< out
        & evapo_pond(:),                  & !< out
        & water_to_soil(:),               & !< out
        & evapo_deficit(:)                & !< out
        & )



      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_10 tile_name=tile%name
        snow_soil_dens(:) = 330._wp
      END IF

    END IF

    IF (tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_11 tile_name=tile%name
      weq_snow(:) = 0._wp
    ELSE IF (tile%is_bare) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_12 tile_name=tile%name
      weq_snow(:) = weq_snow_soil(:)
    ELSE IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_surface_hydrology_expand_13 tile_name=tile%name
      weq_snow(:) = weq_snow_soil(:) + weq_snow_can(:)
    END IF

    !$claw acc end data

  END SUBROUTINE update_surface_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options


    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_skin,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_can,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snowmelt,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, q_snocpymlt,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_glac,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_accum,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_soil, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_skin,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_snow,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_pond,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_to_soil,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_deficit,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_pond,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond_max,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_freeze,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_melt,        weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_hydrology
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_properties(tile, options)

    USE mo_hydro_constants, ONLY: vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matrix_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_porosity_org_below, &
      & hyd_cond_sat_org_below, bclapp_org_below, matrix_pot_org_below, pore_size_index_org_below, &
      & vol_field_cap_org_below, vol_p_wilt_org_below, thresh_org, beta_perc

    USE mo_hydro_util,         ONLY: get_amount_in_rootzone


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & hyd_cond_sat,         &
      & vol_porosity,         &
      & bclapp,               &
      & matrix_pot,           &
      & pore_size_index,      &
      & vol_field_cap,        &
      & vol_p_wilt,           &
      & weq_rootzone_max


    dsl4jsb_Real3D_onChunk :: &
!pdv >
      & soil_depth_sl,        &
      & root_depth_sl,        &
!pdv <
      & fract_org_sl,         &
      & hyd_cond_sat_sl,      &
      & vol_porosity_sl,      &
      & bclapp_sl,            &
      & matrix_pot_sl,        &
      & pore_size_index_sl,   &
      & vol_field_cap_sl,     &
      & vol_p_wilt_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    INTEGER  :: iblk, ics, ice, i
    INTEGER  :: nsoil

    REAL(wp) :: N_perc
    REAL(wp), ALLOCATABLE ::  fract_perc(:,:)          ! fraction of the soil with connected organic pathways
    REAL(wp), ALLOCATABLE ::  fract_uncon(:,:)         ! fraction of the soil with unconnected organic matter
    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_uncon(:,:)  ! saturated hydraulic conductivity of 'unconnected' fraction

    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_porosity_org(:,:)
    REAL(wp), ALLOCATABLE ::  bclapp_org(:,:)
    REAL(wp), ALLOCATABLE ::  matrix_pot_org(:,:)
    REAL(wp), ALLOCATABLE ::  pore_size_index_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_field_cap_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_p_wilt_org(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%is_lake .OR. tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)

!pdv >
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    root_depth_sl)          ! in
!pdv <



    dsl4jsb_Get_var2D_onChunk(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, matrix_pot)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_p_wilt)         ! in

    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
      &  dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)      ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_p_wilt_sl)      ! out

!pdv >
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_rootzone_max)   ! out
!pdv <

    DO i=1,nsoil
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_soil_properties_expand_0 tile_name=tile%name soil=i
      hyd_cond_sat_sl   (:,i) = hyd_cond_sat(:)
      vol_porosity_sl   (:,i) = vol_porosity(:)
      bclapp_sl         (:,i) = bclapp(:)
      matrix_pot_sl     (:,i) = matrix_pot(:)
      pore_size_index_sl(:,i) = pore_size_index(:)
      vol_field_cap_sl  (:,i) = vol_field_cap(:)
      vol_p_wilt_sl     (:,i) = vol_p_wilt(:)
      !$claw end expand
    END DO

    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      ALLOCATE(hyd_cond_sat_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_porosity_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(bclapp_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(matrix_pot_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(pore_size_index_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_field_cap_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_p_wilt_org(SIZE(hyd_cond_sat),nsoil))

      ! Use different parameters for top layer and deeper layers
      hyd_cond_sat_org(:,1) = hyd_cond_sat_org_top
      vol_porosity_org(:,1) = vol_porosity_org_top
      bclapp_org(:,1) = bclapp_org_top
      matrix_pot_org(:,1) = matrix_pot_org_top
      pore_size_index_org(:,1) = pore_size_index_org_top
      vol_field_cap_org(:,1) = vol_field_cap_org_top
      vol_p_wilt_org(:,1) = vol_p_wilt_org_top

      DO i=2,nsoil
        hyd_cond_sat_org(:,i) = hyd_cond_sat_org_below
        vol_porosity_org(:,i) = vol_porosity_org_below
        bclapp_org(:,i) = bclapp_org_below
        matrix_pot_org(:,i) = matrix_pot_org_below
        pore_size_index_org(:,i) = pore_size_index_org_below
        vol_field_cap_org(:,i) = vol_field_cap_org_below
        vol_p_wilt_org(:,i) = vol_p_wilt_org_below
      END DO

      ! Saturated hydraulic conductivity following percolation theory 
      ! see CLM45 Tech Note, section 7.4.1 Hydraulic Properties, p.161 f.
      ! @todo: add CLAW directives

      N_perc = (1._wp - thresh_org)**(-beta_perc)

      ALLOCATE(fract_perc(SIZE(hyd_cond_sat),nsoil))          ! fraction of the soil with connected organic pathways
      ALLOCATE(fract_uncon(SIZE(hyd_cond_sat),nsoil))         ! fraction of the soil with unconnected organic matter
      ALLOCATE(hyd_cond_sat_uncon(SIZE(hyd_cond_sat),nsoil))  ! saturated hydraulic conductivity of 'unconnected' fraction

      WHERE(fract_org_sl(:,:) > thresh_org)
        ! Connected flow pathways consisting of organic material only exist
        fract_perc(:,:)  = MIN(1._wp, N_perc * fract_org_sl(:,:) * (fract_org_sl(:,:)-thresh_org)**beta_perc)
        fract_uncon(:,:) = 1._wp - fract_perc(:,:)
      ELSEWHERE
        ! No connected organic matter pathways exist, and flow passes mineral and organic soil components in series
        fract_perc(:,:)  = 0._wp
        fract_uncon(:,:) = 1._wp
      END WHERE

      WHERE(fract_org_sl(:,:) > 0._wp)
        hyd_cond_sat_uncon(:,:) = fract_uncon(:,:)                                           &
          &                       * ( (1._wp - fract_org_sl(:,:)) / hyd_cond_sat_sl(:,:)     &
          &                            + (fract_org_sl(:,:) - fract_perc(:,:)) / hyd_cond_sat_org(:,:) )**(-1._wp)
        hyd_cond_sat_sl(:,:)    = fract_uncon(:,:) * hyd_cond_sat_uncon(:,:)  &
          &                       + fract_perc(:,:) * hyd_cond_sat_org(:,:)
      END WHERE

      DEALLOCATE(fract_perc)
      DEALLOCATE(fract_uncon)
      DEALLOCATE(hyd_cond_sat_uncon)

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_soil_properties_expand_2 tile_name=tile%name soil=i
      vol_porosity_sl   (:,:) = (1._wp - fract_org_sl(:,:)) * vol_porosity_sl   (:,:) &
        &                       + fract_org_sl(:,:) * vol_porosity_org(:,:)
      bclapp_sl         (:,:) = (1._wp - fract_org_sl(:,:)) * bclapp_sl         (:,:) &
        &                       + fract_org_sl(:,:) * bclapp_org(:,:)
      matrix_pot_sl     (:,:) = (1._wp - fract_org_sl(:,:)) * matrix_pot_sl     (:,:) &
        &                       + fract_org_sl(:,:) * matrix_pot_org(:,:)
      pore_size_index_sl(:,:) = (1._wp - fract_org_sl(:,:)) * pore_size_index_sl(:,:) &
        &                       + fract_org_sl(:,:) * pore_size_index_org(:,:)
      vol_field_cap_sl  (:,:) = (1._wp - fract_org_sl(:,:)) * vol_field_cap_sl  (:,:) &
        &                       + fract_org_sl(:,:) * vol_field_cap_org(:,:)
      vol_p_wilt_sl     (:,:) = (1._wp - fract_org_sl(:,:)) * vol_p_wilt_sl     (:,:) &
        &                       + fract_org_sl(:,:) * vol_p_wilt_org(:,:)

      DEALLOCATE(hyd_cond_sat_org)
      DEALLOCATE(vol_porosity_org)
      DEALLOCATE(bclapp_org)
      DEALLOCATE(matrix_pot_org)
      DEALLOCATE(pore_size_index_org)
      DEALLOCATE(vol_field_cap_org)
      DEALLOCATE(vol_p_wilt_org)

!pdv > 
    CALL get_amount_in_rootzone(vol_field_cap_sl(:,:)*soil_depth_sl(:,:),  &
                             &  soil_depth_sl(:,:), root_depth_sl(:,:), weq_rootzone_max(:))
!pdv <

      !$claw end expand
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, hyd_cond_sat_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_porosity_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, bclapp_sl,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, matrix_pot_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pore_size_index_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_field_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_p_wilt_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_org_sl,       weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_properties
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: calc_soil_hydrology
    USE mo_hydro_constants,        ONLY: oro_var_min, oro_var_max

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory

    dsl4jsb_Real3D_onChunk :: fract_org_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl

    dsl4jsb_Real2D_onChunk :: fract_pond
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_max

    dsl4jsb_Real2D_onChunk :: wtr_pond_max

    ! @todo check if this declaration is even necessary
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pwp_sl
    dsl4jsb_Real3D_onChunk :: root_depth_sl
    dsl4jsb_Real2D_onChunk :: tpe_overflow

    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1
    dsl4jsb_Real3D_onChunk :: latdrain_sl
    dsl4jsb_Real3D_onChunk :: drain_lat_sl
    dsl4jsb_Real3D_onChunk :: drain_perm_sl


    dsl4jsb_Real2D_onChunk :: pdv_1d_out
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: micro_slope
    dsl4jsb_Real2D_onChunk :: macro_slope

    ! Locally allocated vectors

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    LOGICAL  :: ltpe_open
    LOGICAL  :: ltpe_closed
! pdv >
    LOGICAL  :: l_latdrain
    LOGICAL  :: l_frzdrain
    LOGICAL  :: l_pf_soil(options%nc)
! pdv <

    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime, max_slope_tile, &
              & sigma_0, sigma_max


    REAL(wp) :: slope(options%nc)

    INTEGER  :: nsoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    ! Runoff and drainage for glaciers already computed in update_surface_hydrology
    ! Glaciers have no water in soil
    IF (tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)

    ltpe_open = .FALSE.
    ltpe_closed = .FALSE.
    IF (model%config%tpe_scheme == 'open')   ltpe_open = .TRUE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    dsl4jsb_Get_var3D_onChunk(SSE_,    t_soil_sl)       ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)   ! in
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$claw acc enter data create(fract_org_sl)
      !$claw expand parallel update
      fract_org_sl(:,:) = 0._wp
    END IF

    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl      )
    dsl4jsb_Get_var3D_onChunk(HYDRO_, latdrain_sl        )
    dsl4jsb_Get_var3D_onChunk(HYDRO_, drain_lat_sl       )
    dsl4jsb_Get_var3D_onChunk(HYDRO_, drain_perm_sl      )
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pdv_1d_out         )
    dsl4jsb_Get_var2D_onChunk(HYDRO_, oro_stddev         )
    dsl4jsb_Get_var2D_onChunk(HYDRO_, micro_slope        )
    dsl4jsb_Get_var2D_onChunk(HYDRO_, macro_slope        )

    pdv_1d_out = 0._wp

    ! @todo Temporary implementation of a surface pond reservoir to be used as target
    !   for lateral flows from upland areas. Needs to be revised as soon as real ponding
    !   is implemented in the model
    dsl4jsb_var2D_onChunk (HYDRO_, water_to_soil) = dsl4jsb_var2D_onChunk (HYDRO_, water_to_soil) &
      & + dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond)
    dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond) = 0._wp

    !pdv >
    ! For testing only
    ! 1.) Pond depth
    ALLOCATE(wtr_pond_max(nc))

    IF (tile%contains_upland) THEN
      max_slope_tile = dsl4jsb_Config(HYDRO_)%max_slope
    ELSE
      max_slope_tile = 0.01_wp
    END IF

    IF (tile%contains_upland) THEN
      sigma_0     = oro_var_min
      sigma_max   = oro_var_max * 64._wp / REAL(grid%nlat_g, wp) ! @todo: site level
      macro_slope(:) = MIN(max_slope_tile, &
                & MAX(0.01_wp, (oro_stddev(:)-sigma_0) / (oro_stddev(:)+ sigma_max)))
    ELSE
      macro_slope(:) = 0.01_wp
    END IF

    slope(:) = MIN(0.5_wp, macro_slope(:) + micro_slope(:))

    IF (tile%contains_depression) THEN
      wtr_pond_max(:) = MAX(0._wp,dsl4jsb_Config(HYDRO_)%max_pond_depth &
                   &  - dsl4jsb_var2D_onChunk (HYDRO_, ice_pond))
    ELSE
      wtr_pond_max(:) = 0.0_wp
    END IF

    ! Latdrain/ frozen drain flags
    ! For testing purposes assuming no lateral drainage from lowland depression
    IF (tile%contains_lowland .AND. tile%contains_depression) THEN
      l_frzdrain = .FALSE.
      l_latdrain = .FALSE.
    ELSE
      l_frzdrain = .TRUE.
      l_latdrain = .TRUE.
    END IF
    !pdv <
    !pdv >
    ! Simple diagnostic for permafrost reaching into the soil above the berock
    ! impeding fluxes into the "groundwater" 
    dsl4jsb_Get_var2D_onChunk(SSE_,    thaw_depth_max_ym1)       ! in
    DO i = 1, nc
      IF (thaw_depth_max_ym1(i) < SUM(soil_depth_sl(i,:)) .AND. &  ! Thaw depth below soil depth ... 
       &  thaw_depth_max_ym1(i) < SUM(soil_w%dz(:)) - soil_w%dz(nsoil)/2._wp - 0.01_wp) THEN ! ... and not at maximum.
        l_pf_soil(i)=.TRUE.
      ELSE
        l_pf_soil(i)=.FALSE.
      END IF
    END DO
    !pdv <


    !$claw sca forward update savepoint=ser_jsb_hydro_calc_soil_hydrology_1 tile_name=tile%name
    CALL calc_soil_hydrology( &
      ! in
      & nc, &
      & nsoil, &
      & grid%nlat_g, & !< (Effective) number of latitudes for steepness parameter in infiltration computation<
      & dtime, &
      & ltpe_closed, ltpe_open,                                 &
      & dsl4jsb_Config(HYDRO_)%l_infil_subzero,                 &
      & dsl4jsb_Config(HYDRO_)%l_soilhydnew,                    &
      & dsl4jsb_Config(HYDRO_)%l_infiltnew,                     &
      & dsl4jsb_Config(HYDRO_)%l_hydro_q,                       &
      & l_latdrain,                                             &
      & l_frzdrain,                                             &
      & dsl4jsb_Config(HYDRO_)%manning_n,                       &
      & dsl4jsb_Config(HYDRO_)%dt_sat,                          &
      & max_slope_tile,                                         &
      & dsl4jsb_Config(HYDRO_)%i_shm,                           &
      & dsl4jsb_Config(HYDRO_)%i_interpol,                      &
      & dsl4jsb_Config(HYDRO_)%i_infil,                         &
      & dsl4jsb_Config(HYDRO_)%i_testing1,                      &
      & dsl4jsb_Config(HYDRO_)%i_testing2,                      &
      & l_pf_soil(:),                                           &
      & dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl     ),     &
      &                                fract_org_sl  (:,:),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, hyd_cond_sat_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, matrix_pot_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, bclapp_sl         ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, pore_size_index_sl),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_porosity_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_field_cap_sl  ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_p_wilt_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pot_scool_sl ), &
      & slope(:),                                               &
      & oro_stddev(:),                                          &
      & t_soil_sl(:,1),                                         &
      & wtr_pond_max(:),                                        &
      & dsl4jsb_var2D_onChunk (HYDRO_, water_to_soil     ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, evapotrans_soil   ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, transpiration     ),     &
      ! inout
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sl       ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, ice_soil_sl       ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond),               &
      & dsl4jsb_var2D_onChunk (HYDRO_, tpe_overflow      ),     &
!pdv >
      & dsl4jsb_var3D_onChunk (HYDRO_, weq_soilhyd_err   ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, pdv_1d_out        ),     &
!pdv <
      ! out
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sat_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_fc_sl    ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pwp_sl   ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, infilt            ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, runoff            ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, runoff_horton     ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, runoff_dunne      ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, drainage          ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, drain_rock        ),     &
      & drain_lat_sl(:,:),     &
      & drain_perm_sl(:,:)     &
      & )

    !pdv >
    ! For testing only
    ! l_macrotrans_subsurf == .false. -> on macroscale, lateral drainage in unfrozen soils
    ! is assumed to reach bedrock boundary before reaching the lower tile.
    ! In frozen soils lateral drainage can not reach bedrock boundary and 
    ! flows to lower (macroscale) tile, however this is only the case if
    ! the permafrost extends upwards into the soil above bedrock

    pdv_1d_out(:) = macro_slope(:)

    dsl4jsb_Get_var2D_onChunk(SSE_,    thaw_depth_max_ym1)       ! in
    IF (tile%contains_upland .AND. tile%contains_depression) THEN
      IF (.NOT. dsl4jsb_Config(HYDRO_)%l_macrotrans_subsurf) THEN
        DO i = 1, nc
          IF (l_pf_soil(i)) THEN ! PF in soil above bedrock
            latdrain_sl(i,:) = drain_lat_sl(i,:) + drain_perm_sl(i,:)
          ELSE
             latdrain_sl(i,:) = 0._wp
          END IF
        END DO
      ELSE
        latdrain_sl(:,:) = drain_lat_sl(:,:) + drain_perm_sl(:,:)
      END IF
    ELSE
      latdrain_sl(:,:) = drain_lat_sl(:,:) + drain_perm_sl(:,:)
    END IF
    !pdv <

    dsl4jsb_var2D_onChunk (HYDRO_, weq_pond           ) = &
    &  dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond        ) + &
    &  dsl4jsb_var2D_onChunk (HYDRO_, ice_pond        )

    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_pond)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond_max)

    DO i = 1, nc
      IF (weq_pond_max(i) > 1.E-20_wp .AND. weq_pond(i) > 1.E-20_wp) THEN
        fract_pond(i) = MIN(1._wp, (weq_pond(i) / weq_pond_max(i))**0.5_wp)
      ELSE
        fract_pond(i) = 0._wp
      END IF
    END DO

    IF (.NOT. dsl4jsb_Config(HYDRO_)%l_organic) THEN
      !$claw acc exit data delete(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_soil_sl,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sat_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_fc_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pwp_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_horton,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_dunne,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drain_rock,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drain_lat_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drain_perm_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_pond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_pond,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond_max,    weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, pdv_1d_out, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_hydrology
  !
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "canopy_cond_unstressed"
  !! Task "update_canopy_cond_unstressed" calculates canopy conductance without water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_unstressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited

    INTEGER :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      swpar_srf_down)        ! IN
    dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                   ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited) ! OUT

    ! Compute (max) canopy (stomatal) conductance with no water stress
    !$claw expand parallel
    canopy_cond_unlimited(:) = get_canopy_conductance( lai(:), swpar_srf_down(:) )

  END SUBROUTINE update_canopy_cond_unstressed
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_unstressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_unlimited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "water_stress"
  !! Task "update_water_stress" calculates water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_water_stress(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process,      ONLY: get_water_stress_factor
    USE mo_hydro_util,         ONLY: get_amount_in_rootzone

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)

    dsl4jsb_Def_memory(HYDRO_)

    ! HYDRO_
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_depth_sl
    dsl4jsb_Real2D_onChunk :: wtr_rootzone           !< Liquid water content in the root zone
    dsl4jsb_Real2D_onChunk :: ice_rootzone           !< Frozen water content in the root zone
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl            !< Liquid water content in soil column
    dsl4jsb_Real3D_onChunk :: ice_soil_sl            !< Frozen water content in soil column
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl  !< Potentially suppercooled water content in soil column
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_pot !< Potentially suppercooled water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_act !< Actual supercooled water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail     !< Plant available liquid water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail_max !< Maximum plant available liquid water content in root zone
    dsl4jsb_Real2D_onChunk :: weq_rootzone_max       !< Maximum possible amount of water/ice in the rootzone
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_rel       !< Relative root zone soil moisture

    INTEGER :: iblk, ics, ice, ic, nc

    REAL(wp) :: config_w_soil_crit_fract, config_w_soil_wilt_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    !
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_pot_scool_sl)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_rootzone_max)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone)           ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    root_depth_sl)          ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    ice_rootzone)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_pot) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_act) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail_max) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_rel)       ! out

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    ! @todo - do we need these?

    config_w_soil_crit_fract = dsl4jsb_Config(HYDRO_)%w_soil_crit_fract
    config_w_soil_wilt_fract = dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract

    ! Calculate actual plant available water
    !   Meant here is all liquid water (without supercooled water) in the root zone. Water
    !   below the wilting point is not subtracted, so calling it available water might be misleading.

    ! Formerly, the root zone moisture was computed as part of the soil hydrology prior to diffusion
    !   and percolation. This seems wrong as the actual water stress should be the result of the
    !   updated soil state and not of the one at the start of the time step.
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:),  &
                             &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      ! Supercooled water is not available to plants
      CALL get_amount_in_rootzone(wtr_soil_pot_scool_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_pot(:))
      CALL get_amount_in_rootzone(MIN(wtr_soil_pot_scool_sl(:,:), wtr_soil_sl(:,:)),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_act(:))
      wtr_rootzone_avail(:) = MAX(wtr_rootzone(:) - wtr_rootzone_scool_act(:), 0._wp)
    ELSE
      ! Available liquid water in the root zone
      wtr_rootzone_avail(:) = wtr_rootzone(:)
    END IF

    ! Calculate maximum plant available water
    !   As in the the available water calculation above, water below the wilting point is not
    !   subtracted.
    ! pdv: Here also the question is how to treat the potential rootzone soil moisture
    !   Definitly ice should be removed - otherwise plants suffer constant waterstress in permafrost-regions:
    !   If supercooled water is asssumed to be unavailable for plants then this should also be accounted for in wtr_rootzone_avail_max

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      wtr_rootzone_avail_max(:) = MAX(weq_rootzone_max(:) - ice_rootzone(:) - wtr_rootzone_scool_pot(:), 0._wp)
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
    ELSE IF (dsl4jsb_Config(SSE_)%l_freeze .AND. .NOT. dsl4jsb_Config(SSE_)%l_supercool) THEN
      wtr_rootzone_avail_max(:) = MAX(weq_rootzone_max(:) - ice_rootzone(:), 0._wp)
    ELSE
      wtr_rootzone_avail_max(:) = weq_rootzone_max(:)
    END IF

    ! Calculate water stress
    ! @todo: maybe use the distributed permanent wilting point/field cap from above and change get_water_stress_factor function
    water_stress(:) = get_water_stress_factor(wtr_rootzone_avail(:), wtr_rootzone_avail_max(:), &
            &                 config_w_soil_crit_fract, config_w_soil_wilt_fract)

    ! Calculate relative soil moisture (used in LoGro-P Phenology)
    !$claw acc data present(wtr_rootzone_avail, wtr_rootzone_avail_max, soil_depth_sl, root_depth_sl, wtr_rootzone_rel)
    !$claw acc data present(wtr_soil_pot_scool_sl, wtr_rootzone_scool_act, wtr_rootzone_avail)
    !$claw acc parallel loop gang vector
    DO ic= 1, nc
      IF (wtr_rootzone_avail_max(ic) > 0._wp) THEN
        wtr_rootzone_rel(ic) =  MIN(wtr_rootzone_avail(ic) / wtr_rootzone_avail_max(ic), 1._wp)
      ELSE
        ! weq_rootzone_max is zero for glacier tiles
        wtr_rootzone_rel(ic) = 0._wp
      END IF
    END DO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "water_stress"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(ASSIMI_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! HYDRO_
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_rel,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail_max, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_stress,           weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "canopy_cond_stressed"
  !! Task "update_canopy_cond_stressed" calculates canopy conductance under water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_stressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance
    USE mo_phy_schemes,   ONLY: qsat

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: iblk, ics, ice, nc
    ! LOGICAL, DIMENSION(options%nc) :: q_air_gt_qsat_tmp
    LOGICAL, ALLOCATABLE :: q_air_gt_qsat_tmp(:)


    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited) ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_limited) ! out

    !$cla$w acc data create(q_air_gt_qsat_tmp)
    ALLOCATE(q_air_gt_qsat_tmp(nc))
    !$claw acc enter data create(q_air_gt_qsat_tmp)
    !$claw expand parallel update savepoint=ser_jsb_hydro_update_canopy_cond_stressed_expand_0 tile_name=tile%name
    q_air_gt_qsat_tmp(:) = q_air(:) > qsat(t(:),press_srf(:))
    ! Compute (actual) canopy (stomatal) conductance under water stress.
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond_unlimited(:),  & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air_gt_qsat_tmp(:)       & ! in, atmosphere saturated?
                                                    )
    !$claw end expand
    !$cla$w acc end data
    !$claw acc exit data delete(q_air_gt_qsat_tmp)
    DEALLOCATE(q_air_gt_qsat_tmp)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_stressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_stressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_limited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "evaporation"
  !! Task "update_evaporation" calculates evapo(transpi)ration
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_evaporation(tile, options)

    USE mo_phy_schemes,            ONLY: q_effective, heat_transfer_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & q_acoef, &
      & q_bcoef, &
      & q_acoef_wtr, &
      & q_bcoef_wtr, &
      & q_acoef_ice, &
      & q_bcoef_ice, &
      & drag_srf, &
      & drag_wtr, &
      & drag_ice, &
      & qsat_star, &
      & qsat_lwtr, &
      & qsat_lice, &
      & fact_qsat_srf, &
      & fact_qsat_trans_srf, &
      & fact_q_air, &
      & fract_lice, &
      & evapotrans_lnd, &
      & evapo_wtr, &
      & evapo_ice, &
      & evapopot,       &
      & evapotrans, &
      & transpiration

    ! Locally allocated vectors
    !
    ! REAL(wp) ::                  &
    !   & q_air        (options%nc), &  !< Humidity at lowest atmospheric level [kg kg-1]
    !   & q_air_eff    (options%nc), &  !< Effective humidity at lowest atmospheric level
    !   & qsat_srf_eff (options%nc), &  !< Effective surface saturation specific humidity
    !   & heat_tcoef   (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)
    REAL(wp), ALLOCATABLE ::                  &
      & q_air        (:), &  !< Humidity at lowest atmospheric level [kg kg-1]
      & q_air_eff    (:), &  !< Effective humidity at lowest atmospheric level
      & qsat_srf_eff (:), &  !< Effective surface saturation specific humidity
      & heat_tcoef   (:)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_evaporation'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !$cla$w acc data create(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
    ALLOCATE(q_air(nc),q_air_eff(nc),qsat_srf_eff(nc),heat_tcoef(nc))
    !$claw acc enter data create(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans) ! out

    IF (tile%contains_lake) THEN

      dsl4jsb_Get_var2D_onChunk(A2L_,     q_acoef_wtr)   ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,     q_bcoef_wtr)   ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,     drag_wtr)      ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,     qsat_lwtr)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_wtr)     ! out

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_1 tile_name=tile%name
      q_air(:) = q_acoef_wtr(:) * qsat_lwtr(:) + q_bcoef_wtr(:) ! New moisture at lowest atmospheric level by back-substitution
      heat_tcoef(:) = heat_transfer_coef(drag_wtr(:), steplen)  ! Transfer coefficient
      evapo_wtr (:) = heat_tcoef(:) * (q_air(:) - qsat_lwtr(:)) ! Potential evaporation
      !$claw end expand

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(A2L_,   q_acoef_ice)       ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,   q_bcoef_ice)       ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,   drag_ice)          ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   qsat_lice)         ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)        ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)         ! out

        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_2 tile_name=tile%name
        q_air(:) = q_acoef_ice(:) * qsat_lice(:) + q_bcoef_ice(:) ! New moisture at lowest atmospheric level by back-substitution
        heat_tcoef(:) = heat_transfer_coef(drag_ice(:), steplen)  ! Transfer coefficient
        evapo_ice (:) = heat_tcoef(:) * (q_air(:) - qsat_lice(:)) ! Potential evaporation
        evapopot  (:) = (1._wp - fract_lice(:)) * evapo_wtr   (:) + fract_lice(:) * evapo_ice   (:)
        !$claw end expand
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_3 tile_name=tile%name
        evapopot  (:) = evapo_wtr(:)
      END IF
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_4 tile_name=tile%name
      evapotrans(:) = evapopot(:)

    ELSE IF (tile%contains_soil .OR. tile%contains_glacier) THEN

      dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef)             ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef)             ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,  drag_srf)            ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,  qsat_star)           ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_q_air)          ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_srf)       ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_trans_srf) ! in

      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans_lnd)           ! out
      IF (tile%contains_vegetation) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    transpiration)          ! out
      END IF

      !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_5 tile_name=tile%name
      q_air(:) = q_acoef(:) * qsat_star(:) + q_bcoef(:) ! New moisture at lowest atmospheric level by back-substitution
      heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)  ! Transfer coefficient

      ! Compute effective air moisture and surface saturation humidity
      q_air_eff(:) = q_effective( SPREAD(0._wp, NCOPIES=nc, DIM=1), q_air(:), &
                              & SPREAD(1._wp, NCOPIES=nc, DIM=1), SPREAD(0._wp, NCOPIES=nc, DIM=1))
      qsat_srf_eff(:) = q_effective(qsat_star(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))

      evapotrans_lnd(:)          = heat_tcoef(:) * (q_air_eff(:) - qsat_srf_eff(:)) ! Evapotranspiration
      evapotrans(:) = evapotrans_lnd(:)
      evapopot(:)           = heat_tcoef(:) * (q_air    (:) - qsat_star   (:)) ! Potential evaporation
      !$claw end expand
      IF (tile%contains_vegetation) THEN
        !$claw expand parallel update savepoint=ser_jsb_hydro_update_evaporation_expand_6 tile_name=tile%name
        transpiration(:)    = fact_qsat_trans_srf(:) * evapopot(:)             ! Transpiration
      END IF

    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

    !$cla$w acc end data
    !$claw acc exit data delete(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
    DEALLOCATE(q_air,q_air_eff,qsat_srf_eff,heat_tcoef)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_evaporation
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "evaporation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_evaporation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_evaporation'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapopot,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, transpiration,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_lnd,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_wtr,          weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_ice,          weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_evaporation
  !
  !>
  !! @brief Implementation of "update" for task "snow_and_ice_hydrology"
  !! Task "update_snow_hydrology" calculates the new snow depth and fraction
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_ice_hydrology'

    IF (options%nc > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%contains_lake) THEN
      !! Currently done in SEB process
    ELSE IF (tile%contains_land) THEN
!!$      CALL update_snow_and_ice_hydrology_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_and_ice_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_ice_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_ice_hydrology'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_ice_hydrology
  !
  ! ================================================================================================================================
  !
  !! @brief Implementation of "update" for task "snow_and_skin_fraction"
  !! Task "update_snow_and_skin_fraction" calculates the new snow and wet skin fractions
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_skin_fraction( tile, options)

    ! Use declarations
    USE mo_hydro_process,          ONLY: calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare
    ! sollte jetzt in dieses File hier kommen: USE mo_hydro_process,  ONLY: calc_hydro_snow_and_skin_fraction

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    REAL(wp), DIMENSION(options%nc) :: skinres_canopy_max, skinres_max

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_skin_fraction'

    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime, steplen, config_w_skin_max
    !$claw verbatim INTEGER, SAVE :: task_counter = 1

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: fract_wet
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_lice
    dsl4jsb_Real2D_onChunk :: weq_snow_can
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: weq_snow_lice
    dsl4jsb_Real2D_onChunk :: wtr_skin
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air

    ! ---------------------------
    ! Go

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime
    steplen = options%steplen

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(HYDRO_)

    ! Use simple scalar for GPU
    config_w_skin_max = dsl4jsb_Config(HYDRO_)%w_skin_max

    !$claw acc enter data pcreate(skinres_max, skinres_canopy_max)
    !$claw expand parallel update
    skinres_max(:) = config_w_skin_max

    dsl4jsb_Get_memory(HYDRO_)

    ! First handle LAKE_TYPE lct
    IF (tile%is_lake) THEN
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_lice)    ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_lice)  ! out
        !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_lice_1   &
        !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
        CALL calc_wskin_fractions_lice( &
          & weq_snow_lice(:),           & ! in
          & fract_snow_lice(:)          & ! out
          & )
      END IF
      !$claw acc exit data delete(skinres_max, skinres_canopy_max)
      RETURN
    END IF

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    ! Set pointers
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)         ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)        ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                ! in

    IF (tile%is_bare .OR. tile%is_vegetation .OR. tile%is_land) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    oro_stddev)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_wet)         ! out
    END IF
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)               ! in
      dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)    ! out
    END IF

    IF (tile%is_vegetation .OR. tile%is_land .OR. tile%is_bare) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)          ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_soil)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)  ! out

    ! Maximum capacity of skin reservoir (soil + canopy)
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_0 &
      !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      skinres_canopy_max(:) = config_w_skin_max * lai(:) * fract_fpc_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_1 &
      !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      skinres_canopy_max(:) = 0._wp
    END IF
    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_2 &
      !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      skinres_max(:) = config_w_skin_max + skinres_canopy_max(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_3 &
      !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      skinres_max(:) = 0._wp
    END IF

    IF (tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_snow_and_skin_fraction_expand_4 &
      !$claw tile_name=tile%name dtime=dtime steplen=steplen task_counter=task_counter
      fract_snow(:)      = 1._wp
      fract_snow_soil(:) = 1._wp
      !$claw end expand
    END IF

    IF (tile%is_vegetation .OR. tile%is_land) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_can)       ! in
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_veg_1 tile_name=tile%name &
      !$claw dtime=dtime steplen=steplen task_counter=task_counter
      CALL calc_wskin_fractions_veg( &
        & dtime,                     & ! in
        & steplen,                   & ! in
        & skinres_max(:),            & ! in
        & oro_stddev(:),             & ! in
        & t(:),                      & ! in, from the previous time step as long as this is called before the asselin filter
        & press_srf(:),              & ! in
        & drag_srf(:),               & ! in
        & q_air(:),                  & ! in
        & wtr_skin(:),               & ! in
        & weq_snow_soil(:),          & ! in
        & weq_snow_can(:),           & ! in
        & fract_snow_can(:),         & ! out
        & fract_wet(:),              & ! out
        & fract_snow_soil(:)         & ! out
        & )
      !$claw expand parallel update
      fract_snow(:) = fract_snow_soil(:)
    END IF


    IF (tile%is_bare) THEN
      !$claw sca forward update savepoint=ser_jsb_hydro_calc_wskin_fractions_bare_1 tile_name=tile%name &
      !$claw dtime=dtime steplen=steplen  task_counter=task_counter
      CALL calc_wskin_fractions_bare( &
        & dtime,                      & ! in
        & steplen,                    & ! in
        & skinres_max,                & ! in
        & oro_stddev(:),              & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & drag_srf(:),                & ! in
        & q_air(:),                   & ! in
        & wtr_skin(:),                & ! in
        & weq_snow_soil(:),           & ! in
        & fract_wet(:),               & ! out
        & fract_snow_soil(:)          & ! out
        & )
      !$claw expand parallel update
      fract_snow(:)  = fract_snow_soil(:)
   END IF

   !$claw acc exit data delete(skinres_max, skinres_canopy_max)
   !$claw verbatim  task_counter = task_counter + 1
  END SUBROUTINE update_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_skin_fraction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_skin_fraction(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_skin_fraction'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_wet,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_can,   weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_lice,  weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_skin_fraction

  SUBROUTINE update_water_balance( tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(MIX_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rain, &
      & snow, &
      & runoff, &
      & drainage, &
      & evapotrans, &
      & wtr_latflow, &
      & wtr_skin, &
      & weq_snow, &
      & runoff_retention, &
      & wtr_pond, &
      & ice_pond, &
      & weq_pond, &
      & weq_fluxes, &    ! [m3 s-1]
      & weq_land, &      ! [m3]
      & weq_balance_err  ! [m3/(time step)]
    dsl4jsb_Real3D_onChunk :: &
      & wtr_soil_sl, &
      & ice_soil_sl, &
      & latdrain_retention

    REAL(wp), POINTER :: &
      & area(:)

    ! Locally allocated vectors
    !
    REAL(wp) :: &
      & tile_fract (options%nc), &
      & weq_land_new(options%nc)

    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: dtime

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_balance'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    ! IF (tile%lcts() RETURN   ! Check water balance only on root tile
    ! IF (tile%is_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    area => grid%area (ics:ice, iblk)

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)
    IF (TRIM(model%config%usecase) == 'jsbach_qarctic') THEN
      ! MIX_ Process currently only exists for qarctic usecase
      dsl4jsb_Get_memory(MIX_)
    END IF

    dsl4jsb_Get_var2D_onChunk(A2L_,      rain)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      snow)                ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)            ! in
    IF (TRIM(model%config%usecase) == 'jsbach_qarctic') THEN
      dsl4jsb_Get_var2D_onChunk(MIX_,      wtr_latflow)       ! in
    ELSE
      !$claw acc enter data create(wtr_latflow)
      !$claw expand parallel update
      ALLOCATE(wtr_latflow(nc))
      wtr_latflow(:) = 0._wp
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff_retention)       ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, latdrain_retention)     ! in
    IF (tile%contains_soil .AND. .NOT. tile%is_glacier) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow)            ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_pond)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    ice_pond)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans)          ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond)            ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_fluxes)          ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_land)            ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_balance_err)     ! out

    !$claw acc data create(tile, tile%fract, tile_fract, weq_land_new)

    !$claw acc parallel loop gang vector
    DO i=1,nc
      tile_fract(i) = tile%fract(ics-1+i,iblk)
    END DO

    weq_pond(:) = wtr_pond(:) + ice_pond(:)
    !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_0 tile_name=tile%name
    weq_fluxes(:) = rain(:) + snow(:) + evapotrans(:) - runoff(:) - drainage(:) + wtr_latflow(:)
    weq_fluxes(:) = weq_fluxes(:) * tile_fract(:) * area(:) / rhoh2o                    ! kg m-2 s-1 -> m3/s
    !$claw end expand

    IF (tile%contains_soil .AND. .NOT. tile%is_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_1 tile_name=tile%name
      weq_land_new(:) = (  wtr_skin(:) + weq_snow(:) + weq_pond(:) +  runoff_retention(:) &
        &                + SUM(wtr_soil_sl(:,:), DIM=2) + SUM(ice_soil_sl(:,:), DIM=2)  &
        &                + SUM(latdrain_retention(:,:), DIM=2)                          &
        &               ) * tile_fract(:) * area(:)  ! m -> m^3
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_2 tile_name=tile%name
      weq_land_new(:) = 0._wp
    END IF

    IF (.NOT. is_time_experiment_start(options%current_datetime)) THEN
      !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_3 tile_name=tile%name
      weq_balance_err(:) = weq_land(:) + weq_fluxes(:) * dtime - weq_land_new(:)
    END IF

    !$claw expand parallel update savepoint=ser_jsb_hydro_update_water_balance_expand_4 tile_name=tile%name
    weq_land(:) = weq_land_new(:)

    IF (TRIM(model%config%usecase) /= 'jsbach_qarctic') THEN
      !$claw acc exit data delete(wtr_latflow)
      DEALLOCATE(wtr_latflow)
    END IF

    !$claw acc end data

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_balance

  SUBROUTINE aggregate_water_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !dsl4jsb_Def_memory(HYDRO_)

    !CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_balance'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !dsl4jsb_Get_memory(HYDRO_)

    ! weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! Don't aggregate, but explicitely compute water balance on each tile
    CALL update_water_balance(tile, options)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_fluxes,       weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_land,         weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_balance_err,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_balance

  !-----------------------------------------------------------------------------------------------------
  !> Global land mean hydrology output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_hydrology_diagnostics(tile)
#ifndef __ICON__
    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_hydrology_diagnostics'
    IF (debug_on()) CALL message(TRIM(routine), 'Global diagnostics only available with ICON')
#else

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_hydrology_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: transpiration
    dsl4jsb_Real2D_onDomain :: evapotrans
    dsl4jsb_Real2D_onDomain :: weq_land
    dsl4jsb_Real2D_onDomain :: discharge_ocean
    dsl4jsb_Real2D_onDomain :: wtr_rootzone_rel
    dsl4jsb_Real2D_onDomain :: fract_snow
    dsl4jsb_Real2D_onDomain :: weq_snow
    dsl4jsb_Real2D_onDomain :: weq_balance_err

    REAL(wp), POINTER       :: trans_gmean(:)
    REAL(wp), POINTER       :: evapotrans_gmean(:)
    REAL(wp), POINTER       :: weq_land_gsum(:)
    REAL(wp), POINTER       :: discharge_ocean_gsum(:)
    REAL(wp), POINTER       :: wtr_rootzone_rel_gmean(:)
    REAL(wp), POINTER       :: fract_snow_gsum(:)
    REAL(wp), POINTER       :: weq_snow_gsum(:)
    REAL(wp), POINTER       :: weq_balance_err_gsum(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp)               :: global_land_area


    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  transpiration)              ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  evapotrans)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_land)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  discharge_ocean)            ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  wtr_rootzone_rel)           ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  fract_snow)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_snow)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_balance_err)            ! in


    trans_gmean          => HYDRO__mem%trans_gmean%ptr(:)           ! out
    evapotrans_gmean     => HYDRO__mem%evapotrans_gmean%ptr(:)      ! out
    weq_land_gsum        => HYDRO__mem%weq_land_gsum%ptr(:)         ! out
    discharge_ocean_gsum => HYDRO__mem%discharge_ocean_gsum%ptr(:)  ! out
    wtr_rootzone_rel_gmean => HYDRO__mem%wtr_rootzone_rel_gmean%ptr(:)  ! out
    fract_snow_gsum      => HYDRO__mem%fract_snow_gsum%ptr(:)       ! out
    weq_snow_gsum        => HYDRO__mem%weq_snow_gsum%ptr(:)         ! out
    weq_balance_err_gsum => HYDRO__mem%weq_balance_err_gsum%ptr(:)  ! out


    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea


    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')


    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these
    ! cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))

    ! Calculate 1d global land variables, if requested for output
    global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
    scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (HYDRO__mem%trans_gmean%is_in_output)           &
      &  trans_gmean          = global_sum_array(transpiration(:,:)   * scaling(:,:)) / global_land_area
    IF (HYDRO__mem%evapotrans_gmean%is_in_output)      &
      &  evapotrans_gmean     = global_sum_array(evapotrans(:,:)      * scaling(:,:)) / global_land_area
    ! Unit transformation from [m3] to [km3]: 1.e-9
    IF (HYDRO__mem%weq_land_gsum%is_in_output)           &
      &  weq_land_gsum        = global_sum_array(weq_land(:,:)  * notsea(:,:) * in_domain(:,:)) * 1.e-9_wp
    ! Unit transformation from [m3/s] to [Sv] (1 Sv = 1.e6 m3/s)
    IF (HYDRO__mem%discharge_ocean_gsum%is_in_output)  &
      &  discharge_ocean_gsum = global_sum_array(discharge_ocean(:,:) * in_domain(:,:)) * 1.e-6_wp
    IF (HYDRO__mem%wtr_rootzone_rel_gmean%is_in_output)  &
      &  wtr_rootzone_rel_gmean = global_sum_array(wtr_rootzone_rel(:,:) * scaling(:,:)) / global_land_area
    ! Unit transformation from [m2] to [Mio km2]: 1.e-12
    IF (HYDRO__mem%fract_snow_gsum%is_in_output)       &
      &  fract_snow_gsum      = global_sum_array(fract_snow(:,:)      * scaling(:,:)) * 1.e-12_wp
    ! Unit transformation from [m water equivalent](= [t]) to [Gt]: 1.e-9
    IF (HYDRO__mem%weq_snow_gsum%is_in_output)           &
      &  weq_snow_gsum        = global_sum_array(weq_snow(:,:)        * scaling(:,:)) * 1.e-9_wp
    ! No unit transformation [m3]: 1
    IF (HYDRO__mem%weq_balance_err_gsum%is_in_output)     &
      &  weq_balance_err_gsum = global_sum_array(weq_balance_err(:,:) * notsea(:,:) * in_domain(:,:))

    DEALLOCATE (scaling, in_domain)
#endif
  END SUBROUTINE global_hydrology_diagnostics

#endif
END MODULE mo_hydro_interface
