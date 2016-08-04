!>
!! @brief Contains the interfaces to the hydro processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!

MODULE mo_hydro_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, VEG_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
!  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
!  USE mo_util,               ONLY: one_of, int2string

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB, SOIL, VEG, PHENO, ATM2LAND, RAD, HYDRO

  ! Use of process configurations
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(HYDRO)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(PHENO)
  dsl4jsb_Use_memory(RAD)
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks !,t_hydro_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

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

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

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

    return_ptr%name = 'snow_ice_hydrology'
    return_ptr%process_id = HYDRO
    return_ptr%owner_model_id = model_id

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

    return_ptr%name = 'snow_and_skin_fraction'
    return_ptr%process_id = HYDRO
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for hydrology process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of HYDRO process class
  !!
  SUBROUTINE Register_hydro_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    !CALL process%Register_task(tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task(tsk_snow_and_skin_fraction (model_id))

  END SUBROUTINE Register_hydro_tasks
  
  ! ================================================================================================================================
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

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snow_and_ice_hydrology'

    IF (.NOT. tile%Is_process_active(HYDRO)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      !! Currently done in SEB process
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
!!$      CALL update_snow_and_ice_hydrology_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_and_ice_hydrology

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

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snow_and_ice_hydrology'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

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
    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: calc_wskin_fractions_lice, calc_wskin_fractions_veg, calc_wskin_fractions_bare
    ! sollte jetzt in dieses File hier kommen: USE mo_hydro_process,  ONLY: calc_hydro_snow_and_skin_fraction

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    REAL(wp) :: wsr_capacity

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snow_and_skin_fraction'

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime, steplen

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)
    dsl4jsb_Def_memory(SOIL)

    ! Declare pointers for variables in memory
    !dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: fract_water
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_snow_veg
    dsl4jsb_Real2D_onChunk :: fract_water_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_water_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_lice
    dsl4jsb_Real2D_onChunk :: w_snow_can
    dsl4jsb_Real2D_onChunk :: w_snow_soil
    dsl4jsb_Real2D_onChunk :: w_snow_lice
    dsl4jsb_Real2D_onChunk :: w_skin_can
    dsl4jsb_Real2D_onChunk :: w_skin_soil
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
    nc      = options%nc
    dtime   = options%dtime
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(HYDRO)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(SOIL)   
    
    wsr_capacity = dsl4jsb_Config(SOIL)%wsr_capacity

    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(PHENO)
    dsl4jsb_Get_memory(HYDRO)
    dsl4jsb_Get_memory(SOIL) 

    ! Set pointers
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf)         ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)            ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)        ! in

    dsl4jsb_Get_var2D_onChunk(SOIL,     oro_stddev)       ! in
    dsl4jsb_Get_var2D_onChunk(SEB,      t)                ! in
    dsl4jsb_Get_var2D_onChunk(PHENO,    lai)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG,      fract_fpc_max)    ! in

    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_soil)      ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_soil)      ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_soil) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_can)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_veg)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_can)   ! out

    IF (tile%lcts(1)%id == LAKE_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_lice)      ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_lice)  ! out
      CALL calc_wskin_fractions_lice( &
        & w_snow_lice(:),             & ! in
        & fract_snow_lice(:)          & ! out
        & )
    END IF

    IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      ! TBD: On glacier: Wet skin fraction = 0, snow fraction =1, snow depth = 0
      !fract_skin_glacier(:) = 0._wp
      !fract_snow(:)     = 1._wp
    END IF

    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_can)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_can)       ! in
      CALL calc_wskin_fractions_veg( &
        & dtime,                     & ! in
        & steplen,                   & ! in
        & wsr_capacity,              & ! in
        & lai(:),                    & ! in
        & oro_stddev(:),             & ! in
        & fract_fpc_max(:),          & ! in
        & t(:),                      & ! in
        & press_srf(:),              & ! in
        & drag_srf(:),               & ! in
        & q_air(:),                  & ! in
        & w_skin_soil(:),            & ! in
        & w_skin_can(:),             & ! in
        & w_snow_soil(:),            & ! in
        & w_snow_can(:),             & ! in
        & fract_snow_can(:),         & ! inout
        & fract_water_can(:),        & ! out
        & fract_water_soil(:),       & ! out
        & fract_water(:),            & ! out
        & fract_snow_soil(:),        & ! out
        & fract_snow_veg(:)          & ! out
        & )
   END IF

    
    IF (tile%lcts(1)%id == BARE_TYPE) THEN
      CALL calc_wskin_fractions_bare( &
        & dtime,                      & ! in
        & steplen,                    & ! in
        & wsr_capacity,               & ! in
        & lai(:),                     & ! in
        & oro_stddev(:),              & ! in
        & fract_fpc_max(:),           & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & drag_srf(:),                & ! in
        & q_air(:),                   & ! in
        & w_skin_soil(:),             & ! in
        & w_snow_soil(:),             & ! in
        & fract_water_soil(:),        & ! out
        & fract_snow_soil(:)          & ! out
        & )
   END IF




      
  END SUBROUTINE update_snow_and_skin_fraction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_skin_fraction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_skin_fraction(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snow_and_skin_fraction'

    INTEGER  :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO, fract_water,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_water_soil, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_water_can,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_snow_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_snow_can,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, fract_snow_lice,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_skin_fraction

#endif
END MODULE mo_hydro_interface
