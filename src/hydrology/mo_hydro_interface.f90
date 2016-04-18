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
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
!  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
!  USE mo_util,               ONLY: one_of, int2string

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB, SOIL, VEG, ATM2LAND, RAD, HYDRO

  ! Use process configurations
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  ! Use of process memories
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(HYDRO)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks !,t_hydro_process

  dsl4jsb_Def_memory(RAD)
  dsl4jsb_Def_memory(SOIL)
  dsl4jsb_Def_memory(HYDRO)
  dsl4jsb_Def_memory(VEG)
  dsl4jsb_Def_memory(ATM2LAND)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

!!$  !> Type definition for hydrology process
!!$  TYPE, EXTENDS(t_jsb_process) :: t_hydro_process
!!$  CONTAINS
!!$    PROCEDURE :: Register_tasks => Register_hydro_tasks        !< Registers tasks for hydro process
!!$  END TYPE t_hydro_process

  !> Type definition for snow_hydrology task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_hydrology
  END TYPE tsk_snow_hydrology

  !> Constructor interface for snow_hydrology task
  INTERFACE tsk_snow_hydrology
    PROCEDURE Create_task_snow_hydrology
  END INTERFACE tsk_snow_hydrology

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
  !!
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snow_hydrology"
  !!
  FUNCTION Create_task_snow_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_hydrology::return_ptr)

    return_ptr%name = 'snow_hydrology'
    return_ptr%process_id = HYDRO
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_snow_hydrology
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
  SUBROUTINE Register_hydro_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_snow_hydrology         (model_id))
    CALL this%Register_task(tsk_snow_and_skin_fraction (model_id))

  END SUBROUTINE Register_Hydro_Tasks
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "snow_hydrology"
  !! Task "update_snow_hydrology" calculates the new snow depth and fraction
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_hydrology(tile, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snow_hydrology'

    IF (.NOT. tile%Is_process_active(HYDRO)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      CALL update_surface_energy_lake(tile, options)
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_hydrology

  SUBROUTINE aggregate_snow_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snow_hydrology'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_hydrology
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

    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_phy_schemes,            ONLY: heat_transfer_coef, qsat
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: calc_snow_and_skin_fraction
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

    !-----------------------------------------------------------------------
    ! POINTERS TO VARIABLES IN MEMORY

    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)

    ! Pointers to variables in memory
    !dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_snow_srf
    dsl4jsb_Real2D_onChunk :: fract_water_srf
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_water_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_water_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_ice
    dsl4jsb_Real2D_onChunk :: w_snow_can
    dsl4jsb_Real2D_onChunk :: w_snow_soil
    dsl4jsb_Real2D_onChunk :: w_snow_ice
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

    ! Get reference to variables
    dsl4jsb_Get_config(SOIL)
    wsr_capacity = dsl4jsb_Config(SOIL)%wsr_capacity

    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(HYDRO)

    ! Set pointers
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf)         ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)            ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)        ! IN

    dsl4jsb_Get_var2D_onChunk(SOIL,     oro_stddev)       ! IN
    dsl4jsb_Get_var2D_onChunk(SEB,      t)                ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      lai)              ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      fract_fpc_max)    ! IN

    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_can)       ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_can)       ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_soil)      ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_soil)      ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_ice)       ! IN

    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_srf)  ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_soil) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_can)  ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_srf)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)  ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_can)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_ice)   ! OUT

    CALL calc_snow_and_skin_fraction( &
      & dtime,               & ! in
      & steplen,             & ! in
      & wsr_capacity,        & ! in
      & lai(:),              & ! in
      & oro_stddev(:),       & ! in
      & fract_fpc_max(:),    & ! in
      & t(:),                & ! in
      & press_srf(:),        & ! in
      & drag_srf(:),         & ! in
      & q_air(:),            & ! in
      & w_skin_soil(:),      & ! in
      & w_skin_can(:),       & ! in
      & w_snow_soil(:),      & ! in
      & w_snow_can(:),       & ! in
      & w_snow_ice(:),       & ! in
      & fract_snow_can(:),   & ! inout
      & fract_water_can(:),  & ! out
      & fract_water_soil(:), & ! out
      & fract_water_srf(:),  & ! out
      & fract_snow_soil(:),  & ! out
      & fract_snow_srf(:),   & ! out
      & fract_snow_ice(:)  )   ! out

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

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snow_and_skin_fraction'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !X Implementation: Start your process scheme here...

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_skin_fraction

#endif
END MODULE mo_hydro_interface
