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

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_grid,           ONLY: Get_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, VEG_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
!  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
!  USE mo_util,               ONLY: one_of, int2string

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB_, TURB_, SOIL_, VEG_, PHENO_, ATM2LAND_, RAD_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Use_config(RAD_)

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SOIL_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks !,t_hydro_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

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
  !> Constructor for evaporation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "evaporation"
  !!
  FUNCTION Create_task_evaporation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_evaporation::return_ptr)

    return_ptr%name = 'evaporation'
    return_ptr%process_id = HYDRO_
    return_ptr%owner_model_id = model_id

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

    return_ptr%name = 'snow_ice_hydrology'
    return_ptr%process_id = HYDRO_
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
    return_ptr%process_id = HYDRO_
    return_ptr%owner_model_id = model_id

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

    return_ptr%name = 'water_balance'
    return_ptr%process_id = HYDRO_
    return_ptr%owner_model_id = model_id

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

    CALL process%Register_task(tsk_evaporation (model_id))
    !CALL process%Register_task(tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task(tsk_snow_and_skin_fraction (model_id))
    CALL process%Register_task(tsk_water_balance (model_id))

  END SUBROUTINE Register_hydro_tasks
  
  ! ================================================================================================================================
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
    dsl4jsb_Def_memory(ATM2LAND_)
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
      & evapo_lnd, &
      & evapo_wtr, &
      & evapo_ice, &
      & evapopot,       &
      & evapotranspiration, &
      & trans

    ! Locally allocated vectors
    !
    REAL(wp) ::                  &
      & q_air        (options%nc), &  !< Humidity at lowest atmospheric level
      & q_air_eff    (options%nc), &  !< Effective humidity at lowest atmospheric level
      & qsat_srf_eff (options%nc), &  !< Effective surface saturation specific humidity
      & heat_tcoef   (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_evaporation'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotranspiration) ! out

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef_wtr)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef_wtr)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef_ice)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef_ice)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_wtr)          ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_ice)          ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_lwtr)         ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_lice)         ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,      fract_lice)        ! in


      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_wtr)          ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_ice)          ! out

      q_air(:) = q_acoef_wtr(:) * qsat_lwtr(:) + q_bcoef_wtr(:) ! New moisture at lowest atmospheric level by back-substitution

      heat_tcoef(:) = heat_transfer_coef(drag_wtr(:), steplen)  ! Transfer coefficient

      evapo_wtr (:) = heat_tcoef(:) * (q_air(:) - qsat_lwtr(:)) ! Potential evaporation

      q_air(:) = q_acoef_ice(:) * qsat_lice(:) + q_bcoef_ice(:) ! New moisture at lowest atmospheric level by back-substitution

      heat_tcoef(:) = heat_transfer_coef(drag_ice(:), steplen)  ! Transfer coefficient

      evapo_ice (:) = heat_tcoef(:) * (q_air(:) - qsat_lice(:)) ! Potential evaporation

      evapopot          (:) = (1._wp - fract_lice(:)) * evapo_wtr   (:) + fract_lice(:) * evapo_ice   (:)
      evapotranspiration(:) = evapopot(:)

    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE) .OR. Contains_lct(tile%lcts, VEG_TYPE) &
      & .OR. Contains_lct(tile%lcts, BARE_TYPE) .OR. Contains_lct(tile%lcts, GLACIER_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef)             ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef)             ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_srf)            ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_star)           ! in
      dsl4jsb_Get_var2D_onChunk(TURB_,     fact_q_air)          ! in
      dsl4jsb_Get_var2D_onChunk(TURB_,     fact_qsat_srf)       ! in
      dsl4jsb_Get_var2D_onChunk(TURB_,     fact_qsat_trans_srf) ! in

      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_lnd)           ! out
      IF (.NOT. Contains_lct(tile%lcts, GLACIER_TYPE)) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    trans)               ! out
      END IF

      q_air(:) = q_acoef(:) * qsat_star(:) + q_bcoef(:) ! New moisture at lowest atmospheric level by back-substitution

      heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)  ! Transfer coefficient

      ! Compute effective air moisture and surface saturation humidity
      q_air_eff(:) = q_effective( SPREAD(0._wp, NCOPIES=nc, DIM=1), q_air(:), &
                              & SPREAD(1._wp, NCOPIES=nc, DIM=1), SPREAD(0._wp, NCOPIES=nc, DIM=1))
      qsat_srf_eff(:) = q_effective(qsat_star(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))

      evapo_lnd(:)          = heat_tcoef(:) * (q_air_eff(:) - qsat_srf_eff(:)) ! Evapotranspiration
      evapotranspiration(:) = evapo_lnd(:)
      evapopot(:)           = heat_tcoef(:) * (q_air    (:) - qsat_star   (:)) ! Potential evaporation
      IF (.NOT. Contains_lct(tile%lcts, GLACIER_TYPE)) THEN
        trans(:)              = fact_qsat_trans_srf(:) * evapopot(:)             ! Transpiration
      END IF

    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

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

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_evaporation'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract", copy=.TRUE.)
    CALL weighted_by_fract%Set_option('l_all_tiles', .TRUE.)

    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotranspiration, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapopot,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, trans,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_lnd,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_wtr,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_ice,          weighted_by_fract)

    DEALLOCATE(weighted_by_fract)

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

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snow_and_ice_hydrology'

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      !! Currently done in SEB process
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
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
    dsl4jsb_Def_memory(ATM2LAND_)
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(SOIL_)

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

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(SOIL_)   
    
    wsr_capacity = dsl4jsb_Config(SOIL_)%wsr_capacity

    dsl4jsb_Get_memory(HYDRO_)

    ! First handle LAKE_TYPE lct
    IF (tile%lcts(1)%id == LAKE_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow_lice)      ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_lice)  ! out
      CALL calc_wskin_fractions_lice( &
        & w_snow_lice(:),             & ! in
        & fract_snow_lice(:)          & ! out
        & )
      RETURN
    END IF

    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(SEB_)
    IF (tile%lcts(1)%id == BARE_TYPE .OR. tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_memory(SOIL_)
    END IF
    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_memory(VEG_)
      dsl4jsb_Get_memory(PHENO_)
    END IF

    ! Set pointers
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_srf)         ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_air)            ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, press_srf)        ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                ! in

    IF (tile%lcts(1)%id == BARE_TYPE .OR. tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(SOIL_,     oro_stddev)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_water_soil) ! out
    END IF
    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)              ! in
      dsl4jsb_Get_var2D_onChunk(VEG_,      fract_fpc_max)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_water_can)  ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)   ! out
    END IF

    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_skin_soil)      ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow_soil)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_water)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)  ! out

    IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      fract_snow(:)      = 1._wp
      fract_snow_soil(:) = 1._wp
      fract_water(:)     = 0._wp
    END IF

    IF (tile%lcts(1)%id == VEG_TYPE .OR. tile%lcts(1)%id == LAND_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow_can)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_skin_can)       ! in
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
        & fract_snow(:)              & ! out
        & )
   END IF

    
    IF (tile%lcts(1)%id == BARE_TYPE) THEN
      CALL calc_wskin_fractions_bare( &
        & dtime,                      & ! in
        & steplen,                    & ! in
        & wsr_capacity,               & ! in
        & oro_stddev(:),              & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & drag_srf(:),                & ! in
        & q_air(:),                   & ! in
        & w_skin_soil(:),             & ! in
        & w_snow_soil(:),             & ! in
        & fract_water_soil(:),        & ! out
        & fract_snow_soil(:)          & ! out
        & )
      fract_water(:) = fract_water_soil(:)
      fract_snow(:)  = fract_snow_soil(:)
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

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snow_and_skin_fraction'

    INTEGER  :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_water,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_water_soil, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_water_can,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_can,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_lice,  weighted_by_fract)

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
    dsl4jsb_Def_memory(ATM2LAND_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rain, &
      & snow, &
      & runoff, &
      & drainage, &
      & evapotranspiration, &
      & w_soil_column, &
      & w_skin, &
      & w_snow, &
      & water_flux, &    ! [m3 s-1]
      & water_content, & ! [m3]
      & water_budget     ! [m3]

    REAL(wp), POINTER :: &
      & area(:)

    ! Locally allocated vectors
    !
    REAL(wp) :: &
      & tile_fract (options%nc), &
      & water_content_new(options%nc)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_water_balance'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    ! IF (tile%lcts() RETURN   ! Check water balance only on root tile
    IF (tile%lcts(1)%id == LAKE_TYPE) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    area => grid%area (ics:ice, iblk)

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(ATM2LAND_)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, rain)                ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, snow)                ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)            ! in
    IF (tile%lcts(1)%id /= GLACIER_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_soil_column)       ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotranspiration)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_skin)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_flux)          ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_content)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_budget)        ! out

    CALL tile%Get_fraction(ics, ice, iblk, fract=tile_fract(:))

    water_flux(:) = rain(:) + snow(:) + evapotranspiration(:) - runoff(:) - drainage(:)
    water_flux(:) = water_flux(:) * tile_fract(:) * area(:) / rhoh2o                    ! kg m-2 s-1 -> m3/s

    IF (tile%lcts(1)%id /= GLACIER_TYPE) THEN
      water_content_new(:) = (w_skin(:) + w_snow(:) + w_soil_column(:)) * tile_fract(:) * area(:)  ! m -> m^3
    ELSE
      water_content_new(:) = 0._wp
    END IF

    water_budget(:) = (water_content_new(:) - water_content(:)) - water_flux(:) * dtime

    water_content(:) = water_content_new(:)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_balance

  SUBROUTINE aggregate_water_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_water_balance'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, water_flux,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_content,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_budget,     weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_balance

#endif
END MODULE mo_hydro_interface
