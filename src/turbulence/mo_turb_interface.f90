!>
!! @brief Contains the interfaces to the turbulence process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-04-01)
!!
MODULE mo_turb_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, VEG_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes TURB, SEB, HYDRO, SOIL, VEG, ATM2LAND

  ! Use process configurations
  dsl4jsb_Use_config(TURB)

  ! Use process memories
  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(TURB)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(HYDRO)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_turb_tasks !,t_turb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_interface'

  !> Type definition for surface_roughness task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_roughness
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_roughness     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_roughness  !< Aggregates computed task variables
  END TYPE tsk_surface_roughness

  !> Constructor interface for surface_roughness task
  INTERFACE tsk_surface_roughness
    PROCEDURE Create_task_surface_roughness                        !< Constructor function for task
  END INTERFACE tsk_surface_roughness

  !> Type definition for humidity_scaling task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_humidity_scaling
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_humidity_scaling     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_humidity_scaling  !< Aggregates computed task variables
  END TYPE tsk_humidity_scaling

  !> Constructor interface for humidity scaling task
  INTERFACE tsk_humidity_scaling
    PROCEDURE Create_task_humidity_scaling                        !< Constructor function for task
  END INTERFACE tsk_humidity_scaling

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_roughness task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_surface_roughness"
  !!
  FUNCTION Create_task_surface_roughness(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_roughness::return_ptr)

    return_ptr%name = 'surface_roughness'
    return_ptr%process_id = TURB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_roughness

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for humidity_scaling task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_humidity_scaling"
  !!
  FUNCTION Create_task_humidity_scaling(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_humidity_scaling::return_ptr)

    return_ptr%name = 'humidity_scaling'
    return_ptr%process_id = TURB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for turbulence process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of TURB process class
  !!
  SUBROUTINE Register_turb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_roughness(model_id))
    CALL this%Register_task(tsk_humidity_scaling(model_id))

  END SUBROUTINE Register_Turb_Tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface_roughness"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_roughness(tile, options)

    USE mo_turb_constants, ONLY: blending_height, roughness_snow

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare process memories
    dsl4jsb_Def_memory(TURB)
    dsl4jsb_Def_memory(HYDRO)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_lice

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_roughness'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TURB)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process memories
    dsl4jsb_Get_memory(TURB)
    dsl4jsb_Get_memory(HYDRO)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(TURB,  rough_m)         ! IN
    dsl4jsb_Get_var2D_onChunk(TURB,  rough_h)         ! OUT

    ! RMS TBD: This has to be revised!
    IF (tile%Is_process_active(SOIL)) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_soil) ! IN
      rough_h(:) =          fract_snow_soil(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
        &        + (1._wp - fract_snow_soil(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
    ELSE IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_lice) ! IN
      rough_h(:) =          fract_snow_lice(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
        &        + (1._wp - fract_snow_lice(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
    ELSE
      CALL finish(TRIM(routine), 'Don"t know what to do.')
    END IF

    ! Limit roughness length for heat by 1 m (see below Eq.5.6 in ECHAM5 manual)
    rough_h(:) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h(:))))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_roughness

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_roughness"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_roughness(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(TURB)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_roughness'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(TURB)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(TURB, rough_h,        weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_roughness

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "humidity_scaling"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_humidity_scaling(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_soil_constants,         ONLY: ws_wilt_fract
    USE mo_soil_process,           ONLY: relative_humidity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(TURB)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & fract_fpc,            &
      & fract_snow,           &
      & fract_water,          &
      & fact_q_air,           &
      & fact_qsat_srf,        &
      & fact_qsat_trans_srf,  &
      & qsat_star,            &
      & canopy_cond,          &
      & pch,                  &
      & wind_air,             &
      & q_air

    dsl4jsb_Real3D_onChunk :: &
      & w_soil_sl,            &
      & w_soil_fc_sl,         &
      & w_soil_pwp_sl

    ! Locally variables
    REAL(wp) ::                    &
      fact_qsat_veg(options%nc),   &
      fact_qsat_trans(options%nc), &
      fact_qsat_soil(options%nc),  &
      fact_qair_veg(options%nc),   &
      fact_qair_soil(options%nc),  &
      rel_hum(options%nc),         &
      zfract_fpc(options%nc)

    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_humidity_scaling'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(TURB) .OR. Contains_lct(tile%lcts, LAKE_TYPE)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get pointers to process configs and memory
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(TURB)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(HYDRO)
    IF (Contains_lct(tile%lcts, VEG_TYPE) .OR. Contains_lct(tile%lcts, LAND_TYPE)) THEN
      dsl4jsb_Get_memory(VEG)
    END IF

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, pch)                 ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, wind_air)            ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)               ! in

    dsl4jsb_Get_var2D_onChunk(SEB,      qsat_star)           ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow)          ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water)         ! in
    dsl4jsb_Get_var2D_onChunk(TURB,     fact_q_air)          ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB,     fact_qsat_srf)       ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB,     fact_qsat_trans_srf) ! OUT

    IF (Contains_lct(tile%lcts, BARE_TYPE) .OR. Contains_lct(tile%lcts, VEG_TYPE) .OR. Contains_lct(tile%lcts, LAND_TYPE)) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_sl)           ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_fc_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_pwp_sl)       ! in
    END IF
  
    ! ---------------------------
    ! Go
    
    IF (Contains_lct(tile%lcts, VEG_TYPE) .OR. Contains_lct(tile%lcts, LAND_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(VEG,      fract_fpc)           ! in
      dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)         ! in

      WHERE (w_soil_sl(:,1) > w_soil_pwp_sl(:,1))
        fact_qsat_veg(:) = fract_snow(:) + (1._wp - fract_snow(:)) *   &
          & (fract_water(:) + ( 1._wp - fract_water(:)) /                &
          & (1._wp + pch(:) * MAX(1._wp, wind_air(:)) / MAX(1.E-20_wp, canopy_cond(:))))
        fact_qsat_trans(:) = (1._wp - fract_snow(:)) * &
          & (1._wp - fract_water(:)) /                  &
          & (1._wp + pch(:) * MAX(1._wp, wind_air(:)) / MAX(1.E-20_wp, canopy_cond(:)))
      ELSEWHERE
        fact_qsat_veg(:)   = fract_snow(:) + (1._wp - fract_snow(:)) * fract_water(:)
        fact_qsat_trans(:) = 0._wp
      END WHERE
      fact_qair_veg(:) = fact_qsat_veg(:)
      zfract_fpc(:)      = fract_fpc(:)

    ELSE

      fact_qsat_veg(:)   = 0._wp
      fact_qsat_trans(:) = 0._wp
      fact_qair_veg(:)   = 0._wp
      zfract_fpc(:)      = 0._wp
    END IF

    IF (Contains_lct(tile%lcts, BARE_TYPE) .OR. Contains_lct(tile%lcts, VEG_TYPE) .OR. Contains_lct(tile%lcts, LAND_TYPE)) THEN

      rel_hum(:) = relative_humidity(w_soil_sl(:,1), w_soil_fc_sl(:,1))

      WHERE (qsat_star(:) * rel_hum(:) > q_air(:) .AND. rel_hum(:) > 1.e-10_wp)
        fact_qsat_soil(:) = fract_snow(:) + (1._wp - fract_snow(:)) * &
          & (fract_water(:) + (1._wp - fract_water(:)) * rel_hum(:))
        fact_qair_soil(:) = 1._wp
      ELSEWHERE
        fact_qsat_soil(:) = fract_snow(:) + (1._wp - fract_snow(:)) * fract_water(:)
        fact_qair_soil(:) = fact_qsat_soil(:)
      END WHERE
      WHERE (q_air(:) > qsat_star(:))
        fact_qsat_soil(:) = 1._wp
        fact_qair_soil(:) = 1._wp
      END WHERE

    ELSE IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN

      fact_qsat_soil(:) = 1._wp
      fact_qair_soil(:) = 1._wp

    ELSE

      fact_qsat_soil(:) = 0._wp
      fact_qair_soil(:) = 0._wp

    END IF

    fact_qsat_srf(:) = zfract_fpc(:) * fact_qsat_veg(:) + (1._wp - zfract_fpc(:)) * fact_qsat_soil(:)
    fact_q_air(:)    = zfract_fpc(:) * fact_qair_veg(:) + (1._wp - zfract_fpc(:)) * fact_qair_soil(:)
    fact_qsat_trans_srf(:) = zfract_fpc(:) * fact_qsat_trans(:)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_humidity_scaling

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "humidity_scaling"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_humidity_scaling(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(TURB)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_humidity_scaling'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(TURB)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(TURB, fact_q_air,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB, fact_qsat_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB, fact_qsat_trans_srf, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_humidity_scaling

  ! ================================================================================================================================
  !>
  !! @brief Another task for "<PROCESS_NAME_LOWER_CASE>"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !  ....

#endif
END MODULE mo_turb_interface
