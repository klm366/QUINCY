!>
!! @brief Contains the interfaces to the turbulence process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-04-01)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_turb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_control,        ONLY: jsbach_runs_standalone

  ! Use of processes in this module
  dsl4jsb_Use_processes TURB_, SEB_, HYDRO_, PHENO_, ASSIMI_, SSE_, A2L_

  ! Use process configurations
  dsl4jsb_Use_config(TURB_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SEB_)

  ! Use process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)

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
    CALL return_ptr%Construct(name='surface_roughness', process_id=TURB_, owner_model_id=model_id)

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
    CALL return_ptr%Construct(name='humidity_scaling', process_id=TURB_, owner_model_id=model_id)

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

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_config(SEB_)

    ! Declare process memories
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(SEB_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: htc_stable
    dsl4jsb_Real2D_onChunk :: htc_unstable
    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: rough_m_star
    dsl4jsb_Real2D_onChunk :: rough_h_star
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_lice
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2d_onChunk :: fract_fpc_max

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER  :: iblk, ics, ice
    REAL(wp), DIMENSION(options%nc) :: &
      & rough_veg_star, rough_bare_star, rough_snow_star
    REAL(wp) :: &
      & config_rough_bare, config_rough_snow, config_rough_water, config_blending_height, config_roughness_lai_saturation, &
      & config_roughness_momentum_to_heat, lctlib_MinVegRoughness, lctlib_MaxVegRoughness, lctlib_VegRoughness,            &
      & lctlib_HeatTransCoeffStable, lctlib_HeatTransCoeffUnstable

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_roughness'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TURB_)
    dsl4jsb_Get_config(SEB_)

    ! Get process memories
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(TURB_,  htc_stable)      ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  htc_unstable)    ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m)         ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h)         ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m_star)    ! out
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h_star)    ! out

    config_blending_height = dsl4jsb_Config(TURB_)%blending_height
    config_rough_snow                 = dsl4jsb_Config(TURB_)%roughness_snow
    config_rough_water                = dsl4jsb_Config(TURB_)%roughness_water

    lctlib_HeatTransCoeffStable   = dsl4jsb_Lctlib_param(HeatTransCoeffStable)
    lctlib_HeatTransCoeffUnstable = dsl4jsb_Lctlib_param(HeatTransCoeffUnstable)

    htc_stable(:)   = lctlib_HeatTransCoeffStable
    htc_unstable(:) = lctlib_HeatTransCoeffUnstable

    IF (tile%Is_process_active(SSE_)) THEN  ! Should catch all vegetated, bare and glacier tiles

      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil) ! in

      config_rough_bare                 = dsl4jsb_Config(TURB_)%roughness_bare
      config_roughness_momentum_to_heat = dsl4jsb_Config(TURB_)%roughness_momentum_to_heat

      !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_0 tile_name=tile%name
      rough_bare_star(:) = 1._wp  / LOG(config_blending_height / config_rough_bare)**2
      rough_snow_star(:) = 1._wp  / LOG(config_blending_height / config_rough_snow)**2

      IF (tile%contains_vegetation) THEN  ! vegetated tile

        dsl4jsb_Get_memory(PHENO_)

        dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max)   ! in

        IF (dsl4jsb_Config(TURB_)%l_roughness_lai) THEN

          config_roughness_lai_saturation  = dsl4jsb_Config(TURB_)%roughness_lai_saturation
          lctlib_MinVegRoughness           = dsl4jsb_Lctlib_param(MinVegRoughness)
          lctlib_MaxVegRoughness           = dsl4jsb_Lctlib_param(MaxVegRoughness)

          dsl4jsb_Get_var2D_onChunk(PHENO_, lai)              ! in

          !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_1 tile_name=tile%name
          rough_veg_star(:) =                                                                           &
            & 1._wp / LOG(config_blending_height                                                        &
            &             / (lctlib_MinVegRoughness + (lctlib_MaxVegRoughness - lctlib_MinVegRoughness) &
            &                                         * TANH(lai(:) * config_roughness_lai_saturation)) &
            &            )**2

        ELSE
          lctlib_VegRoughness = dsl4jsb_Lctlib_param(VegRoughness)
          !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_2 tile_name=tile%name
          rough_veg_star(:) = 1._wp / LOG(config_blending_height / lctlib_VegRoughness)**2
        ENDIF

        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_3 tile_name=tile%name
        rough_m_star(:) = (1._wp - fract_fpc_max(:)) * rough_bare_star(:) + fract_fpc_max(:) * rough_veg_star(:)

        ! Modify roughness length heat for non-vegetated part by snow cover
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_4 tile_name=tile%name
        rough_bare_star(:) = (1._wp - fract_snow_soil(:)) * rough_bare_star(:) + fract_snow_soil(:) * rough_snow_star(:)
        ! Modify roughness length heat for vegetation by snow cover if not forest or if standalone
        IF (.NOT. dsl4jsb_Lctlib_param(ForestFlag) .OR. jsbach_runs_standalone()) THEN
          !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_5 tile_name=tile%name
          rough_veg_star(:) = (1._wp - fract_snow_soil(:)) * rough_veg_star(:) + fract_snow_soil(:) * rough_snow_star(:)
        END IF

        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_6 tile_name=tile%name
        rough_h_star(:) = (1._wp - fract_fpc_max(:)) * rough_bare_star(:) + fract_fpc_max(:) * rough_veg_star(:)

      ELSE  ! bare tile
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_7 tile_name=tile%name
        rough_m_star(:) = rough_bare_star(:)
        ! Modify roughness length heat for non-vegetated part by snow cover
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_8 tile_name=tile%name
        rough_bare_star(:) = (1._wp - fract_snow_soil(:)) * rough_bare_star(:) + fract_snow_soil(:) * rough_snow_star(:)
        rough_h_star(:) = rough_bare_star(:)
      END IF  ! vegetated or bare tile

      !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_10 tile_name=tile%name
      rough_m(:) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(:)))
      rough_h(:) = config_blending_height * EXP(-1._wp / SQRT(rough_h_star(:))) / config_roughness_momentum_to_heat

    ELSE IF (tile%contains_lake) THEN

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_memory(SEB_)
        dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice) ! IN
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_11 tile_name=tile%name
        rough_m_star(:) =            fract_lice(:)  / (LOG(config_blending_height / config_rough_snow)  ** 2)  &
          &               + (1._wp - fract_lice(:)) / (LOG(config_blending_height / config_rough_water) ** 2)
      ELSE
        !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_12 tile_name=tile%name
        rough_m_star(:) = 1._wp / (LOG(config_blending_height / config_rough_water) ** 2)
      END IF
      rough_h_star(:) = rough_m_star(:)

      !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_13 tile_name=tile%name
      rough_m(:) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(:)))
      rough_h(:) = rough_m(:)

    ELSE
      CALL finish(TRIM(routine), 'Don"t know what to do.')
    END IF

    ! Limit roughness length for heat by 1 m (see below Eq.5.6 in ECHAM5 manual)
    ! !$claw expand parallel update savepoint=ser_jsb_turb_update_roughness_expand_14 tile_name=tile%name
    ! rough_h(:) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h(:))))

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

    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(TURB_)

    dsl4jsb_Real2D_onChunk :: htc_stable
    dsl4jsb_Real2D_onChunk :: htc_unstable
    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: rough_m_star
    dsl4jsb_Real2D_onChunk :: rough_h_star

    REAL(wp) :: config_blending_height, config_roughness_momentum_to_heat

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_roughness'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(TURB_)
    dsl4jsb_Get_memory(TURB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_m_star)
    dsl4jsb_Get_var2D_onChunk(TURB_,  rough_h_star)

    dsl4jsb_Aggregate_onChunk(TURB_, rough_m_star,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, rough_h_star,   weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(TURB_, htc_stable,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, htc_unstable,   weighted_by_fract)

    config_blending_height            = dsl4jsb_Config(TURB_)%blending_height
    config_roughness_momentum_to_heat = dsl4jsb_Config(TURB_)%roughness_momentum_to_heat

    !$claw expand parallel update savepoint=ser_jsb_turb_aggregate_roughness_expand_1 tile_name=tile%name
    rough_m(:) = config_blending_height * EXP(-1._wp / SQRT(rough_m_star(:)))
    rough_h(:) = config_blending_height * EXP(-1._wp / SQRT(rough_h_star(:))) / config_roughness_momentum_to_heat

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

    !USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_sse_process,            ONLY: relative_humidity_soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(HYDRO_)

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(A2L_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & fract_fpc,            &
      & fract_snow,           &
      & fract_water,          &
      & max_moist,            &
      & w_soil_column,        &
      & fact_q_air,           &
      & fact_qsat_srf,        &
      & fact_qsat_trans_srf,  &
      & qsat_star,            &
      & canopy_cond_limited,  &
      & pch,                  &
      & wind_air,             &
      & q_air

    dsl4jsb_Real3D_onChunk :: &
      & w_soil_sl,            &
      & w_soil_fc_sl,         &
      & w_soil_pwp_sl

    ! Local variables
    REAL(wp) ::                    &
      fact_qsat_veg(options%nc),   &
      fact_qsat_trans(options%nc), &
      fact_qsat_soil(options%nc),  &
      fact_qair_veg(options%nc),   &
      fact_qair_soil(options%nc),  &
      rel_hum(options%nc),         & ! relative humidity of uppermost soil layer as compared to its field capacity
      zfract_fpc(options%nc),      &
      w_soil_wilt_fract_tmp
    LOGICAL :: &
      mask(options%nc)

    INTEGER :: iblk, ics, ice, i

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_humidity_scaling'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(TURB_) .OR. tile%contains_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)

    !$claw acc data create(fact_qsat_veg,fact_qsat_trans,fact_qsat_soil, &
    !$claw acc fact_qair_veg,fact_qair_soil,rel_hum,zfract_fpc,mask)

    ! Get pointers to process configs and memory
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
      IF (tile%Is_process_active(ASSIMI_)) THEN
        dsl4jsb_Get_memory(ASSIMI_)
      END IF
    END IF

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   pch)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_air)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)          ! in
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_q_air)          ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_srf)       ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_trans_srf) ! OUT
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,    qsat_star)           ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  fract_water)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  max_moist)           ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,  w_soil_column)       ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  w_soil_sl)           ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  w_soil_fc_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,  w_soil_pwp_sl)       ! in
    END IF

    ! ---------------------------
    ! Go

    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc)           ! in
      IF (tile%Is_process_active(ASSIMI_)) THEN
        dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_limited)         ! in
      ELSE
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_limited)         ! in
      END IF

      IF (model%config%l_compat401) THEN
        !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_0 tile_name=tile%name
        mask(:) = w_soil_sl(:,1) > w_soil_pwp_sl(:,1)
      ELSE
        w_soil_wilt_fract_tmp = dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract
        !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_1 tile_name=tile%name
        mask(:) = w_soil_column(:) > w_soil_wilt_fract_tmp * max_moist(:)
      END IF

      !$claw acc data present(fract_snow,fract_water,pch,wind_air,canopy_cond_limited,fact_qsat_veg,fact_qsat_trans,mask)
      !$claw acc parallel loop gang vector
      DO i = 1, size(mask, dim=1)
        IF (mask(i)) THEN
          fact_qsat_veg(i) = fract_snow(i) + (1._wp - fract_snow(i)) *   &
            & (fract_water(i) + ( 1._wp - fract_water(i)) /                &
            & (1._wp + pch(i) * MAX(1._wp, wind_air(i)) / MAX(1.E-20_wp, canopy_cond_limited(i))))
          fact_qsat_trans(i) = (1._wp - fract_snow(i)) * &
            & (1._wp - fract_water(i)) /                  &
            & (1._wp + pch(i) * MAX(1._wp, wind_air(i)) / MAX(1.E-20_wp, canopy_cond_limited(i)))
        ELSE
          fact_qsat_veg(i)   = fract_snow(i) + (1._wp - fract_snow(i)) * fract_water(i)
          fact_qsat_trans(i) = 0._wp
        END IF
      END DO
      !$claw acc end data

      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_2 tile_name=tile%name
      fact_qair_veg(:) = fact_qsat_veg(:)
      zfract_fpc(:)      = fract_fpc(:)
      !$claw end expand
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_3 tile_name=tile%name
      fact_qsat_veg(:)   = 0._wp
      fact_qsat_trans(:) = 0._wp
      fact_qair_veg(:)   = 0._wp
      zfract_fpc(:)      = 0._wp
      !$claw end expand
    END IF

    IF (tile%contains_soil) THEN
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_4 tile_name=tile%name
      rel_hum(:) = relative_humidity_soil(w_soil_sl(:,1), w_soil_fc_sl(:,1))

      !$claw acc data present(fact_qsat_soil,fact_qair_soil,qsat_star,rel_hum,q_air,fract_snow,fract_water)
      !$claw acc parallel loop gang vector
      DO i = 1, size(qsat_star,dim=1)
        IF (qsat_star(i) * rel_hum(i) > q_air(i) .AND. rel_hum(i) > 1.e-10_wp) THEN
          fact_qsat_soil(i) = fract_snow(i) + (1._wp - fract_snow(i)) * &
            & (fract_water(i) + (1._wp - fract_water(i)) * rel_hum(i))
          fact_qair_soil(i) = 1._wp
        ELSE
          fact_qsat_soil(i) = fract_snow(i) + (1._wp - fract_snow(i)) * fract_water(i)
          fact_qair_soil(i) = fact_qsat_soil(i)
        END IF
        IF (q_air(i) > qsat_star(i)) THEN
          fact_qsat_soil(i) = 1._wp
          fact_qair_soil(i) = 1._wp
        END IF
      END DO
      !$claw acc end data

    ELSE IF (tile%contains_glacier) THEN
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_5 tile_name=tile%name
      fact_qsat_soil(:) = 1._wp
      fact_qair_soil(:) = 1._wp
      !$claw end expand
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_6 tile_name=tile%name
      fact_qsat_soil(:) = 0._wp
      fact_qair_soil(:) = 0._wp
      !$claw end expand
    END IF

    !$claw expand parallel update savepoint=ser_jsb_turb_update_humidity_scaling_expand_7 tile_name=tile%name
    fact_qsat_srf(:) = zfract_fpc(:) * fact_qsat_veg(:) + (1._wp - zfract_fpc(:)) * fact_qsat_soil(:)
    fact_q_air(:)    = zfract_fpc(:) * fact_qair_veg(:) + (1._wp - zfract_fpc(:)) * fact_qair_soil(:)
    fact_qsat_trans_srf(:) = zfract_fpc(:) * fact_qsat_trans(:)
    !$claw end expand

    !$claw acc end data

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

    dsl4jsb_Def_memory(TURB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_humidity_scaling'

    ! Local variables
    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(TURB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(TURB_, fact_q_air,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_trans_srf, weighted_by_fract)

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
