!> Contains the interfaces to the multi-scale intertile exchange process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_mix_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: ON_TILE_
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes MIX_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(MIX_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(MIX_)
  dsl4jsb_Use_memory(HYDRO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_mix_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_interface'

  !> Type definition for reset_mix task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_mix
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_mix     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_mix  !< Aggregates computed task variables
  END TYPE tsk_reset_mix

  !> Constructor interface for reset_mix task
  INTERFACE tsk_reset_mix
    PROCEDURE Create_task_reset_mix                        !< Constructor function for task
  END INTERFACE tsk_reset_mix

  !> Type definition for hydro_mix task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_hydro_mix
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_hydro_mix     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_hydro_mix  !< Aggregates computed task variables
  END TYPE tsk_hydro_mix

  !> Constructor interface for hydro_mix task
  INTERFACE tsk_hydro_mix
    PROCEDURE Create_task_hydro_mix                        !< Constructor function for task
  END INTERFACE tsk_hydro_mix

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !!
  !> Constructor for reset_mix task
  FUNCTION Create_task_reset_mix(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_mix::return_ptr)
    CALL return_ptr%Construct(name='reset_mix', process_id=MIX_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_mix

  !> Constructor for hydro_mix task
  FUNCTION Create_task_hydro_mix(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_hydro_mix::return_ptr)
    CALL return_ptr%Construct(name='hydro_mix', process_id=MIX_, owner_model_id=model_id)

  END FUNCTION Create_task_hydro_mix

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for multi-scale intertile exchange process
  !!
  SUBROUTINE Register_mix_tasks(process, model_id)

    CLASS(t_jsb_process),  INTENT(inout) :: process
    INTEGER,               INTENT(in)    :: model_id

    CALL process%Register_task(tsk_reset_mix(model_id))
    CALL process%Register_task(tsk_hydro_mix(model_id))

  END SUBROUTINE Register_mix_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "reset_mix"
  !!
  SUBROUTINE update_reset_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_mix'

    ! Declare process memories
    dsl4jsb_Def_memory(MIX_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_latflow

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(MIX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_,   wtr_latflow)   ! inout

    ! Reset lateral intertile fluxes prior to update_*_mix
    wtr_latflow(:) = 0._wp

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_reset_mix

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "reset_mix"
  !!
  SUBROUTINE aggregate_reset_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(MIX_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_mix'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(MIX_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latflow,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_mix

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "hydro_mix"
  !!
  SUBROUTINE update_hydro_mix(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_process,          ONLY: get_soilhyd_properties
    USE mo_mix_process,            ONLY: calc_runoff_transport,   &
      &                                  calc_drainage_transport, &
      &                                  calc_intertile_diffusion

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc, ic
    INTEGER  :: nsoil, is              !< Counter variables for soil layers
    LOGICAL  :: l_lag                  !< T: Apply lag time to lateral flow
    REAL(wp) :: dtime                  !< Model time step [s]
    REAL(wp) :: retention_time         !< retention time for lateral water flow, @todo: needs validation
    REAL(wp) :: trans_min_fract        !< min. source or target tile fraction required for lateral transport

    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_vgrid), POINTER          :: soil_w
    CLASS(t_jsb_tile_abstract), POINTER :: trg_tile !< Target tile with MIX_ process running on them

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hydro_mix'

    ! Declare process memories
    dsl4jsb_Def_config(MIX_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(MIX_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Def_memory_tile(MIX_,   trg_tile)
    dsl4jsb_Def_memory_tile(HYDRO_, trg_tile)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_latflow
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: wtr_pond
    dsl4jsb_Real2D_onChunk :: ice_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_max
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: runoff_retention
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real3D_onChunk :: drainage_sl
    dsl4jsb_Real3D_onChunk :: drainage_retention_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl
    dsl4jsb_Real3D_onChunk :: vol_p_wilt_sl
    dsl4jsb_Real3D_onChunk :: vol_wres_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sat_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pwp_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl
    dsl4jsb_Real3D_onChunk :: matric_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl
    dsl4jsb_Real3D_onChunk :: wtr_latdiff_sl

    dsl4jsb_Real2D_onChunk :: wtr_latflow_trg_tile
    dsl4jsb_Real2D_onChunk :: wtr_pond_trg_tile

    dsl4jsb_Real3D_onChunk :: soil_depth_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: vol_p_wilt_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: vol_wres_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_fc_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_sat_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_pwp_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: ice_soil_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: matric_pot_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: bclapp_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl_trg_tile
    dsl4jsb_Real3D_onChunk :: wtr_latdiff_sl_trg_tile

    REAL(wp), DIMENSION(options%nc) :: source_fract    !< Fractions of water source tile []
    REAL(wp), DIMENSION(options%nc) :: target_fract    !< Fractions of water target tile []
    REAL(wp), DIMENSION(options%nc) :: wtr_pond_max    !< Maximum pond water depth       [m]
    REAL(wp), DIMENSION(options%nc) :: wtr_pond_tmp    !< Temporary pond water storage   [m]

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::         &
    & ice_imp_source,                                & !< Soil layer ice impedance factor on source tile []
    & ice_imp_target,                                & !< Soil layer ice impedance factor on target tile []
    & hyd_cond_source,                               & !< Soil layer hyd. conductivity on source tile [m s-1]
    & hyd_cond_target,                               & !< Soil layer hyd. conductivity on target tile [m s-1]
    & mpot_act_source,                               & !< Soil layer matric potential on source tile [?]
    & mpot_act_target,                               & !< Soil layer matric potential on target tile [?]
    & weq_wsat_sl,                                   & !< Maximum water content of soil layer on source tile [m]
    & weq_wsat_sl_trg_tile,                          & !< Maximum water content of soil layer on target tile [m]
    & weq_field_cap_sl,                              & !< Field capacity of soil layer on source tile [m]
    & weq_field_cap_sl_trg_tile,                     & !< Field capacity of soil layer on target tile [m]
    & weq_p_wilt_sl,                                 & !< Wilting point of soil layer on source tile [m]
    & weq_p_wilt_sl_trg_tile,                        & !< Wilting point of soil layer on target tile [m]
    & weq_wres_sl,                                   & !< Water residual content of soil layer on source tile [m]
    & weq_wres_sl_trg_tile                             !< Water residual content of soil layer on target tile [m]

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(MIX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(MIX_)
    dsl4jsb_Get_config(HYDRO_)

    trans_min_fract = dsl4jsb_Config(MIX_)%trans_min_fract ! min. tile fraction for lateral transport

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_,   wtr_latflow)            ! inout
    dsl4jsb_Get_var3D_onChunk(MIX_,   wtr_latdiff_sl)         ! inout
    dsl4jsb_Get_var2D_onChunk(MIX_,   runoff_retention)       ! inout
    dsl4jsb_Get_var3D_onChunk(MIX_,   drainage_retention_sl)  ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)                 ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)               ! inout
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, drainage_sl)          ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)          ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)             ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)             ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)             ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond_max)         ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_Cap_sl)     ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_p_wilt_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_wres_sl)          ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_fc_sl)       ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pwp_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matric_pot_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl)
    END IF

    ! Transport water from between tiles
    ! Note, currently, lateral transport is only considered on tiles that contain soil
    IF (dsl4jsb_Config(MIX_)%l_hydro_mix .AND. tile%contains_soil) THEN
      ! Allocate local variables
      ALLOCATE(ice_imp_source(nc,nsoil))
      ALLOCATE(ice_imp_target(nc,nsoil))
      ALLOCATE(hyd_cond_source(nc,nsoil))
      ALLOCATE(hyd_cond_target(nc,nsoil))
      ALLOCATE(mpot_act_source(nc,nsoil))
      ALLOCATE(mpot_act_target(nc,nsoil))
      ALLOCATE(weq_wsat_sl(nc,nsoil))
      ALLOCATE(weq_wsat_sl_trg_tile(nc,nsoil))
      ALLOCATE(weq_field_cap_sl(nc,nsoil))
      ALLOCATE(weq_field_cap_sl_trg_tile(nc,nsoil))
      ALLOCATE(weq_p_wilt_sl(nc,nsoil))
      ALLOCATE(weq_p_wilt_sl_trg_tile(nc,nsoil))
      ALLOCATE(weq_wres_sl(nc,nsoil))
      ALLOCATE(weq_wres_sl_trg_tile(nc,nsoil))

      ! This check should be moved to init routine as soon as enough stuff is there to justify one
      IF (tile%process_action(HYDRO_) /= ON_TILE_) THEN
        CALL finish(routine, 'MIX_ processes must be activate on the same tiles as HYDRO_')
      END IF

      ! Compute how much water still fits into the pond storage (needed even without lateral transport)
      wtr_pond_max(:) = MAX(0._wp, dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max) - ice_pond(:))

      ! Check whether actual tile has a flow target for lateral surface runoff transport and drainage
      IF (ASSOCIATED(tile%downstream_tile)) THEN
        ! Get downstream tile for lagged flow
        trg_tile => tile%downstream_tile
        l_lag = .TRUE.                 ! Lagged flow due to significant distance 
      ELSE IF (ASSOCIATED(tile%adjacent_tile)) THEN
        ! Get adjacent tile for instant flow
        trg_tile => tile%adjacent_tile
        l_lag = .FALSE.                ! Instant flow due to adjacency
      ELSE
        ! No flow target available for this tile
        trg_tile => NULL()
      END IF

      ! Transport surface runoff and drainage to downstream or adjacent tile
      IF (ASSOCIATED(trg_tile)) THEN

        ! Get area fraction of flow source tile and target tile
        CALL tile%Get_fraction(ics, ice, iblk, fract=source_fract(:))
        CALL trg_tile%Get_fraction(ics, ice, iblk, fract=target_fract(:))

        ! Get memory MIX_ and HYDRO_ memory for target tile needed for
        ! lateral surface runoff and drainage transport
        dsl4jsb_Get_memory_tile(MIX_,   trg_tile)
        dsl4jsb_Get_memory_tile(HYDRO_, trg_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow, trg_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_pond,    trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sl, trg_tile)

        ! Transport surface runoff from source to target tile
        retention_time = 36000._wp   ! Ten hours retention time (
                                     ! outflow = h_re * dt/retention_time)
        CALL calc_runoff_transport(   &
          ! in
          & l_lag,                    & !< T: Apply lag (only for lagged flow)
          & dtime,                    & !< model time step              [s]
          & retention_time,           & !< lagged flow retention time   [s]
          & wtr_pond_max(:),          & !< maximum pond water content   [m]
          & source_fract(:),          & !< area fraction of source tile []
          & target_fract(:),          & !< area fraction of target tile []
          & trans_min_fract,          & !< min. tile fraction           []
          ! inout
          & runoff(:),                & !< surface runoff               [kg m-2 s-1]
          & wtr_latflow(:),           & !< lateral flow at source tile  [kg m-2 s-1]
          & wtr_latflow_trg_tile(:),  & !< lateral flow at target tile  [kg m-2 s-1]
          & wtr_pond(:),              & !< pond storage at source tile  [m]
          & wtr_pond_trg_tile(:),     & !< pond storage at target tile  [m]
          & runoff_retention(:)       & !< runoff retention storage     [m]
          & )

        ! Transport drainage from source to target tile within each soil layer
        retention_time = 360000._wp  ! Hundert hours retention time (
                                     ! outflow = h_re * dt/retention_time)
        DO is = 1, nsoil
          CALL calc_drainage_transport( &
            ! in
            & l_lag,                                 & !< T: Apply lag
            & dtime,                                 & !< model time step                          [s]
            & retention_time,                        & !< lagged flow retention time               [s]
            & wtr_soil_fc_sl(:,is),                  & !< soil layer field capacity at source tile [m]
            & source_fract(:),                       & !< area fraction of source tile             []
            & target_fract(:),                       & !< area fraction of target tile             []
            & trans_min_fract,                       & !< min. tile fraction                       []
            ! inout
            & drainage_sl(:,is),                     & !< lateral drainage at source tile          [kg m-2 s-1]
            & drainage(:),                           & !< total drainage at source tile            [kg m-2 s-1]
            & wtr_latflow(:),                        & !< lateral flow at source tile              [kg m-2 s-1]
            & wtr_latflow_trg_tile(:),               & !< lateral flow at target tile              [kg m-2 s-1]
            & wtr_soil_sl(:,is),                     & !< soil water content at source tile        [m]
            & wtr_soil_sl_trg_tile(:,is),            & !< soil water content at target tile        [m]
            & drainage_retention_sl(:,is)            & !< lateral drainage retention storage       [m]
            & )
        ENDDO
      END IF

      ! Check whether actual tile has an adjacent tile for lateral diffusion
      IF (ASSOCIATED(tile%adjacent_tile)) THEN
        ! Get adjacent tile for lateral diffusion
        trg_tile => tile%adjacent_tile
      ELSE
        ! No flow target available for this tile
        trg_tile => NULL()
      END IF

      ! Enable lateral diffusion between adjacent tiles
      IF (ASSOCIATED(trg_tile)) THEN
        ! Get memory MIX_ and HYDRO_ memory for target tile
        dsl4jsb_Get_memory_tile(MIX_,   trg_tile)
        dsl4jsb_Get_memory_tile(HYDRO_, trg_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow,              trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(MIX_,   wtr_latdiff_sl,           trg_tile)
        dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_pond,                 trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sl,              trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_fc_sl,           trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, soil_depth_sl,            trg_tile)
        dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_porosity_sl,          trg_tile) ! in
        dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_field_cap_sl,         trg_tile) ! in
        dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_p_wilt_sl,            trg_tile) ! in
        dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_wres_sl,              trg_tile) ! in
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sat_sl,          trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_pwp_sl,          trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, ice_soil_sl,              trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_pot_scool_sl,    trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, hyd_cond_sat_sl,          trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, matric_pot_sl,            trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, bclapp_sl,                trg_tile)
        dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, pore_size_index_sl,       trg_tile)

        ! Compute soil ice impedance, diffusivity and matric potential ...
        ! ... for source tile
        weq_wsat_sl     (:,:) = MAX(vol_porosity_sl (:,:) * soil_depth_sl(:,:), 0._wp)
        weq_field_cap_sl(:,:) = MAX(vol_field_cap_sl(:,:) * soil_depth_sl(:,:), 0._wp)
        weq_p_wilt_sl   (:,:) = MAX(vol_p_wilt_sl(:,:)    * soil_depth_sl(:,:), 0._wp)
        weq_wres_sl     (:,:) = MAX(vol_wres_sl(:,:)      * soil_depth_sl(:,:), 0._wp)
        CALL get_soilhyd_properties(       &
          dsl4jsb_Config(HYDRO_)%soilhydmodel,  & !< Soil hydrology model
          1,                                    & !< interpolation methods (not applied in this case)
          nc, nsoil,                            & !< number of cells and soil layers
          soil_depth_sl(:,:),                   & !< soil layer depth   [m]
          wtr_soil_sl(:,:),                     & !< soil water content [m]
          ice_soil_sl(:,:),                     & !< soil ice content   [m]
          wtr_soil_pot_scool_sl(:,:),           & !< potential supercooled water [m]
          weq_wsat_sl(:,:),                     & !< maximum soil moisture storage (water + ice) [m]
          weq_field_cap_sl(:,:),                & !< soil field capacity (water + ice) [m]
          weq_p_wilt_sl(:,:),                   & !< soil wilting point (water + ice) [m]
          weq_wres_sl(:,:),                     & !< soil water residual [m]
          hyd_cond_sat_sl(:,:),                 & !< saturated hydraulic conductivity [m s-1]
          matric_pot_sl(:,:),                   & !< matric potential
          bclapp_sl(:,:),                       & !< clapp hornberger b parameter
          pore_size_index_sl(:,:),              & !< pore size index
          ice_impedance=ice_imp_source(:,:),    & !< ice impedance factor
          K=hyd_cond_source(:,:),               & !< hydraulic conductivity
          mpot_act=mpot_act_source(:,:)         & !< matric potential
        )

        ! ... and target tile
        weq_wsat_sl_trg_tile     (:,:) = MAX(vol_porosity_sl_trg_tile (:,:) * soil_depth_sl_trg_tile(:,:), 0._wp)
        weq_field_cap_sl_trg_tile(:,:) = MAX(vol_field_cap_sl_trg_tile(:,:) * soil_depth_sl_trg_tile(:,:), 0._wp)
        weq_p_wilt_sl_trg_tile   (:,:) = MAX(vol_p_wilt_sl_trg_tile(:,:)    * soil_depth_sl_trg_tile(:,:), 0._wp)
        weq_wres_sl_trg_tile     (:,:) = MAX(vol_wres_sl_trg_tile(:,:)      * soil_depth_sl_trg_tile(:,:), 0._wp)
        CALL get_soilhyd_properties(       &
          dsl4jsb_Config(HYDRO_)%soilhydmodel, & !< Soil hydrology model
          1,                                   & !< interpolation methods (not applied in this case)
          nc, nsoil,                           & !< number of cells and soil layers
          soil_depth_sl_trg_tile(:,:),         & !< soil layer depth   [m]
          wtr_soil_sl_trg_tile(:,:),           & !< soil water content [m]
          ice_soil_sl_trg_tile(:,:),           & !< soil ice content   [m]
          wtr_soil_pot_scool_sl_trg_tile(:,:), & !< potential supercooled water [m]
          weq_wsat_sl_trg_tile(:,:),           & !< maximum soil moisture storage (water + ice) [m]
          weq_field_cap_sl_trg_tile(:,:),      & !< soil field capacity (water + ice) [m]
          weq_p_wilt_sl_trg_tile(:,:),         & !< soil wilting point (water + ice) [m]
          weq_wres_sl_trg_tile(:,:),           & !< soil water residual [m]
          hyd_cond_sat_sl_trg_tile(:,:),       & !< saturated hydraulic conductivity [m s-1]
          matric_pot_sl_trg_tile(:,:),         & !< matric potential
          bclapp_sl_trg_tile(:,:),             & !< clapp hornberger b parameter
          pore_size_index_sl_trg_tile(:,:),    & !< pore size index
          ice_impedance=ice_imp_target(:,:),   & !< ice impedance factor
          K=hyd_cond_target(:,:),              & !< hydraulic conductivity
          mpot_act=mpot_act_target(:,:)        & !< matric potential
        )

        ! Lateral soil water diffusion between adjacent tiles
        DO ic = 1, nc
          IF (source_fract(ic) > trans_min_fract .AND. target_fract(ic) > trans_min_fract .AND. &
              wtr_soil_fc_sl(ic,1) > 1.0e-10_wp  .AND. wtr_soil_fc_sl_trg_tile(ic,1) > 1.0e-10_wp) THEN

            CALL calc_intertile_diffusion(nsoil,                        & !< Number of soil layers
               & dsl4jsb_Config(HYDRO_)%soilhydmodel,                   & !< Soil hydrology model
               & dtime,                                                 & !< time step                   [s]
               & weq_pond_max(ic),                                      & !< height difference btw source and target tile [m]
               & source_fract(ic),      target_fract(ic),               & !< tile fractions              []
               & soil_depth_sl(ic,:),   soil_depth_sl_trg_tile(ic,:),   & !< soil layer depths           [m]
               & wtr_soil_fc_sl(ic,:),  wtr_soil_fc_sl_trg_tile(ic,:),  & !< soil field capacity         [m]
               & wtr_soil_sat_sl(ic,:), wtr_soil_sat_sl_trg_tile(ic,:), & !< soil maximum water content  [m]
               & wtr_soil_pwp_sl(ic,:), wtr_soil_pwp_sl_trg_tile(ic,:), & !< soil wilting point          [m]
               & matric_pot_sl(ic,:),   matric_pot_sl_trg_tile(ic,:),   & !< soil matric potential       [m]
               & hyd_cond_sat_sl(ic,:), hyd_cond_sat_sl_trg_tile(ic,:), & !< soil hydraulic conductivity [m s-1]
               & bclapp_sl(ic,:),       bclapp_sl_trg_tile(ic,:),       & !< soil bclapp parameter       []
               & pore_size_index_sl(ic,:), pore_size_index_sl_trg_tile(ic,:), & !< pore size index       []
               & ice_imp_source(ic,:),  ice_imp_target(ic,:),           & !< ice impedance factor        []
               & hyd_cond_source(ic,:), hyd_cond_target(ic,:),          & !< hyd. conductivity           [m s-1]
               & mpot_act_source(ic,:), mpot_act_target(ic,:),          & !< matric potential            [?]
               & ice_soil_sl(ic,:),     ice_soil_sl_trg_tile(ic,:),     & !< soil ice storage            [m]
               & MIN(wtr_soil_pot_scool_sl(ic,:), wtr_soil_sl(ic,:)),   & !< supercooled water           [m]
               & MIN(wtr_soil_pot_scool_sl_trg_tile(ic,:), wtr_soil_sl_trg_tile(ic,:)), &  !< supercooled water [m]
               & wtr_pond(ic),          wtr_pond_trg_tile(ic),          & !< pond water storage          [m]
               & wtr_latflow(ic),       wtr_latflow_trg_tile(ic),       & !< lateral water flow          [kg m-2 s-1]
               & wtr_soil_sl(ic,:),     wtr_soil_sl_trg_tile(ic,:),     & !< soil water storage          [m]
               & wtr_latdiff_sl(ic,:),  wtr_latdiff_sl_trg_tile(ic,:)   & !< lateral diffusion flux      [kg m-2 s-1]
               & )

          END IF
        ENDDO
      END IF

      ! Runoff on floodplains is only generated if ponds are filled and overflowing
      IF (tile%contains_floodplain) THEN
        WHERE (wtr_pond(:) < wtr_pond_max(:))
            wtr_pond_tmp(:) = wtr_pond(:) + runoff(:) * dtime / rhoh2o
            runoff(:)       = MAX(0._wp, wtr_pond_tmp(:) - wtr_pond_max(:)) * rhoh2o / dtime
            wtr_pond(:)     = MIN(wtr_pond_tmp(:), wtr_pond_max(:))
        END WHERE
      END IF

      ! Pond content which exceeds maximum pond storage capacity goes directly into surface runoff 
      WHERE (wtr_pond(:) > wtr_pond_max(:))
        runoff(:)   = runoff(:) + MAX(0._wp, wtr_pond(:) - wtr_pond_max(:)) * rhoh2o / dtime
        wtr_pond(:) = wtr_pond_max(:)
      END WHERE

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      ! Update total pond storage after lateral water flow
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        weq_pond(ic) = wtr_pond(ic) + ice_pond(ic)
      END DO
      !$ACC END LOOP

      ! Deallocate local fields
      DEALLOCATE(ice_imp_source,ice_imp_target, hyd_cond_source,hyd_cond_target, mpot_act_source,mpot_act_target, &
        & weq_wsat_sl,weq_wsat_sl_trg_tile,weq_field_cap_sl,weq_field_cap_sl_trg_tile, weq_p_wilt_sl,           &
        & weq_p_wilt_sl_trg_tile, weq_wres_sl, weq_wres_sl_trg_tile)

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_hydro_mix

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "hydro_mix"
  !!
  SUBROUTINE aggregate_hydro_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(MIX_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_hydro_mix'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(MIX_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latflow,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   runoff_retention,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,                 weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage_sl,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   drainage_retention_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latdiff_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,            weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_hydro_mix

#endif
END MODULE mo_mix_interface
