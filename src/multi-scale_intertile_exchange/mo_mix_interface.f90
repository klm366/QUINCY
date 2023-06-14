!>
!! @brief Contains the interfaces to the multi-scale intertile exchange process
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-07-22)
!! Tobias Stacke: Revisions of lateral transport for surface runoff and lateral drainage
!!                as well as inter-tile diffusion based on input from Philipp de Vrese, MPI-M (2023-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

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
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: ON_TILE_
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes MIX_, HYDRO_, SSE_

  ! Use of process configurations
  dsl4jsb_Use_config(MIX_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)

  ! Use of process memories
  dsl4jsb_Use_memory(MIX_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SSE_)

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

    USE mo_jsb_math_constants,     ONLY: pi
    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_hydro_constants,        ONLY: k_brock
    USE mo_hydro_process,          ONLY: get_soilhyd_cond_properties
    USE mo_mix_process,            ONLY: calc_runoff_transport,   &
      &                                  calc_drainage_transport, &
      &                                  calc_intertile_diffusion, &
      &                                  calc_snow_redistribution, &
      &                                  diagnose_intertile_heatflux, &
      &                                  get_geometry, calc_slopes_tiles

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    INTEGER  :: nsibl, ic              !< Counter variables for sibling tiles
    INTEGER  :: nsoil, i, j            !< Counter variables for soil layers
    REAL(wp) :: dtime                  !< Model time step [s]
    LOGICAL  :: l_lag                  !< T: Apply lag time to lateral flow
    LOGICAL  :: l_latflow              !< T: Allow for lateral flow between sibling tiles
    REAL(wp) :: retention_time         !< retention time for lateral water flow, @todo: needs validation

    REAL(wp) :: weq_snow_old
    REAL(wp) :: weq_snow_soil_old
    REAL(wp) :: weq_snow_new_soil_old

    REAL(wp) :: weq_snow_old_sbl
    REAL(wp) :: weq_snow_soil_old_sbl
    REAL(wp) :: weq_snow_new_soil_old_sbl
    CHARACTER(len=5000)         :: message_text_long
    !<pdv
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_vgrid), POINTER          :: soil_w
    TYPE(t_jsb_vgrid), POINTER          :: soil_e
    CLASS(t_jsb_tile_abstract), POINTER :: sibl_tile !< Sibling tile with MIX_ process running on them

    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hydro_mix'

    ! Declare process memories
    dsl4jsb_Def_config(MIX_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)

    dsl4jsb_Def_memory(MIX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    dsl4jsb_Def_memory_tile(MIX_,   sibl_tile)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl_tile)
    dsl4jsb_Def_memory_tile(SSE_, sibl_tile)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: thaw_depth_max_ym1
    dsl4jsb_Real2D_onChunk :: wtr_latflow
    dsl4jsb_Real2D_onChunk :: wtr_pond
    dsl4jsb_Real2D_onChunk :: weq_pond

    dsl4jsb_Real2D_onChunk :: snow_rdst
    dsl4jsb_Real2D_onChunk :: weq_snow_new_soil
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: weq_snow
    dsl4jsb_Real2D_onChunk :: snow_soil_dens

    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: runoff_retention
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: drain_rock
    dsl4jsb_Real3D_onChunk :: drain_lat_sl
    dsl4jsb_Real3D_onChunk :: drain_lat_retention_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl
    dsl4jsb_Real3D_onChunk :: vol_p_wilt_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sat_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pwp_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl
    dsl4jsb_Real3D_onChunk :: wtr_diflat_sl

    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl
    dsl4jsb_Real3D_onChunk :: heat_cond_sl
    dsl4jsb_Real3D_onChunk :: heat_flx_lat

    dsl4jsb_Real2D_onChunk :: fract_perm_ref

    dsl4jsb_Real2D_onChunk :: area
    dsl4jsb_Real2D_onChunk :: h_surf
    dsl4jsb_Real2D_onChunk :: oro_stddev

    dsl4jsb_Real2D_onChunk :: steepness
    dsl4jsb_Real2D_onChunk :: micro_slope
    dsl4jsb_Real2D_onChunk :: macro_slope


    dsl4jsb_Real2D_onChunk :: wtr_latflow_sibl_tile
    dsl4jsb_Real2D_onChunk :: wtr_pond_sibl_tile
    dsl4jsb_Real2D_onChunk :: weq_pond_sibl_tile

    dsl4jsb_Real2D_onChunk :: snow_rdst_sibl_tile
    dsl4jsb_Real2D_onChunk :: weq_snow_new_soil_sibl_tile
    dsl4jsb_Real2D_onChunk :: weq_snow_soil_sibl_tile
    dsl4jsb_Real2D_onChunk :: weq_snow_sibl_tile
    dsl4jsb_Real2D_onChunk :: snow_soil_dens_sibl_tile

    dsl4jsb_Real3D_onChunk :: soil_depth_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: vol_p_wilt_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_fc_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_sat_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_pwp_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: ice_soil_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: bclapp_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: wtr_diflat_sl_sibl_tile

    dsl4jsb_Real3D_onChunk :: t_soil_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: heat_cond_sl_sibl_tile
    dsl4jsb_Real3D_onChunk :: heat_flx_lat_sibl_tile

    dsl4jsb_Real2D_onChunk :: area_sibl_tile
    dsl4jsb_Real2D_onChunk :: h_surf_sibl_tile

    REAL(wp), DIMENSION(options%nc) :: fract           !< Fractions of water source tile []
    REAL(wp), DIMENSION(options%nc) :: fract_sibl_tile !< Fractions of water target tile []
    REAL(wp), DIMENSION(options%nc) :: wtr_pond_max    !< Maximum pond water depth       [m]
    REAL(wp), DIMENSION(options%nc) :: wtr_pond_tmp    !< Temporary pond water storage   [m]
    REAL(wp), DIMENSION(options%nc) :: fr_dpr_rel      !< ???
    LOGICAL,  DIMENSION(options%nc) :: l_pf_soil       !< Permafrost mask                []
    LOGICAL,  DIMENSION(options%nc) :: l_pf_soil_pi    !< Permafrost mask
    LOGICAL,  DIMENSION(options%nc) :: l_latdrain      !< T: allow lateral drainage      []

    REAL(wp), ALLOCATABLE, DIMENSION(:) ::           &
    & drain_bot_retention,                           & !< bottom_drainage from retention reservoir
    & h_offset,                                      & !< surface height difference between tiles
    & length_scale,                                  & !< Characteristic length scale of gradients between tiles
    & length_inter                                     !< Length of interface between tiles

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::         &
    & ice_imp,                                       & !< Soil layer ice impedance factor on source tile []
    & ice_imp_sibl_tile,                             & !< Soil layer ice impedance factor on target tile []
    & hyd_cond,                                      & !< Soil layer hyd. conductivity on source tile [m s-1]
    & hyd_cond_sibl_tile,                            & !< Soil layer hyd. conductivity on target tile [m s-1]
    & mpot_act,                                      & !< Soil layer matric potential on source tile [?]
    & mpot_act_sibl_tile,                            & !< Soil layer matric potential on target tile [?]
    & weq_wsat_sl,                                   & !< Maximum water content of soil layer on source tile [m]
    & weq_wsat_sl_sibl_tile,                         & !< Maximum water content of soil layer on target tile [m]
    & weq_field_cap_sl,                              & !< Field capacity of soil layer on source tile [m]
    & weq_field_cap_sl_sibl_tile,                    & !< Field capacity of soil layer on target tile [m]
    & weq_p_wilt_sl,                                 & !< Wilting point of soil layer on source tile [m]
    & weq_p_wilt_sl_sibl_tile                          !< Wilting point of soil layer on target tile [m]

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
    grid   => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    soil_e => Get_vgrid('soil_depth_energy')

    dsl4jsb_Get_config(MIX_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_,   wtr_latflow)            ! inout
    dsl4jsb_Get_var3D_onChunk(MIX_,   wtr_diflat_sl)          ! inout
    dsl4jsb_Get_var2D_onChunk(MIX_,   runoff_retention)       ! inout
    dsl4jsb_Get_var3D_onChunk(MIX_,   drain_lat_retention_sl) ! inout
    dsl4jsb_Get_var2D_onChunk(MIX_,   snow_rdst)              ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)                 ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drain_rock)             ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, drain_lat_sl)           ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)               ! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_new_soil)      ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)          ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)         ! inout

    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)           ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(MIX_, heat_flx_lat)        ! out

    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_perm_ref)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, area          )         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, h_surf        )         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, oro_stddev    )         ! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_, steepness)              ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, micro_slope)            ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, macro_slope)            ! out

    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth_max_ym1)     ! in
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)          ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_Cap_sl)     ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_p_wilt_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_fc_sl)       ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pwp_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)
      dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl)
    END IF

    ! Transport water from between tiles
    ! @todo currently, lateral transport is only considered on tiles that contain soil
    IF (dsl4jsb_Config(MIX_)%l_hydro_mix .AND. tile%contains_soil) THEN
      ! Allocate local variables
      ALLOCATE(drain_bot_retention(nc))
      ALLOCATE(h_offset(nc))
      ALLOCATE(length_scale(nc))
      ALLOCATE(length_inter(nc))
      ALLOCATE(ice_imp(nc,nsoil))
      ALLOCATE(ice_imp_sibl_tile(nc,nsoil))
      ALLOCATE(hyd_cond(nc,nsoil))
      ALLOCATE(hyd_cond_sibl_tile(nc,nsoil))
      ALLOCATE(mpot_act(nc,nsoil))
      ALLOCATE(mpot_act_sibl_tile(nc,nsoil))
      ALLOCATE(weq_wsat_sl(nc,nsoil))
      ALLOCATE(weq_wsat_sl_sibl_tile(nc,nsoil))
      ALLOCATE(weq_field_cap_sl(nc,nsoil))
      ALLOCATE(weq_field_cap_sl_sibl_tile(nc,nsoil))
      ALLOCATE(weq_p_wilt_sl(nc,nsoil))
      ALLOCATE(weq_p_wilt_sl_sibl_tile(nc,nsoil))

      ! This check should be moved to init routine as soon as enough stuff is there to justify one
      IF (tile%process_action(HYDRO_) /= ON_TILE_) THEN
        CALL finish(routine, 'MIX_ processes must be activate on the same tiles as HYDRO_')
      END IF

      !Initiaizing / Updating features of subgrid scale topography
      ! Do not forget lakes maximum pond depth !!!
      dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max) = MAX(0._wp, - h_surf(:))
      dsl4jsb_var2D_onChunk (HYDRO_, fract_pond_max) = dsl4jsb_var2D_onChunk (HYDRO_, fract_pond_max)
      ! Compute steepness parameter for surface runoff and lateral drainage
      ! Part I: First for depressions (only macro slope)
      ! Ridges are calculated below when memory of corresponding sibl tile exists
      IF (tile%contains_depression) THEN
        CALL calc_slopes_tiles(                             &
          ! in
          & nc,                                             & !< Number of gridcells
          & grid%nlat_g,                                    & !< (Effective) number of latitudes for steepness parameter
          & tile%contains_upland,                           & !< T: Upland HydroTile
          & dsl4jsb_Config(HYDRO_)%max_slope,               & !< maximum slope threshold [m/m]
          & oro_stddev(:),                                  & !< Standard deviation of orography [m]
          ! out
          & steepness(:),                                   & !< parameter defining the subgrid slope distribution [m/m]
          & micro_slope(:),                                 & !< slope on microtopography tiles [m/m]
          & macro_slope(:)                                  & !< slope on macrotopography tiles [m/m]
          & )
      END IF

      wtr_pond_max(:) = MAX(0._wp, dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max) &
        &                        - dsl4jsb_var2D_onChunk (HYDRO_, ice_pond))

      DO i = 1, nc
        IF (thaw_depth_max_ym1(i) < SUM(soil_depth_sl(i,:)) .AND. &  ! Thaw depth < soil depth ...
         &  thaw_depth_max_ym1(i) < SUM(soil_w%dz(:)) - soil_w%dz(nsoil)/2._wp - 0.01_wp) THEN ! ... but not at maximum.
          l_pf_soil(i)=.TRUE.
        ELSE
          l_pf_soil(i)=.FALSE.
        END IF
      END DO

      !!!!!!! HERE IT WOULD BE POSSIBLE TO PRESCRIBE A l_pf_soil
      !!!!!!! TO MAINTAIN SOME OF THE KEY EFFECTS OF FROZEN SOILS
      !!!!!!! EVEN AFTER PF HAS THAWED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!! l_pf_soil_pi JUST A DUMMY FOR NOW
      WHERE (fract_perm_ref(:) > 0.5_wp)
        l_pf_soil_pi(:)=.TRUE.
      ELSEWHERE
        l_pf_soil_pi(:)=.FALSE.
      END WHERE

      IF(dsl4jsb_Config(HYDRO_)%l_jsb_wet_test) THEN
        l_pf_soil(:)=l_pf_soil_pi(:)
      END IF

      ! Source tiles are always either upland or ridge tiles
      IF (tile%contains_upland .OR. tile%contains_ridge) THEN

        ! Get area fraction of source tile
        CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))

        ! Get number of tiles in this hierarchy level and loop over then
        nsibl = tile%parent_tile%Get_no_of_children()
        DO ic = 1, nsibl
          IF (ic == 1) THEN
            sibl_tile => tile%parent_tile%Get_first_child_tile()
          ELSE
            sibl_tile => sibl_tile%Get_next_sibling_tile()
          END IF

          ! Identify correct downhill target sibling for runoff
          l_latflow = (tile%contains_upland       .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_upland  .AND. sibl_tile%contains_depression) .OR. & ! upland ridge      --> upland depression
            &         (tile%contains_upland       .AND. tile%contains_depression      .AND. &
            &          sibl_tile%contains_lowland .AND. sibl_tile%contains_depression) .OR. & ! upland depression --> lowland depression
            &         (tile%contains_lowland      .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_lowland .AND. sibl_tile%contains_depression)        ! lowland ridge     --> lowland depression

          IF (l_latflow) THEN
            ! Get target tile area fraction
            CALL sibl_tile%Get_fraction(ics, ice, iblk, fract=fract_sibl_tile(:))
            ! Get memory MIX_ and HYDRO_ memory for target tile
            dsl4jsb_Get_memory_tile(MIX_,   sibl_tile)
            dsl4jsb_Get_memory_tile(HYDRO_, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_pond,    sibl_tile)

            l_lag = tile%contains_upland .AND. sibl_tile%contains_lowland
            retention_time = 36000._wp   ! Ten hours retention time (
                                         ! outflow = h_re * dt/retention_time)

            ! Move surface runoff from source to target tile
            CALL calc_runoff_transport(   &
              ! in
              & l_lag,                    & !< T: Apply lag
              & dtime,                    & !< model time step              [s]
              & retention_time,           & !< lagged flow retention time   [s]
              & wtr_pond_max(:),          & !< maximum pond water content   [m]
              & fract(:),                 & !< area fraction of source tile []
              & fract_sibl_tile(:),       & !< area fraction of target tile []
              ! inout
              & runoff(:),                & !< surface runoff               [kg m-2 s-1]
              & wtr_latflow(:),           & !< lateral flow at source tile  [kg m-2 s-1]
              & wtr_latflow_sibl_tile(:), & !< lateral flow at target tile  [kg m-2 s-1]
              & wtr_pond(:),              & !< pond storage at source tile  [m]
              & wtr_pond_sibl_tile(:),    & !< pond storage at target tile  [m]
              & runoff_retention(:)       & !< runoff retention storage     [m]
              & )
          END IF

          ! Identify correct downhill target sibling for drainage
          l_latflow = (tile%contains_upland       .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_upland  .AND. sibl_tile%contains_depression) .OR. & ! upland ridge      --> upland depression
            ! This condition should only be used if soils are frozen or large scale transport is explicitly required
            &         (tile%contains_upland       .AND. tile%contains_depression      .AND. &
            &          sibl_tile%contains_lowland .AND. sibl_tile%contains_depression) .OR. & ! upland depression --> lowland depression
            !
            &         (tile%contains_lowland      .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_lowland .AND. sibl_tile%contains_depression)        ! lowland ridge     --> lowland depression

          IF (l_latflow) THEN
            ! Get target tile area fraction
            CALL sibl_tile%Get_fraction(ics, ice, iblk, fract=fract_sibl_tile(:))
            ! Get memory MIX_ and HYDRO_ memory for target tile
            dsl4jsb_Get_memory_tile(MIX_,   sibl_tile)
            dsl4jsb_Get_memory_tile(HYDRO_, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow,    sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sl,    sibl_tile)

            l_lag = tile%contains_upland .AND. sibl_tile%contains_lowland
            retention_time = 360000._wp  ! Hundert hours retention time (
                                         ! outflow = h_re * dt/retention_time)

            ! Additional condition: long distance lateral drainage from upland to lowland tiles is allowed
            ! either under permafrost conditions or if l_macrotrans_subsurf is True
            DO i = 1, nc
              IF (tile%contains_upland .AND. tile%contains_depression) THEN    ! Long-distance and
                IF (dsl4jsb_Config(MIX_)%l_macrotrans_subsurf .OR. &           ! Switch true or
                  &  l_pf_soil(i)) THEN                                        ! Permafrost
                  l_latdrain(i) = .TRUE.
                ELSE                                                           ! Long-distance but not allowed
                  l_latdrain(i) = .FALSE.
                END IF
              ELSE                                                             ! Other distances
                l_latdrain(i) = .TRUE.
              END IF
            END DO

            ! Move lateral drainage from source to target tile
            WHERE (l_pf_soil(:))
              drain_bot_retention(:) = 0._wp
            ELSEWHERE
              drain_bot_retention(:) = k_brock * dtime
            END WHERE
            DO i = 1, nsoil
              !First bottom drainage is removed from reservoir
              DO j = 1, nc
                IF (drain_bot_retention(j) > 0._wp) THEN
                  IF (drain_bot_retention(j)   >= drain_lat_retention_sl(j,i)) THEN
                    drain_bot_retention(j)      = drain_bot_retention(j) - drain_lat_retention_sl(j,i)
                    drainage(j)                 = drainage(j)   + drain_lat_retention_sl(j,i) * rhoh2o / dtime
                    drain_rock(j)               = drain_rock(j) + drain_lat_retention_sl(j,i) * rhoh2o / dtime
                    drain_lat_retention_sl(j,i) = 0._wp
                  ELSE
                    drain_lat_retention_sl(j,i) = drain_lat_retention_sl(j,i) - drain_bot_retention(j)
                    drainage(j)                 = drainage(j)   + drain_bot_retention(j) * rhoh2o / dtime
                    drain_rock(j)               = drain_rock(j) + drain_bot_retention(j) * rhoh2o / dtime
                    drain_bot_retention(j)      = 0._wp
                  END IF
                END IF
              END DO
              !! Now water is moved between tiles
              CALL calc_drainage_transport( &
                ! in
                & l_lag,                                 & !< T: Apply lag
                & l_latdrain(:),                         & !< T: allow lateral drainage in gridcell
                & dtime,                                 & !< model time step                          [s]
                & retention_time,                        & !< lagged flow retention time               [s]
                & wtr_soil_fc_sl(:,i),                   & !< soil layer field capacity at source tile [m]
                & fract(:),                              & !< area fraction of source tile             []
                & fract_sibl_tile(:),                    & !< area fraction of target tile             []
                ! inout
                & drain_lat_sl(:,i),                     & !< lateral drainage at source tile          [kg m-2 s-1]
                & drainage(:),                           & !< total drainage at source tile            [kg m-2 s-1]
                & wtr_latflow(:),                        & !< lateral flow at source tile              [kg m-2 s-1]
                & wtr_latflow_sibl_tile(:),              & !< lateral flow at target tile              [kg m-2 s-1]
                & wtr_soil_sl(:,i),                      & !< soil water content at source tile        [m]
                & wtr_soil_sl_sibl_tile(:,i),            & !< soil water content at target tile        [m]
                & drain_lat_retention_sl(:,i)            & !< lateral drainage retention storage       [m]
                & )
            ENDDO
          END IF

          ! Identify correct downhill target sibling for lateral diffusion
          l_latflow = (tile%contains_upland       .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_upland  .AND. sibl_tile%contains_depression) .OR. & ! upland ridge  --> upland depression
            &         (tile%contains_lowland      .AND. tile%contains_ridge           .AND. &
            &          sibl_tile%contains_lowland .AND. sibl_tile%contains_depression)        ! lowland ridge --> lowland depression

          IF (l_latflow) THEN
            ! Get target tile area fraction
            CALL sibl_tile%Get_fraction(ics, ice, iblk, fract=fract_sibl_tile(:))
            ! Get memory MIX_ and HYDRO_ memory for target tile
            dsl4jsb_Get_memory_tile(MIX_,   sibl_tile)
            dsl4jsb_Get_memory_tile(HYDRO_, sibl_tile)
            dsl4jsb_Get_memory_tile(SSE_,   sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow,              sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(MIX_,   wtr_diflat_sl,            sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   snow_rdst,                sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_pond,                 sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, weq_pond,                 sibl_tile) ! in
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, weq_snow_new_soil,        sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, weq_snow_soil,            sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, weq_snow,                 sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, snow_soil_dens,           sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sl,              sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_fc_sl,           sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, soil_depth_sl,            sibl_tile)
            dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_porosity_sl,          sibl_tile) ! in
            dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_field_Cap_sl,         sibl_tile) ! in
            dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_p_wilt_sl,            sibl_tile) ! in
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_sat_sl,          sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_pwp_sl,          sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, ice_soil_sl,              sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, wtr_soil_pot_scool_sl,    sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, hyd_cond_sat_sl,          sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, matrix_pot_sl,            sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, bclapp_sl,                sibl_tile)
            dsl4jsb_Get_var3d_onChunk_tile_name(HYDRO_, pore_size_index_sl,       sibl_tile)

            dsl4jsb_Get_var3d_onChunk_tile_name(SSE_, t_soil_sl,                  sibl_tile)! in
            dsl4jsb_Get_var3d_onChunk_tile_name(SSE_, vol_heat_cap_sl,            sibl_tile)! in
            dsl4jsb_Get_var3d_onChunk_tile_name(SSE_, heat_cond_sl,               sibl_tile)! in
            dsl4jsb_Get_var3d_onChunk_tile_name(MIX_, heat_flx_lat,               sibl_tile)! out

            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, area,                     sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, h_surf,                   sibl_tile)

            DO i = 1, nc
              IF (fract(i) > 1.0e-10_wp .AND. fract_sibl_tile(i) > 1.0e-10_wp) THEN
             ! Initiaizing / Updating features of subgrid scale topography
             ! Part II: For ridges
             ! Get geometric relation between surface features
                CALL get_geometry(fract(i),              &
                  fract_sibl_tile(i),                    &
                  area(i),                               &
                  area_sibl_tile(i),                     &
                  h_surf(i),                             &
                  h_surf_sibl_tile(i),                   &
                  h_offset(i),                           &
                  length_scale(i),                       &
                  length_inter(i)                        &
                  )
              END IF
            END DO

            ! Calculating steepness based on macro an dmicro slopes
            CALL calc_slopes_tiles(                      &
            ! in
              & nc,                                             & !< Number of gridcells
              & grid%nlat_g,                                    & !< (Effective) number of latitudes for steepness parameter
              & tile%contains_upland,                           & !< T: Upland HydroTile
              & dsl4jsb_Config(HYDRO_)%max_slope,               & !< maximum slope threshold [m/m]
              & oro_stddev(:),        & !< Standard deviation of orography [m]
              ! out
              & steepness(:),                                   & !< parameter defining the subgrid slope distribution [m/m]
              & micro_slope(:),                                 & !< slope on microtopography tiles [m/m]
              & macro_slope(:),                                 & !< slope on macrotopography tiles [m/m]
              ! optional
              & h_offset(:),                                    & !< Maximum pond / depression depth [m]
              & length_scale(:),                                & !< Characteristic length scale of surface features (Radius circle (depression + elevation)
              & fract(:), fract_sibl_tile(:)                    & !< Fractions of source and target tile []
              & )

            ! Compute soil ice impedance, diffusivity and matrix potential ...
            ! ... for source tile
            weq_wsat_sl     (:,:) = MAX(vol_porosity_sl (:,:) * soil_depth_sl(:,:), 0._wp)
            weq_field_cap_sl(:,:) = MAX(vol_field_cap_sl(:,:) * soil_depth_sl(:,:), 0._wp)
            weq_p_wilt_sl   (:,:) = MAX(vol_p_wilt_sl(:,:)    * soil_depth_sl(:,:), 0._wp)
            CALL get_soilhyd_cond_properties(       &
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
              hyd_cond_sat_sl(:,:),                 & !< saturated hydraulic conductivity [m s-1]
              matrix_pot_sl(:,:),                   & !< matric potential
              bclapp_sl(:,:),                       & !< clapp hornberger b parameter
              pore_size_index_sl(:,:),              & !< pore size index
              ice_impedance=ice_imp(:,:),           & !< ice impedance factor
              K=hyd_cond(:,:),                      & !< hydraulic conductivity
              mpot_act=mpot_act(:,:)                & !< matrix potential
            )

            ! ... and target tile
            weq_wsat_sl_sibl_tile     (:,:) = MAX(vol_porosity_sl_sibl_tile (:,:) * soil_depth_sl_sibl_tile(:,:), 0._wp)
            weq_field_cap_sl_sibl_tile(:,:) = MAX(vol_field_cap_sl_sibl_tile(:,:) * soil_depth_sl_sibl_tile(:,:), 0._wp)
            weq_p_wilt_sl_sibl_tile   (:,:) = MAX(vol_p_wilt_sl_sibl_tile(:,:)    * soil_depth_sl_sibl_tile(:,:), 0._wp)
            CALL get_soilhyd_cond_properties(       &
              dsl4jsb_Config(HYDRO_)%soilhydmodel,  & !< Soil hydrology model
              1,                                    & !< interpolation methods (not applied in this case)
              nc, nsoil,                            & !< number of cells and soil layers
              soil_depth_sl_sibl_tile(:,:),         & !< soil layer depth   [m]
              wtr_soil_sl_sibl_tile(:,:),           & !< soil water content [m]
              ice_soil_sl_sibl_tile(:,:),           & !< soil ice content   [m]
              wtr_soil_pot_scool_sl_sibl_tile(:,:), & !< potential supercooled water [m]
              weq_wsat_sl_sibl_tile(:,:),           & !< maximum soil moisture storage (water + ice) [m]
              weq_field_cap_sl_sibl_tile(:,:),      & !< soil field capacity (water + ice) [m]
              weq_p_wilt_sl_sibl_tile(:,:),         & !< soil wilting point (water + ice) [m]
              hyd_cond_sat_sl_sibl_tile(:,:),       & !< saturated hydraulic conductivity [m s-1]
              matrix_pot_sl_sibl_tile(:,:),         & !< matric potential
              bclapp_sl_sibl_tile(:,:),             & !< clapp hornberger b parameter
              pore_size_index_sl_sibl_tile(:,:),    & !< pore size index
              ice_impedance=ice_imp_sibl_tile(:,:), & !< ice impedance factor
              K=hyd_cond_sibl_tile(:,:),            & !< hydraulic conductivity
              mpot_act=mpot_act_sibl_tile(:,:)      & !< matrix potential
            )

            ! Lateral soil water diffusion between ridge and depression tiles of the same macrotopography
            DO i = 1, nc
              IF (dsl4jsb_Config(MIX_)%l_itdiff .AND. &
                  fract(i) > 1.0e-10_wp .AND. fract_sibl_tile(i) > 1.0e-10_wp .AND. &
                  wtr_soil_fc_sl(i,1) > 1.0e-10_wp .AND. wtr_soil_fc_sl_sibl_tile(i,1) > 1.0e-10_wp) THEN

                IF (h_offset(i) >= 0._wp) THEN
                  CALL calc_intertile_diffusion(nsoil,                       & !< Number of soil layers
                     & dsl4jsb_Config(HYDRO_)%soilhydmodel,                  & !< Soil hydrology model
                     & dtime,                                                & !< time step                                        [s]
                     & abs(h_offset(i)),                                     & !< height difference between source and target tile [m]
                     & length_scale(i),                                      & !< Characteristic length of lateral gradients
                     & length_inter(i),                                      & !< Length of interface between tile (relative to 1m^2)
                     & fract(i),             fract_sibl_tile(i),             & !< tile fractions              []
                     & soil_depth_sl(i,:),   soil_depth_sl_sibl_tile(i,:),   & !< soil layer depths           [m]
                     & wtr_soil_fc_sl(i,:),  wtr_soil_fc_sl_sibl_tile(i,:),  & !< soil field capacity         [m]
                     & wtr_soil_sat_sl(i,:), wtr_soil_sat_sl_sibl_tile(i,:), & !< soil maximum water content  [m]
                     & wtr_soil_pwp_sl(i,:), wtr_soil_pwp_sl_sibl_tile(i,:), & !< soil wilting point          [m]
                     & matrix_pot_sl(i,:),   matrix_pot_sl_sibl_tile(i,:),   & !< soil matrix potential       [m]
                     & hyd_cond_sat_sl(i,:), hyd_cond_sat_sl_sibl_tile(i,:), & !< soil hydraulic conductivity [m s-1]
                     & bclapp_sl(i,:),       bclapp_sl_sibl_tile(i,:),       & !< soil bclapp parameter       []
                     & pore_size_index_sl(i,:), pore_size_index_sl_sibl_tile(i,:), & !< pore size index       []
                     & ice_imp(i,:),         ice_imp_sibl_tile(i,:),         & !< ice impedance factor        []
                     & hyd_cond(i,:),        hyd_cond_sibl_tile(i,:),        & !< hyd. conductivity           [m s-1]
                     & mpot_act(i,:),        mpot_act_sibl_tile(i,:),        & !< matric potential            [?]
                     & ice_soil_sl(i,:),     ice_soil_sl_sibl_tile(i,:),     & !< soil ice storage            [m]
                     & MIN(wtr_soil_pot_scool_sl(i,:), wtr_soil_sl(i,:)),    & !< supercooled water           [m]
                     & MIN(wtr_soil_pot_scool_sl_sibl_tile(i,:), wtr_soil_sl_sibl_tile(i,:)), &  !< supercooled water [m]
                     & wtr_pond(i),          wtr_pond_sibl_tile(i),          & !< pond water storage          [m]
                     & wtr_latflow(i),       wtr_latflow_sibl_tile(i),       & !< lateral water flow          [kg m-2 s-1]
                     & wtr_soil_sl(i,:),     wtr_soil_sl_sibl_tile(i,:),     & !< soil water storage          [m]
                     & wtr_diflat_sl(i,:),   wtr_diflat_sl_sibl_tile(i,:)    & !< lateral diffusion flux      [kg m-2 s-1]
                     & )
                ELSE
                  CALL calc_intertile_diffusion(nsoil,                       & !< Number of soil layers
                     & dsl4jsb_Config(HYDRO_)%soilhydmodel,                  & !< Soil hydrology model
                     & dtime,                                                & !< time step                                        [s]
                     & abs(h_offset(i)),                                     & !< height difference between source and target tile [m]
                     & length_scale(i),                                      & !< Characteristic length of lateral gradients
                     & length_inter(i),                                      & !< Length of interface between tile (relative to 1m^2)
                     & fract_sibl_tile(i),             fract(i),             & !< tile fractions              []
                     & soil_depth_sl_sibl_tile(i,:),   soil_depth_sl(i,:),   & !< soil layer depths           [m]
                     & wtr_soil_fc_sl_sibl_tile(i,:),  wtr_soil_fc_sl(i,:),  & !< soil field capacity         [m]
                     & wtr_soil_sat_sl_sibl_tile(i,:), wtr_soil_sat_sl(i,:), & !< soil maximum water content  [m]
                     & wtr_soil_pwp_sl_sibl_tile(i,:), wtr_soil_pwp_sl(i,:), & !< soil wilting point          [m]
                     & matrix_pot_sl_sibl_tile(i,:),   matrix_pot_sl(i,:),   & !< soil matrix potential       [m]
                     & hyd_cond_sat_sl_sibl_tile(i,:), hyd_cond_sat_sl(i,:), & !< soil hydraulic conductivity [m s-1]
                     & bclapp_sl_sibl_tile(i,:),       bclapp_sl(i,:),       & !< soil bclapp parameter       []
                     & pore_size_index_sl_sibl_tile(i,:), pore_size_index_sl(i,:), & !< pore size index       []
                     & ice_imp_sibl_tile(i,:),         ice_imp(i,:),         & !< ice impedance factor        []
                     & hyd_cond_sibl_tile(i,:),        hyd_cond(i,:),        & !< hyd. conductivity           [m s-1]
                     & mpot_act_sibl_tile(i,:),        mpot_act(i,:),        & !< matric potential            [?]
                     & ice_soil_sl_sibl_tile(i,:),     ice_soil_sl_sibl_tile(i,:),     & !< soil ice storage            [m]
                     & MIN(wtr_soil_pot_scool_sl_sibl_tile(i,:), wtr_soil_sl_sibl_tile(i,:)),    & !< supercooled water           [m]
                     & MIN(wtr_soil_pot_scool_sl(i,:), wtr_soil_sl(i,:)), &  !< supercooled water [m]
                     & wtr_pond_sibl_tile(i),          wtr_pond(i),          & !< pond water storage          [m]
                     & wtr_latflow_sibl_tile(i),       wtr_latflow(i),       & !< lateral water flow          [kg m-2 s-1]
                     & wtr_soil_sl_sibl_tile(i,:),     wtr_soil_sl(i,:),     & !< soil water storage          [m]
                     & wtr_diflat_sl_sibl_tile(i,:),   wtr_diflat_sl(i,:)    & !< lateral diffusion flux      [kg m-2 s-1]
                     & )
                END IF
              ELSE
                wtr_diflat_sl(i,:)           = 0._wp
                wtr_diflat_sl_sibl_tile(i,:) = 0._wp
              END IF
            ENDDO
            ! Lateral snow redistribution between ridge and depression tiles of the same macrotopography
            DO i = 1, nc
              IF (fract(i) > 1.0e-10_wp .AND. fract_sibl_tile(i) > 1.0e-10_wp) THEN
                IF (h_offset(i) >= 0._wp) THEN
                  CALL calc_snow_redistribution(dsl4jsb_Config(SSE_)%l_dynsnow, & !< snow properties
                     & dtime,                                                & !< time step                   [s]
                     & abs(h_offset(i)),                                     & !< height difference between source and target tile [m]
                     & fract(i),             fract_sibl_tile(i),             & !< tile fractions              []
                     & weq_pond(i),          weq_pond_sibl_tile(i),          & !< pond water storage          [m]
                     & weq_snow_soil(i),     weq_snow_soil_sibl_tile(i),     & !< snow on soil                [m]
                     & weq_snow_new_soil(i), weq_snow_new_soil_sibl_tile(i), & !< new snow on soil            [m]
                     & weq_snow(i),          weq_snow_sibl_tile(i),          & !< snow (soil + can)           [m]
                     & snow_soil_dens(i),    snow_soil_dens_sibl_tile(i),    & !< snow density                [m]
                     & snow_rdst(i),         snow_rdst_sibl_tile(i)          & !< lateral snow redistribution flux  [kg m-2 s-1]
                     & )
                ELSE
                  CALL calc_snow_redistribution(dsl4jsb_Config(SSE_)%l_dynsnow, & !< snow properties
                     & dtime,                                                & !< time step                   [s]
                     & abs(h_offset(i)),           & !< height difference between source and target tile [m]
                     & fract_sibl_tile(i),             fract(i),             & !< tile fractions              []
                     & weq_pond_sibl_tile(i),          weq_pond(i),          & !< pond water storage          [m]
                     & weq_snow_soil_sibl_tile(i),     weq_snow_soil(i),     & !< snow on soil                [m]
                     & weq_snow_new_soil_sibl_tile(i), weq_snow_new_soil(i), & !< new snow on soil            [m]
                     & weq_snow_sibl_tile(i),          weq_snow(i),          & !< snow (soil + can)           [m]
                     & snow_soil_dens_sibl_tile(i),    snow_soil_dens(i),    & !< snow density                [m]
                     & snow_rdst_sibl_tile(i),         snow_rdst(i)          & !< lateral snow redistribution flux  [kg m-2 s-1]
                     & )
                END IF
              END IF
            ENDDO
            DO i = 1, nc
              IF (dsl4jsb_Config(MIX_)%l_itdiff .AND. &
                  fract(i) > 1.0e-10_wp .AND. fract_sibl_tile(i) > 1.0e-10_wp) THEN
                IF (h_offset(i) >= 0._wp) THEN
                  CALL diagnose_intertile_heatflux(soil_e%n_levels, dtime,       &
                     & soil_e%dz(:),         abs(h_offset(i)),                   &
                     & length_scale(i),      length_inter(i),                    &
                     & fract(i),             fract_sibl_tile(i),                 &
                     & t_soil_sl(i,:),       t_soil_sl_sibl_tile(i,:),           &
                     & heat_cond_sl(i,:),    heat_cond_sl_sibl_tile(i,:),        &
                     & vol_heat_cap_sl(i,:), vol_heat_cap_sl_sibl_tile(i,:),     &
                     & heat_flx_lat(i,:),    heat_flx_lat_sibl_tile(i,:)         &
                     & )
                ELSE
                  CALL diagnose_intertile_heatflux(soil_e%n_levels, dtime,       &
                     & soil_e%dz(:),                   abs(h_offset(i)),         &
                     & length_scale(i),                length_inter(i),          &
                     & fract_sibl_tile(i),             fract(i),                 &
                     & t_soil_sl_sibl_tile(i,:),       t_soil_sl(i,:),           &
                     & heat_cond_sl_sibl_tile(i,:),    heat_cond_sl(i,:),        &
                     & vol_heat_cap_sl_sibl_tile(i,:), vol_heat_cap_sl(i,:),     &
                     & heat_flx_lat_sibl_tile(i,:),    heat_flx_lat(i,:)         &
                     & )
                END IF
              ELSE
                heat_flx_lat(i,:)           = 0._wp
                heat_flx_lat_sibl_tile(i,:) = 0._wp
              END IF

            END DO
          END IF
        END DO
      END IF

      ! Runoff on lowland depressions is only generated if pond storage is at maximum
      IF (tile%contains_lowland .AND. tile%contains_depression) THEN
        WHERE (wtr_pond(:) < wtr_pond_max(:))
            wtr_pond_tmp(:) = wtr_pond(:) + runoff(:) * dtime / rhoh2o
            runoff(:)       = MAX(0._wp, wtr_pond_tmp(:) - wtr_pond_max(:)) * rhoh2o / dtime
            wtr_pond(:)     = MIN(wtr_pond_tmp(:), wtr_pond_max(:))
        END WHERE
      END IF

      ! Deallocate local fields
      DEALLOCATE(h_offset, length_scale, length_inter, ice_imp,ice_imp_sibl_tile, hyd_cond,   &
        & hyd_cond_sibl_tile, mpot_act,mpot_act_sibl_tile, weq_wsat_sl,weq_wsat_sl_sibl_tile, &
        & weq_field_cap_sl,weq_field_cap_sl_sibl_tile, weq_p_wilt_sl, weq_p_wilt_sl_sibl_tile)

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
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_hydro_mix'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(MIX_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latflow,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   runoff_retention,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,                 weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drain_rock,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drain_lat_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   drain_lat_retention_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_diflat_sl,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,            weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_new_soil,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   snow_rdst,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MIX_,   heat_flx_lat,           weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_perm_ref,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, area,                   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, h_surf,                 weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, micro_slope,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, macro_slope,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, steepness,              weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond_max,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_pond_max,         weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_hydro_mix

#endif
END MODULE mo_mix_interface
