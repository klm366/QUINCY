!> QUINCY soil-biogeochemistry process interface
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of tasks for the soil-biogeochemistry process incl. update and aggregate routines
!>
MODULE mo_q_sb_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_grid_class,              ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, SB_, HYDRO_, SSE_, VEG_, A2L_

  ! "Integrate" routines
  USE mo_q_sb_rate_modifier,          ONLY: update_sb_rate_modifier_real => update_sb_rate_modifier
  USE mo_q_sb_ssm_main,               ONLY: update_sb_simple_model_real => update_sb_simple_model
  USE mo_q_sb_jsm_main,               ONLY: update_sb_jsm_real => update_sb_jsm
  USE mo_q_sb_update_pools,           ONLY: update_sb_pools_real => update_sb_pools

  dsl4jsb_Use_processes SB_, VEG_, HYDRO_
  dsl4jsb_Use_config(SB_)
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_sb_tasks_quincy
  PUBLIC :: global_q_sb_diagnostics

  ! ======================================================================================================= !
  !>Type definition: reset_sb_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_sb_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_sb_fluxes  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_sb_fluxes    !< Aggregates computed task variables
  END TYPE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Type definition: read_soil_forcing task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_read_soil_forcing
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_read_soil_forcing  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_read_soil_forcing    !< Aggregates computed task variables
  END TYPE tsk_read_soil_forcing

  ! ======================================================================================================= !
  !>Type definition: sb_rate_modifier task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_rate_modifier
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_rate_modifier       !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_rate_modifier    !< Aggregates computed task variables
  END TYPE tsk_sb_rate_modifier

  ! ======================================================================================================= !
  !>Type definition: sb_simple_model task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_simple_model
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_simple_model        !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_simple_model     !< Aggregates computed task variables
  END TYPE tsk_sb_simple_model

  ! ======================================================================================================= !
  !>Type definition: sb_jsm task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_jsm
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_jsm      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_jsm   !< Aggregates computed task variables
  END TYPE tsk_sb_jsm

  ! ======================================================================================================= !
  !>Type definition: sb_pools task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_pools         !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_pools      !< Aggregates computed task variables
  END TYPE tsk_sb_pools

  ! ======================================================================================================= !
  !>Type definition: time_average_soilbiogeochemistry task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_soilbiogeochemistry
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_soilbiogeochemistry    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_soilbiogeochemistry !< Aggregates computed task variables
  END TYPE tsk_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !>Constructor interface: reset_sb_fluxes task
  !>
  INTERFACE tsk_reset_sb_fluxes
    PROCEDURE Create_task_reset_sb_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor interface: read_soil_forcing task
  !>
  INTERFACE tsk_read_soil_forcing
    PROCEDURE Create_task_read_soil_forcing        !< Constructor function for task
  END INTERFACE tsk_read_soil_forcing

  ! ======================================================================================================= !
  !>Constructor interface: sb_rate_modifier task
  !>
  INTERFACE tsk_sb_rate_modifier
    PROCEDURE Create_task_sb_rate_modifier         !< Constructor function for task
  END INTERFACE tsk_sb_rate_modifier

  ! ======================================================================================================= !
  !>Constructor interface: sb_simple_model task
  !>
  INTERFACE tsk_sb_simple_model
    PROCEDURE Create_task_sb_simple_model         !< Constructor function for task
  END INTERFACE tsk_sb_simple_model

  ! ======================================================================================================= !
  !>Constructor interface: update_sb_jsm task
  !>
  INTERFACE tsk_sb_jsm
    PROCEDURE Create_task_sb_jsm         !< Constructor function for task
  END INTERFACE tsk_sb_jsm

  ! ======================================================================================================= !
  !>Constructor interface: sb_pools task
  !>
  INTERFACE tsk_sb_pools
    PROCEDURE Create_task_sb_pools         !< Constructor function for task
  END INTERFACE tsk_sb_pools

  ! ======================================================================================================= !
  !>Constructor interface: time_average_soilbiogeochemistry task
  !>
  INTERFACE tsk_time_average_soilbiogeochemistry
    PROCEDURE Create_task_update_time_average_soilbiogeochemistry         !< Constructor function for task
  END INTERFACE tsk_time_average_soilbiogeochemistry

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_sb_interface'

CONTAINS

  ! ======================================================================================================= !
  !>Register tasks: SB_ process
  !>
  SUBROUTINE Register_sb_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_sb_fluxes(model_id))
    CALL this%Register_task(tsk_read_soil_forcing(model_id))
    CALL this%Register_task(tsk_sb_rate_modifier(model_id))
    CALL this%Register_task(tsk_sb_simple_model(model_id))
    CALL this%Register_task(tsk_sb_jsm(model_id))
    CALL this%Register_task(tsk_sb_pools(model_id))
    CALL this%Register_task(tsk_time_average_soilbiogeochemistry(model_id))
  END SUBROUTINE Register_sb_tasks_quincy

  ! ======================================================================================================= !
  !>Constructor: reset_sb_fluxes task
  !>
  FUNCTION Create_task_reset_sb_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_sb_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_sb_fluxes', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor: read_soil_forcing task
  !>
  FUNCTION Create_task_read_soil_forcing(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_read_soil_forcing::return_ptr)
    CALL return_ptr%Construct(name='read_soil_forcing', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_read_soil_forcing

  ! ======================================================================================================= !
  !>Constructor: sb_rate_modifier task
  !>
  FUNCTION Create_task_sb_rate_modifier(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_rate_modifier::return_ptr)
    CALL return_ptr%Construct(name='sb_rate_modifier', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_rate_modifier

  ! ======================================================================================================= !
  !>Constructor: sb_simple_model task
  !>
  FUNCTION Create_task_sb_simple_model(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_simple_model::return_ptr)
    CALL return_ptr%Construct(name='sb_simple_model', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_simple_model

  ! ======================================================================================================= !
  !>Constructor: sb_jsm task
  !>
  FUNCTION Create_task_sb_jsm(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_jsm::return_ptr)
    CALL return_ptr%Construct(name='sb_jsm', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_jsm

  ! ======================================================================================================= !
  !>Constructor: sb_pools task
  !>
  FUNCTION Create_task_sb_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_pools::return_ptr)
    CALL return_ptr%Construct(name='sb_pools', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_pools

  ! ======================================================================================================= !
  !>Constructor: tsk_time_average_soilbiogeochemistry task
  !>
  FUNCTION Create_task_update_time_average_soilbiogeochemistry(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_soilbiogeochemistry::return_ptr)
    CALL return_ptr%Construct(name='tavrg_soilbiogeochemistry', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_soilbiogeochemistry

  ! ------------------------------------------------------------------------------------------------------- !
  ! update_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "integrate": reset_sb_fluxes task
  !>
  SUBROUTINE update_reset_sb_fluxes(tile, options)

    USE mo_sb_util, ONLY: real_reset_sb_fluxes => reset_sb_fluxes
    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_sb_fluxes'

    iblk = options%iblk
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_reset_sb_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_sb_fluxes

  ! ======================================================================================================= !
  !>calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE update_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_sb_util,             ONLY: calculate_time_average_soilbiogeochemistry
    dsl4jsb_Use_processes SB_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_soilbiogeochemistry(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !>Implementation of "integrate": read_soil_forcing task
  !>
  SUBROUTINE update_read_soil_forcing(tile, options)

    USE mo_jsb_io_netcdf,           ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
    USE mo_lnd_bgcm_idx
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID
    USE mo_jsb_time,                ONLY: t_datetime, get_date_components, get_year_day

    ! dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model            !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store       !< the bgcm store of this tile
    TYPE(t_input_file)              :: input_file       !< file with to be read veg state variables
    TYPE(t_datetime),       POINTER :: current_datetime !< date and time at the beginning of this time step
    REAL(wp)                        :: dtime, year_instant
    INTEGER                         :: iblk    !< Dimensions of the horizontal grid
    INTEGER                         :: veg_pool_idx, litterfall_pool_idx
    INTEGER                         :: year, month, day, hour, minute, second, yearday, &
                                    &  ith_timestep     !< the row number of current timestep in soil forcing
    REAL(wp), POINTER               :: veg_pool_mt_domain(:,:,:,:)      !< dim: compartments, elements, nc, nblks -- on domain!
    REAL(wp), POINTER               :: veg_litterfall_mt_domain(:,:,:,:)      !< dim: compartments, elements, nc, nblks -- on domain!
    REAL(wp), POINTER :: &
      & ptr_1D(:    ),   &                      !< Pointer to 1d variable
      & ptr_2D(:,  :),   &                      !< Pointer to 2d variable
      & ptr_3D(:,:,:),   &                      !< Pointer to 3d variable
      & ptr_4D(:,:,:,:)                         !< Pointer to 4d variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_read_soil_forcing'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(SB_)
    ! A2L_
    dsl4jsb_Real2D_onDomain      :: nhx_deposition
    dsl4jsb_Real2D_onDomain      :: noy_deposition
    dsl4jsb_Real2D_onDomain      :: nhx_n15_deposition
    dsl4jsb_Real2D_onDomain      :: noy_n15_deposition
    dsl4jsb_Real2D_onDomain      :: p_deposition
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: demand_uptake_p_tuptake_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain      :: root_fraction_sl
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain      :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onDomain      :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onDomain      :: wtr_soil_sl
    dsl4jsb_Real3D_onDomain      :: vol_porosity_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain      :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    iblk = options%iblk
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  p_deposition)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var3D_onDomain(VEG_,  root_fraction_sl)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onDomain(SSE_, t_soil_sl)

    !> get the timestep information of current timestep
    current_datetime = options%current_datetime
    dtime   = options%dtime
    CALL get_date_components(current_datetime, year, month, day, hour, minute, second)
    year_instant = get_year_day(current_datetime)
    yearday = INT(year_instant)
    ith_timestep = IMerge_TIME2Sec(yearday, hour, minute, second)/int(dtime) + 1 ! the row number of current timestep in forcing file

    !> 1.0 the first soil forcing file
    ! ----------------------------------------------------------------------------------------------------- !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%soilphysics_forcing_filename), model%grid_id)

    !> read the 3D forcing variables
    !> - Fraction of fineroot
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='veg_root_fraction_sl_'//tile%name,          &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    root_fraction_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Fraction of vertical transport
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_frac_wtr_transp_down_sl_'//tile%name,   &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    frac_wtr_transp_down_sl(:,:,:) = ptr_4D(:,:,:,1)

    !> - Fraction of lateral water loss
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_frac_w_lat_loss_sl_'//tile%name,        &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    frac_w_lat_loss_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Water content
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_wtr_soil_sl_'//tile%name,              &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    wtr_soil_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Volumetric water content
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_vol_porosity_sl_'//tile%name,          &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    vol_porosity_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Soil temperature
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_t_soil_sl_'//tile%name,                &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    t_soil_sl(:,:,:) = ptr_4D(:,:,:,1)

    !> read the 2D forcing variables
    !> - nhx deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_nhx_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    nhx_deposition(:,:) = ptr_3D(:,:,1)

    !> - noy deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_noy_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    noy_deposition(:,:) = ptr_3D(:,:,1)

    !> - nhx n15 deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_nhx_n15_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    nhx_n15_deposition(:,:) = ptr_3D(:,:,1)

    !> - noy n15 deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_noy_n15_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    noy_n15_deposition(:,:) = ptr_3D(:,:,1)

    !> - p deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_p_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    p_deposition(:,:) = ptr_3D(:,:,1)

    !> - plant N uptake demand
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_demand_uptake_n_tuptake_mavg_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    demand_uptake_n_tuptake_mavg(:,:) = ptr_3D(:,:,1)

    !> - plant P uptake demand
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_demand_uptake_p_tuptake_mavg_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    demand_uptake_p_tuptake_mavg(:,:) = ptr_3D(:,:,1)

    CALL input_file%Close()

    !> 2.0 the second soil forcing file
    ! ----------------------------------------------------------------------------------------------------- !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%veglitterfall_forcing_filename), model%grid_id)
    IF (ASSOCIATED(tile%bgcm_store)) THEN
      bgcm_store => tile%bgcm_store
      litterfall_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile%name, routine)
      veg_litterfall_mt_domain => bgcm_store%store_2l_2d_bgcms(bgcm_store%idx_in_store(litterfall_pool_idx))%mt_2l_2d_bgcm(:,:,:,:)
      veg_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_POOL_ID, tile%name, routine)
      veg_pool_mt_domain => bgcm_store%store_2l_2d_bgcms(bgcm_store%idx_in_store(veg_pool_idx))%mt_2l_2d_bgcm(:,:,:,:)
    END IF

    !> - veg fineroot C
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_pool_fine_root_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_pool_mt_domain(ix_fine_root, ixC, :,:) = ptr_3D(:,:,1)

    !> - litterfall leaf
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_leaf_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_leaf, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall fine root
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fine_root_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fine_root, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall coarse root
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_coarse_root_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_coarse_root, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall sap wood
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_sap_wood_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_sap_wood, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall heart wood
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_heart_wood_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_heart_wood, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall fruit
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_fruit_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_fruit, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall labile pool
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_labile_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_labile, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall reserve pool
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_reserve_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_reserve, ixN15, :,:) = ptr_3D(:,:,1)

    !> - litterfall seed bed
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixC, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_N_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixN, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_P_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixP, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_C13_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixC13, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_C14_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixC14, :,:) = ptr_3D(:,:,1)
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_litterfall_seed_bed_N15_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_litterfall_mt_domain(ix_seed_bed, ixN15, :,:) = ptr_3D(:,:,1)

    CALL input_file%Close()

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
    ! ! IN
    ! !$ACC UPDATE DEVICE(root_fraction_sl(:,:,:), num_sl_above_bedrock(:,:), bulk_dens_sl(:,:,:), vol_porosity_sl(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(soil_depth_sl(:,:,:), soil_lay_width_sl(:,:,:), soil_lay_depth_center_sl(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(soil_lay_depth_ubound_sl(:,:,:), silt_sl(:,:,:), clay_sl(:,:,:)) ASYNC(1)

    ! ! IN-OUT
    ! !$ACC UPDATE DEVICE(qmax_org_min_sl(:,:,:), qmax_po4_min_sl(:,:,:), qmax_nh4_min_sl(:,:,:), qmax_po4_om_sl(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(volume_min_sl(:,:,:)) ASYNC(1)

    ! ! SB_
    ! !$ACC UPDATE DEVICE(nh4_assoc(:,:,:), nh4_solute(:,:,:), no3_solute(:,:,:), po4_solute(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(nh4_n15_assoc(:,:,:), nh4_n15_solute(:,:,:), no3_n15_solute(:,:,:), bulk_soil_carbon_sl(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(soil_litter_carbon_sl(:,:,:), bulk_dens_corr_sl(:,:,:), qmax_org(:,:,:), qmax_po4(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(qmax_nh4(:,:,:), qmax_fast_po4(:,:,:), qmax_slow_po4(:,:,:), km_fast_po4(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(km_slow_po4(:,:,:), km_adsorpt_po4_sl(:,:,:), km_adsorpt_nh4_sl(:,:,:), k_bioturb(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(po4_assoc_fast(:,:,:), po4_assoc_slow(:,:,:), po4_occluded(:,:,:), ph_sl(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(Qmax_AlFe_cor(:,:,:), po4_primary(:,:,:), microbial_cue_eff_tmic_mavg(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(microbial_nue_eff_tmic_mavg(:,:,:), microbial_pue_eff_tmic_mavg(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(vmax_weath_mineral_sl(:,:,:), enzyme_frac_poly(:,:,:), enzyme_frac_residue(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(enzyme_frac_poly_c(:,:,:), enzyme_frac_poly_n(:,:,:), enzyme_frac_poly_p(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(dom_cn(:,:,:), dom_cp(:,:,:), fact_n_status_mic_c_growth(:,:,:)) ASYNC(1)
    ! !$ACC UPDATE DEVICE(fact_p_status_mic_c_growth(:,:,:), enzyme_frac_AP(:,:,:)) ASYNC(1)

    ! !$ACC UPDATE DEVICE(sb_pool_mt_domain(:,:,:,:,:)) ASYNC(1)
  END SUBROUTINE update_read_soil_forcing
  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_sb_rate_modifier(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_rate_modifier_real(tile, options)

  END SUBROUTINE update_sb_rate_modifier

  SUBROUTINE update_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_pools_real(tile, options)

  END SUBROUTINE update_sb_pools

  SUBROUTINE update_sb_jsm(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_jsm_real(tile, options)

  END SUBROUTINE update_sb_jsm

  SUBROUTINE update_sb_simple_model(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_simple_model_real(tile, options)

  END SUBROUTINE update_sb_simple_model

  ! ------------------------------------------------------------------------------------------------------- !
  ! aggregate_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_sb_fluxes task
  !>
  SUBROUTINE aggregate_reset_sb_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_sb_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": read_soil_forcing task
  !>
  SUBROUTINE aggregate_read_soil_forcing(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_read_soil_forcing'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_read_soil_forcing

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_rate_modifier task
  !>
  SUBROUTINE aggregate_sb_rate_modifier(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_rate_modifier'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_rate_modifier

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_simple_model task
  !>
  SUBROUTINE aggregate_sb_simple_model(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_simple_model'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_simple_model

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_jsm task
  !>
  SUBROUTINE aggregate_sb_jsm(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_jsm'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_jsm

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_pools task
  !>
  SUBROUTINE aggregate_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)

    ! fluxes
    dsl4jsb_Aggregate_onChunk(SB_, het_respiration              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, ecosystem_total_n_loss       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, fertiliser_nh4               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, fertiliser_no3               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, fertiliser_po4               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, emission_noy                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, emission_n2o                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, emission_n2                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, net_mineralisation_nh4       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, net_mineralisation_no3       , weighted_by_fract)

    ! pools
    ! 2D & 3D mixed
    dsl4jsb_Aggregate_onChunk(SB_, total_soil_n                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_soil_inorg_n           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_n              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_p              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c13            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c14            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_n15            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_c      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_n      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_p      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_n    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_p    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c13  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c14  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_n15  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_ag_woody_litter_c   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_ag_woody_litter_n   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_ag_woody_litter_p   , weighted_by_fract)

    ! VEG_ fluxes
    dsl4jsb_Aggregate_onChunk(VEG_,   net_biosphere_production      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   net_biosphere_production_c13  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   net_biosphere_production_c14  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   biological_n_fixation         , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_pools

  ! ======================================================================================================= !
  !>Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_soilbiogeochemistry(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)

    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_c        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_n        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_p        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_c13      , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_c14      , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_bg_soil_n15      , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_n         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_p         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c13       , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_c14       , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_ag_litter_n15       , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !>
  !> calculations of soil-biogeochemical diagnostic global sums
  !>
  SUBROUTINE global_q_sb_diagnostics(tile)

    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N
    USE mo_jsb_math_constants,      ONLY: one_day, one_year
    USE mo_sync,                    ONLY: global_sum_array
    USE mo_jsb_grid,                ONLY: Get_grid, Get_vgrid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid
    TYPE(t_jsb_vgrid), POINTER      :: soil_w

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp), ALLOCATABLE  :: sb_pool_total_c_2d(:,:)
    REAL(wp), ALLOCATABLE  :: total_soil_inorg_n_2d(:,:)
    REAL(wp), ALLOCATABLE  :: het_respiration_2d(:,:)
    REAL(wp), POINTER      :: total_som_c_gsum(:)
    REAL(wp), POINTER      :: total_inorg_n_gsum(:)
    REAL(wp), POINTER      :: rh_gsum(:)
    REAL(wp), POINTER      :: nbp_gsum(:)
    REAL(wp), POINTER      :: fnloss_gsum(:)
    REAL(wp), POINTER      :: fn2o_gsum(:)
    INTEGER                :: nsoil, is

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_q_sb_diagnostics'

    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real3D_onDomain :: sb_pool_total_c
    dsl4jsb_Real3D_onDomain :: total_soil_inorg_n
    dsl4jsb_Real3D_onDomain :: het_respiration
    dsl4jsb_Real2D_onDomain :: emission_n2o
    dsl4jsb_Real2D_onDomain :: ecosystem_total_n_loss
    dsl4jsb_Real2D_onDomain :: net_biosphere_production
    dsl4jsb_Real3D_onDomain :: soil_depth_sl

    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_w%n_levels
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')
    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onDomain(SB_,     sb_pool_total_c)              ! in
    dsl4jsb_Get_var3D_onDomain(SB_,     total_soil_inorg_n)           ! in
    dsl4jsb_Get_var3D_onDomain(SB_,     het_respiration)              ! in
    dsl4jsb_Get_var2D_onDomain(SB_,     emission_n2o)                 ! in
    dsl4jsb_Get_var2D_onDomain(SB_,     ecosystem_total_n_loss)       ! in
    dsl4jsb_Get_var2D_onDomain(VEG_,    net_biosphere_production)     ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_,  soil_depth_sl)                ! in

    total_som_c_gsum        => SB__mem%total_som_c_gsum%ptr(:)       ! out
    total_inorg_n_gsum      => SB__mem%total_inorg_n_gsum%ptr(:)     ! out
    rh_gsum                 => SB__mem%rh_gsum%ptr(:)                ! out
    nbp_gsum                => SB__mem%nbp_gsum%ptr(:)               ! out
    fn2o_gsum               => SB__mem%fn2o_gsum%ptr(:)              ! out
    fnloss_gsum             => SB__mem%fnloss_gsum%ptr(:)            ! out

    ! Domain Mask - to mask all halo cells for global sums (otherwise these cells are counted twice)
    ALLOCATE(in_domain(grid%nproma,grid%nblks))
    ALLOCATE(scaling(grid%nproma,grid%nblks))
    ALLOCATE(sb_pool_total_c_2d(grid%nproma, grid%nblks))
    ALLOCATE(total_soil_inorg_n_2d(grid%nproma, grid%nblks))
    ALLOCATE(het_respiration_2d(grid%nproma, grid%nblks))

    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ! Sum soil layers
    sb_pool_total_c_2d(:,:) = 0._wp
    total_soil_inorg_n_2d(:,:) = 0._wp
    het_respiration_2d(:,:) = 0._wp
    DO is = 1, nsoil
      sb_pool_total_c_2d(:,:)    = sb_pool_total_c_2d(:,:) + sb_pool_total_c(:,is,:)*soil_depth_sl(:,is,:)
      total_soil_inorg_n_2d(:,:) = total_soil_inorg_n_2d(:,:) + total_soil_inorg_n(:,is,:)*soil_depth_sl(:,is,:)
      het_respiration_2d(:,:)    = het_respiration_2d(:,:) + het_respiration(:,is,:)*soil_depth_sl(:,is,:)
    END DO

    ! Calculate global carbon inventories, if requested for output
    !  => Conversion from [mol(C)/m^2] to [PgC]
    !     1 mol C = molar_mass_C g C   => 1 Pg C = molar_mass_C * e-15
    scaling(:,:) = molar_mass_C * 1.e-15_wp * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (SB__mem%total_som_c_gsum%is_in_output) &
          &  total_som_c_gsum = global_sum_array(sb_pool_total_c_2d(:,:) * scaling(:,:))
    ! Calculate global carbon inventories, if requested for output
    !  => Conversion from [mol(N)/m^2] to [TgN]
    !     1 mol C = molar_mass_N g N   => 1 Pg N = molar_mass_N * e-12
    scaling(:,:) = molar_mass_N * 1.e-12_wp * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (SB__mem%total_som_c_gsum%is_in_output) &
          &  total_inorg_n_gsum = global_sum_array(total_soil_inorg_n_2d(:,:) * scaling(:,:))

    ! Calculate global fluxes, if requested for output
    !  => Conversion from [mumol(C)/m^2/s] to [PgC/yr]
    !     1 mol C = molar_mass_C g C   => 1 Pg C = molar_mass_C * e-15 * e-6
    !     1 s times one_day (60*60*24) time one_year (365)
    scaling(:,:) = molar_mass_C * 1.e-21_wp * one_day * one_year * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (SB__mem%rh_gsum%is_in_output)   &
          &  rh_gsum          = global_sum_array(het_respiration_2d(:,:) * scaling(:,:))
    IF (SB__mem%nbp_gsum%is_in_output)   &
          &  nbp_gsum         = global_sum_array(net_biosphere_production(:,:) * scaling(:,:))
    !  => Conversion from [mumol(N)/m^2/s] to [TgN/yr]
    !     1 mol C = molar_mass_N g N   => 1 Tg N = molar_mass_N * e-12 * e-6
    !     1 s times one_day (60*60*24) time one_year (365)
    scaling(:,:) = molar_mass_C * 1.e-18_wp * one_day * one_year * notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (SB__mem%fnloss_gsum%is_in_output) &
          &  fnloss_gsum      = global_sum_array(ecosystem_total_n_loss(:,:) * scaling(:,:))
    IF (SB__mem%fn2o_gsum%is_in_output) &
          &  fn2o_gsum        = global_sum_array(emission_n2o(:,:) * scaling(:,:))

    DEALLOCATE (scaling, in_domain, sb_pool_total_c_2d, het_respiration_2d)

  END SUBROUTINE global_q_sb_diagnostics

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  FUNCTION IMerge_TIME2Sec (kd, kh, km, ks) RESULT (ix)
    !+
    !
    ! IMerge_TIME2Sec [function, integer]
    !    merge day, hour, minute and second to internal time
    !    (
    !    day    [integer] input (day of year, 0...364)
    !    hour   [integer] input (hour of day, 0...23)
    !    minute [integer] input (minute of hour, 0...59)
    !    second [integer] input (second of minute, 0...59)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: kd, kh, km, ks
    INTEGER             :: ix

    ix = kd*86400 + kh*3600 + km*60 + ks

  END FUNCTION IMerge_TIME2Sec

#endif
END MODULE mo_q_sb_interface
