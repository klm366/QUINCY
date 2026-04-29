!> Contains the interfaces to the lex process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains interfaces to the lex processes
!>
!> The subroutine [[reset_lex_fluxes]] resets the lateral flux fields after these
!> have been used in HYDRO_ and SSE_ modules. The subroutine [[update_lex]] determines
!> the transport of water and heat between connected tiles. Note, that in the current
!> version water is transported directly in this routine, while heat transport is only
!> diagnosed and later applied in the SSE_ module. This should be harmonized later.
!>

MODULE mo_lex_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on, acc_stream
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes LEX_, HYDRO_, SSE_

  ! Use of process configurations
  dsl4jsb_Use_config(LEX_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_config(SSE_)

  ! Use of process memories
  dsl4jsb_Use_memory(LEX_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SSE_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_lex_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_interface'

  !> Type definition for reset_lex_fluxes task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_lex_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_lex_fluxes     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_lex_fluxes  !< Aggregates computed task variables
  END TYPE tsk_reset_lex_fluxes

  !> Constructor interface for reset_lex_fluxes task
  INTERFACE tsk_reset_lex_fluxes
    PROCEDURE Create_task_reset_lex_fluxes                        !< Constructor function for task
  END INTERFACE tsk_reset_lex_fluxes

  !> Type definition for lex task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_lex
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_lex                  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_lex               !< Aggregates computed task variables
  END TYPE tsk_lex

  !> Constructor interface for lex task
  INTERFACE tsk_lex
    PROCEDURE Create_task_lex                                     !< Constructor function for task
  END INTERFACE tsk_lex

CONTAINS

  ! =======================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for reset_lex_fluxes task
  !>
  FUNCTION Create_task_reset_lex_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "reset_lex_fluxes"

    ALLOCATE(tsk_reset_lex_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_lex_fluxes', process_id=LEX_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_lex_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for lex task
  !>
  FUNCTION Create_task_lex(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id      !< Model ID
    CLASS(t_jsb_process_task), POINTER    :: return_ptr    !< Instance of process task "lex"

    ALLOCATE(tsk_lex::return_ptr)
    CALL return_ptr%Construct(name='lex', process_id=LEX_, owner_model_id=model_id)

  END FUNCTION Create_task_lex

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for lex process
  !>
  SUBROUTINE Register_lex_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process         !< Process
    INTEGER,              INTENT(in)    :: model_id        !< Model ID

    CALL process%Register_task(tsk_reset_lex_fluxes(model_id))
    CALL process%Register_task(tsk_lex(model_id))

  END SUBROUTINE Register_lex_tasks

  ! ======================================================================================================= !
  !>
  !>#### Reset lateral exchange fluxes
  !>
  !> This subroutine resets the lateral flux fields after these have been used in the HYDRO_ and SSE_
  !> routines. It thereby ensures that fluxes are not accumulated over multiple time steps but reset
  !> properly at the begin of each time step.
  !>
  SUBROUTINE update_reset_lex_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile          !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options       !< Options of the current block

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_memory(LEX_)       !< Memory of the lateral exchange process
    dsl4jsb_Def_memory(HYDRO_)     !< Memory of the hydrology process
    dsl4jsb_Def_memory(SSE_)       !< Memory of the soil and snow energy process

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_runoff_in         !< Surface runoff received from adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: weq_snowblow_in       !< Snow blown received from adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: wtr_runoff_out        !< Surface runoff to adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: wtr_infil_out         !< Infiltration from ponds to soil in adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: weq_snowblow_out      !< Snow blown off to adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: water_to_soil_lex     !< Infiltratable water from connected tile [m]
    dsl4jsb_Real2D_onChunk :: wtr_lex_net_srf       !< Surface water exchange with connected tiles [kg m-2 s-1]

    dsl4jsb_Real3D_onChunk :: wtr_drainage_in_sl    !< Drainage received from adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_diffusion_in_sl   !< Water diffusion received from adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_drainage_out_sl   !< Drainage to adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_diffusion_out_sl  !< Water diffusion to adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: wtr_infil_in_sl       !< Infiltration into ponds received from adjacent tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: nrg_diffusion_in_sl   !< Energy diffusion received from adjacent tiles [W m-2]
    dsl4jsb_Real3D_onChunk :: nrg_convection_in_sl  !< Energy advection received from adjacent tiles [W m-2]
    dsl4jsb_Real3D_onChunk :: nrg_diffusion_out_sl  !< Energy diffusion to adjacent tiles [W m-2]
    dsl4jsb_Real3D_onChunk :: nrg_convection_out_sl !< Energy advection to adjacent tiles [W m-2]
    dsl4jsb_Real3D_onChunk :: wtr_lex_net_sl        !< Below ground water exchange with connected tiles [kg m-2 s-1]
    dsl4jsb_Real3D_onChunk :: nrg_lex_net_sl        !< Below ground energy exchange with connected tiles [W m-2]

    ! Local variables
    INTEGER  :: &
      & iblk,   &  !< Current block index
      & ics,    &  !< Index of first cell of block
      & ice        !< Index of last cell of block

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_lex_fluxes'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! Return, if process is not active on this tile.
    IF (.NOT. tile%Is_process_calculated(LEX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(LEX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(LEX_,   wtr_runoff_in)         ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   weq_snowblow_in)       ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   wtr_runoff_out)        ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   wtr_infil_out)         ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   weq_snowblow_out)      ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, water_to_soil_lex)     ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_lex_net_srf)       ! inout

    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_drainage_in_sl)    ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_diffusion_in_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_drainage_out_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_diffusion_out_sl)  ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_infil_in_sl)       ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_diffusion_in_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_convection_in_sl)  ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_diffusion_out_sl)  ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_convection_out_sl) ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_lex_net_sl)        ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_,   nrg_lex_net_sl)        ! inout

    !> Reset all lateral exchange fluxes to zero
    !>
    wtr_runoff_in        = 0._wp
    weq_snowblow_in      = 0._wp
    wtr_runoff_out       = 0._wp
    wtr_infil_out        = 0._wp
    weq_snowblow_out     = 0._wp
    wtr_drainage_in_sl   = 0._wp
    wtr_diffusion_in_sl  = 0._wp
    wtr_drainage_out_sl  = 0._wp
    wtr_diffusion_out_sl = 0._wp
    wtr_infil_in_sl      = 0._wp
    wtr_lex_net_srf      = 0._wp
    wtr_lex_net_sl       = 0._wp
    water_to_soil_lex    = 0._wp
    nrg_diffusion_in_sl  = 0._wp
    nrg_convection_in_sl = 0._wp
    nrg_diffusion_out_sl = 0._wp
    nrg_convection_out_sl= 0._wp
    nrg_lex_net_sl       = 0._wp

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_reset_lex_fluxes

  ! ======================================================================================================= !
  !>
  !>#### Aggregation of variables from reset_lateral_exchange_fluxes
  !>
  !> Lateral exchange as well as some hydrology and soil energy variables are aggregated in this routine.
  !>
  SUBROUTINE aggregate_reset_lex_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    ! Declare pointers to process memory
    dsl4jsb_Def_memory(LEX_)    !< Memory of the lateral exchange process
    dsl4jsb_Def_memory(HYDRO_)  !< Memory of the hydrology process
    dsl4jsb_Def_memory(SSE_)    !< Memory of the soil and snow energy process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions

    ! Local variables
    INTEGER :: &
      & iblk, &                 !< Number of current block
      & ics, &                  !< Index of first cell of current block
      & ice                     !< Index of last cell of current block

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_lex_fluxes'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LEX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_runoff_in,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   weq_snowblow_in,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_runoff_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_infil_out,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   weq_snowblow_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_infil_in_sl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_drainage_in_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_drainage_out_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_convection_in_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_convection_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   slope_eff,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fact_lag_srf,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fact_lag_blg,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   h_offset_ct,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   length_scale_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   length_iface_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fract_latflow_srf_ct,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fract_latflow_blg_ct,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   nrg_lex_net_sl,        weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_lex_fluxes

  ! ======================================================================================================= !
  !>
  !>#### Update routine for lateral exchange
  !>
  !> This routine determines the transport of surface runoff, lateral drainage, lateral ingress from ponds,
  !> the diffusive transport of water and energy and the redistribution of snow between connected tiles.
  !>
  SUBROUTINE update_lex(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_lex_process,            ONLY: calc_runoff_transport,           &
      &                                  calc_drainage_transport,         &
      &                                  calc_lat_fluxes_pond,            &
      &                                  calc_lat_fluxes_soil,            &
      &                                  calc_convection_heat,            &
      &                                  calc_snow_redistribution,        &
      &                                  shift_fluxes_sl,                 &
      &                                  remove_excess_water

    USE mo_hydro_process, ONLY: get_soilhyd_properties

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile           !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options        !< Options of current block

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(LEX_)              !< Configurable parameters of the lateral exchange process
    dsl4jsb_Def_config(HYDRO_)            !< Configurable parameters of the hydrology process
    dsl4jsb_Def_config(SSE_)              !< Configurable parameters of the soil and snow energy process

    dsl4jsb_Def_memory(LEX_)              !< Memory of the lateral exchange process
    dsl4jsb_Def_memory(HYDRO_)            !< Memory of the hydrology process
    dsl4jsb_Def_memory(SSE_)              !< Memory of the soil and snow energy process

    dsl4jsb_Def_memory_tile(LEX_, sibl)   !< Memory of the lateral exchange process for sibling tile
    dsl4jsb_Def_memory_tile(HYDRO_, sibl) !< Memory of the hydrology process for sibling tile
    dsl4jsb_Def_memory_tile(SSE_, sibl)   !< Memory of the soil and snow energy process for sibling tile

    ! Declare pointers to variables in memory for the current tile
    dsl4jsb_Real2D_onChunk :: &
      & wtr_latflow_srf,      &     !< Outflow of intermediate surface runoff storage on source tile [m]
      & wtr_lex_net_srf,      &     !< Surface water exchange with connected tiles [kg m-2 s-1]
      & weq_snowblow_in,      &     !< Snow blown received from adjacent tiles [kg m-2 s-1]
      & wtr_runoff_out,       &     !< Surface runoff to adjacent tiles [kg m-2 s-1]
      & weq_snowblow_out,     &     !< Snow blown off to adjacent tiles [kg m-2 s-1]
      & weq_pond,             &     !< Surface pond water and ice content [weq m]
      & wtr_pond,             &     !< Surface pond water content [m]
      & weq_snow,             &     !< Snow water equivalent on surface [m]
      & weq_snow_soil,        &     !< Snow water equivalent on soil [m]
      & new_snow_soil,        &     !< Snow water equivalent on soil after redistribution [m]
      & snow_soil_dens,       &     !< Snow density on soil [kg m-3]
      & runoff,               &     !< Surface runoff generated on tile [kg m-2 s-1]
      & drainage                    !< Subsurface drainage generated on tile [kg m-2 s-1]

    dsl4jsb_Real3D_onChunk ::    &
      &  microcon_ctc_ct,        &  !< micro-scale connectivity factor for direct contact []
      &  microcon_hfp_blg_ct,    &  !< micro-scale connectivity factor for below-ground hydrological flow paths []
      &  fract_latflow_srf_ct,   &  !< fraction of lateral surface flow to adjacent tile []
      &  fract_latflow_blg_ct,   &  !< fraction of lateral below-ground flow to adjacent tile []
      &  wtr_latflow_sl,         &  !< Outflow of intermediate subsurface lateral flow storage on source tile [m]
      &  soil_depth_sl,          &  !< Soil depth until bedrock within each layer [m]
      &  wtr_lex_net_sl,         &  !< Below ground water exchange with connected tiles [kg m-2 s-1]
      &  wtr_infil_in_sl,        &  !< Infiltration into ponds received from adjacent tiles [kg m-2 s-1]
      &  wtr_diffusion_in_sl,    &  !< Water diffusion received from adjacent tiles [kg m-2 s-1]
      &  drainage_sl,            &  !< Subsurface drainage generated on tile [kg m-2 s-1]
      &  wtr_drainage_out_sl,    &  !< Drainage to adjacent tiles [kg m-2 s-1]
      &  wtr_diffusion_out_sl,   &  !< Water diffusion to adjacent tiles [kg m-2 s-1]
      &  wtr_soil_sl,            &  !< Soil water content [m]
      &  ice_soil_sl,            &  !< Soil ice content [m]
      &  wtr_soil_pot_scool_sl,  &  !< Potentially supercooled soil water content [m]
      &  wtr_soil_sat_sl,        &  !< Soil moisture at saturation [m]
      &  wtr_soil_fc_sl,         &  !< Soil moisture at field capacity [m]
      &  wtr_soil_res_sl,        &  !< Soil moisture at residual water content [m]
      &  hyd_cond_sat_sl,        &  !< Hydraulic conductivity of soil layers [m s-1]
      &  vol_porosity_sl,        &  !< Volumetric soil porosity [m/m]
      &  vol_wres_sl,            &  !< Volumetric residual water content [m/m]
      &  matric_pot_sl,          &  !< Matric potential of soil layers at saturation [m]
      &  bclapp_sl,              &  !< Clapp and Hornberger exponent b of soil layers []
      &  pore_size_index_sl,     &  !< Pore size index of soil layers []
      &  t_soil_sl,              &  !< Soil layer temperature [K]
      &  heat_cond_sl,           &  !< Thermal conductivity of the soil layer [J m-1 s-1 K-1]
      &  vol_heat_cap_sl,        &  !< Volumetric heat capacity of the soil layer [J m-3 K-1]
      &  nrg_diffusion_in_sl,    &  !< Energy diffusion received from adjacent tiles [W m-2]
      &  nrg_convection_in_sl,   &  !< Energy advection received from adjacent tiles [W m-2]
      &  nrg_diffusion_out_sl,   &  !< Energy diffusion to adjacent tiles [W m-2]
      &  nrg_convection_out_sl,  &  !< Energy advection to adjacent tiles [W m-2]
      &  nrg_lex_net_sl,         &  !< Below ground energy exchange with connected tiles [W m-2]
      &  h_offset_ct,            &  !< Difference in surface elevation between adjacent tiles [m]
      &  length_scale_ct,        &  !< char. length scale of surface heterogeneity [m]
      &  length_iface_ct            !< char. interface length between adjacent tiles [m m-2]

    ! Pointers to variables in memory for sibling tile
    dsl4jsb_Real2D_onChunk ::    &
      &  wtr_lex_net_srf_sibl,   &  !< Surface water exchange with connected tiles on sibling tile [kg m-2 s-1]
      &  wtr_runoff_in_sibl,     &  !< Surface runoff received from adjacent tiles on sibling tile [kg m-2 s-1]
      &  weq_snowblow_in_sibl,   &  !< Snow blown received from adjacent tiles on sibling tile [kg m-2 s-1]
      &  wtr_infil_out_sibl,     &  !< Infiltration from ponds to soil in adjacent tiles on sibling tile [kg m-2 s-1]
      &  weq_snowblow_out_sibl,  &  !< Snow blown off to adjacent tiles on sibling tile [kg m-2 s-1]
      &  weq_pond_sibl,          &  !< Surface pond water and ice content on sibling tile [weq m]
      &  wtr_pond_sibl,          &  !< Surface pond water content on sibling tile [m]
      &  ice_pond_sibl,          &  !< Surface pond ice content on sibling tile [m]
      &  fract_pond_max_sibl,    &  !< Maximum pond fraction on sibling tile [/]
      &  weq_pond_max_sibl,      &  !< Maximum pond water and ice content on sibling tile [weq m]
      &  weq_snow_sibl,          &  !< Snow water equivalent on surface on sibling tile [m]
      &  weq_snow_soil_sibl,     &  !< Snow water equivalent on soil on sibling tile [m]
      &  new_snow_soil_sibl,     &  !< Snow water equivalent on soil after redistribution on sibling tile [m]
      &  snow_soil_dens_sibl,    &  !< Snow density on soil on sibling tile [kg m-3]
      &  water_to_soil_lex_sibl     !< Infiltratable water from connected tile on sibling tile [m]

    dsl4jsb_Real3D_onChunk ::         &
      &  soil_depth_sl_sibl,          &  !< Soil depth until bedrock within each layer on sibling tile [m]
      &  wtr_lex_net_sl_sibl,         &  !< Subsurface water exchange with connected tiles on sibling tile [kg m-2 s-1]
      &  wtr_drainage_in_sl_sibl,     &  !< Drainage received from adjacent tiles on sibling tile [kg m-2 s-1]
      &  wtr_diffusion_in_sl_sibl,    &  !< Water diffusion received from adjacent tiles on sibling tile [kg m-2 s-1]
      &  wtr_diffusion_out_sl_sibl,   &  !< Water diffusion to adjacent tiles on sibling tile [kg m-2 s-1]
      &  wtr_soil_sl_sibl,            &  !< Soil water content on sibling tile [m]
      &  ice_soil_sl_sibl,            &  !< Soil ice content on sibling tile [m]
      &  wtr_soil_pot_scool_sl_sibl,  &  !< Potentially supercooled soil water content on sibling tile [m]
      &  wtr_soil_sat_sl_sibl,        &  !< Soil moisture at saturation on sibling tile [m]
      &  wtr_soil_fc_sl_sibl,         &  !< Soil moisture at field capacity on sibling tile [m]
      &  wtr_soil_res_sl_sibl,        &  !< Soil moisture at residual water content on sibling tile [m]
      &  hyd_cond_sat_sl_sibl,        &  !< Hydraulic conductivity of soil layers on sibling tile [m s-1]
      &  vol_porosity_sl_sibl,        &  !< Volumetric soil porosity on sibling tile [m/m]
      &  vol_wres_sl_sibl,            &  !< Volumetric residual water content on sibling tile [m/m]
      &  matric_pot_sl_sibl,          &  !< Matric potential of soil layers at saturation on sibling tile [m]
      &  bclapp_sl_sibl,              &  !< Clapp and Hornberger exponent b of soil layers on sibling tile []
      &  pore_size_index_sl_sibl,     &  !< Pore size index of soil layers on sibling tile []
      &  t_soil_sl_sibl,              &  !< Soil layer temperature on sibling tile [K]
      &  heat_cond_sl_sibl,           &  !< Thermal conductivity of the soil layer on sibling tile [J m-1 s-1 K-1]
      &  vol_heat_cap_sl_sibl,        &  !< Volumetric heat capacity of the soil layer on sibling tile [J m-3 K-1]
      &  nrg_diffusion_in_sl_sibl,    &  !< Energy diffusion received from adjacent tiles on sibling tile [W m-2]
      &  nrg_convection_in_sl_sibl,   &  !< Energy advection received from adjacent tiles on sibling tile [W m-2]
      &  nrg_diffusion_out_sl_sibl,   &  !< Energy diffusion to adjacent tiles on sibling tile [W m-2]
      &  nrg_convection_out_sl_sibl,  &  !< Energy advection to adjacent tiles on sibling tile [W m-2]
      &  nrg_lex_net_sl_sibl             !< Subsurface energy exchange with connected tiles on sibling tile [W m-2]

    CLASS(t_jsb_tile_abstract), POINTER :: sibl  !< Sibling to current tile

    TYPE(t_jsb_model), POINTER    :: model       !< Current instance of the model
    TYPE(t_jsb_vgrid), POINTER    :: soil_e      !< Vertical grid for soil heat transfer

    INTEGER  :: &
      &  iblk,                               &  !< Current block index
      &  ics,                                &  !< First cell in block
      &  ice,                                &  !< Last cell in block
      &  nc,                                 &  !< Number of cells in block
      &  nsoil,                              &  !< Number of vertical soil layers (grid dimension)
      &  nsibl,                              &  !< Number of siblings on tile level
      &  nrgrd,                              &  !< Number of layers on helper grid determined by layer overlaps
      &  ic, is, it,                         &  !< counters for loops over cells (ic), soil layers (is), tiles (it)
      &  lex_idt,                            &  !< lex id of the tile
      &  lex_ids,                            &  !< lex id of the sibling tile
      &  lex_idc,                            &  !< lex id of the current tile
      &  lex_id_first,                       &  !< lex id of the first tile on the present tile level
      &  lex_id_last                            !< lex id of the last tile on the present tile level

    REAL(wp) :: dtime                           !< Timestep length [s]

    REAL(wp), PARAMETER :: fos = 0.99_wp        !< Factor of safety to avoid negative fluxes due to roundoff errors [-]

    REAL(wp), DIMENSION(options%nc) ::  &
      &  fract,                              &  !< Cover fractions of tile [-]
      &  fract_sibl,                         &  !< Cover fractions of sibling tile [-]
      &  snow_soil_tle_tm1,                  &  !< Snow on soil on tile at end of previous timestep [weq m]
      &  snow_new_tle_tm1,                   &  !< New snow on soil on tile at end of previous timestep [weq m]
      &  snow_dens_tle_tm1,                  &  !< Snow density on tile at end of previous timestep [kg m-3]
      &  wtr_incomming_srf,                  &  !< Lateral flux (surface) from sibling to tile [kg m-2 s-1]
      &  wtr_incomming_srf_sibl                 !< Lateral flux (surface) from tile to sibling [kg m-2 s-1]

    REAL(wp), DIMENSION(options%nc, options%nsoil_w) ::  &
      &  wtr_latflow_sl_used,    &  !< Outflow from latflow res. with water in bedrock-layers on sibling
                                    !< moved to the layer above [m timestep-1]
      &  hyd_cond_sl,            &  !< Saturation hydraulic conductivity for tile [m s-1]
      &  hyd_cond_sl_sibl,       &  !< Saturation hydraulic conductivity for sibling [m s-1]
      &  diffus_sl,              &  !< diffusivity for tile [m2 s-1]
      &  diffus_sl_sibl,         &  !< diffusivity for sibling [m2 s-1]
      &  mpot_act,               &  !< matric potential for tile  [m]
      &  mpot_act_sibl,          &  !< matric potential for sibling [m]
      &  ice_impedance_sl,       &  !< ice impedence factor for tile [/]
      &  ice_impedance_sl_sibl,  &  !< ice impedence factor for sibling [/]
      &  weq_wsat_sl,            &  !< Soil moisture at saturation for tile [m]
      &  weq_wsat_sl_sibl,       &  !< Soil moisture at saturation for sibling [m]
      &  weq_wres_sl,            &  !< moisture residual for tile [m]
      &  weq_wres_sl_sibl,       &  !< moisture residual for sibling [m]
      &  wtr_incomming_sl,       &  !< Lateral flux (soil) from sibling to tile [kg m-2 s-1]
      &  wtr_incomming_sl_sibl      !< Lateral flux (soil) from tile to sibling [kg m-2 s-1]

    REAL(wp) :: hlp1                !< helper variable

    LOGICAL  :: l_dynsnow           !< Switch from SSE to use the dynamical snow heat properties computation

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_lex'


    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsoil = options%nsoil_w
    dtime = options%dtime

    ! Get local variables from tile information
    model  => Get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')
    nsibl  = tile%parent_tile%Get_no_of_children()

    ! If process is not active on this tile, do nothing and return
    IF (.NOT. tile%Is_process_calculated(LEX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! First get first and last lex_id on current tile level
    lex_id_first = 0
    lex_id_last  = 0

    ! Loop over sibling tiles to get lex ids
    DO it = 1, nsibl
      IF (it == 1) THEN
        sibl => tile%parent_tile%Get_first_child_tile()
        lex_idc = sibl%lcts(1)%lex_id
      ELSE
        sibl => sibl%Get_next_sibling_tile()
        lex_idc = sibl%lcts(1)%lex_id
      END IF

      IF (lex_id_first == 0 .AND. lex_idc > 0 .AND. lex_idc == lex_idc) THEN
        lex_id_first = lex_idc
      END IF

      IF (lex_idc > 0 .AND. lex_idc == lex_idc) THEN
        lex_id_last = lex_idc
      END IF
    END DO

    ! Set pointers to the types containing configurable parameters of the respective processes
    dsl4jsb_Get_config(LEX_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    dsl4jsb_Get_memory(LEX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    dsl4jsb_Get_var3D_onChunk(LEX_,   microcon_ctc_ct)        ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   microcon_hfp_blg_ct)    ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   fract_latflow_srf_ct)   ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   fract_latflow_blg_ct)   ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   h_offset_ct)            ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   length_scale_ct)        ! in
    dsl4jsb_Get_var3D_onChunk(LEX_,   length_iface_ct)        ! in
    dsl4jsb_Get_var2D_onChunk(LEX_,   weq_snowblow_in)        ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   wtr_runoff_out)         ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,   weq_snowblow_out)       ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_infil_in_sl)        ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_diffusion_in_sl)    ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_drainage_out_sl)    ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   wtr_diffusion_out_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_diffusion_in_sl)    ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_convection_in_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_diffusion_out_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,   nrg_convection_out_sl)  ! inout

    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_latflow_srf)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_latflow_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sat_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_fc_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_res_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_wres_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, matric_pot_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_lex_net_srf)        ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)                 ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)               ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)          ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, new_snow_soil)          ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)         ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_lex_net_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)            ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_, drainage_sl)            ! inout

    dsl4jsb_Get_var3D_onChunk(SSE_,   heat_cond_sl)           ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   vol_heat_cap_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   t_soil_sl)              ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   nrg_lex_net_sl)         ! inout

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(fract(1:nc), fract_sibl(1:nc), fract_sibl(1:nc), snow_soil_tle_tm1(1:nc), snow_new_tle_tm1(1:nc)) &
    !$ACC   CREATE(snow_dens_tle_tm1(1:nc), wtr_incomming_srf(1:nc), wtr_incomming_srf_sibl(1:nc))                   &
    !$ACC   CREATE(wtr_latflow_sl_used(1:nc,1:nsoil), hyd_cond_sl(1:nc,1:nsoil), hyd_cond_sl_sibl(1:nc,1:nsoil))     &
    !$ACC   CREATE(diffus_sl(1:nc,1:nsoil), diffus_sl_sibl(1:nc,1:nsoil), mpot_act(1:nc,1:nsoil))                    &
    !$ACC   CREATE(mpot_act_sibl(1:nc,1:nsoil), ice_impedance_sl(1:nc,1:nsoil), ice_impedance_sl_sibl(1:nc,1:nsoil)) &
    !$ACC   CREATE(weq_wsat_sl(1:nc,1:nsoil), weq_wsat_sl_sibl(1:nc,1:nsoil), weq_wres_sl(1:nc,1:nsoil))             &
    !$ACC   CREATE(weq_wres_sl_sibl(1:nc,1:nsoil), wtr_incomming_sl(1:nc,1:nsoil), wtr_incomming_sl_sibl(1:nc,1:nsoil))

    l_dynsnow = dsl4jsb_Config(SSE_)%l_dynsnow
    ! Get lex_id of tile and its cover fraction
    lex_idt = tile%lcts(1)%lex_id
    CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))

    ! Initialze those variables, where the values at the end of ast time step,
    ! i.e. unmodified by interactions with siblings are needed.
    DO ic = 1, nc
      snow_soil_tle_tm1(ic) = weq_snow_soil(ic)
      snow_new_tle_tm1(ic)  = new_snow_soil(ic)
      snow_dens_tle_tm1(ic) = snow_soil_dens(ic)
    END DO

    ! Maximum number of intersecting soil layers for lateral exchange
    nrgrd   = (nsoil + 2) * 2

    !> Loop over all tiles on the current tile level to check if it should exchange water and heat
    !> with any connected tiles, compute the fluxes and update the respective variables.
    !>
    ! This includes tiles that may not participate in the lateral exchange, such as glaciers.
    DO it = 1, nsibl

      ! Get lex_id of sibling
      IF (it == 1) THEN
        sibl => tile%parent_tile%Get_first_child_tile()
        lex_ids = sibl%lcts(1)%lex_id
      ELSE
        sibl => sibl%Get_next_sibling_tile()
        lex_ids = sibl%lcts(1)%lex_id
      END IF

      !> Only calculate lateral exchange on tiles for which connections have been defined
      !> and skip those tiles that are not participating in lateral exchange
      !>
      ! i.e. lex_id > 0 and lex_id /= nan
      IF (lex_idt == 0 .OR. lex_ids == 0 .OR. &
        & lex_idt /= lex_idt .OR. lex_ids /= lex_ids) CYCLE

      ! Get process memories for sibling
      dsl4jsb_Get_memory_tile(LEX_, sibl)
      dsl4jsb_Get_memory_tile(HYDRO_, sibl)
      dsl4jsb_Get_memory_tile(SSE_, sibl)

      dsl4jsb_Get_var2D_onChunk_tile_name(LEX_,   wtr_runoff_in, sibl)          ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(LEX_,   weq_snowblow_in, sibl)        ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(LEX_,   wtr_infil_out, sibl)          ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(LEX_,   weq_snowblow_out, sibl)       ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   wtr_drainage_in_sl, sibl)     ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   wtr_diffusion_in_sl, sibl)    ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   wtr_diffusion_out_sl, sibl)   ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   nrg_diffusion_in_sl, sibl)    ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   nrg_convection_in_sl, sibl)   ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   nrg_diffusion_out_sl, sibl)   ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(LEX_,   nrg_convection_out_sl, sibl)  ! inout

      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, ice_pond, sibl)               ! in
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, fract_pond_max, sibl)         ! in
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, weq_pond_max, sibl)           ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, soil_depth_sl, sibl)          ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_soil_pot_scool_sl, sibl)  ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_soil_sat_sl, sibl)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_soil_fc_sl, sibl)         ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_soil_res_sl, sibl)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_porosity_sl, sibl)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, vol_wres_sl, sibl)            ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, matric_pot_sl, sibl)          ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, bclapp_sl, sibl)              ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, pore_size_index_sl, sibl)     ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, hyd_cond_sat_sl, sibl)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, ice_soil_sl, sibl)            ! in
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, wtr_lex_net_srf, sibl)        ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, weq_pond, sibl)               ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, wtr_pond, sibl)               ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, water_to_soil_lex, sibl)      ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, weq_snow, sibl)               ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, weq_snow_soil, sibl)          ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, new_snow_soil, sibl)          ! inout
      dsl4jsb_Get_var2D_onChunk_tile_name(HYDRO_, snow_soil_dens, sibl)         ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_lex_net_sl, sibl)         ! inout
      dsl4jsb_Get_var3D_onChunk_tile_name(HYDRO_, wtr_soil_sl, sibl)            ! inout

      dsl4jsb_Get_var3D_onChunk_tile_name(SSE_,   heat_cond_sl, sibl)           ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(SSE_,   vol_heat_cap_sl, sibl)        ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(SSE_,   t_soil_sl, sibl)              ! in
      dsl4jsb_Get_var3D_onChunk_tile_name(SSE_,   nrg_lex_net_sl, sibl)         ! inout

      ! Get cover fraction of sibling tile
      CALL sibl%Get_fraction(ics, ice, iblk, fract=fract_sibl(:))

      ! Check if lateral flux fields have been reset
      IF (lex_idt == lex_id_first) THEN
        IF (ANY(wtr_lex_net_srf_sibl(:) /= 0._wp)) THEN
          WRITE (message_text,*) 'Reset of lateral flux fields did not work', &
            &  MAXVAL(wtr_lex_net_srf_sibl(:))
          CALL finish(TRIM(routine), message_text)
        END IF
      END IF

      ! Start calculating exchange processes between different tiles
      IF (lex_idt /= lex_ids) THEN

        wtr_incomming_srf(:)       = 0._wp
        wtr_incomming_srf_sibl(:)  = 0._wp
        wtr_incomming_sl(:,:)      = 0._wp
        wtr_incomming_sl_sibl(:,:) = 0._wp

        ! Incomming later subsurface runoff is adjusted according to differences in soil depth (until bedrock),
        ! and for micro-scale connections, as well as for differences in surface elevation
        wtr_latflow_sl_used(:,:) = wtr_latflow_sl(:,:)

        CALL shift_fluxes_sl( nc, nsoil, nrgrd,    &
          &  microcon_hfp_blg_ct(:,lex_ids),       &
          &  h_offset_ct(:,lex_ids),               &
          &  soil_e%ubounds(:),                    &
          &  soil_depth_sl(:,:),                   &
          &  soil_depth_sl_sibl(:,:),              &
          &  wtr_latflow_sl_used(:,:)              &
          )

        !> Compute the lateral exchange of surface runoff. The outflow from each tile is already calculated in HYDRO_.
        !> Here only the transfer of surface runoff onto siblings is done
        !>
        IF (dsl4jsb_Config(LEX_)%l_runoff_trns) THEN
          DO ic = 1, nc

            hlp1 = wtr_latflow_srf(ic) * fos
            wtr_incomming_srf_sibl(ic) = wtr_incomming_srf_sibl(ic) - wtr_runoff_in_sibl(ic)

            CALL calc_runoff_transport( dtime,              &
              & fract_latflow_srf_ct(ic,lex_ids),           &
              & fract(ic),                                  &
              & fract_sibl(ic),                             &
              & hlp1,                                       &
              & wtr_runoff_in_sibl(ic),                     &
              & wtr_runoff_out(ic),                         &
              & wtr_lex_net_srf(ic),                        &
              & wtr_lex_net_srf_sibl(ic),                   &
              & wtr_pond_sibl(ic)                           &
              )

            wtr_incomming_srf_sibl(ic) = wtr_incomming_srf_sibl(ic) + wtr_runoff_in_sibl(ic)
          END DO
        END IF

        ! Lateral exchange of subsurface runoff and soil moisture between tiles only allowed
        ! if both tiles actually contain soil
        IF (tile%contains_soil .AND. sibl%contains_soil) THEN
          !> Compute the lateral exchange of subsurface runoff. The outflow from each tile is already
          !> calculated in HYDRO_. Here only the transfer of surface runoff onto siblings is done
          !>
          ! As we already updated the vertical distribution of fluxes in `shift_fluxes_sl`, we can loop directly
          ! over the soil layers and move the fluxes without any further vertical interpolation.
          IF (dsl4jsb_Config(LEX_)%l_drainage_trns) THEN
            DO is = 1, nsoil
              DO ic = 1, nc

                hlp1 = wtr_latflow_sl_used(ic,is) * fos
                wtr_incomming_sl_sibl(ic,is) = wtr_incomming_sl_sibl(ic,is) &
                  &                          - wtr_drainage_in_sl_sibl(ic,is)
                CALL calc_drainage_transport( dtime,        &
                  & fract_latflow_blg_ct(ic,lex_ids),       &
                  & fract(ic),                              &
                  & fract_sibl(ic),                         &
                  & hlp1,                                   &
                  & wtr_drainage_in_sl_sibl(ic,is),         &
                  & wtr_drainage_out_sl(ic,is),             &
                  & wtr_lex_net_sl(ic,is),                  &
                  & wtr_lex_net_sl_sibl(ic,is),             &
                  & wtr_soil_sl_sibl(ic,is)                 &
                  )
                wtr_incomming_sl_sibl(ic,is) = wtr_incomming_sl_sibl(ic,is) &
                  &                          + wtr_drainage_in_sl_sibl(ic,is)

              END DO
            END DO

            ! Sanity checke: remove excess water from soil and surface
            CALL remove_excess_water(nc,                   &
              & nsoil,                                     &
              & weq_pond_max_sibl(:),                      &
              & ice_pond_sibl(:),                          &
              & wtr_pond_sibl(:),                          &
              & water_to_soil_lex_sibl(:),                 &
              & wtr_soil_sat_sl=wtr_soil_sat_sl_sibl(:,:), &
              & wtr_soil_sl=wtr_soil_sl_sibl(:,:)          &
              )

          END IF

          ! Calculate soil moisture content at saturation, field capacity and moisture residual for the tile
          ! without substracting the soil ice content
          DO is = 1, nsoil
            DO ic = 1, nc
              weq_wsat_sl(ic,is) = vol_porosity_sl (ic,is) * soil_depth_sl(ic,is)
              weq_wres_sl(ic,is) = vol_wres_sl(ic,is)      * soil_depth_sl(ic,is)
            END DO
          END DO

          !> Compute the hydrological water flow properties on the tile and sibling tile,
          !> needed for lateral diffusion of soil moisture between connected tiles.
          !> This is done after the soil water content is updated with the inflow of lateral subsurface drainage.
          !>
          ! Compute hydrological water flow properties on tile (with updated soil water content)
          CALL get_soilhyd_properties(              &
            & dsl4jsb_Config(HYDRO_)%soilhydmodel,  &
            & dsl4jsb_Config(HYDRO_)%interpol_mean, &
            & nc, nsoil,                            &
            & soil_depth_sl(:,:),                   &
            & wtr_soil_sl(:,:),                     &
            & ice_soil_sl(:,:),                     &
            & wtr_soil_pot_scool_sl(:,:),           &
            & weq_wsat_sl(:,:),                     &
            & weq_wres_sl(:,:),                     &
            & hyd_cond_sat_sl(:,:),                 &
            & matric_pot_sl(:,:),                   &
            & bclapp_sl(:,:),                       &
            & pore_size_index_sl(:,:),              &
            & ice_impedance=ice_impedance_sl(:,:),  &
            & K=hyd_cond_sl(:,:),                   &
            & D=diffus_sl(:,:),                     &
            & mpot_act=mpot_act(:,:)                &
            )

          ! Calculate soil moisture content at saturation, field capacity and moisture residual for the sibling tile
          ! without substracting the soil ice content
          DO is = 1, nsoil
            DO  ic = 1, nc
              weq_wsat_sl_sibl     (ic,is) = vol_porosity_sl_sibl (ic,is) &
                &                          * soil_depth_sl_sibl(ic,is)
              weq_wres_sl_sibl     (ic,is) = vol_wres_sl_sibl(ic,is)      &
                &                          * soil_depth_sl_sibl(ic,is)
            END DO
          END DO

          ! Compute hydrological water flow properties on sibling tile (with updated soil water content)
          CALL get_soilhyd_properties(                  &
            & dsl4jsb_Config(HYDRO_)%soilhydmodel,      &
            & dsl4jsb_Config(HYDRO_)%interpol_mean,     &
            & nc, nsoil,                                &
            & soil_depth_sl_sibl(:,:),                  &
            & wtr_soil_sl_sibl(:,:),                    &
            & ice_soil_sl_sibl(:,:),                    &
            & wtr_soil_pot_scool_sl_sibl(:,:),          &
            & weq_wsat_sl_sibl(:,:),                    &
            & weq_wres_sl_sibl(:,:),                    &
            & hyd_cond_sat_sl_sibl(:,:),                &
            & matric_pot_sl_sibl(:,:),                  &
            & bclapp_sl_sibl(:,:),                      &
            & pore_size_index_sl_sibl(:,:),             &
            & ice_impedance=ice_impedance_sl_sibl(:,:), &
            & K=hyd_cond_sl_sibl(:,:),                  &
            & D=diffus_sl_sibl(:,:),                    &
            & mpot_act=mpot_act_sibl(:,:)               &
            )

          !> Compute the lateral ingress of water and heat from ponds into the soil of a connected tile.
          !>
          ! Since an update of the actual soil temperatures outside of the SSE routines is not intended,
          ! only the lateral heat flux is calculated and passed to the SSE routine as an energy source/sink.
          ! However, the soil hydrological state is updated here.
          CALL calc_lat_fluxes_pond(               &
            &  nc, nsoil, nrgrd, dtime,            &
            &  dsl4jsb_Config(LEX_)%l_infil_trns,  &
            &  dsl4jsb_Config(LEX_)%l_dffsn_nrg,   &
            &  microcon_ctc_ct(:,lex_ids),         &
            &  h_offset_ct(:,lex_ids),             &
            &  length_scale_ct(:,lex_ids),         &
            &  length_iface_ct(:,lex_ids),         &
            &  fract(:),                           &
            &  fract_sibl(:),                      &
            &  fract_pond_max_sibl(:),             &
            !  Hydrological properties
            &  soil_depth_sl(:,:),                 &
            &  wtr_soil_fc_sl(:,:),                &
            &  hyd_cond_sl(:,:),                   &
            &  mpot_act(:,:),                      &
            !  Thermal properties
            &  soil_e%ubounds(:),                  &
            &  t_soil_sl(:,:),                     &
            &  t_soil_sl_sibl(:,:),                &
            &  heat_cond_sl(:,:),                  &
            &  vol_heat_cap_sl(:,:),               &
            &  vol_heat_cap_sl_sibl(:,:),          &
            &  ice_pond_sibl(:),                   &
            !  Water states and fluxes between pond and soil
            &  wtr_infil_out_sibl(:),              &
            &  wtr_infil_in_sl(:,:),               &
            &  wtr_incomming_sl(:,:),              &
            &  wtr_lex_net_srf_sibl(:),            &
            &  wtr_lex_net_sl(:,:),                &
            &  wtr_pond_sibl(:),                   &
            &  wtr_soil_sl(:,:),                   &
            !  Heat fluxes between pond and soil
            &  nrg_diffusion_in_sl_sibl(:,1),      &
            &  nrg_diffusion_out_sl_sibl(:,1),     &
            &  nrg_lex_net_sl_sibl(:,1),           &
            &  nrg_diffusion_in_sl(:,:),           &
            &  nrg_diffusion_out_sl(:,:),          &
            &  nrg_lex_net_sl(:,:)                 &
            &  )

          !> Compute the lateral diffusion of water and heat between the soil layer of the current tile and
          !> a connected sibling tile.
          !>
          CALL calc_lat_fluxes_soil(                    &
            &  nc, nsoil, nrgrd, dtime,                 &
            &  dsl4jsb_Config(LEX_)%l_dffsn_wtr,        &
            &  dsl4jsb_Config(LEX_)%l_dffsn_nrg,        &
            &  microcon_ctc_ct(:,lex_ids),              &
            &  h_offset_ct(:,lex_ids),                  &
            &  length_scale_ct(:,lex_ids),              &
            &  length_iface_ct(:,lex_ids),              &
            &  fract(:),                                &
            &  fract_sibl(:),                           &
            ! Soil hydrological properties
            &  soil_depth_sl(:,:),                      &
            &  soil_depth_sl_sibl(:,:),                 &
            &  wtr_soil_fc_sl(:,:),                     &
            &  wtr_soil_fc_sl_sibl(:,:),                &
            &  wtr_soil_sat_sl(:,:),                    &
            &  wtr_soil_sat_sl_sibl(:,:),               &
            &  wtr_soil_res_sl(:,:),                    &
            &  wtr_soil_res_sl_sibl(:,:),               &
            &  hyd_cond_sat_sl(:,:),                    &
            &  hyd_cond_sat_sl_sibl(:,:),               &
            &  hyd_cond_sl(:,:),                        &
            &  hyd_cond_sl_sibl(:,:),                   &
            &  mpot_act(:,:),                           &
            &  mpot_act_sibl(:,:),                      &
            ! Soil thermal properties
            &  soil_e%ubounds(:),                       &
            &  t_soil_sl(:,:),                          &
            &  t_soil_sl_sibl(:,:),                     &
            &  heat_cond_sl(:,:),                       &
            &  heat_cond_sl_sibl(:,:),                  &
            &  vol_heat_cap_sl(:,:),                    &
            &  vol_heat_cap_sl_sibl(:,:),               &
            ! Soil water fluxes
            &  wtr_diffusion_in_sl(:,:),                &
            &  wtr_diffusion_in_sl_sibl(:,:),           &
            &  wtr_diffusion_out_sl(:,:),               &
            &  wtr_diffusion_out_sl_sibl(:,:),          &
            &  wtr_incomming_sl(:,:),                   &
            &  wtr_incomming_sl_sibl(:,:),              &
            &  wtr_lex_net_sl(:,:),                     &
            &  wtr_lex_net_sl_sibl(:,:),                &
            &  wtr_soil_sl(:,:),                        &
            &  wtr_soil_sl_sibl(:,:),                   &
            ! Soil heat fluxes
            &  nrg_diffusion_in_sl(:,:),                &
            &  nrg_diffusion_in_sl_sibl(:,:),           &
            &  nrg_diffusion_out_sl(:,:),               &
            &  nrg_diffusion_out_sl_sibl(:,:),          &
            &  nrg_lex_net_sl(:,:),                     &
            &  nrg_lex_net_sl_sibl(:,:)                 &
            &  )

          IF (dsl4jsb_Config(LEX_)%l_cnvctn_nrg) THEN
            !> Compute the advection of heat due to the lateral transport of water between the
            !> current tile and a connected sibling tile. This is done for surface and subsurface lateral fluxes.
            !>
            ! Since an update of the actual soil temperatures outside of the SSE routines is not intended,
            ! only the lateral heat flux is calculated and passed to the SSE routine as an energy source/sink.
            CALL calc_convection_heat(                   &
              &  nc, nsoil, nrgrd, dtime,                &
              &  soil_e%dz(:),                           &
              &  h_offset_ct(:,lex_ids),                 &
              &  fract(:),                               &
              &  fract_sibl(:),                          &
              &  wtr_incomming_srf(:),                   &
              &  wtr_incomming_srf_sibl(:),              &
              &  t_soil_sl(:,:),                         &
              &  t_soil_sl_sibl(:,:),                    &
              &  wtr_pond(:),                            &
              &  wtr_pond_sibl(:),                       &
              &  wtr_incomming_sl(:,:),                  &
              &  wtr_incomming_sl_sibl(:,:),             &
              &  soil_depth_sl(:,:),                     &
              &  soil_depth_sl_sibl(:,:),                &
              &  vol_heat_cap_sl(:,:),                   &
              &  vol_heat_cap_sl_sibl(:,:),              &
              &  nrg_convection_in_sl(:,:),              &
              &  nrg_convection_in_sl_sibl(:,:),         &
              &  nrg_convection_out_sl(:,:),             &
              &  nrg_convection_out_sl_sibl(:,:),        &
              &  nrg_lex_net_sl(:,:),                    &
              &  nrg_lex_net_sl_sibl(:,:)                &
              &  )
          END IF

          ! Sanity check: remove excess water from soil and surface
          CALL remove_excess_water(nc,                   &
            & nsoil,                                     &
            & weq_pond_max_sibl(:),                      &
            & ice_pond_sibl(:),                          &
            & wtr_pond_sibl(:),                          &
            & water_to_soil_lex_sibl(:),                 &
            & wtr_soil_sat_sl=wtr_soil_sat_sl_sibl(:,:), &
            & wtr_soil_sl=wtr_soil_sl_sibl(:,:)          &
            )
        ELSE
          ! Sanity check: remove excess water from surface
          CALL remove_excess_water(nc, 0,                &
            & weq_pond_max_sibl(:),                      &
            & ice_pond_sibl(:),                          &
            & wtr_pond_sibl(:),                          &
            & water_to_soil_lex_sibl(:)                  &
            )
        END IF

        ! Update total pond water and ice content on the sibling tile.
        ! This is required because the water balance does not use wtr_pond and ice_pond directly.
        DO ic = 1,nc
          weq_pond_sibl(ic) = wtr_pond_sibl(ic) + ice_pond_sibl(ic)
        END DO

        IF (dsl4jsb_Config(LEX_)%l_snow_trns) THEN
          !> Estimate lateral snow redistribution between connected tiles based solely on the assumption
          !> that snow is redistributed from ridges to fill the depressions on the landscape. Factors like
          !> wind direction/speed, vegetation effects or surface geometry are not considered here.
          !>
          DO ic = 1,nc
            IF (microcon_ctc_ct(ic,lex_ids) > 0._wp) THEN
              hlp1 = MIN(1._wp, microcon_ctc_ct(ic,lex_ids) / SUM(microcon_ctc_ct(ic,:)))
              ! Snowblow between connected tiles
              CALL calc_snow_redistribution(                  &
                & l_dynsnow,                                  &
                & dtime,                                      &
                & hlp1,                                       &
                & h_offset_ct(ic,lex_ids),                    &
                & fract(ic),                                  &
                & fract_sibl(ic),                             &
                & weq_pond(ic),                               &
                & weq_pond_sibl(ic),                          &
                & snow_soil_tle_tm1(ic),                      &
                & snow_new_tle_tm1(ic),                       &
                & snow_dens_tle_tm1(ic),                      &
                & weq_snowblow_in(ic),                        &
                & weq_snowblow_in_sibl(ic),                   &
                & weq_snowblow_out(ic),                       &
                & weq_snowblow_out_sibl(ic),                  &
                & wtr_lex_net_srf(ic),                        &
                & wtr_lex_net_srf_sibl(ic),                   &
                & weq_snow_soil(ic),                          &
                & weq_snow_soil_sibl(ic),                     &
                & new_snow_soil(ic),                          &
                & new_snow_soil_sibl(ic),                     &
                & weq_snow(ic),                               &
                & weq_snow_sibl(ic),                          &
                & snow_soil_dens(ic),                         &
                & snow_soil_dens_sibl(ic)                     &
                & )
            END IF
          END DO
        END IF
      END IF
    END DO

    ! In case of tiles that are not connected to downstream tiles, the outflow is given to runoff and drainage
    DO ic = 1,nc
      runoff(ic) = runoff(ic) &
        &        + wtr_latflow_srf(ic) / dtime * rhoh2o - wtr_runoff_out(ic)
    END DO

    ! Update the total and layered drainage to account for the water given to the connected sibling tile.
    IF (tile%contains_soil) THEN
      DO ic = 1,nc
        drainage(ic) = drainage(ic) &
          & + SUM(wtr_latflow_sl(ic,:)) / dtime * rhoh2o - SUM(wtr_drainage_out_sl(ic,:))

        drainage_sl(ic,:) = drainage_sl(ic,:) &
          & + wtr_latflow_sl(ic,:) / dtime * rhoh2o - wtr_drainage_out_sl(ic,:)
      END DO
    END IF

    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_lex

  ! ======================================================================================================= !
  !>
  !>#### Aggregation of variables from lateral exchange process
  !>
  !> Lateral exchange as well as some hydrology and soil energy variables are aggregated in this routine.
  !>
  SUBROUTINE aggregate_lex(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    dsl4jsb_Def_memory(LEX_)    !> Memory definition for lateral exchange process
    dsl4jsb_Def_memory(HYDRO_)  !> Memory definition for hydrology process
    dsl4jsb_Def_memory(SSE_)    !> Memory definition for soil energy process

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions

    INTEGER :: &
      & iblk,  &                !< Number of current block
      & ics,   &                !< Index of first cell of current block
      & ice                     !< Index of last cell of current block

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_lex'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LEX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(SSE_,   nrg_lex_net_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_srf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_to_soil_lex,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,                weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, new_snow_soil,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_runoff_in,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   weq_snowblow_in,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_infil_in_sl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_drainage_in_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_runoff_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_infil_out,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   weq_snowblow_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_drainage_out_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   wtr_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_convection_in_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   nrg_convection_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   microcon_ctc_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   microcon_hfp_srf_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   microcon_hfp_blg_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   macrocon_hfp_srf_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   macrocon_hfp_blg_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   slope_eff,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fact_lag_srf,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fact_lag_blg,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   h_offset_ct,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   length_scale_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   length_iface_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fract_latflow_srf_ct,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_,   fract_latflow_blg_ct,  weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lex


  ! ================================================================================================================================

#endif
END MODULE mo_lex_interface
