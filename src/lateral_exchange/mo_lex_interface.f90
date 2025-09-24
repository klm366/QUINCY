!> Contains the interfaces to the lex process
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
!>#### Contains interfaces to lex processes and reset of lateral flux fields after these have been used in HYDRO_ and SSE_ modules
!>
MODULE mo_lex_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on, acc_stream
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

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
    PROCEDURE, NOPASS :: Integrate => update_lex     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_lex  !< Aggregates computed task variables
  END TYPE tsk_lex

  !> Constructor interface for lex task
  INTERFACE tsk_lex
    PROCEDURE Create_task_lex                        !< Constructor function for task
  END INTERFACE tsk_lex

  !X Repeat here for each additional task of the process the two blocks above ...

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for reset_lex_fluxes task
  !!
  FUNCTION Create_task_reset_lex_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_lex_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_lex_fluxes', process_id=LEX_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_lex_fluxes


  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for lex task
  !!
  FUNCTION Create_task_lex(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_lex::return_ptr)
    CALL return_ptr%Construct(name='lex', process_id=LEX_, owner_model_id=model_id)

  END FUNCTION Create_task_lex

  !X Repeat here for each additional task of the process the function above ...

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for lex process
  !!
  SUBROUTINE Register_lex_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,               INTENT(in)    :: model_id

    CALL process%Register_task(tsk_reset_lex_fluxes(model_id))
    CALL process%Register_task(tsk_lex(model_id))

  END SUBROUTINE Register_lex_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "reset_lex_fluxes"
  !! Task "reset_lex_fluxes" resets the lateral flux fields after these have been used in the HYDRO_ and SSE_ routines.
  !!
  SUBROUTINE update_reset_lex_fluxes(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice

    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_lex_fluxes'


    ! Declare process memories
    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_runoff_in
    dsl4jsb_Real2D_onChunk :: weq_snowblow_in
    dsl4jsb_Real2D_onChunk :: wtr_runoff_out
    dsl4jsb_Real2D_onChunk :: wtr_infil_out
    dsl4jsb_Real2D_onChunk :: weq_snowblow_out

    dsl4jsb_Real3D_onChunk :: wtr_drainage_in_sl
    dsl4jsb_Real3D_onChunk :: wtr_diffusion_in_sl
    dsl4jsb_Real3D_onChunk :: wtr_drainage_out_sl
    dsl4jsb_Real3D_onChunk :: wtr_diffusion_out_sl

    dsl4jsb_Real3D_onChunk :: wtr_infil_in_sl
    dsl4jsb_Real3D_onChunk :: nrg_diffusion_in_sl
    dsl4jsb_Real3D_onChunk :: nrg_convection_in_sl
    dsl4jsb_Real3D_onChunk :: nrg_diffusion_out_sl
    dsl4jsb_Real3D_onChunk :: nrg_convection_out_sl

    dsl4jsb_Real2D_onChunk :: water_to_soil_lex
    dsl4jsb_Real2D_onChunk :: wtr_lex_net_srf

    dsl4jsb_Real3D_onChunk :: wtr_lex_net_sl

    dsl4jsb_Real3D_onChunk :: nrg_lex_net_sl

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(LEX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(LEX_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)

    ! Get process variables

    dsl4jsb_Get_var2D_onChunk(LEX_,    wtr_runoff_in)        ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,    weq_snowblow_in)      ! inout

    dsl4jsb_Get_var2D_onChunk(LEX_,    wtr_runoff_out)       ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,    wtr_infil_out)        ! inout
    dsl4jsb_Get_var2D_onChunk(LEX_,    weq_snowblow_out)     ! inout

    dsl4jsb_Get_var3D_onChunk(LEX_,    wtr_drainage_in_sl)   ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,    wtr_diffusion_in_sl)  ! inout

    dsl4jsb_Get_var3D_onChunk(LEX_,    wtr_drainage_out_sl)  ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,    wtr_diffusion_out_sl) ! inout

    dsl4jsb_Get_var3D_onChunk(LEX_,    wtr_infil_in_sl)      ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,    nrg_diffusion_in_sl)  ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,    nrg_convection_in_sl) ! inout

    dsl4jsb_Get_var3D_onChunk(LEX_,    nrg_diffusion_out_sl) ! inout
    dsl4jsb_Get_var3D_onChunk(LEX_,    nrg_convection_out_sl)! inout

    dsl4jsb_Get_var2D_onChunk(HYDRO_,  water_to_soil_lex)    ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_,  wtr_lex_net_srf)      ! inout
    dsl4jsb_Get_var3D_onChunk(HYDRO_,  wtr_lex_net_sl)       ! inout

    dsl4jsb_Get_var3D_onChunk(SSE_,    nrg_lex_net_sl)       ! inout

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
    nrg_diffusion_in_sl  = 0._wp
    nrg_convection_in_sl = 0._wp
    nrg_diffusion_out_sl = 0._wp
    nrg_convection_out_sl= 0._wp

    water_to_soil_lex    = 0._wp
    wtr_lex_net_srf      = 0._wp
    wtr_lex_net_sl       = 0._wp


    nrg_lex_net_sl       = 0._wp


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_reset_lex_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "reset_lex_fluxes"
  !!
  SUBROUTINE aggregate_reset_lex_fluxes(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice

    !X if necessary: REAL(wp) :: dtime, steplen
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

    dsl4jsb_Aggregate_onChunk(LEX_, wtr_runoff_in,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, weq_snowblow_in,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_runoff_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_infil_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, weq_snowblow_out,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_infil_in_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_drainage_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_diffusion_in_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_drainage_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_diffusion_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_diffusion_in_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_convection_in_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_diffusion_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_convection_out_sl,weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, slope_eff,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fact_lag_srf,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fact_lag_blg,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, h_offset_ct,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, length_scale_ct,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, length_iface_ct,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fract_latflow_srf_ct, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fract_latflow_blg_ct, weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_srf,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_sl,     weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(SSE_, nrg_lex_net_sl,       weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_lex_fluxes

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "lex"
  !! Task "lex" determines the transport of surface runoff, lateral drainage, lateral ingress from ponds,
  !! the diffusive transport of water and energy and the redistribution of snow.
  !!
  SUBROUTINE update_lex(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_lex_process,            ONLY: calc_runoff_transport,           &
      &                                  calc_drainage_transport,         &
      &                                  calc_pond_soil_flux_water,       &
      &                                  calc_pond_soil_flux_heat,        &
      &                                  calc_diffusion_water,            &
      &                                  calc_diffusion_heat,             &
      &                                  calc_convection_heat,            &
      &                                  calc_snow_redistribution,        &
      &                                  shift_fluxes_sl,                 &
      &                                  remove_excess_water

    USE mo_hydro_process, ONLY: get_soilhyd_properties

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile           !< Sibling tile with LEX_ process running on them
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables

    CLASS(t_jsb_tile_abstract), POINTER :: sibl                 !< Sibling tile with LEX_ process running on them

    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_vgrid), POINTER    :: soil_e

    INTEGER  :: iblk,                                        &  !< options%iblk
      &         ics,                                         &  !< options%ics
      &         ice,                                         &  !< options%ice
      &         nc,                                          &  !< options%nc
      &         nsoil,                                       &  !< options%nsoil_w
      &         nsibl,                                       &  !< number of siblings on tile level
      &         ic, is, it,                                  &  !< counters for loops over cells (ic), soil layers (is), tiles (it)
      &         lex_idt,                                     &  !< lex id of tile
      &         lex_ids,                                     &  !< lex id of sibling tile
      &         lex_idc, lex_id_first, lex_id_last              !< lex id of current tile (lex_idc) and of the first (lex_id_first) and last tile (lex_id_last) on the present tile level

    REAL(wp) :: dtime                               !< Timestep length [s]

    REAL(wp), PARAMETER :: fos = 0.99_wp            !< Factor of safety to avoid negative fluxes due to roundoff errors [/]

    REAL(wp), DIMENSION(options%nc) ::  &
      &  fract,                                  &  !< Cover fractions of tile [/]
      &  fract_sibl,                             &  !< Cover fractions of sibling tile [/]
      &  snow_soil_tle_tm1,                      &  !< Snow soil (WEQ) on tile at end of previous timestep [m]
      &  snow_new_tle_tm1,                       &  !< New snow soil (WEQ) on tile at end of previous timestep [m]
      &  snow_dens_tle_tm1,                      &  !< Snow density on tile at end of previous timestep [kg m-3]
      &  wtr_incomming_srf,                      &  !< Lateral flux (surface) from sibling to tile [kg m-2 s-1]
      &  wtr_incomming_srf_sibl                     !< Lateral flux (surface) from tile to sibling [kg m-2 s-1]

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
    CHARACTER(len=5000)         :: message_text_long

    ! Declare process memories
    dsl4jsb_Def_config(LEX_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    dsl4jsb_Def_memory_tile(LEX_, sibl)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl)
    dsl4jsb_Def_memory_tile(SSE_, sibl)

    ! Declare pointers to variables in memory
    ! Pointers to variables in memory for tile
    dsl4jsb_Real2D_onChunk :: &
      & wtr_latflow_srf,      &
      & wtr_lex_net_srf,      &
      & weq_snowblow_in,      &
      & wtr_runoff_out,       &
      & weq_snowblow_out,     &
      & weq_pond,             &
      & weq_snow,             &
      & weq_snow_soil,        &
      & new_snow_soil,        &
      & snow_soil_dens,       &
      & runoff,               &
      & drainage

    dsl4jsb_Real3D_onChunk ::    &
      &  microcon_ctc_ct,        &
      &  microcon_hfp_blg_ct,    &
      &  fract_latflow_srf_ct,   &
      &  fract_latflow_blg_ct,   &
      &  wtr_latflow_sl,         &
      &  soil_depth_sl,          &
      &  wtr_lex_net_sl,         &
      &  wtr_infil_in_sl,        &
      &  wtr_diffusion_in_sl,    &
      &  drainage_sl,            &
      &  wtr_drainage_out_sl,    &
      &  wtr_diffusion_out_sl,   &
      &  wtr_soil_sl,            &
      &  ice_soil_sl,            &
      &  wtr_soil_pot_scool_sl,  &
      &  wtr_soil_sat_sl,        &
      &  wtr_soil_fc_sl,         &
      &  wtr_soil_res_sl,        &
      &  hyd_cond_sat_sl,        &
      &  vol_porosity_sl,        &
      &  vol_wres_sl,            &
      &  matric_pot_sl,          &
      &  bclapp_sl,              &
      &  pore_size_index_sl,     &
      &  t_soil_sl,              &
      &  heat_cond_sl,           &
      &  vol_heat_cap_sl,        &
      &  nrg_diffusion_in_sl,    &
      &  nrg_convection_in_sl,   &
      &  nrg_diffusion_out_sl,   &
      &  nrg_convection_out_sl,  &
      &  nrg_lex_net_sl,         &
      &  h_offset_ct,            &
      &  length_scale_ct,        &
      &  length_iface_ct

    ! Pointers to variables in memory for sibl
    dsl4jsb_Real2D_onChunk ::    &
      &  wtr_lex_net_srf_sibl,   &
      &  wtr_runoff_in_sibl,     &
      &  weq_snowblow_in_sibl,   &
      &  wtr_infil_out_sibl,     &
      &  weq_snowblow_out_sibl,  &
      &  wtr_pond_sibl,          &
      &  weq_pond_sibl,          &
      &  ice_pond_sibl,          &
      &  fract_pond_max_sibl,    &
      &  weq_pond_max_sibl,      &
      &  weq_snow_sibl,          &
      &  weq_snow_soil_sibl,     &
      &  new_snow_soil_sibl,     &
      &  snow_soil_dens_sibl,    &
      &  water_to_soil_lex_sibl

    dsl4jsb_Real3D_onChunk ::         &
      &  soil_depth_sl_sibl,          &
      &  wtr_lex_net_sl_sibl,         &
      &  wtr_drainage_in_sl_sibl,     &
      &  wtr_diffusion_in_sl_sibl,    &
      &  wtr_diffusion_out_sl_sibl,   &
      &  wtr_soil_sl_sibl,            &
      &  ice_soil_sl_sibl,            &
      &  wtr_soil_pot_scool_sl_sibl,  &
      &  wtr_soil_sat_sl_sibl,        &
      &  wtr_soil_fc_sl_sibl,         &
      &  wtr_soil_res_sl_sibl,        &
      &  hyd_cond_sat_sl_sibl,        &
      &  vol_porosity_sl_sibl,        &
      &  vol_wres_sl_sibl,            &
      &  matric_pot_sl_sibl,          &
      &  bclapp_sl_sibl,              &
      &  pore_size_index_sl_sibl,     &
      &  t_soil_sl_sibl,              &
      &  heat_cond_sl_sibl,           &
      &  vol_heat_cap_sl_sibl,        &
      &  nrg_diffusion_in_sl_sibl,    &
      &  nrg_convection_in_sl_sibl,   &
      &  nrg_diffusion_out_sl_sibl,   &
      &  nrg_convection_out_sl_sibl,  &
      &  nrg_lex_net_sl_sibl

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    nsoil   = options%nsoil_w
    dtime   = options%dtime

    ! Get local variables from tile information
    model  => Get_model(tile%owner_model_id)
    soil_e => Get_vgrid('soil_depth_energy')
    nsibl  = tile%parent_tile%Get_no_of_children()

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(LEX_)) RETURN

    ! First get first and last lex_id on current tile level
    lex_id_first = 0
    lex_id_last  = 0
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

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process configs
    dsl4jsb_Get_config(LEX_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get process memories for tile
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
    ! Get lex_id of tile and the cover fraction
    lex_idt = tile%lcts(1)%lex_id
    CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))

    ! Initialze those variables, where the values at the end of ast time step,
    ! i.e. unmodified by interactions with siblings are needed.
    DO ic = 1, nc
      snow_soil_tle_tm1(ic) = weq_snow_soil(ic)
      snow_new_tle_tm1(ic)  = new_snow_soil(ic)
      snow_dens_tle_tm1(ic) = snow_soil_dens(ic)
    END DO

    ! Loop over all tiles on the current tile level
    ! This includes tiles that ma not participate in the lateral exchange, such as glaciers.
    DO it = 1, nsibl

      ! Get lex_id of sibling
      IF (it == 1) THEN
        sibl => tile%parent_tile%Get_first_child_tile()
        lex_ids = sibl%lcts(1)%lex_id
      ELSE
        sibl => sibl%Get_next_sibling_tile()
        lex_ids = sibl%lcts(1)%lex_id
      END IF

      ! Only calculate processes on tiles for which connections have been defined
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

      CALL sibl%Get_fraction(ics, ice, iblk, fract=fract_sibl(:))

      ! Check if lateral flux fields have been reset
      IF (lex_idt == lex_id_first) THEN
        IF (ANY(wtr_lex_net_srf_sibl(:) /= 0._wp)) THEN
          WRITE (message_text_long,*) 'Reset of lateral flux fields did not work', &
            &                          MAXVAL(wtr_lex_net_srf_sibl(:))
          CALL finish(TRIM(routine), message_text_long)
        END IF
      END IF

      ! Start calculating exchange processes between different tiles
      IF (lex_idt /= lex_ids) THEN

        DO ic = 1, nc
          wtr_incomming_srf(ic)       = 0._wp
          wtr_incomming_srf_sibl(ic)  = 0._wp
        END DO

        DO is = 1, nsoil
          DO ic = 1, nc
            wtr_incomming_sl(ic,is)      = 0._wp
            wtr_incomming_sl_sibl(ic,is) = 0._wp
            wtr_latflow_sl_used(ic,is)   = wtr_latflow_sl(ic,is)
          END DO
        END DO

        ! Incomming later subsurface runoff is adjusted according to differences in soil depth (until bedrock),
        ! and for micro-scale connections, also for differences in surface elevation
        CALL shift_fluxes_sl(nc, nsoil,     &
          & microcon_hfp_blg_ct(:,lex_ids), &
          & h_offset_ct(:,lex_ids),         &
          & soil_depth_sl(:,:),             &
          & soil_depth_sl_sibl(:,:),        &
          & wtr_latflow_sl_used(:,:)        &
          )

        ! Lateral fluxes of surface runoff.
        ! The outflow from tile is actually calculated in HYDRO_ and here only the distributions on siblings is done
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

        IF (tile%contains_soil .AND. sibl%contains_soil) THEN
          ! Lateral fluxes of subsurface runoff.
          ! The outflow from tile is actually calculated in HYDRO_ and here only the distributions on siblings is done
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

            ! Remove excess water from soil and surface
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

          ! Calculating wsat, field capacity and mositure residual for tile disregarding soil ice content
          DO is = 1, nsoil
            DO ic = 1, nc
              weq_wsat_sl(ic,is) = vol_porosity_sl (ic,is) * soil_depth_sl(ic,is)
              weq_wres_sl(ic,is) = vol_wres_sl(ic,is)      * soil_depth_sl(ic,is)
            END DO
          END DO

          ! Get hydrological properties on tile
          CALL get_soilhyd_properties(                      &
            & dsl4jsb_Config(HYDRO_)%soilhydmodel,          &
            & dsl4jsb_Config(HYDRO_)%interpol_mean,         &
            & nc, nsoil,                                    &
            & soil_depth_sl(:,:),                           &
            & wtr_soil_sl(:,:),                             &
            & ice_soil_sl(:,:),                             &
            & wtr_soil_pot_scool_sl(:,:),                   &
            & weq_wsat_sl(:,:),                             &
            & weq_wres_sl(:,:),                             &
            & hyd_cond_sat_sl(:,:),                         &
            & matric_pot_sl(:,:),                           &
            & bclapp_sl(:,:),                               &
            & pore_size_index_sl(:,:),                      &
            & ice_impedance=ice_impedance_sl(:,:),          &
            & K=hyd_cond_sl(:,:),                           &
            & D=diffus_sl(:,:),                             &
            & mpot_act=mpot_act(:,:)                        &
            )

          ! Calculating wsat, field capacity and mositure residual for sibling disregarding soil ice content
          DO is = 1, nsoil
            DO  ic = 1, nc
              weq_wsat_sl_sibl     (ic,is) = vol_porosity_sl_sibl (ic,is) &
                &                          * soil_depth_sl_sibl(ic,is)
              weq_wres_sl_sibl     (ic,is) = vol_wres_sl_sibl(ic,is)      &
                &                          * soil_depth_sl_sibl(ic,is)
            END DO
          END DO

          ! Get hydrological properties on sibling (with updated soil water content)
          CALL get_soilhyd_properties(                      &
            & dsl4jsb_Config(HYDRO_)%soilhydmodel,          &
            & dsl4jsb_Config(HYDRO_)%interpol_mean,         &
            & nc, nsoil,                                    &
            & soil_depth_sl_sibl(:,:),                      &
            & wtr_soil_sl_sibl(:,:),                        &
            & ice_soil_sl_sibl(:,:),                        &
            & wtr_soil_pot_scool_sl_sibl(:,:),              &
            & weq_wsat_sl_sibl(:,:),                        &
            & weq_wres_sl_sibl(:,:),                        &
            & hyd_cond_sat_sl_sibl(:,:),                    &
            & matric_pot_sl_sibl(:,:),                      &
            & bclapp_sl_sibl(:,:),                          &
            & pore_size_index_sl_sibl(:,:),                 &
            & ice_impedance=ice_impedance_sl_sibl(:,:),     &
            & K=hyd_cond_sl_sibl(:,:),                      &
            & D=diffus_sl_sibl(:,:),                        &
            & mpot_act=mpot_act_sibl(:,:)                   &
            )

          IF (dsl4jsb_Config(LEX_)%l_infil_trns) THEN

            DO is = 1, nsoil
              DO  ic = 1, nc
                wtr_incomming_sl(ic,is) = wtr_incomming_sl(ic,is) - wtr_infil_in_sl(ic,is)
              END DO
            END DO

            ! Lateral ingress from pond to connected tile
            CALL calc_pond_soil_flux_water(              &
              &  nc, nsoil, dtime,                       &
              &  microcon_ctc_ct(:,lex_ids),             &
              &  h_offset_ct(:,lex_ids),                 &
              &  length_scale_ct(:,lex_ids),             &
              &  length_iface_ct(:,lex_ids),             &
              &  fract_pond_max_sibl(:),                 &
              &  fract(:),                               &
              &  fract_sibl(:),                          &
              &  soil_depth_sl(:,:),                     &
              &  wtr_soil_fc_sl(:,:),                    &
              &  hyd_cond_sl(:,:),                       &
              &  mpot_act(:,:),                          &
              &  wtr_infil_out_sibl(:),                  &
              &  wtr_infil_in_sl(:,:),                   &
              &  wtr_lex_net_srf_sibl(:),                &
              &  wtr_lex_net_sl(:,:),                    &
              &  wtr_pond_sibl(:),                       &
              &  wtr_soil_sl(:,:)                        &
              &  )

            DO is = 1, nsoil
              DO  ic = 1, nc
                wtr_incomming_sl(ic,is) = wtr_incomming_sl(ic,is) + wtr_infil_in_sl(ic,is)
              END DO
            END DO

          END IF

          IF (dsl4jsb_Config(LEX_)%l_dffsn_wtr) THEN
            DO ic = 1, nc
              IF (microcon_ctc_ct(ic,lex_ids) > 0._wp) THEN

                wtr_incomming_sl(ic,:)      = wtr_incomming_sl(ic,:)      &
                  &                         - wtr_diffusion_in_sl(ic,:)
                wtr_incomming_sl_sibl(ic,:) = wtr_incomming_sl_sibl(ic,:) &
                  &                         - wtr_diffusion_in_sl_sibl(ic,:)

                ! Diffusion between connected tiles (water)
                CALL calc_diffusion_water(nsoil, dtime,       &
                  &  h_offset_ct(ic,lex_ids),                 &
                  &  length_scale_ct(ic,lex_ids),             &
                  &  length_iface_ct(ic,lex_ids),             &
                  &  fract(ic),                               &
                  &  fract_sibl(ic),                          &
                  &  soil_depth_sl(ic,:),                     &
                  &  soil_depth_sl_sibl(ic,:),                &
                  &  wtr_soil_fc_sl(ic,:),                    &
                  &  wtr_soil_fc_sl_sibl(ic,:),               &
                  &  wtr_soil_sat_sl(ic,:),                   &
                  &  wtr_soil_sat_sl_sibl(ic,:),              &
                  &  wtr_soil_res_sl(ic,:),                   &
                  &  wtr_soil_res_sl_sibl(ic,:),              &
                  &  hyd_cond_sat_sl(ic,:),                   &
                  &  hyd_cond_sat_sl_sibl(ic,:),              &
                  &  hyd_cond_sl(ic,:),                       &
                  &  hyd_cond_sl_sibl(ic,:),                  &
                  &  mpot_act(ic,:),                          &
                  &  mpot_act_sibl(ic,:),                     &
                  &  wtr_diffusion_in_sl(ic,:),               &
                  &  wtr_diffusion_in_sl_sibl(ic,:),          &
                  &  wtr_diffusion_out_sl(ic,:),              &
                  &  wtr_diffusion_out_sl_sibl(ic,:),         &
                  &  wtr_lex_net_sl(ic,:),                    &
                  &  wtr_lex_net_sl_sibl(ic,:),               &
                  &  wtr_soil_sl(ic,:),                       &
                  &  wtr_soil_sl_sibl(ic,:)                   &
                  &  )

                wtr_incomming_sl(ic,:)      = wtr_incomming_sl(ic,:)      &
                  &                         + wtr_diffusion_in_sl(ic,:)
                wtr_incomming_sl_sibl(ic,:) = wtr_incomming_sl_sibl(ic,:) &
                  &                         + wtr_diffusion_in_sl_sibl(ic,:)

              END IF
            END DO
          END IF

          IF (dsl4jsb_Config(LEX_)%l_dffsn_nrg) THEN
            DO ic = 1, nc
              IF (microcon_ctc_ct(ic,lex_ids) > 0._wp) THEN

                ! Diffusion between connected tiles (energy).
                ! Since an update of the actual soil temperatures outside of the SSE routines is not intended,
                ! only the lateral heat flux is calculated and passed to the SSE routine as an energy source/sink.

                IF (h_offset_ct(ic,lex_ids) > 0._wp) THEN
                  hlp1 = length_iface_ct(ic,lex_ids) * fract_pond_max_sibl(ic)
                  CALL calc_pond_soil_flux_heat(              &
                    &  nsoil, dtime,                          &
                    &  soil_e%dz(:),                          &
                    &  h_offset_ct(ic,lex_ids),               &
                    &  length_scale_ct(ic,lex_ids),           &
                    &  hlp1,                                  &
                    &  fract(ic),                             &
                    &  fract_sibl(ic),                        &
                    &  t_soil_sl_sibl(ic,1),                  &
                    &  vol_heat_cap_sl_sibl(ic,1),            &
                    &  wtr_pond_sibl(ic), ice_pond_sibl(ic),  &
                    &  t_soil_sl(ic,:),                       &
                    &  heat_cond_sl(ic,:),                    &
                    &  vol_heat_cap_sl_sibl(ic,:),            &
                    &  nrg_diffusion_in_sl_sibl(ic,1),        &
                    &  nrg_diffusion_out_sl_sibl(ic,1),       &
                    &  nrg_lex_net_sl_sibl(ic,1),             &
                    &  nrg_diffusion_in_sl(ic,:),             &
                    &  nrg_diffusion_out_sl(ic,:),            &
                    &  nrg_lex_net_sl(ic,:)                   &
                    &  )
                END IF

                CALL calc_diffusion_heat(nsoil, dtime,        &
                  &  soil_e%dz(:),                            &
                  &  h_offset_ct(ic,lex_ids),                 &
                  &  length_scale_ct(ic,lex_ids),             &
                  &  length_iface_ct(ic,lex_ids),             &
                  &  fract(ic),                               &
                  &  fract_sibl(ic),                          &
                  &  t_soil_sl(ic,:),                         &
                  &  t_soil_sl_sibl(ic,:),                    &
                  &  heat_cond_sl(ic,:),                      &
                  &  heat_cond_sl_sibl(ic,:),                 &
                  &  vol_heat_cap_sl(ic,:),                   &
                  &  vol_heat_cap_sl_sibl(ic,:),              &
                  &  nrg_diffusion_in_sl(ic,:),               &
                  &  nrg_diffusion_in_sl_sibl(ic,:),          &
                  &  nrg_diffusion_out_sl(ic,:),              &
                  &  nrg_diffusion_out_sl_sibl(ic,:),         &
                  &  nrg_lex_net_sl(ic,:),                    &
                  &  nrg_lex_net_sl_sibl(ic,:)                &
                  &  )
              END IF
            END DO
          END IF

          IF (dsl4jsb_Config(LEX_)%l_cnvctn_nrg) THEN
            ! Fluxes reulting from the heat contained in the lateral fluxes of water
            ! Since an update of the actual soil temperatures outside of the SSE routines is not intended,
            ! only the lateral heat flux is calculated and passed to the SSE routine as an energy source/sink.
            CALL calc_convection_heat(nc, nsoil, dtime,  &
              &  soil_e%dz(:),                           &
              &  h_offset_ct(:,lex_ids),                 &
              &  fract(:),                               &
              &  fract_sibl(:),                          &
              &  wtr_incomming_srf(:),                   &
              &  wtr_incomming_srf_sibl(:),              &
              &  t_soil_sl(:,:),                         &
              &  t_soil_sl_sibl(:,:),                    &
              &  wtr_incomming_sl(:,:),                  &
              &  wtr_incomming_sl_sibl(:,:),             &
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

          ! Remove excess water from soil and surface
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
          ! Remove excess water from surface only
          CALL remove_excess_water(nc, 0,                &
            & weq_pond_max_sibl(:),                      &
            & ice_pond_sibl(:),                          &
            & wtr_pond_sibl(:),                          &
            & water_to_soil_lex_sibl(:)                  &
            )
        END IF

        ! Update total weq_pond on sibling.
        ! This is required because the water balance does not use wtr_pond and ice_pond directly
        DO ic = 1,nc
          weq_pond_sibl(ic) = wtr_pond_sibl(ic) + ice_pond_sibl(ic)
        END DO

        IF (dsl4jsb_Config(LEX_)%l_snow_trns) THEN
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

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "lex"
  !!
  SUBROUTINE aggregate_lex(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice

    !X if necessary: REAL(wp) :: dtime, steplen
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


    dsl4jsb_Aggregate_onChunk(SSE_, nrg_lex_net_sl,    weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_srf,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_srf,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_latflow_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_lex_net_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_to_soil_lex,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, new_snow_soil,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,      weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(LEX_, wtr_runoff_in,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, weq_snowblow_in,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_infil_in_sl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_drainage_in_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_runoff_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_infil_out,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, weq_snowblow_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_drainage_out_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, wtr_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_diffusion_in_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_diffusion_out_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_convection_in_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, nrg_convection_out_sl, weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(LEX_, microcon_ctc_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, microcon_hfp_srf_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, microcon_hfp_blg_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, macrocon_hfp_srf_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, macrocon_hfp_blg_ct,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, slope_eff,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fact_lag_srf,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fact_lag_blg,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, h_offset_ct,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, length_scale_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, length_iface_ct,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fract_latflow_srf_ct,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LEX_, fract_latflow_blg_ct,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lex


  ! ================================================================================================================================

#endif
END MODULE mo_lex_interface
