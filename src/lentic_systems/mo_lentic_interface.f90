! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format, e.g. "2016-04-30"                             ===
! === ... LENTIC with your process name, e.g. "ALBEDO"                                                     ===
! === ... lentic with your process name, e.g. "albedo"                                                     ===
! === ... <TASK>                    with the name of a task of the process, e.g. "snow_albedo"                                ===
! === ... <EXPLANATIION_OF_PROCESS> with a short explanation line, e.g.                                                       ===
! ===     "computes the surface temperature from latent and sensible heat flux, net radiation and below ground sensible heat" ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementations points (Implementation:) that you    ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================



!>
!! @brief Contains the interfaces to the processes related to lentic systems
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (2022-08-17)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lentic_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes LENTIC_, A2L_

  ! Use of process configurations
  dsl4jsb_Use_config(LENTIC_)

  ! Use of process memories
  dsl4jsb_Use_memory(LENTIC_)
  dsl4jsb_Use_memory(A2L_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_lentic_tasks !, t_lentic_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lentic_interface'
  
  !> Type definition for lentic process
!   TYPE, EXTENDS(t_jsb_process) :: t_lentic_process
!   CONTAINS
!     PROCEDURE :: Register_tasks => Register_lentic_tasks        !< Registers tasks for lentic process
!   END TYPE t_lentic_process
  
  !> Type definition for lentic task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_lentic_wbphysics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_waterbody_physics     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_waterbody_physics  !< Aggregates computed task variables
  END TYPE tsk_lentic_wbphysics
  
  !> Constructor interface for waterbody_physics
  INTERFACE tsk_lentic_wbphysics
    PROCEDURE Create_task_waterbody_physics                        !< Constructor function for task
  END INTERFACE tsk_lentic_wbphysics

  !X Repeat here for each additional task of the process the two blocks above ...
  
CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for lentic task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "lentic"
  !!
  FUNCTION Create_task_lentic_wbphysics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_lentic_wbphysics::return_ptr)
    CALL return_ptr%Construct(name='lentic_wbphysics', process_id=LENTIC_, owner_model_id=model_id)

  END FUNCTION Create_task_lentic_wbphysics
  

  !X Repeat here for each additional task of the process the function above ...
    
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for lentic process
  !!
  !! @param[in,out] this      Instance of lentic process class  
  !! @param[in]     model_id  Model id  
  !!
  SUBROUTINE Register_lentic_tasks(this, model_id)

    CLASS(t_lentic_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in)    :: model_id

    CALL this%Register_task(tsk_lentic_wbphysics(model_id))
    !X Repeat the line above here for each additional task of the process ...
    
  END SUBROUTINE Register_lentic_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "lentic_wbphysics"
  !! Task "lentic_wbphysics" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_lentic_wbphysics(tile, options)

  !!  USE mo_rad_process, ONLY: calc_radiation_surface_net

!!#ifdef SERIALIZE
    !$claw verbatim USE mo_ser_common, ONLY: init
!!#endif

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_lentic_wbphysics'


    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(LENTIC_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(LENTIC_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
    ! atmospheric forcing
      & wind_air,             &
      & t_air,                &
      & q_air,                &
      & press_srf,            &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
    ! topography and location
      & depth_w,              &
      & fetch,                &
      & par_Coriolis,         &
      & pondtype,             &
    ! FLake inout
      & t_snow,               &
      & t_lice,               &
      & t_mnw,                &
      & t_wml,                &
      & t_bot,                &
      & c_t,                  &
      & h_ice,                &
      & h_ml,                 &
      & t_sfc

    REAL(wp) ::               &
      & w_star(options%nc),   &
      & u_star(options%nc),   &
      & I_bot(options%nc),    &
      & Q_watvap(options%nc), &
      & steplen
    
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    steplen = options%steplen
    
    TYPE(t_jsb_model), POINTER :: model

    

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(LENTIC_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    ! Get process memories
    dsl4jsb_Get_memory(LENTIC_)
    dsl4jsb_Get_memory(A2L_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)     ! in

    dsl4jsb_Get_var2D_onChunk(LENTIC_, depth_w)         ! in
    dsl4jsb_Get_var2D_onChunk(LENTIC_, fetch)           ! in
    dsl4jsb_Get_var2D_onChunk(LENTIC_, par_Coriolis)    ! in
    dsl4jsb_Get_var2D_onChunk(LENTIC_, pondtype)        ! in

    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_lice)          ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_mnw)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_wml)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_bot)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, c_t)             ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, h_ice)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, h_ml)            ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_sfc)           ! inout

    USE flake_derivedtypes         
    USE flake_parameters , ONLY :   &
      tpl_T_f                     , & 
      tpl_rho_w_r                 , & 
      h_Snow_min_flk              , & 
      h_Ice_min_flk                   
    USE flake_paramoptic_ref       
    USE flake_albedo_ref           
    USE flake           , ONLY :    &
      flake_driver                , & 
      flake_radflux               , & 
      T_snow_p_flk, T_snow_n_flk  , &
      T_ice_p_flk, T_ice_n_flk    , & 
      T_mnw_p_flk, T_mnw_n_flk    , & 
      T_wML_p_flk, T_wML_n_flk    , & 
      T_bot_p_flk, T_bot_n_flk    , & 
      T_B1_p_flk, T_B1_n_flk      , & 
      C_T_p_flk, C_T_n_flk        , & 
      h_snow_p_flk, h_snow_n_flk  , & 
      h_ice_p_flk, h_ice_n_flk    , & 
      h_ML_p_flk, h_ML_n_flk      , & 
      H_B1_p_flk, H_B1_n_flk      , & 
      Q_snow_flk                  , & 
      Q_ice_flk                   , & 
      Q_w_flk                     , & 
      Q_bot_flk                   , & 
      I_atm_flk                   , & 
      I_snow_flk                  , & 
      I_ice_flk                   , & 
      I_w_flk                     , & 
      I_h_flk                     , & 
      I_bot_flk                   , & 
      I_intm_0_h_flk              , & 
      I_intm_h_D_flk              , & 
      Q_star_flk                  , & 
      u_star_w_flk                , & 
      w_star_sfc_flk              , & 
      dMsnowdt_flk                , &
      H_Ice_max
    USE SfcFlx          , ONLY :    &
      SfcFlx_lwradwsfc            , & 
      SfcFlx_momsenlat  

    ! ---------------------------
    ! Go    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_lentic_wbphysics
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "lentic"
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.  
  !! @param[in]     options Additional run-time parameters.  
  !!
  SUBROUTINE aggregate_lentic(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_lentic'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)    
    !X Implementation: Start your process scheme here...
    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lentic

  
  ! ================================================================================================================================
  !>
  !! @brief Another task for "lentic"
  !!
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !  ....
  
#endif
END MODULE mo_lentic_interface
