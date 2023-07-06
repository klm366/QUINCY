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
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
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
    PROCEDURE, NOPASS :: Integrate => update_lentic_wbphysics    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_lentic_wbphysics !< Aggregates computed task variables
  END TYPE tsk_lentic_wbphysics
  
  !> Constructor interface for lentic_wbphysics
  INTERFACE tsk_lentic_wbphysics
    PROCEDURE Create_task_lentic_wbphysics                        !< Constructor function for task
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
  !! @param[in,out] process      Instance of lentic process class  
  !! @param[in]     model_id  Model id  
  !!
  SUBROUTINE Register_lentic_tasks(process, model_id)

    CLASS(t_jsb_process),    INTENT(inout) :: process
    INTEGER,                 INTENT(in)    :: model_id

    CALL process%Register_task(tsk_lentic_wbphysics(model_id))
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

    USE mo_lentic_process, ONLY: calc_surface_temperature, calc_radflux

    USE mo_lentic_constants , ONLY : & 
      albedo_water_ref, &
      albedo_blueice_ref  , & ! Blue ice
      albedo_whiteice_ref , & ! White ice
      c_albice_MR         , & ! Constant in the interpolation formula for the ice albedo (Mironov and Ritter 2004)
      tpl_T_f           ! Fresh water freezing point [K]

!!#ifdef SERIALIZE
    !$claw verbatim USE mo_ser_common, ONLY: init
!!#endif

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen!, dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_lentic_wbphysics'

    TYPE(t_jsb_model), POINTER :: model

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(LENTIC_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(LENTIC_) 



    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
    ! ! atmospheric forcing
    !   & wind_air,             &
    !   & t_air,                &
    !   & q_air,                &
    !   & press_srf,            &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
    ! ! topography and location
    !   & depth_w,              &
    !   & fetch,                &
    !   & par_Coriolis,         &
    !   & pondtype,             &
    ! ! FLake inout
    !   & t_snow,               &
    !   & t_lice,               &
    !   & t_mnw,                &
    !   & t_wml,                &
    !   & t_bot,                &
    !   & c_t,                  &
      & h_ice,                &
      & h_snow,               &
      & h_ml,                 &
      & t_sfc

     REAL(wp) ::                 &
        & albedo_water,   & 
        & albedo_ice(options%nc),     & 
        & albedo_snow(options%nc),    & 
    !   & w_star(options%nc),     &
    !   & u_star(options%nc),     &
        & I_atm(options%nc),      &
        & I_bot(options%nc),      &
        & I_h(options%nc),        &
        & I_ice(options%nc),      &
        & I_intm_0_h(options%nc), &
        & I_intm_h_D(options%nc), &
        & I_snow(options%nc),     & ! might not be needed
        & I_w(options%nc)!,        & 
    !   & Q_watvap(options%nc),   &
    !   & steplen
    
    
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    steplen = options%steplen
    

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(LENTIC_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    ! Get process memories
    dsl4jsb_Get_memory(LENTIC_)
    dsl4jsb_Get_memory(A2L_)

    ! Get process variables
    ! dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)           ! in
    ! dsl4jsb_Get_var2D_onChunk(A2L_, t_air)              ! in
    ! dsl4jsb_Get_var2D_onChunk(A2L_, q_air)              ! in
    ! dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)     ! in

    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, depth_w)         ! in
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, fetch)           ! in
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, par_Coriolis)    ! in
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, pondtype)        ! in

    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, t_snow)          ! inout
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, t_lice)          ! inout
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, t_mnw)           ! inout
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, t_wml)           ! inout
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, t_bot)           ! inout
    ! dsl4jsb_Get_var2D_onChunk(LENTIC_, c_t)             ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, h_ice)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, h_snow)           ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, h_ml)            ! inout
    dsl4jsb_Get_var2D_onChunk(LENTIC_, t_sfc)           ! inout

    ! USE mo_phy_schemes,   ONLY: lwnet_from_lwdown
   
    ! ---------------------------
    ! Go    

    !  Set albedos of the lake water, lake ice and snow

    albedo_water = albedo_water_ref
    ! Use empirical formulation proposed by Mironov and Ritter (2004) for GME 
    !_nu albedo_ice   = albedo_whiteice_ref 
    albedo_ice   = EXP(-c_albice_MR*(tpl_T_f-T_sfc)/tpl_T_f)
    albedo_ice   = albedo_whiteice_ref*(1._wp-albedo_ice) + albedo_blueice_ref*albedo_ice
    ! Snow is not considered
    albedo_snow  = albedo_ice
    

    !  Set initial values
    I_atm = lw_srf_down + swvis_srf_down + swnir_srf_down
    ! T_snow_p_flk = t_snow
    ! T_ice_p_flk  = t_lice
    ! T_mnw_p_flk  = t_mnw
    ! T_wML_p_flk  = t_wml
    ! T_bot_p_flk  = t_bot
    ! C_T_p_flk    = c_t
!    h_ML_n_flk   = h_ml


    !  Compute solar radiation fluxes (positive downward)    (needed?)
    ! I_atm_flk = I_atm_in
    ! Q_bot_flk = Q_bot_in

    call calc_radflux ( &
                        !& depth_w,        & ! in
                        & albedo_water,   & ! in
                        & albedo_ice,     & ! in
                        & albedo_snow,    & ! in
                        & h_ice,          & ! in 
                        & h_ML,           & ! in 
                        & h_snow,         & ! in
                        & I_atm,          & ! in
                        & I_bot,          & ! out
                        & I_h,            & ! out
                        & I_ice,          & ! out
                        & I_intm_0_h,     & ! out
                        & I_intm_h_D,     & ! out
                        & I_snow,         & ! out (but maybe not needed) 
                        & I_w             & ! out
                        )

    ! TODO: COMPARE THE CONSTANTS!
    !  Compute long-wave radiation fluxes (positive downward)
    ! Q_w_flk = Q_atm_lw_in                          ! Radiation of the atmosphere 

    ! lwnet_from_lwdown = zemiss_def * (lwdown - stbo * t_srf**4._wp)
    ! SfcFlx_lwradwsfc = c_lwrad_emis*tpsf_C_StefBoltz*T**4_iintegers
    ! Q_w_flk = Q_w_flk - SfcFlx_lwradwsfc(T_sfc_p)  ! Radiation of the surface (notice the sign)
    
    ! ! JSBACH VERSION!  
    ! ! Compute LW net radiation from incoming and the thermal radiation from the surface
    ! zlw_net = lwnet_from_lwdown(lw_down, t)

    ! TEST
    !call calc_surface_temperature(T_sfc )
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_lentic_wbphysics
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "lentic"
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.  
  !! @param[in]     options Additional run-time parameters.  
  !!
  SUBROUTINE aggregate_lentic_wbphysics(tile, options)

    ! Arguments
     TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(LENTIC_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_lentic_wbphysics'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LENTIC_)

    model             => Get_model(tile%owner_model_id)    
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(LENTIC_, h_ice,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LENTIC_, h_ml,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LENTIC_, h_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LENTIC_, T_sfc,        weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lentic_wbphysics

  
  ! ================================================================================================================================
  !>
  !! @brief Another task for "lentic"
  !!
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !  ....
  
#endif
END MODULE mo_lentic_interface
