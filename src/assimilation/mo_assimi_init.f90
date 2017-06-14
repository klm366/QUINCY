!>
!! @brief Initialization of the the assimilation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_assimi_init

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assimi_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize assimilation process (after memory has been set up)
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE assimi_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! USE PROC_memory
    dsl4jsb_Use_memory(ASSIMI_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':assimi_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: beta_air
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime_dacc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tphen_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, t_jmax_opt)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    beta_air(:)                     = 1.0_wp
    beta_soil_ps(:)                 = 1.0_wp
    beta_soil_gs(:)                 = 1.0_wp
    
    beta_air_daytime(:)             = 1.0_wp
    beta_air_daytime_dacc(:)        = 1.0_wp
    beta_soil_ps_daytime(:)         = 1.0_wp
    beta_soil_ps_daytime_dacc(:)    = 1.0_wp
    beta_soil_gs_daytime(:)         = 1.0_wp
    beta_soil_gs_daytime_dacc(:)    = 1.0_wp
    
    beta_air_tfrac_mavg(:)          = 1.0_wp
    beta_air_tcnl_mavg(:)           = 1.0_wp
    beta_soil_ps_tfrac_mavg(:)      = 1.0_wp
    beta_soil_ps_tcnl_mavg(:)       = 1.0_wp
    beta_soil_gs_tphen_mavg(:)      = 1.0_wp
    beta_soil_gs_tfrac_mavg(:)      = 1.0_wp
    beta_soil_gs_tcnl_mavg(:)       = 1.0_wp

    t_jmax_opt(:)                   = 1.0_wp

  END SUBROUTINE assimi_init

END MODULE mo_assimi_init
