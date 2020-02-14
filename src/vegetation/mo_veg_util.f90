!>
!! @brief diverse subroutines and functions for the vegetation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2019-11-04)
!!
MODULE mo_veg_util
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reset_veg_fluxes

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> init/reset vegetation fluxes (with/to zero)
  !! 
  !! called prior to all Tasks (pheno_interface)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_veg_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_
    USE mo_jsb_math_constants,    ONLY: zero
    
    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_veg_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_) 
    ! Declare pointers to variables in memory
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: veg_growth_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
    !
    ! Get process memories
    dsl4jsb_Get_memory(VEG_) 
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_leaf_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_leaf_carbon)       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)            ! out



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    veg_growth_leaf_carbon(:)       = zero
    veg_litterfall_leaf_carbon(:)   = zero
    
    maint_respiration_pot(:)        = zero


  END SUBROUTINE reset_veg_fluxes

#endif
END MODULE mo_veg_util
