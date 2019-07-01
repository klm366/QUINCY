!>
!! @brief Initialization of the the radiation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_rad_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rad_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize radiation process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE rad_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_rad_constants,         ONLY: rfr_ratio_toc

    ! USE PROC_memory
    dsl4jsb_Use_memory(RAD_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rad_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(RAD_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc_tvegdyn_mavg
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tcnl_mavg_cl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc    
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(RAD_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
!    dsl4jsb_Get_config(RAD_)
    ! Get process memories
    dsl4jsb_Get_memory(RAD_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(RAD_, rfr_ratio_boc)
    dsl4jsb_Get_var2D_onDomain(RAD_, rfr_ratio_boc_tvegdyn_mavg)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_shaded_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(RAD_, ppfd_shaded_tcnl_mavg_cl)
    

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    rfr_ratio_boc(:,:)                  = rfr_ratio_toc
    rfr_ratio_boc_tvegdyn_mavg(:,:)     = rfr_ratio_boc(:,:) 

    ppfd_sunlit_tfrac_mavg_cl(:,:,:)    = 150.0_wp
    ppfd_sunlit_tcnl_mavg_cl(:,:,:)     = 150.0_wp
    
    ppfd_shaded_tfrac_mavg_cl(:,:,:)    = 10.0_wp
    ppfd_shaded_tcnl_mavg_cl(:,:,:)     = 10.0_wp


  END SUBROUTINE rad_init

#endif
END MODULE mo_rad_init
