!>
!! @brief Initialization of the the assimilation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_pheno_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp, l_true, l_false
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize assimilation process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE pheno_init(tile, options)

    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, PHENO_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_pheno_constants

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(PHENO_)

    ! Use of process memories
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,      ONLY: t_pheno_memory
 
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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)                          !   TYPE(t_pheno_memory), POINTER      :: pheno__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk    :: l_growing_season
    dsl4jsb_Real2D_onChunk    :: lai_max
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(PHENO_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(PHENO_, l_growing_season)
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai_max)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    ! set growing season flag depending on whether this is a evergreen plant (TRUE) or not (FALSE)
    IF(lctlib%phenology_type == ievergreen) THEN
      l_growing_season(:) = l_true
    ELSE
      l_growing_season(:) = l_false
    ENDIF

    ! set lai_max from the site specific LAI from all_site_list.dat (saved in pheno_config:lai_max) and check for sensible values
    IF (dsl4jsb_Config(PHENO_)%lai_max > 0.5_wp .AND. dsl4jsb_Config(PHENO_)%lai_max < 8.0_wp) THEN
      lai_max(:) = dsl4jsb_Config(PHENO_)%lai_max
    ELSE
      lai_max(:) = 6.0_wp
    END IF

  END SUBROUTINE pheno_init

#endif
END MODULE mo_pheno_init
