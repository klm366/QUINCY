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

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message

  USE mo_jsb_class,               ONLY: Get_model
  USE mo_base_memory,             ONLY: t_slm_tile, t_model
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, PHENO_, HD_, RAD_, ASSIMI_, PHENO_, SB_

  ! Use of process memories
  dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,      ONLY: t_pheno_memory

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
  SUBROUTINE pheno_init(tile)

    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_pheno_constants
 
    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)       :: tile       !< data structure for dynamic phenoetation and soil variables
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(PHENO_)                          !   TYPE(t_pheno_memory), POINTER      :: pheno__mem
    TYPE(t_pheno_memory), POINTER         :: pheno      ! local pointer to pheno__mem for SLM usability
    dsl4jsb_Get_memory(PHENO_)

    pheno  => pheno__mem
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 

    ! set growing season flag depending on whether this is a evergreen plant (TRUE) or not (FALSE)
    IF(lctlib%phenology_type == ievergreen) THEN
      pheno%growing_season%ptr=.TRUE.
    ELSE
      pheno%growing_season%ptr=.FALSE.
    ENDIF

  END SUBROUTINE pheno_init

END MODULE mo_pheno_init
