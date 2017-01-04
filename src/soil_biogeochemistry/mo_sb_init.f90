!>
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sb_init

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_sb_config_class,       ONLY: sb_config
    USE mo_sb_constants

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)


  IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)        :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)
    TYPE(t_sb_memory),    POINTER                :: sb         ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(SB_)
    sb      => sb__mem 
    model   => Get_model(tile%owner_model_id)
    lctlib  => model%lctlib(tile%lcts%lib_id)


! add code here
! 
! 




  END SUBROUTINE sb_init

END MODULE mo_sb_init
