!>
!! @brief Initialization of the the anthropogenic2land memory.
!!
!! @author
!!  Rackhun Son, MPI-BGC Jena
!!
!! @par Revision History
!!  First version by Rackhun Son (20221229)
!!
MODULE mo_anthro2l_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: anthro2l_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_anthro2l_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Run anthropogenic2land init
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE anthro2l_init(tile)

    USE mo_jsb_class,                        ONLY: Get_model
    USE mo_jsb_tile_class,                   ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,                  ONLY: t_jsb_model
    USE mo_jsb_process_class,                ONLY: ANTHRO2L_

    ! Use of process configurations (t_PROC_config)
    !dsl4jsb_Use_config(ANTHRO2L_)    
    
    ! USE PROC_memory
    dsl4jsb_Use_memory(ANTHRO2L_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':anthro2l_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(ANTHRO2L_)
    dsl4jsb_Def_memory(ANTHRO2L_)
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (ASSOCIATED(tile%parent_tile)) RETURN ! should only run on the most upper tile of the hierarchy

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of anthro2l memory (DL-fire) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    !dsl4jsb_Get_config(ANTHRO2L_)
    ! Get process memories
    dsl4jsb_Get_memory(ANTHRO2L_)
    ! Get process variables (Set pointers to variables in memory)
    !dsl4jsb_Get_var2D_onDomain(ANTHRO2L_, swpar_fdiffuse)


  END SUBROUTINE anthro2l_init

#endif
END MODULE mo_anthro2l_init
