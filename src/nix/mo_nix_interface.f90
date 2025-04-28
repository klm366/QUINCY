!> Contains the interfaces to the NIX processes
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/radiation/mo_rad_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_nix_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options


  ! Use of processes in this module
  dsl4jsb_Use_processes NIX_

  ! Use process configurations
  dsl4jsb_Use_config(NIX_)

  ! Use of process memories
  dsl4jsb_Use_memory(NIX_)


  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_nix_tasks



  !> Type definition for nix layering task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_nix_layering
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_nix_layering
    PROCEDURE, NOPASS :: Aggregate => update_nix_layering
  END TYPE tsk_nix_layering

  INTERFACE tsk_nix_layering
    PROCEDURE Create_task_nix_layering
  END INTERFACE tsk_nix_layering


  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_interface'


CONTAINS



  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for nix layering task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "layering"
  !!
  FUNCTION Create_task_nix_layering(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_nix_layering::return_ptr)
    CALL return_ptr%Construct(name='nix_layering', process_id=NIX_, owner_model_id=model_id)

  END FUNCTION Create_task_nix_layering


  ! ======================================================================================================= !
  !> Register tasks for nix  process
  !>
  SUBROUTINE Register_nix_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_nix_layering                  (model_id))

  END SUBROUTINE Register_nix_tasks


  ! ================================================================================================================================
  !>
  !> Implementation to update the snow cover layering
  !! Task "nix_layering" calculates changes in snow cover layering
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_nix_layering(tile, options)

    USE mo_nix_process,             ONLY: calc_nix_layering

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options



    ! Declare pointers to process configs and memory
    dsl4jsb_Def_memory(NIX_)

   ! Declare pointers for variables in memory - 3D
    dsl4jsb_Real3D_onChunk :: &
      & t_sn

    dsl4jsb_Int2D_onChunk :: &
      & test_int2d


    ! Declare local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: iblk , ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_nix_layering'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    model => Get_model(tile%owner_model_id)

    ! Get memory
    dsl4jsb_Get_memory(NIX_)

    ! Get fields from memory
    dsl4jsb_Get_var3D_onChunk(NIX_,      t_sn)  !< inout
    dsl4jsb_Get_var2D_onChunk(NIX_,      test_int2d)  !< inout

    ! ---------------------------
    ! Go



    CALL calc_nix_layering( &
      & nc,                 &
      & t_sn(:,:)           &
      & )

    ! Clean up and finish
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_nix_layering





#endif

END MODULE mo_nix_interface





