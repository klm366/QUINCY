!> Contains the interfaces to the land ice (lice) process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Tasks of the land ice (lice) process
!>
!> The lice process couples the ICON-Land glacier tile to external land ice models
!> (currently the EBFM firn model) via YAC. The actual data exchange (YAC put/get)
!> is done outside of the block loop in [[mo_interface_lice_ebfm]]. The task defined
!> here runs on the glacier tile and derives diagnostics from the received fields,
!> which are then aggregated to the parent tiles.
!>
!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
MODULE mo_lice_interface
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: lp
  USE mo_exception,         ONLY: message, finish
  USE mo_jsb_control,       ONLY: debug_on, acc_stream
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class, ONLY: t_jsb_process
  USE mo_jsb_task_class,    ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes LICE_

  ! Use of process configurations
  dsl4jsb_Use_config(LICE_)

  ! Use of process memories
  dsl4jsb_Use_memory(LICE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_lice_tasks

  !> Type definition for the lice_diag task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_lice_diag
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_lice_diag
    PROCEDURE, NOPASS :: Aggregate => aggregate_lice_diag
  END TYPE tsk_lice_diag

  !> Constructor interface for the lice_diag task
  INTERFACE tsk_lice_diag
    PROCEDURE Create_task_lice_diag
  END INTERFACE tsk_lice_diag

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lice_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Register tasks of the lice process
  !>
  SUBROUTINE Register_lice_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_lice_diag(model_id))

  END SUBROUTINE Register_lice_tasks

  ! ======================================================================================================= !
  !> Constructor for the lice_diag task
  !>
  FUNCTION Create_task_lice_diag(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_lice_diag::return_ptr)
    CALL return_ptr%Construct(name='lice_diag', process_id=LICE_, owner_model_id=model_id)

  END FUNCTION Create_task_lice_diag

  ! ======================================================================================================= !
  !> Diagnostics of the fields received from the external land ice model (runs on the glacier tile)
  !>
  !> - icefract_valid:  1 where a value has been received (icefract >= 0), 0 elsewhere
  !> - glac_fract_diff: icefract minus the glacier fraction of the grid cell in ICON-Land
  !>                    (the glacier tile fraction is the absolute glacier fraction of the grid cell)
  !>
  SUBROUTINE update_lice_diag(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of current block

    dsl4jsb_Def_config(LICE_)
    dsl4jsb_Def_memory(LICE_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: icefract
    dsl4jsb_Real2D_onChunk :: icefract_valid
    dsl4jsb_Real2D_onChunk :: glac_fract_diff
    dsl4jsb_Real2D_onChunk :: albedo
    dsl4jsb_Real2D_onChunk :: t_sub
    dsl4jsb_Real2D_onChunk :: ghf_cond

    REAL(lp), POINTER :: tile_fract(:)

    TYPE(t_jsb_model), POINTER :: model

    INTEGER :: iblk, ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_lice_diag'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(LICE_)) RETURN

    model => Get_model(tile%owner_model_id)
    IF (.NOT. model%Is_process_enabled(LICE_)) RETURN

    dsl4jsb_Get_config(LICE_)
    IF (.NOT. dsl4jsb_Config(LICE_)%active) RETURN

    IF (.NOT. tile%is_glacier) &
      & CALL finish(TRIM(routine), 'The lice process must run on the glacier tile (ON_TILE_), not on '//TRIM(tile%name))

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LICE_)

    dsl4jsb_Get_var2D_onChunk(LICE_, icefract)         ! in
    dsl4jsb_Get_var2D_onChunk(LICE_, icefract_valid)   ! out
    dsl4jsb_Get_var2D_onChunk(LICE_, glac_fract_diff)  ! out
    dsl4jsb_Get_var2D_onChunk(LICE_, albedo)           ! inout
    dsl4jsb_Get_var2D_onChunk(LICE_, t_sub)            ! inout
    dsl4jsb_Get_var2D_onChunk(LICE_, ghf_cond)         ! inout

    ! Make sure the state variables are undefined (-1) where EBFM does not provide data

    tile_fract => tile%fract(ics:ice, iblk)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      IF (icefract(ic) >= 0._lp) THEN
        icefract_valid(ic)  = 1._lp
        glac_fract_diff(ic) = icefract(ic) - tile_fract(ic)
      ELSE
        icefract_valid(ic)  = 0._lp
        glac_fract_diff(ic) = 0._lp
        albedo(ic)          = -1._lp
        t_sub(ic)           = -1._lp
        ghf_cond(ic)        = -1._lp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_lice_diag

  ! ======================================================================================================= !
  !> Aggregate lice variables from the glacier tile to the parent tiles (land, box)
  !>
  SUBROUTINE aggregate_lice_diag(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of current block

    dsl4jsb_Def_memory(LICE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    INTEGER :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_lice_diag'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LICE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()
    dsl4jsb_Aggregate_onChunk(LICE_, icefract,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, icefract_valid,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, glac_fract_diff, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, albedo,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, t_sub,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, ghf_cond,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, runoff,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, smb,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LICE_, snowmass,        weighted_by_fract)
    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lice_diag

#endif
END MODULE mo_lice_interface
