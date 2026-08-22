!> Contains the YAC coupling interface between the ICON-Land lice process and the EBFM firn model
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
!>#### YAC coupling of ICON-Land to the EBFM firn model
!>
!> ICON-Land is registered as its own YAC component (see `land_comp_name` in ICON's
!> mo_coupling_config), living on the atmosphere processes and using the cell points
!> of the atmosphere grid. The fields are defined from ICON's construct_atmo_coupling
!> (before the YAC definitions are synchronized), the data exchange is done from
!> jsbach_start_timestep, i.e. outside of the block loop.
!>
!> Fields exchanged (phase 1):
!>   - icefract (EBFM -> ICON-Land): ice-covered fraction of the grid cell on the EBFM
!>     mesh, conservatively remapped to the ICON grid; -1 where EBFM has no data
!>   - albedo (EBFM -> ICON-Land): surface albedo of the EBFM cells, area-weighted over
!>     the covered part of the ICON cell; -1 where EBFM has no data
!>
!> The received fields are written into the lice memory of the glacier tile;
!> aggregation to the parent tiles is done by the lice_diag task of the lice process.
!>
!> Follows the pattern of mo_interface_hd_ocean.f90 (HD -> ocean coupling).
!>
MODULE mo_interface_lice_ebfm
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: lp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_control,         ONLY: acc_stream
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_ebfm_land, land_comp_name
  USE mo_coupling_utils,      ONLY: cpl_def_field, cpl_def_field_metadata, cpl_get_field
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  dsl4jsb_Use_processes LICE_
  dsl4jsb_Use_config(LICE_)
  dsl4jsb_Use_memory(LICE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsb_fdef_lice_fields
  PUBLIC :: interface_lice_get
  PUBLIC :: interface_lice_put

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_lice_ebfm'

  !> Names of the coupling fields (must match coupling.yaml and the EBFM side)
  CHARACTER(len=*), PARAMETER :: fname_icefract = 'icefract'
  CHARACTER(len=*), PARAMETER :: fname_albedo   = 'albedo'

  ! YAC field ids (one domain only)
  INTEGER :: field_id_icefract = -1
  INTEGER :: field_id_albedo   = -1

  LOGICAL :: fields_defined = .FALSE.

  !> Cached pointer to the glacier tile (set on first use)
  CLASS(t_jsb_tile_abstract), POINTER :: glacier_tile => NULL()

CONTAINS

  ! ======================================================================================================= !
  !> Define the YAC fields of the ICON-Land/EBFM coupling
  !>
  !> Called from ICON's construct_atmo_coupling before the YAC definitions are synchronized.
  !>
  SUBROUTINE jsb_fdef_lice_fields(land_comp_id, cell_point_id, grid_id, grid_name)

    INTEGER,          INTENT(in) :: land_comp_id  !< yac component id of the land component
    INTEGER,          INTENT(in) :: cell_point_id !< yac point id (cell points of the atmosphere grid)
    INTEGER,          INTENT(in) :: grid_id       !< yac grid id (atmosphere grid)
    CHARACTER(len=*), INTENT(in) :: grid_name     !< yac grid name (for field metadata)

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(LICE_)

    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_fdef_lice_fields'

    IF (grid_id > 0) CONTINUE  ! avoid compiler warning about dummy argument not being used

    ! It is assumed that only one domain is active
    model => Get_model(1)

    dsl4jsb_Get_config(LICE_)
    IF (.NOT. dsl4jsb_Config(LICE_)%active) THEN
      CALL message(routine, 'WARNING: coupled_to_ebfm_land is set, but the ICON-Land lice process is not active ' // &
        &                   '(jsb_lice_nml: active=.FALSE.); no fields defined.')
      RETURN
    END IF

    IF (land_comp_id < 0) CALL finish(routine, 'Invalid YAC component id for the land component')

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ! All fields are defined here unconditionally; which of them are actually exchanged
    ! is determined by the coupling configuration (coupling.yaml). A defined field that is
    ! not coupled is simply never received (YAC_ACTION_NONE).

    ! ---------------------------------------------------------------------
    ! Fields received from EBFM (target fields)
    ! ---------------------------------------------------------------------
    ! No YAC mask is used: YAC writes an interpolated value or the fixed fallback value
    ! (-1, see coupling configuration) to every cell, so that -1 unambiguously means
    ! "no data from EBFM". Land/sea handling remains with ICON-Land.
    CALL cpl_def_field(                                          &
      & land_comp_id, cell_point_id, modelTimeStep,              &
      & fname_icefract, 1, field_id_icefract)
    CALL cpl_def_field_metadata(                                 &
      & routine, land_comp_name, grid_name, fname_icefract,      &
      & 'Ice-covered fraction of grid cell from EBFM (1); -1 where EBFM provides no data')

    CALL cpl_def_field(                                          &
      & land_comp_id, cell_point_id, modelTimeStep,              &
      & fname_albedo, 1, field_id_albedo)
    CALL cpl_def_field_metadata(                                 &
      & routine, land_comp_name, grid_name, fname_albedo,        &
      & 'Surface albedo from EBFM (1); -1 where EBFM provides no data')

    ! ---------------------------------------------------------------------
    ! Fields sent to EBFM (source fields) - not implemented yet
    ! ---------------------------------------------------------------------
    ! CALL cpl_def_field(land_comp_id, cell_point_id, modelTimeStep, '<name>', 1, field_id_<name>)

    fields_defined = .TRUE.

    CALL message(routine, 'Defined YAC fields of the ICON-Land/EBFM coupling on component ' // land_comp_name)

  END SUBROUTINE jsb_fdef_lice_fields

  ! ======================================================================================================= !
  !> Receive fields from EBFM and store them in the lice memory of the glacier tile
  !>
  !> Called from jsbach_start_timestep (outside of the block loop).
  !>
  SUBROUTINE interface_lice_get(model)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model

    dsl4jsb_Def_config(LICE_)
    dsl4jsb_Def_memory(LICE_)

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),           POINTER :: grid

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: icefract
    dsl4jsb_Real2D_onDomain :: albedo

    LOGICAL :: received_data

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_lice_get'

    IF (.NOT. is_coupled_run() .OR. .NOT. is_coupled_to_ebfm_land()) RETURN
    IF (.NOT. fields_defined) RETURN

    dsl4jsb_Get_config(LICE_)
    IF (.NOT. dsl4jsb_Config(LICE_)%active) RETURN

    tile => get_glacier_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No glacier tile found - the lice process requires use_glacier=.TRUE.')
    IF (.NOT. tile%Is_process_calculated(LICE_)) &
      & CALL finish(routine, 'The lice process is not active on the glacier tile')

    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory of the glacier tile
    dsl4jsb_Get_memory(LICE_)
    dsl4jsb_Get_var2D_onDomain(LICE_, icefract)   ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, albedo)     ! OUT

    ! The memory is only updated at coupling time steps (received_data=.TRUE.); in between,
    ! the last received (or restarted) values are kept.

    ! Receive ice fraction
    ! --------------------
    CALL cpl_get_field(                                           &
      & routine, field_id_icefract, fname_icefract, grid%ntotal,  &
      & icefract, first_get=.TRUE., received_data=received_data)
    IF (received_data) THEN
      !$ACC UPDATE DEVICE(icefract) ASYNC(acc_stream)
      CALL message(routine, 'Received ' // fname_icefract // ' from EBFM')
    END IF

    ! Receive surface albedo
    ! ----------------------
    CALL cpl_get_field(                                           &
      & routine, field_id_albedo, fname_albedo, grid%ntotal,      &
      & albedo, first_get=.FALSE., received_data=received_data)
    IF (received_data) THEN
      !$ACC UPDATE DEVICE(albedo) ASYNC(acc_stream)
      CALL message(routine, 'Received ' // fname_albedo // ' from EBFM')
    END IF

  END SUBROUTINE interface_lice_get

  ! ======================================================================================================= !
  !> Send fields to EBFM (not implemented yet)
  !>
  !> To be called from jsbach_finish_timestep (outside of the block loop), using the
  !> aggregated values of the box tile.
  !>
  SUBROUTINE interface_lice_put(model)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model

    dsl4jsb_Def_config(LICE_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_lice_put'

    IF (.NOT. is_coupled_run() .OR. .NOT. is_coupled_to_ebfm_land()) RETURN
    IF (.NOT. fields_defined) RETURN

    dsl4jsb_Get_config(LICE_)
    IF (.NOT. dsl4jsb_Config(LICE_)%active) RETURN

    ! Sending of fields to EBFM is not implemented yet. Note: before cpl_put_field the data must be
    ! brought to the host: !$ACC UPDATE HOST(<var>) ASYNC(acc_stream) followed by !$ACC WAIT(acc_stream), e.g.
    !   CALL cpl_put_field(routine, field_id_<name>, '<name>', grid%ntotal, <var>, write_restart=write_coupler_restart)

  END SUBROUTINE interface_lice_put

  ! ======================================================================================================= !
  !> Find the glacier tile in the tile tree of the model (cached after the first call)
  !>
  FUNCTION get_glacier_tile(model) RESULT(return_ptr)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model
    CLASS(t_jsb_tile_abstract), POINTER      :: return_ptr

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    IF (ASSOCIATED(glacier_tile)) THEN
      return_ptr => glacier_tile
      RETURN
    END IF

    return_ptr => NULL()

    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    DO WHILE (ASSOCIATED(tile))
      IF (tile%is_glacier) THEN
        return_ptr => tile
        EXIT
      END IF
      CALL model%Goto_next_tile(tile)
    END DO
    CALL model%Reset_tiles()

    glacier_tile => return_ptr

  END FUNCTION get_glacier_tile

#endif
END MODULE mo_interface_lice_ebfm
