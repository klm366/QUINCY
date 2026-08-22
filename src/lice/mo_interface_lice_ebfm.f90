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
!> Fields received from EBFM (state at the coupling instant):
!>   - icefract: ice-covered fraction of the grid cell on the EBFM mesh, conservatively
!>     remapped to the ICON grid; -1 where EBFM has no data
!>   - albedo: surface albedo of the EBFM cells, area-weighted over the covered part of
!>     the ICON cell; -1 where EBFM has no data
!>   - t_sub [K], ghf_cond [W m-2 K-1]: temperature of the first firn layer and conductance
!>     between surface and that layer (for the ground heat flux in the SEB); -1 where no data
!>   - runoff, smb [kg m-2 s-1], snowmass [kg m-2]: per unit ice area; 0 where no data
!>
!> Fields sent to EBFM (glacier tile values, put every time step, time-averaged by YAC
!> over the coupling period; JSBACH computes the surface energy balance, EBFM the firn):
!>   - t_srf [K], melt [kg m-2 s-1], evapotrans [kg m-2 s-1, negative upward]
!>   - diagnostics: hfss, hfls, rsns, rlns, ghf [W m-2]
!> The source fields are masked to cells with a glacier fraction > 0, so that the dummy
!> values of the glacier tile on ice-free cells are never mapped to the EBFM mesh.
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
  USE mo_coupling_utils,      ONLY: cpl_def_field, cpl_def_field_metadata, cpl_def_cell_field_mask, &
    &                               cpl_get_field, cpl_put_field
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  dsl4jsb_Use_processes LICE_, SEB_, HYDRO_, RAD_, SSE_
  dsl4jsb_Use_config(LICE_)
  dsl4jsb_Use_memory(LICE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SSE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsb_fdef_lice_fields
  PUBLIC :: interface_lice_get
  PUBLIC :: interface_lice_put

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_lice_ebfm'

  !> Names of the coupling fields (must match coupling.yaml and the EBFM side)
  ! received from EBFM
  CHARACTER(len=*), PARAMETER :: fname_icefract   = 'icefract'
  CHARACTER(len=*), PARAMETER :: fname_albedo     = 'albedo'
  CHARACTER(len=*), PARAMETER :: fname_t_sub      = 't_sub'
  CHARACTER(len=*), PARAMETER :: fname_ghf_cond   = 'ghf_cond'
  CHARACTER(len=*), PARAMETER :: fname_runoff     = 'runoff'
  CHARACTER(len=*), PARAMETER :: fname_smb        = 'smb'
  CHARACTER(len=*), PARAMETER :: fname_snowmass   = 'snowmass'
  ! sent to EBFM
  CHARACTER(len=*), PARAMETER :: fname_t_srf      = 't_srf'
  CHARACTER(len=*), PARAMETER :: fname_melt       = 'melt'
  CHARACTER(len=*), PARAMETER :: fname_evapotrans = 'evapotrans'
  CHARACTER(len=*), PARAMETER :: fname_hfss       = 'hfss'
  CHARACTER(len=*), PARAMETER :: fname_hfls       = 'hfls'
  CHARACTER(len=*), PARAMETER :: fname_rsns       = 'rsns'
  CHARACTER(len=*), PARAMETER :: fname_rlns       = 'rlns'
  CHARACTER(len=*), PARAMETER :: fname_ghf        = 'ghf'

  ! YAC field ids (one domain only)
  INTEGER :: field_id_icefract   = -1
  INTEGER :: field_id_albedo     = -1
  INTEGER :: field_id_t_sub      = -1
  INTEGER :: field_id_ghf_cond   = -1
  INTEGER :: field_id_runoff     = -1
  INTEGER :: field_id_smb        = -1
  INTEGER :: field_id_snowmass   = -1
  INTEGER :: field_id_t_srf      = -1
  INTEGER :: field_id_melt       = -1
  INTEGER :: field_id_evapotrans = -1
  INTEGER :: field_id_hfss       = -1
  INTEGER :: field_id_hfls       = -1
  INTEGER :: field_id_rsns       = -1
  INTEGER :: field_id_rlns       = -1
  INTEGER :: field_id_ghf        = -1

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

    TYPE(t_jsb_model),          POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),           POINTER :: grid

    dsl4jsb_Def_config(LICE_)

    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    LOGICAL, ALLOCATABLE :: is_valid(:)
    INTEGER :: nproma, nblks, jb, jc, cell_mask_id

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_fdef_lice_fields'

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

    CALL def_target_field(fname_t_sub,    'Temperature of the first firn layer of EBFM (K); -1 where no data', field_id_t_sub)
    CALL def_target_field(fname_ghf_cond, 'Conductance surface - first firn layer of EBFM (W m-2 K-1); -1 where no data', &
      &                                                                                                    field_id_ghf_cond)
    CALL def_target_field(fname_runoff,   'Runoff from the EBFM firn column (kg m-2 s-1); 0 where no data', field_id_runoff)
    CALL def_target_field(fname_smb,      'Climatic surface mass balance of EBFM (kg m-2 s-1); 0 where no data', field_id_smb)
    CALL def_target_field(fname_snowmass, 'Snow mass on top of the ice in EBFM (kg m-2); 0 where no data', field_id_snowmass)

    ! ---------------------------------------------------------------------
    ! Fields sent to EBFM (source fields)
    ! ---------------------------------------------------------------------
    ! The values are taken from the glacier tile, which carries dummy values on cells
    ! without glacier. A YAC source mask restricts the fields to cells with a glacier
    ! fraction > 0 (the tile fractions have been initialised in jsbach_init).
    tile => get_glacier_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No glacier tile found - the lice process requires use_glacier=.TRUE.')

    grid => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE(is_valid(nproma*nblks))
    is_valid(:) = .FALSE.
    DO jb = 1, nblks
      DO jc = 1, nproma
        IF (tile%fract(jc, jb) > 0._lp) is_valid((jb-1)*nproma+jc) = .TRUE.
      END DO
    END DO
    CALL cpl_def_cell_field_mask(routine, grid_id, is_valid, cell_mask_id)
    DEALLOCATE(is_valid)

    CALL def_source_field(fname_t_srf,      'Surface temperature of the glacier tile (K)',                  field_id_t_srf)
    CALL def_source_field(fname_melt,       'Snow/ice melt on the glacier tile (kg m-2 s-1)',               field_id_melt)
    CALL def_source_field(fname_evapotrans, 'Evapotranspiration incl. sublimation, negative upward (kg m-2 s-1)', &
      &                                                                                                  field_id_evapotrans)
    CALL def_source_field(fname_hfss,       'Sensible heat flux at the surface (W m-2)',                    field_id_hfss)
    CALL def_source_field(fname_hfls,       'Latent heat flux at the surface (W m-2)',                      field_id_hfls)
    CALL def_source_field(fname_rsns,       'Net surface shortwave radiation (W m-2)',                      field_id_rsns)
    CALL def_source_field(fname_rlns,       'Net surface longwave radiation (W m-2)',                       field_id_rlns)
    CALL def_source_field(fname_ghf,        'Ground heat flux (W m-2)',                                     field_id_ghf)

    fields_defined = .TRUE.

    CALL message(routine, 'Defined YAC fields of the ICON-Land/EBFM coupling on component ' // land_comp_name)

  CONTAINS

    SUBROUTINE def_target_field(name, metadata, field_id)
      CHARACTER(len=*), INTENT(in)  :: name, metadata
      INTEGER,          INTENT(out) :: field_id

      CALL cpl_def_field(                                          &
        & land_comp_id, cell_point_id, modelTimeStep,              &
        & name, 1, field_id)
      CALL cpl_def_field_metadata(routine, land_comp_name, grid_name, name, metadata)
    END SUBROUTINE def_target_field

    SUBROUTINE def_source_field(name, metadata, field_id)
      CHARACTER(len=*), INTENT(in)  :: name, metadata
      INTEGER,          INTENT(out) :: field_id

      CALL cpl_def_field(                                          &
        & land_comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
        & name, 1, field_id)
      CALL cpl_def_field_metadata(routine, land_comp_name, grid_name, name, metadata)
    END SUBROUTINE def_source_field

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
    dsl4jsb_Real2D_onDomain :: t_sub
    dsl4jsb_Real2D_onDomain :: ghf_cond
    dsl4jsb_Real2D_onDomain :: runoff
    dsl4jsb_Real2D_onDomain :: smb
    dsl4jsb_Real2D_onDomain :: snowmass

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
    dsl4jsb_Get_var2D_onDomain(LICE_, t_sub)      ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, ghf_cond)   ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, runoff)     ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, smb)        ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, snowmass)   ! OUT

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

    ! Receive firn state and mass balance
    ! -----------------------------------
    CALL get_target_field(field_id_t_sub,    fname_t_sub,    t_sub)
    CALL get_target_field(field_id_ghf_cond, fname_ghf_cond, ghf_cond)
    CALL get_target_field(field_id_runoff,   fname_runoff,   runoff)
    CALL get_target_field(field_id_smb,      fname_smb,      smb)
    CALL get_target_field(field_id_snowmass, fname_snowmass, snowmass)

  CONTAINS

    SUBROUTINE get_target_field(field_id, name, var)
      INTEGER,          INTENT(in)    :: field_id
      CHARACTER(len=*), INTENT(in)    :: name
      REAL(lp), CONTIGUOUS, TARGET, INTENT(inout) :: var(:,:)

      LOGICAL :: lreceived

      CALL cpl_get_field(routine, field_id, name, grid%ntotal, var, first_get=.FALSE., received_data=lreceived)
      IF (lreceived) THEN
        !$ACC UPDATE DEVICE(var) ASYNC(acc_stream)
        CALL message(routine, 'Received ' // name // ' from EBFM')
      END IF
    END SUBROUTINE get_target_field

  END SUBROUTINE interface_lice_get

  ! ======================================================================================================= !
  !> Send the glacier tile surface energy balance results to EBFM
  !>
  !> Called from jsbach_finish_timestep (outside of the block loop) every time step; YAC
  !> averages the fields over the coupling period.
  !>
  SUBROUTINE interface_lice_put(model)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model

    dsl4jsb_Def_config(LICE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),           POINTER :: grid

    ! Pointers to variables in memory of the glacier tile
    dsl4jsb_Real2D_onDomain :: t
    dsl4jsb_Real2D_onDomain :: snowmelt
    dsl4jsb_Real2D_onDomain :: evapotrans
    dsl4jsb_Real2D_onDomain :: sensible_hflx
    dsl4jsb_Real2D_onDomain :: latent_hflx
    dsl4jsb_Real2D_onDomain :: sw_srf_net
    dsl4jsb_Real2D_onDomain :: lw_srf_net
    dsl4jsb_Real2D_onDomain :: grnd_hflx

    LOGICAL :: write_coupler_restart

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_lice_put'

    IF (.NOT. is_coupled_run() .OR. .NOT. is_coupled_to_ebfm_land()) RETURN
    IF (.NOT. fields_defined) RETURN

    dsl4jsb_Get_config(LICE_)
    IF (.NOT. dsl4jsb_Config(LICE_)%active) RETURN

    tile => get_glacier_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No glacier tile found - the lice process requires use_glacier=.TRUE.')

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SSE_)

    dsl4jsb_Get_var2D_onDomain(SEB_,   t)              ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, snowmelt)       ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, evapotrans)     ! IN
    dsl4jsb_Get_var2D_onDomain(SEB_,   sensible_hflx)  ! IN
    dsl4jsb_Get_var2D_onDomain(SEB_,   latent_hflx)    ! IN
    dsl4jsb_Get_var2D_onDomain(RAD_,   sw_srf_net)     ! IN
    dsl4jsb_Get_var2D_onDomain(RAD_,   lw_srf_net)     ! IN
    dsl4jsb_Get_var2D_onDomain(SSE_,   grnd_hflx)      ! IN

    ! The fields are computed on the device; bring them to the host for the YAC put
    !$ACC UPDATE HOST(t, snowmelt, evapotrans, sensible_hflx, latent_hflx, sw_srf_net, lw_srf_net, grnd_hflx) &
    !$ACC   ASYNC(acc_stream)
    !$ACC WAIT(acc_stream)

    write_coupler_restart = .FALSE.
    CALL cpl_put_field(routine, field_id_t_srf,      fname_t_srf,      grid%ntotal, t,             &
      &                write_restart=write_coupler_restart)
    CALL cpl_put_field(routine, field_id_melt,       fname_melt,       grid%ntotal, snowmelt)
    CALL cpl_put_field(routine, field_id_evapotrans, fname_evapotrans, grid%ntotal, evapotrans)
    CALL cpl_put_field(routine, field_id_hfss,       fname_hfss,       grid%ntotal, sensible_hflx)
    CALL cpl_put_field(routine, field_id_hfls,       fname_hfls,       grid%ntotal, latent_hflx)
    CALL cpl_put_field(routine, field_id_rsns,       fname_rsns,       grid%ntotal, sw_srf_net)
    CALL cpl_put_field(routine, field_id_rlns,       fname_rlns,       grid%ntotal, lw_srf_net)
    CALL cpl_put_field(routine, field_id_ghf,        fname_ghf,        grid%ntotal, grnd_hflx)

    IF (write_coupler_restart) &
      & CALL message(routine, 'YAC says it is put for restart')

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
