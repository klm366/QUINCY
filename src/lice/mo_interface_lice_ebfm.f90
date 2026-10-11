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
!> ICON-Land is its own YAC component (`land_comp_name` in mo_jsb_coupling, predefined
!> by jsb_construct_coupling before ICON defines its own components), living on the
!> atmosphere processes and using the cell points of the atmosphere grid. The fields are
!> defined from ICON's construct_atmo_coupling (before the YAC definitions are
!> synchronized), the data exchange is done from jsbach_start_timestep, i.e. outside of
!> the block loop.
!>
!> Fields received from EBFM (state at the coupling instant):
!>   - icefract: ice-covered fraction of the grid cell on the EBFM mesh, conservatively
!>     remapped to the ICON grid
!>   - albedo: surface albedo of the EBFM cells, area-weighted over the covered part of
!>     the ICON cell
!>   - t_sub [K], ghf_cond [W m-2 K-1], hcap_sub [J m-2 K-1]: temperature of the first firn
!>     layer, conductance between surface and that layer and heat capacity of the surface layer
!>     (for the ground heat flux term in the SEB of the lice tile)
!>   - runoff, smb [kg m-2 s-1], snowmass [kg m-2]: per unit ice area
!>   Where EBFM has no data, the state fields carry the CDI missing value (mo_jsb_io::missval,
!>   set as fixed value in the coupling configuration) and the mass fields are zero.
!>
!> Fields sent to EBFM (lice tile values, put every time step, time-averaged by YAC
!> over the coupling period; JSBACH computes the surface energy balance, EBFM the firn):
!>   - t_srf [K], melt [kg m-2 s-1], evapotrans [kg m-2 s-1, negative upward]
!> The source fields are masked to cells with a lice tile fraction > 0, so that the dummy
!> values of the lice tile on other cells are never mapped to the EBFM mesh.
!>
!> The received fields are written into the lice memory of the lice tile (the land ice tile
!> coupled to EBFM, a second glacier tile next to 'glac'); aggregation to the parent tiles is
!> done by the lice_diag task of the lice process. The fields sent to EBFM are the glacier
!> physics results of the lice tile.
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
  USE mo_coupling_config,     ONLY: is_coupled_run, is_coupled_to_ebfm_land
  USE mo_jsb_coupling,        ONLY: land_comp_name, land_comp_id
  USE mo_coupling_utils,      ONLY: cpl_def_field, cpl_def_field_metadata, cpl_def_cell_field_mask, &
    &                               cpl_get_field, cpl_put_field
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  dsl4jsb_Use_processes LICE_, SEB_, HYDRO_
  dsl4jsb_Use_config(LICE_)
  dsl4jsb_Use_memory(LICE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

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
  CHARACTER(len=*), PARAMETER :: fname_hcap_sub   = 'hcap_sub'
  CHARACTER(len=*), PARAMETER :: fname_runoff     = 'runoff'
  CHARACTER(len=*), PARAMETER :: fname_smb        = 'smb'
  CHARACTER(len=*), PARAMETER :: fname_snowmass   = 'snowmass'
  ! sent to EBFM
  CHARACTER(len=*), PARAMETER :: fname_t_srf      = 't_srf'
  CHARACTER(len=*), PARAMETER :: fname_melt       = 'melt'
  CHARACTER(len=*), PARAMETER :: fname_evapotrans = 'evapotrans'

  ! YAC field ids (one domain only)
  INTEGER :: field_id_icefract   = -1
  INTEGER :: field_id_albedo     = -1
  INTEGER :: field_id_t_sub      = -1
  INTEGER :: field_id_ghf_cond   = -1
  INTEGER :: field_id_hcap_sub   = -1
  INTEGER :: field_id_runoff     = -1
  INTEGER :: field_id_smb        = -1
  INTEGER :: field_id_snowmass   = -1
  INTEGER :: field_id_t_srf      = -1
  INTEGER :: field_id_melt       = -1
  INTEGER :: field_id_evapotrans = -1

  LOGICAL :: fields_defined = .FALSE.
  LOGICAL :: coverage_checked = .FALSE.   !< mesh coverage of the lice tile checked (once per run)

  !> Cached pointer to the lice tile (the tile on which the lice process runs; set on first use)
  CLASS(t_jsb_tile_abstract), POINTER :: lice_tile => NULL()

CONTAINS

  ! ======================================================================================================= !
  !> Define the YAC fields of the ICON-Land/EBFM coupling
  !>
  !> Called from ICON's construct_atmo_coupling before the YAC definitions are synchronized.
  !>
  SUBROUTINE jsb_fdef_lice_fields(cell_point_id, grid_id, grid_name)

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
    IF (.NOT. dsl4jsb_Config(LICE_)%active) CALL finish(routine, &
      & 'coupled_to_ebfm_land is set, but the ICON-Land lice process is not active (jsb_lice_nml: active=.FALSE.)')

    IF (land_comp_id < 0) CALL finish(routine, &
      & 'The YAC component ' // land_comp_name // ' has not been defined (jsb_construct_coupling)')

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ! All fields are defined here unconditionally; which of them are actually exchanged
    ! is determined by the coupling configuration (coupling.yaml). A defined field that is
    ! not coupled is simply never received (YAC_ACTION_NONE).

    ! ---------------------------------------------------------------------
    ! Fields received from EBFM (target fields)
    ! ---------------------------------------------------------------------
    ! No YAC mask is used: YAC writes an interpolated value or the fixed fallback value
    ! (the missing value, see coupling configuration) to every cell, so that the missing
    ! value unambiguously means "no data from EBFM". Land/sea handling remains with ICON-Land.
    CALL def_target_field(fname_icefract, 'Ice-covered fraction of grid cell from EBFM (1); missing value where EBFM provides no data', &
      &                                                                                                  field_id_icefract)
    CALL def_target_field(fname_albedo,   'Surface albedo from EBFM (1); missing value where EBFM provides no data', field_id_albedo)
    CALL def_target_field(fname_t_sub,    'Temperature of the first firn layer of EBFM (K); missing value where no data', field_id_t_sub)
    CALL def_target_field(fname_ghf_cond, 'Conductance surface - first firn layer of EBFM (W m-2 K-1); missing value where no data', &
      &                                                                                                    field_id_ghf_cond)
    CALL def_target_field(fname_hcap_sub, 'Heat capacity of the firn surface layer of EBFM (J m-2 K-1); missing value where no data', &
      &                                                                                                    field_id_hcap_sub)
    CALL def_target_field(fname_runoff,   'Runoff from the EBFM firn column (kg m-2 s-1); 0 where no data', field_id_runoff)
    CALL def_target_field(fname_smb,      'Climatic surface mass balance of EBFM (kg m-2 s-1); 0 where no data', field_id_smb)
    CALL def_target_field(fname_snowmass, 'Snow mass on top of the ice in EBFM (kg m-2); 0 where no data', field_id_snowmass)

    ! ---------------------------------------------------------------------
    ! Fields sent to EBFM (source fields)
    ! ---------------------------------------------------------------------
    ! The values are taken from the lice tile, which carries dummy values on cells without
    ! land ice coupled to EBFM. A YAC source mask restricts the fields to cells with a lice
    ! tile fraction > 0 (the tile fractions have been initialised in jsbach_init).
    tile => get_lice_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No lice tile found - the lice process must be ON_TILE_ on the lice tile')

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

    CALL def_source_field(fname_t_srf,      'Surface temperature of the lice tile (K)',                  field_id_t_srf)
    CALL def_source_field(fname_melt,       'Snow/ice melt on the lice tile (kg m-2 s-1)',               field_id_melt)
    CALL def_source_field(fname_evapotrans, 'Evapotranspiration incl. sublimation, negative upward (kg m-2 s-1)', &
      &                                                                                                  field_id_evapotrans)

    fields_defined = .TRUE.

    CALL message(routine, 'Defined YAC fields of the ICON-Land/EBFM coupling on component ' // land_comp_name)

  CONTAINS

    !> Field received from EBFM: defined without a mask, every cell of the ICON grid receives
    !> (cells without EBFM data get the missing value of the coupling configuration)
    SUBROUTINE def_target_field(name, metadata, field_id)
      CHARACTER(len=*), INTENT(in)  :: name, metadata
      INTEGER,          INTENT(out) :: field_id

      CALL cpl_def_field(                                          &
        & land_comp_id, cell_point_id, modelTimeStep,              &
        & name, 1, field_id)
      CALL cpl_def_field_metadata(routine, field_id, metadata)
    END SUBROUTINE def_target_field

    !> Field sent to EBFM: defined with the cell mask of the lice tile fraction (cell_mask_id),
    !> so that only cells with a lice tile contribute to the interpolation to the EBFM mesh
    SUBROUTINE def_source_field(name, metadata, field_id)
      CHARACTER(len=*), INTENT(in)  :: name, metadata
      INTEGER,          INTENT(out) :: field_id

      CALL cpl_def_field(                                          &
        & land_comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
        & name, 1, field_id)
      CALL cpl_def_field_metadata(routine, field_id, metadata)
    END SUBROUTINE def_source_field

  END SUBROUTINE jsb_fdef_lice_fields

  ! ======================================================================================================= !
  !> Receive fields from EBFM and store them in the lice memory of the lice tile
  !>
  !> Called from jsbach_start_timestep (outside of the block loop).
  !>
  SUBROUTINE interface_lice_get(model, l_experiment_start)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model
    LOGICAL,           INTENT(in)            :: l_experiment_start !< True: very first time step of the
                                                                   !< experiment (not of a restart leg)

    dsl4jsb_Def_memory(LICE_)
    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),           POINTER :: grid

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: icefract
    dsl4jsb_Real2D_onDomain :: albedo
    dsl4jsb_Real2D_onDomain :: t_sub
    dsl4jsb_Real2D_onDomain :: ghf_cond
    dsl4jsb_Real2D_onDomain :: hcap_sub
    dsl4jsb_Real2D_onDomain :: runoff
    dsl4jsb_Real2D_onDomain :: smb
    dsl4jsb_Real2D_onDomain :: snowmass

    ! Surface temperature state of the SEB process on the lice tile (experiment-start init)
    REAL(lp), POINTER :: seb_t(:,:), seb_t_unfilt(:,:), seb_t_filt(:,:), seb_t_srf(:,:)

    LOGICAL :: received_icefract, received_t_sub

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_lice_get'

    IF (.NOT. is_coupled_run() .OR. .NOT. is_coupled_to_ebfm_land()) RETURN
    ! The fields are defined by jsb_fdef_lice_fields (called from the ICON coupling frame when
    ! coupled_to_ebfm_land is set), which aborts if the lice process is not active
    IF (.NOT. fields_defined) CALL finish(routine, &
      & 'coupled_to_ebfm_land is set, but the YAC fields of the ICON-Land/EBFM coupling are not defined')

    tile => get_lice_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No lice tile found - the lice process must be ON_TILE_ on the lice tile')

    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory of the lice tile
    dsl4jsb_Get_memory(LICE_)
    dsl4jsb_Get_var2D_onDomain(LICE_, icefract)   ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, albedo)     ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, t_sub)      ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, ghf_cond)   ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, hcap_sub)   ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, runoff)     ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, smb)        ! OUT
    dsl4jsb_Get_var2D_onDomain(LICE_, snowmass)   ! OUT

    ! The memory is only updated at coupling time steps (received_data=.TRUE.); in between,
    ! the last received (or restarted) values are kept.

    ! Receive ice fraction, albedo, firn state and mass balance
    ! ----------------------------------------------------------
    ! (the first get of the time step is flagged for ICON's coupling timers)
    CALL get_target_field(field_id_icefract, fname_icefract, icefract, first_get=.TRUE., received=received_icefract)
    IF (received_icefract .AND. .NOT. coverage_checked) THEN
      ! Once per run: every cell with a lice tile fraction must be covered by the EBFM mesh
      CALL check_mesh_coverage()
      coverage_checked = .TRUE.
    END IF
    CALL get_target_field(field_id_albedo,   fname_albedo,   albedo)
    CALL get_target_field(field_id_t_sub,    fname_t_sub,    t_sub, received=received_t_sub)
    CALL get_target_field(field_id_ghf_cond, fname_ghf_cond, ghf_cond)
    CALL get_target_field(field_id_hcap_sub, fname_hcap_sub, hcap_sub)
    CALL get_target_field(field_id_runoff,   fname_runoff,   runoff)
    CALL get_target_field(field_id_smb,      fname_smb,      smb)
    CALL get_target_field(field_id_snowmass, fname_snowmass, snowmass)

    ! Initialize the surface temperature of the lice tile from the firn state
    ! ------------------------------------------------------------------------
    ! At the very first time step of an experiment (cold start; restart legs get the temperature
    ! from the restart file) the surface temperature of the lice tile is set to the received first
    ! firn layer temperature t_sub. This replaces the constant initval (280 K) of the SEB memory,
    ! which is above the melting point and would otherwise melt several mm w.e. on the ice sheet
    ! while the surface relaxes during the first coupling window.
    IF (l_experiment_start) THEN
      ! The initial firn state must have arrived in this very first get; otherwise the coupling
      ! configuration is wrong (e.g. no ebfm2land couple), and the run would silently continue
      ! with the SEB initval on the lice tile.
      IF (.NOT. received_t_sub) CALL finish(routine, &
        & 'No ' // fname_t_sub // ' received from EBFM at experiment start: check the coupling configuration')
      dsl4jsb_Get_memory(SEB_)
      seb_t        => SEB__mem%t%ptr(:,:)
      seb_t_unfilt => SEB__mem%t_unfilt%ptr(:,:)
      WHERE (ghf_cond(:,:) >= 0._lp)   ! only where EBFM provides data (no missing value)
        seb_t(:,:)        = t_sub(:,:)
        seb_t_unfilt(:,:) = t_sub(:,:)
      END WHERE
      !$ACC UPDATE DEVICE(seb_t, seb_t_unfilt) ASYNC(acc_stream)
      ! t_filt and t_srf are allocated with TMX (or on LAND_TYPE tiles); check to be safe
      IF (ASSOCIATED(SEB__mem%t_filt%ptr)) THEN
        seb_t_filt => SEB__mem%t_filt%ptr(:,:)
        WHERE (ghf_cond(:,:) >= 0._lp) seb_t_filt(:,:) = t_sub(:,:)
        !$ACC UPDATE DEVICE(seb_t_filt) ASYNC(acc_stream)
      END IF
      IF (ASSOCIATED(SEB__mem%t_srf%ptr)) THEN
        seb_t_srf => SEB__mem%t_srf%ptr(:,:)
        WHERE (ghf_cond(:,:) >= 0._lp) seb_t_srf(:,:) = t_sub(:,:)
        !$ACC UPDATE DEVICE(seb_t_srf) ASYNC(acc_stream)
      END IF
      CALL message(routine, 'Initialized the lice tile surface temperature from the EBFM firn state (t_sub)')
    END IF

  CONTAINS

    ! Warn if grid cells have a lice tile fraction but receive no data from EBFM, i.e. are not
    ! covered by the EBFM mesh (icefract carries the fixed value < 0 there). On such cells the lice
    ! tile runs like a JSBACH glacier tile (runoff P - E, no firn storage), which is consistent
    ! but not what the tile is meant for: the fraction file (fract_lice) and the mesh disagree.
    SUBROUTINE check_mesh_coverage()

      USE mo_sync, ONLY: global_sum_array

      LOGICAL,  POINTER     :: is_owned(:,:)
      REAL(lp), ALLOCATABLE :: uncovered(:,:)
      REAL(lp)              :: n_uncovered
      INTEGER               :: jc, jb
      CHARACTER(len=32)     :: n_str

      is_owned => grid%patch%cells%decomp_info%owner_mask(:,:)
      ALLOCATE(uncovered(SIZE(icefract,1), SIZE(icefract,2)))
      uncovered(:,:) = 0._lp
      DO jb = 1, SIZE(icefract,2)
        DO jc = 1, SIZE(icefract,1)
          IF (is_owned(jc,jb) .AND. tile%fract(jc,jb) > 0._lp .AND. icefract(jc,jb) < 0._lp) uncovered(jc,jb) = 1._lp
        END DO
      END DO
      n_uncovered = global_sum_array(uncovered)
      DEALLOCATE(uncovered)

      WRITE(n_str, '(I0)') NINT(n_uncovered)
      IF (n_uncovered > 0.5_lp) THEN
        CALL message(routine, 'WARNING: ' // TRIM(n_str) // ' grid cells have a lice tile fraction > 0 but ' // &
          & 'are not covered by the EBFM mesh (no data received); the lice tile runs there like a JSBACH ' // &
          & 'glacier tile (runoff P - E). Check fract_lice in bc_land_frac against the EBFM/Elmer mesh.')
      ELSE
        CALL message(routine, 'Mesh coverage check: all grid cells with a lice tile fraction receive EBFM data')
      END IF

    END SUBROUTINE check_mesh_coverage

    SUBROUTINE get_target_field(field_id, name, var, first_get, received)
      INTEGER,          INTENT(in)    :: field_id
      CHARACTER(len=*), INTENT(in)    :: name
      REAL(lp), CONTIGUOUS, TARGET, INTENT(inout) :: var(:,:)
      LOGICAL, OPTIONAL, INTENT(in)   :: first_get   !< first get of the time step (coupling timers)
      LOGICAL, OPTIONAL, INTENT(out)  :: received    !< data received by this call

      LOGICAL :: lreceived

      CALL cpl_get_field(routine, field_id, name, grid%ntotal, var, first_get=first_get, received_data=lreceived)
      IF (lreceived) THEN
        !$ACC UPDATE DEVICE(var) ASYNC(acc_stream)
        CALL message(routine, 'Received ' // name // ' from EBFM')
      END IF
      IF (PRESENT(received)) received = lreceived
    END SUBROUTINE get_target_field

  END SUBROUTINE interface_lice_get

  ! ======================================================================================================= !
  !> Send the lice tile surface energy balance results to EBFM
  !>
  !> Called from jsbach_finish_timestep (outside of the block loop) every time step; YAC
  !> averages the fields over the coupling period.
  !>
  SUBROUTINE interface_lice_put(model)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_grid),           POINTER :: grid

    ! Pointers to variables in memory of the lice tile
    dsl4jsb_Real2D_onDomain :: t
    dsl4jsb_Real2D_onDomain :: snowmelt
    dsl4jsb_Real2D_onDomain :: evapotrans

    LOGICAL :: write_coupler_restart

    CHARACTER(len=*), PARAMETER :: routine = modname//':interface_lice_put'

    IF (.NOT. is_coupled_run() .OR. .NOT. is_coupled_to_ebfm_land()) RETURN
    ! The fields are defined by jsb_fdef_lice_fields (called from the ICON coupling frame when
    ! coupled_to_ebfm_land is set), which aborts if the lice process is not active
    IF (.NOT. fields_defined) CALL finish(routine, &
      & 'coupled_to_ebfm_land is set, but the YAC fields of the ICON-Land/EBFM coupling are not defined')

    tile => get_lice_tile(model)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(routine, 'No lice tile found - the lice process must be ON_TILE_ on the lice tile')

    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(SEB_,   t)              ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, snowmelt)       ! IN
    dsl4jsb_Get_var2D_onDomain(HYDRO_, evapotrans)     ! IN

    ! The fields are computed on the device; bring them to the host for the YAC put
    !$ACC UPDATE HOST(t, snowmelt, evapotrans) ASYNC(acc_stream)
    !$ACC WAIT(acc_stream)

    write_coupler_restart = .FALSE.
    CALL cpl_put_field(routine, field_id_t_srf,      fname_t_srf,      grid%ntotal, t,             &
      &                write_restart=write_coupler_restart)
    CALL cpl_put_field(routine, field_id_melt,       fname_melt,       grid%ntotal, snowmelt)
    CALL cpl_put_field(routine, field_id_evapotrans, fname_evapotrans, grid%ntotal, evapotrans)

    IF (write_coupler_restart) &
      & CALL message(routine, 'YAC says it is put for restart')

  END SUBROUTINE interface_lice_put

  ! ======================================================================================================= !
  !> Find the lice tile, i.e. the tile on which the lice process is calculated (ON_TILE_), in the
  !> tile tree of the model (cached after the first call)
  !>
  FUNCTION get_lice_tile(model) RESULT(return_ptr)

    TYPE(t_jsb_model), INTENT(inout), TARGET :: model
    CLASS(t_jsb_tile_abstract), POINTER      :: return_ptr

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER :: itile

    IF (ASSOCIATED(lice_tile)) THEN
      return_ptr => lice_tile
      RETURN
    END IF

    return_ptr => NULL()

    ! model%tiles(:) is the precomputed pre-order list of the tile tree (parents first)
    DO itile = 1, SIZE(model%tiles)
      tile => model%tiles(itile)%p
      IF (tile%is_glacier .AND. tile%Is_process_calculated(LICE_)) THEN
        return_ptr => tile
        EXIT
      END IF
    END DO

    lice_tile => return_ptr

  END FUNCTION get_lice_tile

#endif
END MODULE mo_interface_lice_ebfm
