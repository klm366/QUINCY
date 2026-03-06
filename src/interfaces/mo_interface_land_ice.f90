!> Interface between JSBACH glacier tile and an external land-ice model (e.g. Elmer/Ice), through YAC coupler
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
!>#### Defines YAC coupling fields and exchanges data with the ice-sheet model each coupling timestep.
!>#### Follows the pattern of mo_interface_hd_ocean.f90 (HD -> ocean coupling).
!>
MODULE mo_interface_land_ice
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_LAND_ICE__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_coupling_utils,      ONLY: cpl_def_cell_field_mask, &
                                    cpl_def_field, &
                                    cpl_put_field, &
                                    cpl_get_field
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

  USE mo_jsb_task_class,      ONLY: t_jsb_task_options
  USE mo_jsb_control,         ONLY: acc_stream

  dsl4jsb_Use_processes LAND_ICE_, HYDRO_
  dsl4jsb_Use_config(LAND_ICE_)
  dsl4jsb_Use_memory(LAND_ICE_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsb_fdef_land_ice_fields
  PUBLIC :: interface_land_ice

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_land_ice'

  ! YAC field IDs for coupling with ice-sheet model
  INTEGER :: field_id_surface_temp
  INTEGER :: field_id_smb
  INTEGER :: field_id_ground_hflx
  INTEGER :: field_id_rainfall
  INTEGER :: field_id_snowfall
  INTEGER :: field_id_ice_thickness
  INTEGER :: field_id_ice_elevation
  INTEGER :: field_id_glacier_fract
  INTEGER :: field_id_basal_melt
  INTEGER :: field_id_calving

CONTAINS

  !> Define YAC coupling fields for the land-ice interface.
  !! Called during ICON initialization to register the coupling fields.
  !!
  SUBROUTINE jsb_fdef_land_ice_fields(comp_id, cell_point_id, grid_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid
    USE mo_jsb_grid,           ONLY: Get_grid

    INTEGER, INTENT(in) :: comp_id       ! yac component id
    INTEGER, INTENT(in) :: cell_point_id ! yac point id
    INTEGER, INTENT(in) :: grid_id       ! yac grid id

    ! Local variables
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    INTEGER :: cell_mask_id
    INTEGER :: model_id
    INTEGER :: jb, jc
    INTEGER :: nproma, nblks
    LOGICAL,  ALLOCATABLE :: is_valid(:)
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    CHARACTER(LEN=*), PARAMETER :: routine = modname // ':jsb_fdef_land_ice_fields'

    dsl4jsb_Def_config(LAND_ICE_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onDomain :: fract_glac

    ! It is assumed that only one domain is active (jg=1), i.e. no nesting.
    model_id = 1

    model => Get_model(model_id)

    dsl4jsb_Get_config(LAND_ICE_)
    IF (.NOT. dsl4jsb_Config(LAND_ICE_)%active) RETURN
    IF (.NOT. model%config%use_land_ice) RETURN

    CALL model%Get_top_tile(tile)
    dsl4jsb_Get_memory(HYDRO_)

    grid  => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ! Define the glacier mask for YAC.
    ! Coupling shall only occur on grid cells that have a glacier fraction.
    ALLOCATE(is_valid(nproma*nblks))
    is_valid(:) = .FALSE.

    dsl4jsb_Get_var2D_onDomain(HYDRO_, fract_glac)

    DO jb = 1, nblks
      DO jc = 1, nproma
        IF ( fract_glac(jc, jb) > 0.0_wp ) THEN
          is_valid((jb-1)*nproma+jc) = .TRUE.
        END IF
      END DO
    END DO

    CALL cpl_def_cell_field_mask(routine, grid_id, is_valid, cell_mask_id)

    DEALLOCATE(is_valid)

    ! --- Fields sent to the ice-sheet model (JSBACH -> Elmer/Ice) ---

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_surface_temp', 1, field_id_surface_temp)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_smb', 1, field_id_smb)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_ground_hflx', 1, field_id_ground_hflx)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_rainfall', 1, field_id_rainfall)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_snowfall', 1, field_id_snowfall)

    ! --- Fields received from the ice-sheet model (Elmer/Ice -> JSBACH) ---

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_ice_thickness', 1, field_id_ice_thickness)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_ice_elevation', 1, field_id_ice_elevation)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_glacier_fract', 1, field_id_glacier_fract)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_basal_melt', 1, field_id_basal_melt)

    CALL cpl_def_field( &
      & comp_id, cell_point_id, cell_mask_id, modelTimeStep, &
      & 'lice_calving', 1, field_id_calving)

  END SUBROUTINE jsb_fdef_land_ice_fields

  !> Exchange data with the ice-sheet model.
  !! Called from jsbach_finish_timestep to send surface fields
  !! and receive ice-sheet state updates via YAC.
  !!
  SUBROUTINE interface_land_ice(tile, options)

    USE mo_fortran_tools, ONLY: init

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(LAND_ICE_)
    dsl4jsb_Def_memory(LAND_ICE_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':interface_land_ice'

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: surface_temp
    dsl4jsb_Real2D_onDomain :: surface_mass_balance
    dsl4jsb_Real2D_onDomain :: ground_heat_flux
    dsl4jsb_Real2D_onDomain :: rainfall_rate
    dsl4jsb_Real2D_onDomain :: snowfall_rate
    dsl4jsb_Real2D_onDomain :: ice_thickness
    dsl4jsb_Real2D_onDomain :: ice_surface_elevation
    dsl4jsb_Real2D_onDomain :: glacier_fraction
    dsl4jsb_Real2D_onDomain :: basal_melt_rate
    dsl4jsb_Real2D_onDomain :: calving_flux

    ! Local variables
    LOGICAL :: write_coupler_restart
    INTEGER :: n

    n = options%nc  ! only to avoid compiler warnings

    IF ( .NOT. is_coupled_run() ) RETURN

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(LAND_ICE_)
    IF (.NOT. dsl4jsb_Config(LAND_ICE_)%active) RETURN
    IF (.NOT. model%config%use_land_ice) RETURN

    IF (ASSOCIATED(tile%parent)) &
      & CALL finish(routine, 'Land-ice coupling works on root tile only!')

    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(LAND_ICE_)

    ! --- Send fields to the ice-sheet model ---

    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, surface_temp)           ! IN
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, surface_mass_balance)   ! IN
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, ground_heat_flux)       ! IN
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, rainfall_rate)          ! IN
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, snowfall_rate)          ! IN

    !$ACC UPDATE HOST(surface_temp, surface_mass_balance, ground_heat_flux, rainfall_rate, snowfall_rate)

    CALL cpl_put_field( &
      & routine, field_id_surface_temp, 'lice_surface_temp', grid%ntotal, &
      & surface_temp, write_restart=write_coupler_restart)

    CALL cpl_put_field( &
      & routine, field_id_smb, 'lice_smb', grid%ntotal, &
      & surface_mass_balance, write_restart=write_coupler_restart)

    CALL cpl_put_field( &
      & routine, field_id_ground_hflx, 'lice_ground_hflx', grid%ntotal, &
      & ground_heat_flux, write_restart=write_coupler_restart)

    CALL cpl_put_field( &
      & routine, field_id_rainfall, 'lice_rainfall', grid%ntotal, &
      & rainfall_rate, write_restart=write_coupler_restart)

    CALL cpl_put_field( &
      & routine, field_id_snowfall, 'lice_snowfall', grid%ntotal, &
      & snowfall_rate, write_restart=write_coupler_restart)

    ! --- Receive fields from the ice-sheet model ---

    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, ice_thickness)          ! OUT
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, ice_surface_elevation)  ! OUT
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, glacier_fraction)       ! OUT
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, basal_melt_rate)        ! OUT
    dsl4jsb_Get_var2D_onDomain(LAND_ICE_, calving_flux)           ! OUT

    CALL cpl_get_field( &
      & routine, field_id_ice_thickness, 'lice_ice_thickness', grid%ntotal, &
      & ice_thickness)

    CALL cpl_get_field( &
      & routine, field_id_ice_elevation, 'lice_ice_elevation', grid%ntotal, &
      & ice_surface_elevation)

    CALL cpl_get_field( &
      & routine, field_id_glacier_fract, 'lice_glacier_fract', grid%ntotal, &
      & glacier_fraction)

    CALL cpl_get_field( &
      & routine, field_id_basal_melt, 'lice_basal_melt', grid%ntotal, &
      & basal_melt_rate)

    CALL cpl_get_field( &
      & routine, field_id_calving, 'lice_calving', grid%ntotal, &
      & calving_flux)

    !$ACC UPDATE DEVICE(ice_thickness, ice_surface_elevation, glacier_fraction, basal_melt_rate, calving_flux)

    IF ( write_coupler_restart ) &
      & CALL message(routine, 'YAC says it is put for restart')

  END SUBROUTINE interface_land_ice
#endif

END MODULE mo_interface_land_ice
