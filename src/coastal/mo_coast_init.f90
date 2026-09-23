! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Initialisation of the coastal-constraint process
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
!>#### Reads the coastal boundary conditions
!>
!> The boundary-condition file is prepared offline on the ICON grid of the simulation.
!> It carries the long-term mean coastal sea surface salinity and the complex harmonic
!> constants of eight tidal constituents. The ICON grid is not checked when the file is
!> read: two R02B04 grids with the same number of cells would be accepted without an
!> error, so the global attribute uuidOfHGrid of the file has to match the grid of the
!> simulation.
!>
MODULE mo_coast_init
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_control,      ONLY: debug_on
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_grid_class,   ONLY: t_jsb_grid
  USE mo_jsb_grid,         ONLY: Get_grid
  USE mo_jsb_tile_class,   ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,    ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_coast_constants,  ONLY: mangrove_lct_id

  dsl4jsb_Use_processes COAST_
  dsl4jsb_Use_config(COAST_)
  dsl4jsb_Use_memory(COAST_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: coast_init


  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_init'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Initialise the coastal-constraint process on a tile
  !>
  !> Called from jsbach_init in mo_jsb_model_init for every tile on which COAST_ is
  !> active. The coastal process has boundary conditions only; there are no initial
  !> conditions to set. Because only the mangrove tile uses the boundary conditions, the
  !> file is read directly for that tile, without the intermediate init_vars structure
  !> that processes shared between many tiles use.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE coast_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CALL coast_init_bc(tile)

  END SUBROUTINE coast_init

  ! ======================================================================================================== !
  !>
  !>#### Read the coastal boundary conditions on the mangrove tile
  !>
  !> We here read the boundary conditions on the mangrove tile. The mask and the salinity
  !> are required whenever the process is active. The tidal constants are required only
  !> when the tidal-inundation factor is applied. Without them the tide height stays at
  !> zero; with the default ventilation weights, which add up to one, the aboveground roots
  !> then ventilate all inundated root and the inundation factor is one.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE coast_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_input_file)         :: input_file

    REAL(wp), ALLOCATABLE, TARGET :: buffer(:,:)
    REAL(wp), POINTER             :: ptr_2d(:,:)

    INTEGER :: lib_id, ivar
    LOGICAL :: have_tides

    CHARACTER(len=11), PARAMETER :: tide_names(16) = [character(len=11) :: &
      & 'tide_hre_m2', 'tide_him_m2', 'tide_hre_s2', 'tide_him_s2',            &
      & 'tide_hre_n2', 'tide_him_n2', 'tide_hre_k2', 'tide_him_k2',            &
      & 'tide_hre_k1', 'tide_him_k1', 'tide_hre_o1', 'tide_him_o1',            &
      & 'tide_hre_p1', 'tide_him_p1', 'tide_hre_q1', 'tide_him_q1' ]

    dsl4jsb_Def_config(COAST_)
    dsl4jsb_Def_memory(COAST_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':coast_init_bc'
    ! ------------------------------------------------------------------------------------------------------ !

    ! The boundary conditions are needed on the mangrove tile only
    IF (.NOT. ALLOCATED(tile%lcts)) RETURN
    lib_id = tile%lcts(1)%lib_id
    IF (lib_id /= mangrove_lct_id) RETURN

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)

    dsl4jsb_Get_config(COAST_)
    dsl4jsb_Get_memory(COAST_)

    CALL message(TRIM(routine), 'Reading coastal boundary conditions from ' &
      &                         //TRIM(dsl4jsb_Config(COAST_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(COAST_)%bc_filename), model%grid_id)

    ALLOCATE(buffer(grid%nproma, grid%nblks))

    ! Mask of the cells that carry mangrove. Outside the mask the update routine sets the
    ! control factors to one without evaluating the equations.
    IF (.NOT. input_file%Has_var('coast_mask')) &
      & CALL finish(TRIM(routine), 'Variable coast_mask not found in '//TRIM(dsl4jsb_Config(COAST_)%bc_filename))
    ptr_2d => input_file%Read_2d(variable_name='coast_mask', fill_array=buffer)
    dsl4jsb_var2D_onDomain(COAST_, coast_mask) = MERGE(1._wp, 0._wp, buffer > 0.5_wp)
    !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, coast_mask)) ASYNC(1)

    ! Salinity. Without it the salinity factor would be evaluated at zero salinity, which
    ! gives the minimum everywhere, so a missing variable is an error.
    IF (.NOT. input_file%Has_var('salinity')) &
      & CALL finish(TRIM(routine), 'Variable salinity not found in '//TRIM(dsl4jsb_Config(COAST_)%bc_filename))
    ptr_2d => input_file%Read_2d(variable_name='salinity', fill_array=buffer)
    dsl4jsb_var2D_onDomain(COAST_, salinity) = MERGE(buffer, 0._wp, buffer >= 0._wp)
    !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, salinity)) ASYNC(1)

    ! Tidal constants
    have_tides = .TRUE.
    DO ivar = 1, SIZE(tide_names)
      IF (.NOT. input_file%Has_var(TRIM(tide_names(ivar)))) have_tides = .FALSE.
    END DO

    IF (have_tides) THEN
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_m2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_m2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_m2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_m2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_m2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_m2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_s2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_s2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_s2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_s2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_s2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_s2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_n2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_n2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_n2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_n2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_n2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_n2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_k2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_k2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_k2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_k2', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_k2) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_k2)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_k1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_k1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_k1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_k1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_k1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_k1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_o1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_o1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_o1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_o1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_o1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_o1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_p1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_p1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_p1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_p1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_p1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_p1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_hre_q1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_hre_q1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_hre_q1)) ASYNC(1)
      ptr_2d => input_file%Read_2d(variable_name='tide_him_q1', fill_array=buffer)
      dsl4jsb_var2D_onDomain(COAST_, tide_him_q1) = buffer
      !$ACC UPDATE DEVICE(dsl4jsb_var2D_onDomain(COAST_, tide_him_q1)) ASYNC(1)
    ELSE IF (dsl4jsb_Config(COAST_)%l_apply_tidal_stress) THEN
      CALL finish(TRIM(routine), 'Tidal constants not found in '//TRIM(dsl4jsb_Config(COAST_)%bc_filename) &
        &                        //' but l_apply_tidal_stress is true')
    ELSE
      CALL message(TRIM(routine), 'No tidal constants in the boundary-condition file; tide height set to zero')
    END IF

    DEALLOCATE(buffer)
    CALL input_file%Close()

    IF (debug_on()) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE coast_init_bc

#endif
END MODULE mo_coast_init
