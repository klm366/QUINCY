!> Initialization of the surface energy balance memory.
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
!>#### Initialization of the surface energy balance process
!>
!> [[seb_init]] is called from [[mo_jsb_model_init:jsbach_init]] for each tile on which the SEB_
!> process is active, after the memory has been set up. It sets the boundary conditions (the skin
!> conductivity, if `l_skin_temp`) and the factors needed for the pseudo soil temperature. Initial
!> surface temperatures are not set here but in [[mo_sse_init:sse_init_ic]] together with the soil
!> temperatures.
!>
MODULE mo_seb_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: lp
  USE mo_exception,         ONLY: message, finish
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid

  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init

  !> Fields read from the bc file on the root tile and used to initialize all SEB_ tiles;
  !> deallocated after the last SEB_ tile has been initialized.
  TYPE t_seb_init_vars
    REAL(lp), ALLOCATABLE :: skin_conductivity(:,:) !< Skin conductivity [W m-2 K-1]
  END TYPE

  TYPE(t_seb_init_vars) :: seb_init_vars !< Fields for the initialization of all SEB_ tiles

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  !
  !> Initialize the SEB_ process on a tile (after memory has been set up)
  !>
  !> Reads the bc fields (on the root tile), sets the boundary conditions and computes the
  !> factors of the pseudo soil temperature (see [[mo_seb_land:calc_pseudo_soil_temp]]):
  !> \( F = \exp(-\Delta t / \tau) \) (`F_pseudo_soil_temp`) and \( N = 1 / (1 - F) \)
  !> (`N_pseudo_soil_temp`), with time step \(\Delta t\) and memory time \(\tau\) in days;
  !> \(\tau\) is the phenology parameter `tau_pseudo_soil` (10 days in the QUINCY standalone model).
  !
  SUBROUTINE seb_init(tile)
    USE mo_jsb_time,          ONLY: timestep_in_days
#ifdef __QUINCY_STANDALONE__
    USE mo_qs_forcing,          ONLY: forcing_options
    USE mo_jsb_math_constants,  ONLY: eps8
#else
    USE mo_pheno_parameters,  ONLY: pheno_param_jsbach
#endif
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Tile to be initialized

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init'

    dsl4jsb_Real2D_onDomain :: &
      & F_pseudo_soil_temp,    & !< Exponential factor for the pseudo soil temperature []
      & N_pseudo_soil_temp       !< Normalisation for the pseudo soil temperature []

    dsl4jsb_Def_memory(SEB_)

    model => get_model(tile%owner_model_id)

#ifdef __QUINCY_STANDALONE__
    CALL qs_seb_init_bc(tile)
#else
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL seb_read_init_vars(tile)
    END IF

    CALL seb_init_bc(tile)
    !CALL seb_init_ic(tile)

    IF (tile%Is_last_process_tile(SEB_)) THEN
      CALL seb_finalize_init_vars()
    END IF
#endif

    ! Initial values for the calculation of pseudo_soil_temp in SUBROUTINE calc_pseudo_soil_temp
    dsl4jsb_Get_memory(SEB_)


    IF (tile%contains_land .OR. model%config%use_tmx) THEN
      dsl4jsb_Get_var2D_onDomain(SEB_, N_pseudo_soil_temp) ! OUT
      dsl4jsb_Get_var2D_onDomain(SEB_, F_pseudo_soil_temp) ! OUT
#ifdef __QUINCY_STANDALONE__
      F_pseudo_soil_temp = EXP(- MAX(eps8, timestep_in_days(forcing_options%forcing_hour)) / 10._lp)
#else
      F_pseudo_soil_temp = EXP(- timestep_in_days(model_id=tile%owner_model_id)/pheno_param_jsbach%EG_SG%tau_pseudo_soil)
#endif
      N_pseudo_soil_temp = 1._lp / (1._lp - F_pseudo_soil_temp)

      !$ACC UPDATE DEVICE(F_pseudo_soil_temp, N_pseudo_soil_temp) ASYNC(1)
    END IF
  END SUBROUTINE seb_init

#ifdef __QUINCY_STANDALONE__
#else
  !> Read the bc fields of the SEB_ process into [[seb_init_vars]]
  !>
  !> Only `skin_conductivity` is read from `bc_filename`, if `l_skin_temp`; the model stops if
  !> it is missing in the file.
  !>
  SUBROUTINE seb_read_init_vars(tile)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Root tile

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER :: grid  !< Horizontal grid

    TYPE(t_input_file) :: input_file !< bc file
    INTEGER :: &
      & nproma, &                    !< Block length
      & nblks                        !< Number of blocks

    REAL(lp), POINTER :: ptr(:,:) !< Pointer to the field read from file

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_read_init_vars'

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      ALLOCATE(seb_init_vars%skin_conductivity(nproma, nblks))

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SEB_)%bc_filename), model%grid_id)

      IF (input_file%Has_var('skin_conductivity')) THEN
        ptr => input_file%Read_2d('skin_conductivity', fill_array=seb_init_vars%skin_conductivity)
      ELSE
        CALL finish(TRIM(routine), '*** Error: BC file does not contain skin_conductivity variable.')
      END IF

      CALL input_file%Close()
    END IF

  END SUBROUTINE seb_read_init_vars
#endif

  !> Deallocate the fields in [[seb_init_vars]]
  SUBROUTINE seb_finalize_init_vars
    IF (ALLOCATED(seb_init_vars%skin_conductivity)) DEALLOCATE(seb_init_vars%skin_conductivity)
  END SUBROUTINE seb_finalize_init_vars

  ! ======================================================================================================= !
  !>
  !> SEB_ init bc - quincy standalone
  !>
  !> Sets a constant skin conductivity of 40 W m-2 K-1 (estimated from the ICON-Land input
  !> file), if `l_skin_temp`.
  !>
  SUBROUTINE qs_seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Tile to be initialized

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':qs_seb_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      IF (tile%is_land .OR. model%config%use_tmx) THEN
        dsl4jsb_var2D_onDomain(SEB_, skin_conductivity) = 40.0_lp   ! estimated from icon-land input file
      END IF
    END IF

  END SUBROUTINE qs_seb_init_bc

  ! ======================================================================================================= !
  !>
  !> SEB_ init bc
  !>
  !> Sets the skin conductivity from [[seb_init_vars]] on land tiles (all tiles with
  !> `use_tmx`), if `l_skin_temp`.
  !>
  SUBROUTINE seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Tile to be initialized

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
      IF (tile%is_land .OR. model%config%use_tmx) THEN
        dsl4jsb_var2D_onDomain(SEB_, skin_conductivity) = seb_init_vars%skin_conductivity
        !$ACC UPDATE ASYNC(1) DEVICE(dsl4jsb_var2D_onDomain(SEB_, skin_conductivity))
      END IF
    END IF

  END SUBROUTINE seb_init_bc

  !> SEB_ init ic - currently empty and not called
  !>
  !> The initial surface temperature is set in [[mo_sse_init:sse_init_ic]] together with the soil
  !> temperatures.
  !>
  SUBROUTINE seb_init_ic(tile)

    !USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< Tile to be initialized

    ! TYPE(t_jsb_model),       POINTER :: model

    ! dsl4jsb_Def_config(SEB_)
    !dsl4jsb_Def_memory(SEB_)

    !REAL(lp), POINTER :: &
    !  & ptr_2D(:,  :)       !< temporary pointer

    !TYPE(t_input_file) :: input_file

    !INTEGER  :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_ic'

    !model => get_model(tile%owner_model_id)

    ! Get seb config
    ! dsl4jsb_Get_config(SEB_)

    ! IF (.NOT. model%Is_process_enabled(SEB_)) RETURN

    ! Get seb memory of the tile
    !dsl4jsb_Get_memory(SEB_)

    ! Initialize parameters
    !
    ! IF (debug_on()) CALL message(TRIM(routine), 'Setting initial state of seb memory for tile '// &
      ! &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SEB_)%ic_filename))

    ! Initial surface temperature
    ! Is initialized in sse_init_ic together with soil temperatures

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
