!>
!! @brief Contains methods for JSBACH initializations.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

  USE mo_jsb_base,         ONLY: jsbach
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_model

  INTERFACE init_model
    MODULE PROCEDURE init_model_standalone
    MODULE PROCEDURE init_model_patch
    MODULE PROCEDURE init_model_decomposition
  END INTERFACE init_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS

  SUBROUTINE init_model_standalone(model_id, model_name)

    USE mo_jsb_grid_class,  ONLY: t_jsb_grid, new_grid
    USE mo_jsb_model_class, ONLY: register_grid, print_model

    INTEGER,                    INTENT(in) :: model_id
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name
!!$#ifdef HAVE_F2003
!!$    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self
!!$#else
!!$    TYPE(t_jsb_model), POINTER                :: self
!!$#endif

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_standalone'

    model => jsbach%models(model_id)%m
    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    ! TBD: loop through already existing grids and use it if already exists with same name
    grid => new_grid(model%config%grid_filename, model%config%grid_name)
    CALL register_grid(model, grid)

    CALL print_model(model)

    NULLIFY(model)

  END SUBROUTINE init_model_standalone

!!$  SUBROUTINE init_model_patch( model_id, patch, &
!!$                             & nsoil_w, dzsoil_w, nsoil_e, dzsoil_e)
  SUBROUTINE init_model_patch( model_id, patch)

    USE mo_jsb_domain_iface, ONLY: t_patch
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_model_class,  ONLY: register_grid

    INTEGER,                    INTENT(in) :: model_id
    TYPE(t_patch),              INTENT(in) :: patch

!!$    ! The soil layer structure has be sent back to ICON because this information is needed to setup the CDI output
!!$    INTEGER,           OPTIONAL, INTENT(out) :: &
!!$      & nsoil_w,                                & !< Number of soil layers (water)
!!$      & nsoil_e                                   !< Number of soil layers (energy)
!!$    REAL(wp), POINTER, OPTIONAL              :: &
!!$      & dzsoil_w(:),                            & !< Soil layer depths (water)
!!$      & dzsoil_e(:)                               !< Soil layer depths (energy)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_patch'

    IF (patch%id < 0) &
      CALL finish(TRIM(routine), 'Improper call with invalid patch')

    model => jsbach%models(model_id)%m

    grid => new_grid(patch=patch)
    CALL register_grid(model, grid)

    CALL init_model_common(model)

!!$    IF (PRESENT(nsoil_w)) nsoil_w = model%config%soil%nsoil_water
!!$    IF (PRESENT(nsoil_e)) nsoil_e = model%config%soil%nsoil_energy
!!$    IF (PRESENT(dzsoil_w)) THEN
!!$      ALLOCATE(dzsoil_w(nsoil_w))
!!$      dzsoil_w(1:nsoil_w) = model%config%soil%dz_water(1:nsoil_w)
!!$    END IF
!!$    IF (PRESENT(dzsoil_e)) THEN
!!$      ALLOCATE(dzsoil_e(nsoil_e))
!!$      dzsoil_w(1:nsoil_e) = model%config%soil%dz_energy(1:nsoil_e)
!!$    END IF

    NULLIFY(model)

  END SUBROUTINE init_model_patch

  SUBROUTINE init_model_decomposition( model_id, decomposition)

    USE mo_jsb_domain_iface, ONLY: pe_decomposed
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_model_class,  ONLY: register_grid

    INTEGER,               INTENT(in) :: model_id
    TYPE(pe_decomposed),   INTENT(in) :: decomposition

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_decomposition'

    IF (decomposition%nlon < 0) &
      CALL finish(TRIM(routine), 'Improper call with invalid decomposition')

    model => jsbach%models(model_id)%m

    grid => new_grid(decomp=decomposition)
    CALL register_grid(model, grid)

    CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE init_model_decomposition

  SUBROUTINE init_model_common(model)

    USE mo_jsb_model_class,    ONLY: max_no_of_var_lists, print_model
    USE mo_jsb_model_usecases, ONLY: init_usecase
    USE mo_srf_init,           ONLY: srf_init_ic, srf_init_bc
    USE mo_soil_init,          ONLY: soil_init
    USE mo_veg_init,           ONLY: veg_init_ic, veg_init_bc
    USE mo_jsb_base,           ONLY: get_tile

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_common'

    ALLOCATE(model%var_list_names(max_no_of_var_lists))
    model%var_list_names(:) = ''

    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
    ELSE
      ! Explicit model configuration through namelist
      ! TBD
    END IF

    CALL print_model(model)

    ! At this point, the model configuration, tile and memory structure is set up.

    ! We continue with initializing the different model components.

    IF (model%config%srf%active) THEN
      CALL srf_init_bc(model%land)
      CALL srf_init_ic(model%land)
    END IF

    IF (model%config%soil%active) THEN
      CALL soil_init(model)
    END IF

    ! Vegetation must be initialized after soil!
    IF (model%config%veg%active) THEN
      CALL veg_init_bc(model%land)
      CALL veg_init_ic(model%land)
    END IF

    ! Initialize tile fractions for all tiles
    ! This is done after the initialization of the initial and boundary conditions of the model components
    ! because e.g. the fraction of the vegetated tile depends on LAI.
    ! TBD: later read veg_ratio_max into 'veg' tile, and set 'bare' tile to 1-veg_ratio_max

  END SUBROUTINE init_model_common

#endif
END MODULE mo_jsb_model_init
