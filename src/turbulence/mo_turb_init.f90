!>
!! @brief Initialization of the the turbulence memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_turb_init
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_process_class, ONLY: ON_TILE_, AGGREGATE_
  USE mo_jsb_task_class,    ONLY: t_jsb_task_options

  dsl4jsb_Use_processes TURB_
  dsl4jsb_Use_config(TURB_)
  dsl4jsb_Use_memory(TURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: turb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_init'

CONTAINS

  !
  !> Intialize turbulence process (after memory has been set up)
  !
  SUBROUTINE turb_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init'

    IF (tile%process_action(TURB_) == ON_TILE_) THEN
      ! CALL turb_init_bc(tile)
      CALL turb_init_ic(tile)
    ELSE IF (tile%process_action(TURB_) == AGGREGATE_) THEN
      ! CALL turb_aggregate_bc(tile)
      CALL turb_aggregate_ic(tile)
    END IF

  END SUBROUTINE turb_init

  SUBROUTINE turb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init_bc'

    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings
    !model => get_model(tile%owner_model_id)

  END SUBROUTINE turb_init_bc

  RECURSIVE SUBROUTINE turb_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(TURB_)

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get turb config
    dsl4jsb_Get_config(TURB_)

    ! Get turb memory of the tile
    dsl4jsb_Get_memory(TURB_)

    ! Initialize parameters
    !
    IF (debug_on()) CALL message(TRIM(routine), &
         'Initializing turb memory for tile '//TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(TURB_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TURB_)%bc_filename), model%grid_id)

    ! Surface roughness length
    ptr_2D => input_file%Read_2d(         &
      & variable_name='roughness_length',         &
      & fill_array = dsl4jsb_var_ptr(TURB_,rough_m))
    dsl4jsb_var2D_onDomain(TURB_,rough_m) = &
      MERGE(dsl4jsb_var2D_onDomain(TURB_,rough_m), 1._wp, dsl4jsb_var2D_onDomain(TURB_,rough_m) > 0._wp)

    !$acc update device(dsl4jsb_var2D_onDomain(TURB_,rough_m))

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(ptr_2D)) CONTINUE
    NULLIFY(ptr_2D)

  END SUBROUTINE turb_init_ic

  SUBROUTINE turb_aggregate_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_memory(TURB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_jsb_grid),        POINTER :: grid
    TYPE(t_jsb_task_options)         :: options

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_aggregate_ic'

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    options%ics = 1
    options%ice = grid%nproma

    ics = options%ics
    ice = options%ice
    nc = ice - ics + 1

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(TURB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    DO iblk = 1, grid%nblks

      CALL model%Set_options(iblk=iblk)
      options%iblk = iblk
      IF (iblk == grid%nblks) THEN
        options%ice = grid%npromz
        ice = grid%npromz
        ! TODO: ics and ice should be retrieved for each iblk (are not necessarily constant for each block in ICON)
      END IF

      dsl4jsb_Aggregate_onChunk(TURB_, rough_m,    weighted_by_fract)

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE turb_aggregate_ic

#endif
END MODULE mo_turb_init
