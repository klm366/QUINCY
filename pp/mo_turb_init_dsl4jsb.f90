!
! This file was automatically generated from:
!    /users/rschnur/icon-les-dev-staging.pgi_gpu/externals/jsbach/src/turbulence/mo_turb_init.f90
!    with md5 hash: 6fa58f7110a6e4627b4cf70d250ea791
!
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
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_process_class, ONLY: turb_
  USE mo_turb_config_class, ONLY: t_turb_config
  USE mo_turb_memory_class, ONLY: t_turb_memory

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

    CALL turb_init_bc(tile)
    CALL turb_init_ic(tile)

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
    TYPE(t_turb_config), POINTER :: turb__conf
    TYPE(t_turb_memory), POINTER :: turb__mem

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get turb config
    SELECT TYPE (cc => model%processes(turb_)%p%config); TYPE IS (t_turb_config) ; turb__conf => cc ; END SELECT

    IF (.NOT. turb__conf%active) RETURN

    ! Get turb memory of the tile
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT

    ! Initialize parameters
    !
    IF (debug_on()) CALL message(TRIM(routine), &
         'Initializing turb memory for tile '//TRIM(tile%name)//' from '//TRIM(turb__conf%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(turb__conf%bc_filename), model%grid_id)

    ! Surface roughness length
    ptr_2D => input_file%Read_2d(         &
      & variable_name='roughness_length',         &
      & fill_array = turb__mem%rough_m%ptr)
    turb__mem%rough_m%ptr(:,:) = &
      MERGE(turb__mem%rough_m%ptr(:,:), 1._wp, turb__mem%rough_m%ptr(:,:) > 0._wp)

    !$acc update device(turb__mem%rough_m%ptr(:,:))

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(ptr_2D)) CONTINUE
    NULLIFY(ptr_2D)

  END SUBROUTINE turb_init_ic

#endif
END MODULE mo_turb_init
