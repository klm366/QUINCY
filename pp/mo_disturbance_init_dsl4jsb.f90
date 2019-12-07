!
! This file was automatically generated from:
!    /users/rschnur/icon-les-dev-staging.pgi_gpu/externals/jsbach/src/disturbance/mo_disturbance_init.f90
!    with md5 hash: cc5877191d91b5b2d15aaa399c3d7d40
!
!>
!! @brief Initialization of the the disturbance memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-06-20)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_init
#ifndef __NO_JSBACH__

  !USE mo_kind,              ONLY: wp
  !USE mo_exception,         ONLY: message

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  !USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input
  !USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  !USE mo_jsb_process_class, ONLY: disturb_
  !USE mo_disturb_config_class, ONLY: t_disturb_config
  !USE mo_disturb_memory_class, ONLY: t_disturb_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: disturb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_init'

CONTAINS

  !
  !> Intialize disturb process (after memory has been set up)
  !
  SUBROUTINE disturb_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init'

    IF (tile%level > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    !TYPE(t_disturb_memory), POINTER :: disturb__mem
    !SELECT TYPE (mm => tile%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT

    !CALL disturb_init_ic(tile)

    !CALL disturb_init_bc(tile)

  END SUBROUTINE disturb_init
!
!   SUBROUTINE disturb_init_bc(tile)
!     !-----------------------------------------------------------------------
!     !  ARGUMENTS
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!
!     !-----------------------------------------------------------------------
!     !  LOCAL VARIABLES
!     TYPE(t_jsb_model),       POINTER :: model
!     TYPE(t_disturb_memory), POINTER :: disturb__mem
!
!     CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init_bc'
!
!     !-----------------------------------------------------------------------
!     ! CONTENT
!
!     model => get_model(tile%owner_model_id)
!
!     ! Get disturb memory of the tile
!     SELECT TYPE (mm => tile%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT
!
!     IF (.NOT. ASSOCIATED(disturb__mem)) RETURN
!
!     CALL message(TRIM(routine), 'Setting  disturb boundary conditions for tile '//TRIM(tile%name))
!
!     ! R: nothing to do in the moment ...
!
!   END SUBROUTINE disturb_init_bc
!
!   RECURSIVE SUBROUTINE disturb_init_ic(tile)
!
!     !USE mo_disturb_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial
!
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!
!     TYPE(t_jsb_model),         POINTER :: model
!     TYPE(t_disturb_config), POINTER :: disturb__conf
!     TYPE(t_disturb_memory), POINTER :: disturb__mem
!
!     !INTEGER :: idisturb
!     REAL(wp), POINTER :: &
!       & return_pointer  (:,  :) !, & !< temporary pointer
!       ! & return_pointer3d(:,:,:)
!     !LOGICAL :: &
!     !  read_cpools
!
!     TYPE(t_input_file) :: input_file
!     !INTEGER :: i, n1, n2
!
!     CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init_ic'
!
!     model => get_model(tile%owner_model_id)
!
!     ! Get disturb memory and config  of the tile
!     SELECT TYPE (mm => tile%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT
!     SELECT TYPE (cc => model%processes(disturb_)%p%config); TYPE IS (t_disturb_config) ; disturb__conf => cc ; END SELECT
!
!     IF (.NOT. disturb__conf%active) RETURN
!
!     CALL message(TRIM(routine), 'Reading disturb pools in memory for tile '//TRIM(tile%name))
!
!     ! Initialize physical state
!     !
!
!     input_file = jsb_netcdf_open_input(TRIM(disturb__conf%ic_filename), model%grid_id)
!
!     ! TBD: adapt names of variables in input files to names in memory
!
!     ! R: These variables should be 3D as they are tile dependent. However for THIS tile they would be 2D but in the input file
!     !    they are 3D!
!     ! R: jedenfalls wurden 3D variablen bisher nirgendwo über einlesen vom File initialisiert!
!     return_pointer => input_file%Read_3d(         &
!       & variable_name='c_green',           &
!       & fill_array = disturb__mem%c_green%ptr)
!     disturb__mem%c_green%ptr(:,:,:) = &
!       MERGE(disturb__mem%c_green%ptr(:,:,:), 0._wp, disturb__mem%c_green%ptr(:,:,:) >= 0._wp)
!
!
!     CALL input_file%Close()
!
!   END SUBROUTINE disturb_init_ic
!
#endif
END MODULE mo_disturb_init
