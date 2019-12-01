!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/carbon/mo_carbon_init.f90
!    with md5 hash: 59a7fc85cb0a09da42206666e3ccd422
!
!>
!! @brief Initialization of the the carbon memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-12)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_carbon_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_time_dt
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_process_class, ONLY: carbon_
  USE mo_carbon_config_class, ONLY: t_carbon_config
  USE mo_carbon_memory_class, ONLY: t_carbon_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: carbon_init, carbon_read_cpools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_init'

CONTAINS

  !
  !> Intialize carbon process (after memory has been set up)
  !
  SUBROUTINE carbon_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    !TYPE(t_jsb_model), POINTER                :: model
    REAL(wp)                                  :: dtime

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init'

    TYPE(t_carbon_memory), POINTER :: carbon__mem
    !TYPE(t_carbon_config), POINTER :: carbon__conf

    REAL(wp), POINTER, DIMENSION(:,:) :: f_pseudo_temp
    REAL(wp), POINTER, DIMENSION(:,:) :: n_pseudo_temp
    REAL(wp), POINTER, DIMENSION(:,:) :: f_pseudo_precip
    REAL(wp), POINTER, DIMENSION(:,:) :: n_pseudo_precip

    !model => get_model(tile%owner_model_id)

    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
    !SELECT TYPE (cc => model%processes(carbon_)%p%config); TYPE IS (t_carbon_config) ; carbon__conf => cc ; END SELECT

    dtime = get_time_dt()

    !CALL carbon_init_ic(tile)
    !CALL carbon_init_bc(tile)


    ! R: This should stand here because these variables are more technical constants (only depending on dtime)
    !    than a part of the carbon process parametrization. Furthermore they should be calculated only once when a simulation
    !    is started or restarted and not for each tile and for each time step again!
    !    It is not necessary to write and read them to/from the restart files.
    !    It is not necessary to write them to the output file.
    !    However, we could think about to kick them also out of the memory...
    !    AND we should make them 1D, later, as soon as we do not want them in the output for debugging!

     !f_pseudo_temp => carbon__mem%f_pseudo_temp%ptr(ics:ice,iblk)    ! for this we would need iblk, but we want the whole field!
     F_pseudo_temp   => carbon__mem%F_pseudo_temp%ptr(:,:)   ! out
     N_pseudo_temp   => carbon__mem%N_pseudo_temp%ptr(:,:)   ! out
     F_pseudo_precip => carbon__mem%F_pseudo_precip%ptr(:,:) ! out
     N_pseudo_precip => carbon__mem%N_pseudo_precip%ptr(:,:) ! out

    ! R: Note, I changed the "writing style" as compared to JSBACH3.
    ! F_pseudo_temp and N_pseudo_temp are constants used for updating the pseudo-15-day mean temperature.
    ! These constants define the contributions of t_air and pseudo_temp to the new pseudo_temp.
    F_pseudo_temp = EXP(-dtime / (15._wp * 86400._wp))   ! F_pseudo_temp decreases expotentially from 1 with higher
                                                         ! (dtime/15DaysInSeconds). For (dtime/15DaysInSeconds)
                                                         !  =1 it would be 0.368 and =0.5 it would be 0.607.
    N_pseudo_temp = 1._wp - F_pseudo_temp                ! Rises with the same amout as F_pseudo_temp falls.

    F_pseudo_precip = EXP(-dtime / (15._wp * 86400._wp))
    N_pseudo_precip = 1._wp - F_pseudo_precip

  END SUBROUTINE carbon_init

  SUBROUTINE carbon_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !TYPE(t_carbon_memory), POINTER :: carbon__mem

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    !SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  carbon boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE carbon_init_bc

  SUBROUTINE carbon_init_ic(tile)

    !USE mo_carbon_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),         POINTER :: model
    !TYPE(t_carbon_config), POINTER :: carbon__conf
    !TYPE(t_carbon_memory), POINTER :: carbon__mem

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init_ic'

   !-----------------------------------------------------------------------
   ! CONTENT

   IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

   IF (debug_on()) CALL message(TRIM(routine), 'Initializing carbon memory for tile '//TRIM(tile%name))

   !model => Get_model(tile%owner_model_id)

   ! Get carbon memory of the tile
   !SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT

   ! R: nothing to do in the moment ...

  END SUBROUTINE carbon_init_ic

  SUBROUTINE carbon_read_cpools(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    TYPE(t_carbon_config), POINTER :: carbon__conf
    TYPE(t_carbon_memory), POINTER :: carbon__mem

    !INTEGER :: icarbon
    REAL(wp), POINTER ::          &
      & return_pointer  (:,  :) , &  !< temporary pointer just used to "call" the reading
      ! & return_pointer3d(:,:,:)
      & memory_pointer  (:,  :)      !< pointer the memory of the current variable of the current tile
    !LOGICAL :: &
    !  read_cpools

    TYPE(t_input_file) :: input_file

    INTEGER, PARAMETER :: npools = 21
    INTEGER            :: i
    CHARACTER(len=40)            :: invar_name
    CHARACTER(len=20), PARAMETER :: Cpool_name(npools) = [character(len=20) ::           &
      & 'c_green',  'c_woods',  'c_reserve',                                             &
      & 'c_acid_ag1',  'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',               &
      & 'c_acid_bg1',  'c_water_bg1', 'c_ethanol_bg1', 'c_nonsoluble_bg1', 'c_humus_1',  &
      & 'c_acid_ag2',  'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2',               &
      & 'c_acid_bg2',  'c_water_bg2', 'c_ethanol_bg2', 'c_nonsoluble_bg2', 'c_humus_2' ]

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_read_cpools'

    model => get_model(tile%owner_model_id)

    ! Get carbon memory and config  of the tile
    SELECT TYPE (mm => tile%mem(carbon_)%p); TYPE IS (t_carbon_memory); carbon__mem => mm ; END SELECT
    SELECT TYPE (cc => model%processes(carbon_)%p%config); TYPE IS (t_carbon_config) ; carbon__conf => cc ; END SELECT

    IF (.NOT. carbon__conf%active) RETURN

    ! Initialize the C Pools


    CALL message(TRIM(routine), 'Reading carbon pools for tile '//TRIM(tile%name)// &
      &                         ' from '//TRIM(carbon__conf%ic_filename))

    input_file = jsb_netcdf_open_input(TRIM(carbon__conf%ic_filename), model%grid_id)

    ! TBD: adapt names of variables in input files to names in memory

    ! Note, this is done to avoid 6 Lines of source code for each of the 21 C-pools and
    !       eventually in addition for each of the existing tile (6x21x14=5166 lines).

    ! If we do not want to read carbon pools for each tile we would have to add this
    ! case construct to choose the behaviour for each tile:
    !SELECT CASE(tile%name)
    !CASE (box)
    !CASE (land)
    !CASE (veg)
    !CASE (pft01)
    ! ...

    DO i=1,npools
      invar_name = 'carbon_'//TRIM(Cpool_name(i))//'_'//TRIM(tile%name)

      SELECT CASE(Cpool_name(i))
      CASE ("c_green")
        memory_pointer => CARBON__mem%c_green%ptr ! CARBON__mem= tile%mem(carbon_)%p
      CASE ("c_woods")
        memory_pointer => CARBON__mem%c_woods%ptr
      CASE ("c_reserve")
        memory_pointer => CARBON__mem%c_reserve%ptr
      CASE ("c_acid_ag1")
        memory_pointer => CARBON__mem%c_acid_ag1%ptr
      CASE ("c_water_ag1")
        memory_pointer => CARBON__mem%c_water_ag1%ptr
      CASE ("c_ethanol_ag1")
        memory_pointer => CARBON__mem%c_ethanol_ag1%ptr
      CASE ("c_nonsoluble_ag1")
        memory_pointer => CARBON__mem%c_nonsoluble_ag1%ptr
      CASE ("c_acid_bg1")
        memory_pointer => CARBON__mem%c_acid_bg1%ptr
      CASE ("c_water_bg1")
        memory_pointer => CARBON__mem%c_water_bg1%ptr
      CASE ("c_ethanol_bg1")
        memory_pointer => CARBON__mem%c_ethanol_bg1%ptr
      CASE ("c_nonsoluble_bg1")
        memory_pointer => CARBON__mem%c_nonsoluble_bg1%ptr
      CASE ("c_humus_1")
        memory_pointer => CARBON__mem%c_humus_1%ptr
      CASE ("c_acid_ag2")
        memory_pointer => CARBON__mem%c_acid_ag2%ptr
      CASE ("c_water_ag2")
        memory_pointer => CARBON__mem%c_water_ag2%ptr
      CASE ("c_ethanol_ag2")
        memory_pointer => CARBON__mem%c_ethanol_ag2%ptr
      CASE ("c_nonsoluble_ag2")
        memory_pointer => CARBON__mem%c_nonsoluble_ag2%ptr
      CASE ("c_acid_bg2")
        memory_pointer => CARBON__mem%c_acid_bg2%ptr
      CASE ("c_water_bg2")
        memory_pointer => CARBON__mem%c_water_bg2%ptr
      CASE ("c_ethanol_bg2")
        memory_pointer => CARBON__mem%c_ethanol_bg2%ptr
      CASE ("c_nonsoluble_bg2")
        memory_pointer => CARBON__mem%c_nonsoluble_bg2%ptr
      CASE ("c_humus_2")
        memory_pointer => CARBON__mem%c_humus_2%ptr
      END SELECT

      return_pointer => input_file%Read_2d(        &
        & variable_name = invar_name,              &
        & fill_array = memory_pointer )

      memory_pointer(:,:) = &
        & MERGE( memory_pointer(:,:), 0._wp,  memory_pointer(:,:) >= 0._wp)

    END DO

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(return_pointer)) CONTINUE
    NULLIFY(return_pointer)

  END SUBROUTINE carbon_read_cpools

#endif
END MODULE mo_carbon_init
