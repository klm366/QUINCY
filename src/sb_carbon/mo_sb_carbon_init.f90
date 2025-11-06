!> Initialization of the the carbon memory.
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
MODULE mo_sb_carbon_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_control,       ONLY: debug_on

  dsl4jsb_Use_processes SB_CARBON_
  dsl4jsb_Use_config(SB_CARBON_)
  dsl4jsb_Use_memory(SB_CARBON_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sb_carbon_init, sb_carbon_read_cpools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_carbon_init'

CONTAINS

  SUBROUTINE sb_carbon_init(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !dsl4jsb_Def_memory(sb_CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_carbon_init'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get sb_carbon memory of the tile
    !dsl4jsb_Get_memory(SB_CARBON_)

    IF (.NOT. tile%Is_process_active(SB_CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  sb_carbon boundary conditions for tile '//TRIM(tile%name))

    CALL sb_carbon_init_ic(tile)
    !CALL carbon_init_bc(tile)

  END SUBROUTINE sb_carbon_init

  SUBROUTINE sb_carbon_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !dsl4jsb_Def_memory(SB_CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_carbon_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    !dsl4jsb_Get_memory(SB_CARBON_)

    IF (.NOT. tile%Is_process_active(SB_CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  sb_carbon boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE sb_carbon_init_bc

  SUBROUTINE sb_carbon_init_ic(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(SB_CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_carbon_init_ic'

   !-----------------------------------------------------------------------
   ! CONTENT

   IF (.NOT. tile%Is_process_active(SB_CARBON_)) RETURN

   IF (debug_on()) CALL message(TRIM(routine), 'Initializing sb_carbon memory for tile '//TRIM(tile%name))

   model => get_model(tile%owner_model_id)

   dsl4jsb_Get_config(SB_CARBON_)
   IF (dsl4jsb_Config(SB_CARBON_)%read_sbcpools) CALL sb_carbon_read_cpools(tile)

  END SUBROUTINE sb_carbon_init_ic

  SUBROUTINE sb_carbon_read_cpools(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(SB_CARBON_)
    dsl4jsb_Def_memory(SB_CARBON_)

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
                                                                          ]

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_read_cpools'

    model => get_model(tile%owner_model_id)

    ! Get carbon memory and config  of the tile
    dsl4jsb_Get_memory(SB_CARBON_)
    dsl4jsb_Get_config(SB_CARBON_)

    IF (.NOT. model%Is_process_enabled(SB_CARBON_)) RETURN

    ! Initialize the C Pools


    CALL message(TRIM(routine), 'Reading sb_carbon pools for tile '//TRIM(tile%name)// &
      &                         ' from '//TRIM(dsl4jsb_Config(SB_CARBON_)%ic_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_CARBON_)%ic_filename), model%grid_id)

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

  END SUBROUTINE sb_carbon_read_cpools

#endif
END MODULE mo_sb_carbon_init
