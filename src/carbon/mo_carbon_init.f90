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
MODULE mo_carbon_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_control,       ONLY: debug_on

  dsl4jsb_Use_processes CARBON_
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_memory(CARBON_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: carbon_init, carbon_read_cpools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_init'

CONTAINS

  SUBROUTINE carbon_init(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !dsl4jsb_Def_memory(CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    !dsl4jsb_Get_memory(CARBON_)

    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  carbon boundary conditions for tile '//TRIM(tile%name))

    CALL carbon_init_ic(tile)
    !CALL carbon_init_bc(tile)

  END SUBROUTINE carbon_init

  SUBROUTINE carbon_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !dsl4jsb_Def_memory(CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    !model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    !dsl4jsb_Get_memory(CARBON_)

    IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  carbon boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE carbon_init_bc

  SUBROUTINE carbon_init_ic(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_init_ic'

   !-----------------------------------------------------------------------
   ! CONTENT

   IF (.NOT. tile%Is_process_active(CARBON_)) RETURN

   IF (debug_on()) CALL message(TRIM(routine), 'Initializing carbon memory for tile '//TRIM(tile%name))

   model => get_model(tile%owner_model_id)

   dsl4jsb_Get_config(CARBON_)
   IF (dsl4jsb_Config(CARBON_)%read_cpools) CALL carbon_read_cpools(tile)

  END SUBROUTINE carbon_init_ic

  SUBROUTINE carbon_read_cpools(tile)

    USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_vgrid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)

    !INTEGER :: icarbon
    REAL(wp), POINTER ::          &
      & return_pointer  (:,  :) , &  !< temporary pointer just used to "call" the reading
      & return_pointer3d(:,:,:) , &
      & memory_pointer  (:,  :) , &  !< pointer the memory of the current variable of the current tile
      & memory_pointer3d(:,:,:)      !< pointer the memory of the current variable of the current tile

    !LOGICAL :: &
    !  read_cpools

    TYPE(t_input_file) :: input_file

    INTEGER, PARAMETER :: npools = 21
    INTEGER            :: i
    CHARACTER(len=40)  :: invar_name
    LOGICAL            :: invar_3D

    INTEGER            :: is, nsoil
    TYPE(t_jsb_vgrid), POINTER   :: soil_c

    CHARACTER(len=20), PARAMETER :: Cpool_name_flat(npools) = [character(len=20) ::      &
      & 'c_green',  'c_woods',  'c_reserve',                                             &
      & 'c_acid_ag1',  'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',               &
      & 'c_acid_bg1',  'c_water_bg1', 'c_ethanol_bg1', 'c_nonsoluble_bg1', 'c_humus_1',  &
      & 'c_acid_ag2',  'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2',               &
      & 'c_acid_bg2',  'c_water_bg2', 'c_ethanol_bg2', 'c_nonsoluble_bg2', 'c_humus_2' ]

    CHARACTER(len=20), PARAMETER :: Cpool_name_deep(npools) = [character(len=20) ::      &
      & 'c_green',  'c_woods',  'c_reserve',                                             &
      & 'c_acid_ag1',  'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',               &
      & 'c_acid_bg1_sl',  'c_water_bg1_sl', 'c_ethanol_bg1_sl', 'c_nonsoluble_bg1_sl', 'c_humus_1_sl',  &
      & 'c_acid_ag2',  'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2',               &
      & 'c_acid_bg2_sl',  'c_water_bg2_sl', 'c_ethanol_bg2_sl', 'c_nonsoluble_bg2_sl', 'c_humus_2_sl' ]

    CHARACTER(len=*), PARAMETER :: routine = modname//':carbon_read_cpools'

    model => get_model(tile%owner_model_id)

    ! Get carbon memory and config  of the tile
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_config(CARBON_)

    IF (.NOT. model%Is_process_enabled(CARBON_)) RETURN

    ! Initialize the C Pools
    invar_3D = .false.

    IF (dsl4jsb_Config(CARBON_)%on_soillev) THEN
      soil_c => Get_vgrid('soil_depth_carbon')
      nsoil = soil_c%n_levels
    END IF

    CALL message(TRIM(routine), 'Reading carbon pools for tile '//TRIM(tile%name)// &
      &                         ' from '//TRIM(dsl4jsb_Config(CARBON_)%ic_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(CARBON_)%ic_filename), model%grid_id)

    IF (.NOT. dsl4jsb_Config(CARBON_)%on_soillev) THEN
      DO i=1,npools
        invar_name = 'carbon_'//TRIM(Cpool_name_flat(i))//'_'//TRIM(tile%name)

        SELECT CASE(Cpool_name_flat(i))
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
    ELSE
      DO i=1,npools
        invar_name = 'carbon_'//TRIM(Cpool_name_deep(i))//'_'//TRIM(tile%name)
        invar_3D = .false.

        SELECT CASE(Cpool_name_deep(i))
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
        CASE ("c_acid_bg1_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_acid_bg1_sl%ptr
        CASE ("c_water_bg1_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_water_bg1_sl%ptr
        CASE ("c_ethanol_bg1_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_ethanol_bg1_sl%ptr
        CASE ("c_nonsoluble_bg1_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_nonsoluble_bg1_sl%ptr
        CASE ("c_humus_1_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_humus_1_sl%ptr
        CASE ("c_acid_ag2")
          memory_pointer => CARBON__mem%c_acid_ag2%ptr
        CASE ("c_water_ag2")
          memory_pointer => CARBON__mem%c_water_ag2%ptr
        CASE ("c_ethanol_ag2")
          memory_pointer => CARBON__mem%c_ethanol_ag2%ptr
        CASE ("c_nonsoluble_ag2")
          memory_pointer => CARBON__mem%c_nonsoluble_ag2%ptr
        CASE ("c_acid_bg2_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_acid_bg2_sl%ptr
        CASE ("c_water_bg2_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_water_bg2_sl%ptr
        CASE ("c_ethanol_bg2_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_ethanol_bg2_sl%ptr
        CASE ("c_nonsoluble_bg2_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_nonsoluble_bg2_sl%ptr
        CASE ("c_humus_2_sl")
          invar_3D = .true.
          memory_pointer3d => CARBON__mem%c_humus_2_sl%ptr
        END SELECT

        IF (invar_3D == .false.) THEN
          return_pointer => input_file%Read_2d(        &
            & variable_name = invar_name,              &
            & fill_array = memory_pointer )

          memory_pointer(:,:) = &
            & MERGE( memory_pointer(:,:), 0._wp,  memory_pointer(:,:) >= 0._wp)
        ELSE
          DO is=1,nsoil
            return_pointer3d => input_file%Read_2d_extdim( &
              & variable_name = invar_name,     &
              & start_extdim = is, end_extdim = is, extdim_name='soillev')
            memory_pointer3d(:,is,:) = MERGE(return_pointer3d(:,:,1), 0._wp, return_pointer3d(:,:,1) >= 0._wp)
          END DO
!           return_pointer3d => input_file%Read_2d_extdim(        &
!             & variable_name = invar_name,              &
!             & fill_array = memory_pointer3d )
!
!           memory_pointer3d(:,:,:) = &
!             & MERGE( memory_pointer3d(:,:,:), 0._wp,  memory_pointer3d(:,:,:) >= 0._wp)

        END IF

      END DO
    END IF

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(return_pointer3d)) CONTINUE
    NULLIFY(return_pointer3d)
    IF (ASSOCIATED(return_pointer)) CONTINUE
    NULLIFY(return_pointer)

  END SUBROUTINE carbon_read_cpools

#endif
END MODULE mo_carbon_init
