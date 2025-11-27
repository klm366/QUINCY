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

    ! Use declarations
    USE mo_sb_carbon_process, ONLY: convert_jsb_lct_to_quincy
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    !TYPE(t_jsb_model),       POINTER :: model
    !dsl4jsb_Def_memory(sb_CARBON_)
    TYPE(t_jsb_model),              POINTER :: model                !< the model
    TYPE(t_jsb_grid),               POINTER :: hgrid                !< Horizontal grid
    INTEGER                                 :: iblk, ic, is         !< loop dimensions
    INTEGER                                 :: nproma, nblks        !< dimensions
    INTEGER                                 :: JSBACH_LctNumber, QUINCY_LctNumber
    REAL(wp)                                :: k_root_dist
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_carbon_init'

    !-----------------------------------------------------------------------
    ! CONTENT
    dsl4jsb_Def_memory(SB_CARBON_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain    :: num_sl_above_bedrock
    dsl4jsb_Real2D_onDomain    :: root_depth
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    ! SB_CARBON_ 3D
    dsl4jsb_Real3D_onDomain    :: root_fraction_sl_fixed

    !model => get_model(tile%owner_model_id)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    hgrid  => Get_grid(model%grid_id)
    nblks  =  hgrid%nblks
    nproma =  hgrid%nproma



    ! Get sb_carbon memory of the tile
    !dsl4jsb_Get_memory(SB_CARBON_)
    dsl4jsb_Get_memory(SB_CARBON_)
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   num_sl_above_bedrock)        ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,   root_depth)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_,   soil_lay_depth_ubound_sl)
    ! SB_CARBON_ 3D
    dsl4jsb_Get_var3D_onDomain(SB_CARBON_,   root_fraction_sl_fixed)
    
    IF (.NOT. tile%Is_process_active(SB_CARBON_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting  sb_carbon boundary conditions for tile '//TRIM(tile%name))
    
    JSBACH_LctNumber               = dsl4jsb_Lctlib_param(LctNumber)
    CALL convert_jsb_lct_to_quincy(JSBACH_LctNumber, QUINCY_LctNumber, &
      &   lctlib_k_root_dist= k_root_dist) 
    
    !>     - fixed root fractions (no changes after init) depending on root_depth (SPQ variable)
    !>
    ! root fractions assuming an exponential distribution of roots \n
    !   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
    ! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
    !   thus for each layer
    !
    ! NOTE: upper bound of a layer is the larger value compared to lower bound !
    !
    DO iblk = 1,nblks
      DO ic = 1,nproma
        ! root depth larger than thickness of first soil layer
        IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, 1, iblk)) THEN
          root_fraction_sl_fixed(ic, 1, iblk) = (1._wp - EXP(-k_root_dist * soil_lay_depth_ubound_sl(ic, 1, iblk))) &
            &                             / (1._wp - EXP(-k_root_dist * root_depth(ic, iblk)))
          ! soil layers 2 to x
          !
          DO is = 2,INT(num_sl_above_bedrock(ic, iblk))
            IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, is-1, iblk)) THEN
              IF (root_depth(ic, iblk) > soil_lay_depth_ubound_sl(ic, is, iblk)) THEN
                root_fraction_sl_fixed(ic, is, iblk) = (EXP(-k_root_dist * soil_lay_depth_ubound_sl(ic, is-1, iblk)) &
                  &                                - EXP(-k_root_dist * soil_lay_depth_ubound_sl(ic, is, iblk))) &
                  &                                / (1._wp - EXP(-k_root_dist * root_depth(ic, iblk)))
              ELSE
                root_fraction_sl_fixed(ic, is, iblk) = (EXP(-k_root_dist * soil_lay_depth_ubound_sl(ic, is-1, iblk)) &
                  &                                - EXP(-k_root_dist * root_depth(ic, iblk))) &
                  &                                / (1._wp - EXP(-k_root_dist * root_depth(ic, iblk)))
              END IF
            ELSE
              root_fraction_sl_fixed(ic, is, iblk) = 0.0_wp
            END IF
          END DO
        ! root depth lower than / equal thickness of first soil layer
        ! all roots are in the first soil layer
        ELSE
          root_fraction_sl_fixed(ic, :, iblk) = 0.0_wp
          root_fraction_sl_fixed(ic, 1, iblk) = 1.0_wp
        END IF
      END DO
    END DO    

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
