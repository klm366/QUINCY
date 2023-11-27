!>
!! @brief Initialization of the methane production and transport module.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-17
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_methane_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, logical2string

  dsl4jsb_Use_processes HYDRO_, CH4_
  dsl4jsb_Use_config(CH4_)
  dsl4jsb_Use_memory(CH4_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: methane_init

  TYPE t_ch4_init_vars
!    REAL(wp), POINTER ::                 &
  END TYPE t_ch4_init_vars

  TYPE(t_ch4_init_vars) :: ch4_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_init'

CONTAINS

  !
  !> Intialize methane production and transport process (after memory has been set up)
  !
  SUBROUTINE methane_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL methane_read_init_vars(tile)
    END IF

    CALL methane_init_bc(tile)
    CALL methane_init_ic(tile)

    IF (tile%Is_last_process_tile(CH4_)) THEN
      CALL methane_finalize_init_vars()
    END IF

  END SUBROUTINE methane_init

  SUBROUTINE methane_finalize_init_vars

!    DEALLOCATE(                          &
!      & top_init_vars%cti_basin,         &
!      & )

!    IF (ASSOCIATED(ch4_init_vars%frac_inun_clim)) DEALLOCATE(ch4_init_vars%frac_inun_clim)

  END SUBROUTINE methane_finalize_init_vars

  SUBROUTINE methane_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(CH4_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get methane config
    dsl4jsb_Get_config(CH4_)

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Reading methane boundary conditions for tile '//TRIM(tile%name))

  END SUBROUTINE methane_read_init_vars

  SUBROUTINE methane_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get methane config
    dsl4jsb_Get_config(CH4_)

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Setting methane boundary conditions for tile '//TRIM(tile%name))

  END SUBROUTINE methane_init_bc

  SUBROUTINE methane_init_ic(tile)

    USE mo_methane_util, ONLY: find_bedrock_layer, find_layer
    USE mo_methane_process, ONLY: determine_root_dist
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: layer

    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)
    dsl4jsb_Def_memory(HYDRO_)

    INTEGER :: nblks, i, j, nproma

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_init_ic'

    dsl4jsb_Real2D_onDomain :: soil_depth
    dsl4jsb_Real2D_onDomain :: mxlay_above_rock
    dsl4jsb_Real3D_onDomain :: root_dist_sl

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    model => get_model(tile%owner_model_id)
    layer => Get_vgrid('soil_depth_water')
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(CH4_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get soil config
    dsl4jsb_Get_config(CH4_)

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Initializing methane memory for tile '//TRIM(tile%name))

    ! Determine dz for soil layer midpoints
    dsl4jsb_memory(CH4_)%layer_midpoint_dz(1) = layer%dz(1) * 0.5_wp
    DO i=2,layer%n_levels
      dsl4jsb_memory(CH4_)%layer_midpoint_dz(i) = (layer%dz(i) + layer%dz(i-1)) * 0.5_wp
    ENDDO
    ! TODO: StW: some $acc stuff missing here

    ! Determine (constant) layer for in which bedrock is located 
    dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)
    dsl4jsb_Get_var2D_onDomain(CH4_, mxlay_above_rock)
    dsl4jsb_Get_var3D_onDomain(CH4_, root_dist_sl)
    DO j=1,nblks
      DO i=1,nproma
        mxlay_above_rock(i,j) = REAL(find_bedrock_layer(layer%n_levels,layer%ubounds(:),soil_depth(i,j)),wp)
      ENDDO
    ENDDO
    ! TODO: StW: some $acc stuff missing here

    ! Determine distribution of roots in soil layers.
    ! NOTE: In case of dynamical roots, this should be moved into the time loop in update_methane 
    !       and be changed from working on domain to chunk
    IF (dsl4jsb_Lctlib_param(NaturalVegFlag) .OR. &
        dsl4jsb_Lctlib_param(PastureFlag)    .OR. &
        dsl4jsb_Lctlib_param(CropFlag)) THEN
      DO j=1,nblks
        DO i=1,nproma
           CALL determine_root_dist(INT(mxlay_above_rock(i,j)+0.5_wp),layer%dz(:),layer%ubounds(:), &
             &                      dsl4jsb_Lctlib_param(root_dist_beta),                           &
             &                      root_dist_sl(i,:,j))
        ENDDO
      ENDDO
    ELSE
      root_dist_sl(:,:,:) = 0._wp
    ENDIF
    ! TODO: StW: some $acc stuff missing here

!      !$acc data present(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
!      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_basin), &
!      !$acc update device(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
!      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_basin), &
!      !$acc end data

  END SUBROUTINE methane_init_ic

#endif

END MODULE mo_methane_init
