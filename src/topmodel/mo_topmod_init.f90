!>
!! @brief Initialization of the the soil and snow memory.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_init
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

  dsl4jsb_Use_processes TOP_, SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_memory(TOP_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: topmodel_init !, init_soil_properties

  TYPE t_top_init_vars
    REAL(wp), POINTER ::                 &
      & cti_mean      (:,:)   => NULL(), &
      & cti_mu        (:,:)   => NULL(), &
      & cti_phi       (:,:)   => NULL(), &
      & cti_chi       (:,:)   => NULL(), &
      & cti_basin     (:,:)   => NULL(), &
      & cti_openwater (:,:)   => NULL(), &
      & frac_inun_clim(:,:,:) => NULL()
  END TYPE t_top_init_vars

  TYPE(t_top_init_vars) :: top_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_topmodel_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE topmodel_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL topmodel_read_init_vars(tile)
    END IF

    CALL topmodel_init_bc(tile)
    CALL topmodel_init_ic(tile)

    IF (tile%Is_last_process_tile(SSE_)) THEN
      CALL topmodel_finalize_init_vars()
    END IF

  END SUBROUTINE topmodel_init

  SUBROUTINE topmodel_finalize_init_vars

    DEALLOCATE(                          &
      & top_init_vars%cti_basin,         &
      & top_init_vars%cti_openwater,     &
      & top_init_vars%cti_mean,          &
      & top_init_vars%cti_mu,            &
      & top_init_vars%cti_phi,           &
      & top_init_vars%cti_chi            &
      & )

    IF (ASSOCIATED(top_init_vars%frac_inun_clim)) DEALLOCATE(top_init_vars%frac_inun_clim)

  END SUBROUTINE topmodel_finalize_init_vars

  SUBROUTINE topmodel_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get topmodel config
    dsl4jsb_Get_config(TOP_)

    ALLOCATE(                                           &
      & top_init_vars%cti_basin    (nproma, nblks    ), &
      & top_init_vars%cti_openwater(nproma, nblks    ), &
      & top_init_vars%cti_mean     (nproma, nblks    ), &
      & top_init_vars%cti_mu       (nproma, nblks    ), &
      & top_init_vars%cti_phi      (nproma, nblks    ), &
      & top_init_vars%cti_chi      (nproma, nblks    )  &
      & )
    IF (dsl4jsb_Config(TOP_)%inun_prescribed) ALLOCATE(top_init_vars%frac_inun_clim(nproma, nblks, 12))

    top_init_vars%cti_mean      (:,:) = missval
    top_init_vars%cti_mu        (:,:) = missval
    top_init_vars%cti_phi       (:,:) = missval
    top_init_vars%cti_chi       (:,:) = missval
    top_init_vars%cti_basin     (:,:) = missval
    top_init_vars%cti_openwater (:,:) = missval
    IF (dsl4jsb_Config(TOP_)%inun_prescribed) top_init_vars%frac_inun_clim(:,:,:) = missval

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_topmodel) &
      CALL message(TRIM(routine), 'Reading topmodel boundary conditions for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%ic_filename), model%grid_id)

    ptr_2D => input_file%Read_2d(variable_name='cti_mean'     , fill_array=top_init_vars%cti_mean)
    ptr_2D => input_file%Read_2d(variable_name='cti_mu'       , fill_array=top_init_vars%cti_mu)
    ptr_2D => input_file%Read_2d(variable_name='cti_phi'      , fill_array=top_init_vars%cti_phi)
    ptr_2D => input_file%Read_2d(variable_name='cti_chi'      , fill_array=top_init_vars%cti_chi)
    ptr_2D => input_file%Read_2d(variable_name='cti_basin'    , fill_array=top_init_vars%cti_basin)
    ptr_2D => input_file%Read_2d(variable_name='cti_openwater', fill_array=top_init_vars%cti_openwater)

    CALL input_file%Close()

    IF (dsl4jsb_Config(TOP_)%inun_prescribed) THEN
      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%bc_filename), model%grid_id)

      ptr_3D => input_file%Read_2d_time(variable_name='frac_inun', fill_array=top_init_vars%frac_inun_clim)

      CALL input_file%Close()
    ENDIF

  END SUBROUTINE topmodel_read_init_vars

  SUBROUTINE topmodel_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_topmodel) &
      CALL message(TRIM(routine), 'Setting top-model boundary conditions for tile '//TRIM(tile%name))

    IF (dsl4jsb_Config(TOP_)%inun_prescribed) THEN

      ! Get soil memory of the tile
      dsl4jsb_Get_memory(TOP_)

      !
      dsl4jsb_memory(TOP_)%frac_inun_clim(:,:,:) = top_init_vars%frac_inun_clim(:,:,:) ! StW: Should dims 2 and 3 be swapped?

      !$acc data present(dsl4jsb_var3D_onDomain(TOP_,frac_inun_clim))
      !$acc update device(dsl4jsb_var3D_onDomain(TOP_,frac_inun_clim))
      !$acc end data

    ENDIF

  END SUBROUTINE topmodel_init_bc

  SUBROUTINE topmodel_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(TOP_)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_topmodel) &
      CALL message(TRIM(routine), 'Initializing top-model memory for tile '//TRIM(tile%name))

    IF (.NOT. dsl4jsb_Config(TOP_)%inun_prescribed) THEN
      dsl4jsb_memory(TOP_)%cti_mean     (:,:) = top_init_vars%cti_mean(:,:)
      dsl4jsb_memory(TOP_)%cti_mu       (:,:) = top_init_vars%cti_mu(:,:)
      dsl4jsb_memory(TOP_)%cti_phi      (:,:) = top_init_vars%cti_phi(:,:)
      dsl4jsb_memory(TOP_)%cti_chi      (:,:) = top_init_vars%cti_chi(:,:)
      dsl4jsb_memory(TOP_)%cti_basin    (:,:) = top_init_vars%cti_basin(:,:)
      dsl4jsb_memory(TOP_)%cti_openwater(:,:) = top_init_vars%cti_openwater(:,:)

      !$acc data present(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_basin), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_openwater), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_mu), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_phi), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_chi))
      !$acc update device(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_basin), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_openwater), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_mu), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_phi), &
      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_chi))
      !$acc end data
    ENDIF

  END SUBROUTINE topmodel_init_ic

#endif
END MODULE mo_top_init
