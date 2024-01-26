!>
!! @brief Initialization of the megan scheme.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-24
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_megan_init
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
  USE mo_util,               ONLY: i2s => int2string, l2s => logical2string

  dsl4jsb_Use_processes MEGAN_
  dsl4jsb_Use_config(MEGAN_)
  dsl4jsb_Use_memory(MEGAN_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: megan_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_megan_init'

CONTAINS

  !
  !> Intialize topmodel wetland process (after memory has been set up)
  !
  SUBROUTINE megan_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL megan_read_init_vars(tile)
    END IF

    CALL megan_init_bc(tile)
    CALL megan_init_ic(tile)

    IF (tile%Is_last_process_tile(MEGAN_)) THEN
      CALL megan_finalize_init_vars()
    END IF

  END SUBROUTINE megan_init

  SUBROUTINE megan_finalize_init_vars

  END SUBROUTINE megan_finalize_init_vars

  SUBROUTINE megan_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(MEGAN_)

!    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

!    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get megan config
    dsl4jsb_Get_config(MEGAN_)

    IF (.NOT. tile%Is_process_active(MEGAN_)) RETURN

!    IF (debug_on() .OR. dsl4jsb_Config(MEGAN_)%dbg_process) &
!      CALL message(TRIM(routine), 'Reading topmodel boundary conditions for tile '//TRIM(tile%name))

!   input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%ic_filename), model%grid_id)

!    ptr_2D => input_file%Read_2d(variable_name='cti_openwater', fill_array=top_init_vars%cti_openwater)

!    CALL input_file%Close()

!    ENDIF

  END SUBROUTINE megan_read_init_vars

  SUBROUTINE megan_init_bc(tile)
    USE mo_megan_memory_class, ONLY: t_megan_specie

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(MEGAN_)
    dsl4jsb_Def_memory(MEGAN_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(MEGAN_)
    dsl4jsb_Get_memory(MEGAN_)

    IF (.NOT. tile%Is_process_active(MEGAN_)) RETURN

    ! Add constant speice properties to memory structure
    ! --- isoprene ---
    dsl4jsb_memory(MEGAN_)%iso%longname  = 'isoprene'
    dsl4jsb_memory(MEGAN_)%iso%shortname = 'C5H8'
    dsl4jsb_memory(MEGAN_)%iso%class    = 1
    dsl4jsb_memory(MEGAN_)%iso%mweight  = 68.12_wp
    dsl4jsb_memory(MEGAN_)%iso%cnumber  = 5
    dsl4jsb_memory(MEGAN_)%iso%cafac    = 2
    dsl4jsb_memory(MEGAN_)%iso%beta     = 0._wp
    dsl4jsb_memory(MEGAN_)%iso%ldf      = 1._wp
    dsl4jsb_memory(MEGAN_)%iso%scale    = dsl4jsb_Config(MEGAN_)%scale_iso
    dsl4jsb_memory(MEGAN_)%iso%ef       = 0._wp ! TODO: lctlib-something
    
    ! --- carbon monoxide ---
    dsl4jsb_memory(MEGAN_)%iso%longname  = 'carbon monoxide'
    dsl4jsb_memory(MEGAN_)%iso%shortname = 'CO'
    dsl4jsb_memory(MEGAN_)%iso%class    = 16
    dsl4jsb_memory(MEGAN_)%iso%mweight  = 28.01_wp
    dsl4jsb_memory(MEGAN_)%iso%cnumber  = 1
    dsl4jsb_memory(MEGAN_)%iso%cafac    = 1
    dsl4jsb_memory(MEGAN_)%iso%beta     = dsl4jsb_Config(MEGAN_)%beta_CO
    dsl4jsb_memory(MEGAN_)%iso%ldf      = 1._wp
    dsl4jsb_memory(MEGAN_)%iso%scale    = dsl4jsb_Config(MEGAN_)%scale_CO
    dsl4jsb_memory(MEGAN_)%iso%ef       = 0._wp ! TODO: lctlib-something
    
    ! --- nitricoxide ---
    dsl4jsb_memory(MEGAN_)%iso%longname  = 'nitricoxide'
    dsl4jsb_memory(MEGAN_)%iso%shortname = 'NO'
    dsl4jsb_memory(MEGAN_)%iso%class    = 17
    dsl4jsb_memory(MEGAN_)%iso%mweight  = 30.01_wp
    dsl4jsb_memory(MEGAN_)%iso%cnumber  = 0
    dsl4jsb_memory(MEGAN_)%iso%cafac    = 1
    dsl4jsb_memory(MEGAN_)%iso%beta     = dsl4jsb_Config(MEGAN_)%beta_NO
    dsl4jsb_memory(MEGAN_)%iso%ldf      = 0._wp
    dsl4jsb_memory(MEGAN_)%iso%scale    = dsl4jsb_Config(MEGAN_)%scale_NO
    dsl4jsb_memory(MEGAN_)%iso%ef       = 0._wp ! TODO: lctlib-something
    
!    IF (debug_on() .OR. dsl4jsb_Config(MEGAN_)%dbg_process) &
!      CALL message(TRIM(routine), 'Setting megan boundary conditions for tile '//TRIM(tile%name))

!    IF (dsl4jsb_Config(TOP_)%inun_prescribed) THEN

      ! Get soil memory of the tile
!      dsl4jsb_Get_memory(MEGAN_)

      !
!      dsl4jsb_memory(TOP_)%frac_inun_clim(:,:,:) = top_init_vars%frac_inun_clim(:,:,:) ! StW: Should dims 2 and 3 be swapped?

!      !$acc data present(dsl4jsb_var3D_onDomain(TOP_,frac_inun_clim))
!      !$acc update device(dsl4jsb_var3D_onDomain(TOP_,frac_inun_clim))
!      !$acc end data

!    ENDIF

  END SUBROUTINE megan_init_bc

  SUBROUTINE megan_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(MEGAN_)
    dsl4jsb_Def_memory(MEGAN_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':megan_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(MEGAN_)

    ! Get soil config
    dsl4jsb_Get_config(MEGAN_)

    IF (.NOT. tile%Is_process_active(MEGAN_)) RETURN

!    IF (debug_on() .OR. dsl4jsb_Config(MEGAN_)%dbg_process) &
!      CALL message(TRIM(routine), 'Initializing megan memory for tile '//TRIM(tile%name))

!    IF (.NOT. dsl4jsb_Config(TOP_)%inun_prescribed) THEN
!      dsl4jsb_memory(TOP_)%cti_openwater(:,:) = top_init_vars%cti_openwater(:,:)

!      !$acc data present(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
!      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_chi))
!      !$acc update device(dsl4jsb_var2D_onDomain(TOP_,cti_mean), &
!      !$acc dsl4jsb_var2D_onDomain(TOP_,cti_chi))
!      !$acc end data
!    ENDIF

  END SUBROUTINE megan_init_ic

#endif
END MODULE mo_megan_init
