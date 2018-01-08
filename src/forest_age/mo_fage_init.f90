!>
!! @brief Initialization of the forest-age memory.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!

MODULE mo_fage_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, t_jsb_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes FAGE_
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_memory(FAGE_)
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: fage_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_init'

CONTAINS

  SUBROUTINE fage_init(tile)
    
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fage_init'
    
    CALL fage_init_bc(tile)
    CALL fage_init_ic(tile)

  END SUBROUTINE fage_init

  SUBROUTINE fage_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    REAL(wp), POINTER :: &
      & ptr_3D(:,  :)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: forest_age

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fage_init_bc'

    IF (.NOT. tile%Is_process_active(FAGE_)) RETURN

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    !TODO-JN any bcs??
!    dsl4jsb_var3D_onDomain(FM_, rotation_cycle) = REAL(dsl4jsb_Config(FM_)%global_rotation_cycle, wp)
!    IF (dsl4jsb_Config(FM_)%scheme == 'rotation_local' ) THEN 
!      IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of fm memory (rotation_cycle) for tile '// &
!        &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(FM_)%bc_filename))
!
!      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(FM_)%bc_filename), model%grid_id)
!      ptr_2D => input_file%Read_2d(        &
!        & variable_name='rotation_cycle',           &
!        & fill_array = dsl4jsb_var_ptr(FM_, rotation_cycle))
!    ENDIF

  END SUBROUTINE fage_init_bc


  SUBROUTINE fage_init_ic(tile)

    USE mo_fage_util,      ONLY: get_mean_age

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    REAL(wp), POINTER :: ptr_2D(:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: forest_age

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fage_init_ic'

    !-----------------------------------------------------------------------
    dsl4jsb_Real3D_onDomain :: fract_per_age


    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onDomain(FAGE_, fract_per_age)

    !TODO-JN-COHORTS: for now all fracts to youngest age, later initialise from file
    dsl4jsb_var3D_onDomain(FAGE_, fract_per_age) = 0._wp
    dsl4jsb_var_ptr(FAGE_, fract_per_age) (:,1,:) = 1._wp

    ! input_file = jsb_netcdf_open_input(dsl4jsb_Config(FAGE_)%ic_filename, model%grid_id)
!    ptr_2D(:,:) = tile%fract(:,:)

!----------------------------------------------------------------------------------------
!TODO-JN seuftz, folgendes wäre nur im absolute Flächen Fall nötig gewesen, aber den haben wir ja noch gar nicht!
!    ALLOCATE(ptr_2D(grid%nproma, grid%nblks))
!    CALL tile%Get_fraction(fract=ptr_2D(:,:))
!    dsl4jsb_var_ptr(FAGE_, fract_per_age) (:,1,:) = ptr_2D(:,:) !TODO-JN check: arrays seem to not be 0 based?!
!    DEALLOCATE(ptr_2D)
!----------------------------------------------------------------------------------------


    IF (tile%Has_children()) THEN

!      DO isoil=1,nsoil
!        ptr_3D => input_file%Read_2d_extdim(  &
!          & variable_name='layer_moist',                &
!          & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
!        dsl4jsb_var_ptr(HYDRO_,w_soil_sl) (:,isoil,:) = &
!          & MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
!      END DO

      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = get_mean_age(dsl4jsb_var3D_onDomain(FAGE_, fract_per_age))
      !avgAge = util function - age weighted with area

    ELSE
      !Only one age -> for now initialise age with zero later e.g. from file
      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = 0._wp
    ENDIF

  END SUBROUTINE fage_init_ic
  
#endif
END MODULE mo_fage_init
