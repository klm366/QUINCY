!>
!! @brief Initialization of the the radiation memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
MODULE mo_radiation_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d, jsb_netcdf_read_2d_time, jsb_netcdf_read_2d_extdim

  dsl4jsb_Use_processes RAD
  dsl4jsb_Use_config(RAD)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: radiation_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_radiation_init'

CONTAINS

  !
  !> Intialize radiation process (after memory has been set up)
  !
  SUBROUTINE radiation_init(tile)

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':radiation_init'

    CALL radiation_init_bc(tile)
    CALL radiation_init_ic(tile)

  END SUBROUTINE radiation_init

  
  SUBROUTINE radiation_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS 
    TYPE(t_jsb_tile), INTENT(inout) :: tile
    
    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(RAD)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':radiation_init_bc'
    
    !-----------------------------------------------------------------------
    ! CONTENT
    
    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD)

    IF (.NOT. ASSOCIATED(rad__mem)) RETURN

    CALL message(TRIM(routine), 'Setting  radiation boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...
    
  END SUBROUTINE radiation_init_bc

  
  RECURSIVE SUBROUTINE radiation_init_ic(tile)

!    USE mo_radiation_process, ONLY: init_radiation_temperature
!    USE mo_radiation_util,    ONLY: radiation_depth_to_layers_2d

    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)

    INTEGER :: iradiation
    REAL(wp), POINTER ::                         &
                         return_pointer (:,  :), & !< temporary pointer
                         return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file
    INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':radiation_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD)

    ! Get radiation config
    dsl4jsb_Get_config(RAD)

    IF (.NOT. dsl4jsb_Config(RAD)%active) RETURN

    CALL message(TRIM(routine), 'Initializing radiation memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(RAD)%bc_filename), model%id)

    ! background albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='alb_background',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD,alb_background))
    dsl4jsb_var2D_onDomain(RAD,alb_background) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_background), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_background) >= 0._wp)

    dsl4jsb_var2D_onDomain(RAD,alb_vis_srf) = dsl4jsb_var2D_onDomain(RAD,alb_background)
    dsl4jsb_var2D_onDomain(RAD,alb_nir_srf) = dsl4jsb_var2D_onDomain(RAD,alb_background)
    
    
    ! soil albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='alb_vis_soil',          &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_vis_soil))
    dsl4jsb_var2D_onDomain(RAD, alb_vis_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_vis_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_vis_soil) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='alb_nir_soil',          &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_nir_soil))
    dsl4jsb_var2D_onDomain(RAD,alb_nir_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_nir_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_nir_soil) >= 0._wp)
    
    ! vegetation albedo
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='alb_vis_can',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_vis_can))
    dsl4jsb_var2D_onDomain(RAD, alb_vis_can) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD, alb_vis_can), 0._wp, dsl4jsb_var2D_onDomain(RAD, alb_vis_can) >= 0._wp)

    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='alb_nir_can',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_nir_can))
    dsl4jsb_var2D_onDomain(RAD, alb_nir_can) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD, alb_nir_can), 0._wp, dsl4jsb_var2D_onDomain(RAD, alb_nir_can) >= 0._wp)

    CALL jsb_netcdf_close_file(input_file)    
   
    
    
    ! Initialize variables
    !

!     dsl4jsb_var2D_onDomain(RAD,sw_srf_net)      = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,swvis_srf_net)   = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,swnir_srf_net)   = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,lw_srf_net)      = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,rad_srf_net)     = 0._wp

    ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
    ! TBD: sub-tiles should not read the netcdf files again
!!$    IF (tile%no_of_sub_tiles > 0) THEN
!!$      DO i=1,tile%no_of_sub_tiles
!!$        CALL radiation_init_ic(tile%sub_tiles(i))
!!$      END DO
!!$    END IF

  END SUBROUTINE radiation_init_ic

#endif
END MODULE mo_radiation_init
