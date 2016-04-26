!>
!! @brief Initialization of the the radiation memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-09)
!!
MODULE mo_rad_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes RAD
  dsl4jsb_Use_config(RAD)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: rad_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_init'

CONTAINS

  !
  !> Intialize radiation process (after memory has been set up)
  !
  SUBROUTINE rad_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rad_init'

    CALL rad_init_bc(tile)
!!$    CALL rad_init_ic(tile)

  END SUBROUTINE rad_init

  SUBROUTINE rad_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(RAD)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rad_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD)

!!$    IF (.NOT. ASSOCIATED(rad__mem)) RETURN

    CALL message(TRIM(routine), 'Setting  radiation boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE rad_init_bc

  RECURSIVE SUBROUTINE rad_init_ic(tile)

!    USE mo_rad_process, ONLY: init_radiation_temperature
!    USE mo_rad_util,    ONLY: radiation_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)

    INTEGER :: iradiation
    REAL(wp), POINTER ::                         &
                         return_pointer (:,  :), & !< temporary pointer
                         return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file
    INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rad_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get radiation memory of the tile
    dsl4jsb_Get_memory(RAD)

    ! Get radiation config
    dsl4jsb_Get_config(RAD)

    IF (.NOT. dsl4jsb_Config(RAD)%active) RETURN

    CALL message(TRIM(routine), 'Initializing radiation memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(RAD)%bc_filename), model%grid_id)

    ! background albedo
    return_pointer => input_file%Read_2d(         &
      & variable_name='alb_background',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD,alb_background))
    dsl4jsb_var2D_onDomain(RAD,alb_background) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_background), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_background) >= 0._wp)

    dsl4jsb_var2D_onDomain(RAD,alb_vis) = dsl4jsb_var2D_onDomain(RAD,alb_background)
    dsl4jsb_var2D_onDomain(RAD,alb_nir) = dsl4jsb_var2D_onDomain(RAD,alb_background)

    ! soil albedo
    return_pointer => input_file%Read_2d(         &
      & variable_name='alb_vis_soil',          &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_vis_soil))
    dsl4jsb_var2D_onDomain(RAD, alb_vis_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_vis_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_vis_soil) >= 0._wp)

    return_pointer => input_file%Read_2d(         &
      & variable_name='alb_nir_soil',          &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_nir_soil))
    dsl4jsb_var2D_onDomain(RAD,alb_nir_soil) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD,alb_nir_soil), 0._wp, dsl4jsb_var2D_onDomain(RAD,alb_nir_soil) >= 0._wp)

    ! vegetation albedo
    return_pointer => input_file%Read_2d(         &
      & variable_name='alb_vis_can',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_vis_can))
    dsl4jsb_var2D_onDomain(RAD, alb_vis_can) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD, alb_vis_can), 0._wp, dsl4jsb_var2D_onDomain(RAD, alb_vis_can) >= 0._wp)

    return_pointer => input_file%Read_2d(         &
      & variable_name='alb_nir_can',           &
      & fill_array = dsl4jsb_var2D_ptr(RAD, alb_nir_can))
    dsl4jsb_var2D_onDomain(RAD, alb_nir_can) = &
      MERGE(dsl4jsb_var2D_onDomain(RAD, alb_nir_can), 0._wp, dsl4jsb_var2D_onDomain(RAD, alb_nir_can) >= 0._wp)

    CALL input_file%Close()

    ! Initialize variables
    !

!     dsl4jsb_var2D_onDomain(RAD,sw_srf_net)      = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,swvis_srf_net)   = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,swnir_srf_net)   = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,lw_srf_net)      = 0._wp
!     dsl4jsb_var2D_onDomain(RAD,rad_srf_net)     = 0._wp

  END SUBROUTINE rad_init_ic

#endif
END MODULE mo_rad_init
