!>
!! @brief Initialization of the the vegetation memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes VEG, PHENO, SOIL, HYDRO
  dsl4jsb_Use_config(VEG)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(HYDRO)

  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(PHENO)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(HYDRO)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: veg_init_bc, veg_init_ic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  SUBROUTINE veg_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(VEG)

    REAL(wp), POINTER :: lai_mon(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    dsl4jsb_Get_memory(VEG)

    ! ...

  END SUBROUTINE veg_init_bc

  SUBROUTINE veg_init_ic(tile)

    USE mo_veg_process, ONLY: get_foliage_projected_cover
    !USE mo_veg_util,    ONLY: get_water_in_root_zone
    USE mo_soil_util,   ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract),  INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    dsl4jsb_Def_config(VEG)
    dsl4jsb_Def_config(HYDRO)
    dsl4jsb_Def_config(SOIL)

    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(HYDRO)

    INTEGER :: i

    REAL(wp), POINTER ::        &
      & lai_mon(:,:,:),         &
      & return_pointer (:,:)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get veg memory of the tile
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(PHENO)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(HYDRO)


    ! Get veg config
    dsl4jsb_Get_config(VEG)
    dsl4jsb_Get_config(SOIL)
    dsl4jsb_Get_config(HYDRO)

    IF (.NOT. dsl4jsb_Config(VEG)%active) RETURN

    ! Initialize tile fractions
    ! TBD

    CALL message(TRIM(routine), 'Initializing veg memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(VEG)%bc_filename), model%grid_id)

    ! Forest fraction
    return_pointer => input_file%Read_2d(        &
      & variable_name='forest_fract',            &
      & fill_array = dsl4jsb_var2D_ptr(VEG, fract_forest))
    dsl4jsb_var2D_onDomain(VEG, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(VEG, fract_forest) >= 0._wp)

! R    dsl4jsb_var2D_onDomain(VEG, canopy_cond) = 0._wp


    CALL input_file%Close()

    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    return_pointer => input_file%Read_2d(        &
      & variable_name='veg_ratio_max',           &
      & fill_array = dsl4jsb_var2D_ptr(VEG, fract_fpc_max))
    dsl4jsb_var2D_onDomain(VEG, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(VEG, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(VEG, fract_fpc_max) >= 0._wp)

!!$    veg_mem%fract_fpc(:,:) = get_foliage_projected_cover(veg_mem%fract_fpc_max%ptr(:,:), veg_mem%lai%ptr(:,:), veg_mem%clump%ptr(:,:))
    dsl4jsb_var2D_onDomain(VEG, fract_fpc) = &
      get_foliage_projected_cover(dsl4jsb_var2D_onDomain(VEG, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO, lai), &
        & SPREAD(SPREAD(2._wp, DIM=1, NCOPIES=SIZE(dsl4jsb_var2D_ptr(PHENO, lai),1)), &
        &        DIM=2, NCOPIES=SIZE(dsl4jsb_var2D_ptr(PHENO, lai),2))) ! Temporary

    CALL input_file%Close()

  END SUBROUTINE veg_init_ic

#endif
END MODULE mo_veg_init
