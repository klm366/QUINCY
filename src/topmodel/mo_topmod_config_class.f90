!>
!! @brief Contains structures and methods for the topmodel for calculating innundation fractions
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

MODULE mo_top_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_top_config

  TYPE, EXTENDS(t_jsb_config) :: t_top_config
    LOGICAL :: dbg_topmodel
  CONTAINS
    PROCEDURE :: Init => Init_topmodel_config
  END TYPE t_top_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_topconfig_class'

CONTAINS

  SUBROUTINE Init_topmodel_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
!    USE mo_jsb_grid,           ONLY: Register_vgrid
!    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input

    CLASS(t_top_config), INTENT(inout) :: config

    LOGICAL :: active, dbg_topmodel
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_topmodel_nml/                      &
      & active, dbg_topmodel,                        &
      & ic_filename, bc_filename

    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_topmodel_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting topmodel configuration')

    ! Set defaults
    active          = .FALSE.
    ic_filename      = 'ic_land_soil.nc'
    bc_filename      = 'bc_land_soil.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_topmodel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_topmodel_nml)

    CALL close_nml(nml_handler)

    config%active        = active
    config%dbg_topmodel  = dbg_topmodel
    config%ic_filename   = ic_filename
    config%bc_filename   = bc_filename

    write(message_text,*) active,dbg_topmodel
    CALL message(routine,'Topmodel settings: '//trim(message_text))

    !! -----------------------------------------------------------------------------------------------------------------------------
    IF (.NOT. active) RETURN

    input_file = jsb_netcdf_open_input(ic_filename)

    ! @todo: At the moment, the soil layers for the energy calculations are the same as for the hydrology
!    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    CALL input_file%Close()

!    vgrid_soil_sb => new_vgrid('soil_layer_sb', ZAXIS_GENERIC, nsoil, &
!                      longname='Soil (water) layers from SSE_ QUINCY', &
!                      units='m', &
!                      ! levels=  , &      ! is calculated from ubounds and dz
!                      ! lbounds=  , &     ! is calculated from ubounds and dz
!                      !ubounds=ubounds_soil_lay_water(:), & 
!                      ubounds=soil_vgrid%ubounds(:), & 
!                      !dz=dz_water(:))
!                      dz=soil_vgrid%dz(:))
!    CALL register_vgrid(vgrid_soil_sb)

  END SUBROUTINE Init_topmodel_config

#endif
END MODULE mo_top_config_class
