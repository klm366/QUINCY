!>
!! @brief Contains the interface to io 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_io_iface
#ifndef __NO_JSBACH__

#ifdef __ICON__
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_grib2,         ONLY: t_grib2_var, grib2_var
  USE mo_cdi,           ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24, ZAXIS_SURFACE, &
                            & ZAXIS_DEPTH_BELOW_LAND, FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, TSTEP_CONSTANT, &
                            & zaxisCreate, zaxisDestroy, zaxisInqName, zaxisInqSize
  USE mo_cdi_constants, ONLY: GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL, ZA_SURFACE, ZA_DEPTH_BELOW_LAND
#else
  USE mo_linked_list,   ONLY: NETCDF2, NETCDF4, GRIB
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var, grib2_var
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
            ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND,                   &
            FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
            GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL, &
            ZA_SURFACE, ZA_DEPTH_BELOW_LAND, TSTEP_CONSTANT, &
            zaxisCreate, zaxisDestroy, zaxisInqName, zaxisInqSize

#ifndef __ICON__
  TYPE t_cf_var
    CHARACTER(len=128) :: standard_name = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=128) :: long_name     = ''
    INTEGER            :: datatype      = -1
  END TYPE t_cf_var

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
  END TYPE t_grib2_var

  INCLUDE 'cdi.inc'

  ! From ICON mo_cdi_constants.f90
  INTEGER, PARAMETER :: ZA_SURFACE             =  1
  INTEGER, PARAMETER :: ZA_DEPTH_BELOW_LAND    = 11
  INTEGER, PARAMETER :: GRID_CELL              =  1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL =  1
  INTEGER, PARAMETER :: GRID_REFERENCE         =  9
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_iface'

#ifndef __ICON__
CONTAINS

  ! constructor for GRIB2 derived data type
  FUNCTION grib2_var(discipline, category, number, bits, gridtype, subgridtype)
    INTEGER, INTENT(IN) :: discipline
    INTEGER, INTENT(IN) :: category
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: bits
    INTEGER, INTENT(IN) :: gridtype
    INTEGER, INTENT(IN) :: subgridtype
    TYPE(t_grib2_var) :: grib2_var

    grib2_var%discipline  = discipline
    grib2_var%category    = category
    grib2_var%number      = number
    grib2_var%bits        = bits
    grib2_var%gridtype    = gridtype
    grib2_var%subgridtype = subgridtype

!!$    grib2_var%additional_keys%nint_keys = 0
  END FUNCTION grib2_var
#endif
  
#endif
END MODULE mo_jsb_io_iface
