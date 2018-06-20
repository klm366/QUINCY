!>
!! @brief Contains the interface to io 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                 by Reiner Schnur (2013-04-04)
!! IO namelist subroutine for jsb4 standalone    by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_io_iface
#ifndef __NO_JSBACH__

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Get_netcdf_precision
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24

  ! The following parameter definition is copied from cdi.inc - Steffen 01.09.2017
  !INCLUDE 'cdi.inc'
  INTEGER, PARAMETER :: DATATYPE_FLT32  = 132
  INTEGER, PARAMETER :: DATATYPE_FLT64  = 164
  INTEGER, PARAMETER :: DATATYPE_PACK16 = 16
  INTEGER, PARAMETER :: DATATYPE_PACK24 = 24

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_iface'

CONTAINS

  FUNCTION Get_netcdf_precision() RESULT(return_value)

    INTEGER :: return_value

#ifdef __ICON__
    IF (lnetcdf_flt64_output) THEN
      return_value = DATATYPE_FLT64
    ELSE
      return_value = DATATYPE_FLT32
    END IF
#else
    ! This might be overwritten by ECHAM in mo_output.f90 based on the value of grib precision.
    return_value = DATATYPE_FLT32
#endif

  END FUNCTION Get_netcdf_precision

#endif
END MODULE mo_jsb_io_iface
