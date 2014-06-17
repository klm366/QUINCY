MODULE mo_jsb_io
#ifndef __NO_JSBACH__

  USE mo_jsb_io_iface, ONLY:                                           &
    DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
    ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND,                             &
    FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2,           &
    GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL,                 &
    TSTEP_CONSTANT,                                                    &
    zaxisCreate, zaxisInqName, zaxisInqSize, zaxisDestroy
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: grib_bits     = DATATYPE_PACK16, &
                        grib_extbits  = DATATYPE_PACK24, &
                        nc_bits       = DATATYPE_FLT32,  &
                        nc_extbits    = DATATYPE_FLT64

  TYPE t_cf
    CHARACTER(len=128) :: standard_name = ''
    CHARACTER(len=128) :: units         = ''
    CHARACTER(len=128) :: long_name     = ''
    INTEGER            :: datatype      = -1
  END TYPE t_cf

  TYPE t_grib1
    INTEGER :: table
    INTEGER :: parameter
    INTEGER :: bits
  END TYPE t_grib1

  TYPE t_grib2
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
  END TYPE t_grib2

  INTEGER, PARAMETER :: max_tables = 5
  INTEGER, PARAMETER :: max_cdi_zaxis_ids = 20

  INTEGER :: tables(max_tables) = (/186, 187, 188, 189, 190/)

  INTEGER, SAVE :: cdi_zaxis_ids(max_cdi_zaxis_ids) = 0  ! Collects CDI zaxis IDs of all JSBACH vgrids

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io'

!CONTAINS

#endif
END MODULE mo_jsb_io
