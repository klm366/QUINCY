!>
!! @brief Contains the interface to io 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! JN (2016-02-20) io namelist subroutine for jsb4 standalone
!!
MODULE mo_jsb_io_iface
#ifndef __NO_JSBACH__
  USE mo_kind,               ONLY: wp
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel
  USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, close_nml, position_nml

#ifdef __ICON__
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_grib2,         ONLY: t_grib2_var, grib2_var
  USE mo_cdi,           ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24, ZAXIS_SURFACE, &
                            & ZAXIS_DEPTH_BELOW_LAND, FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, TSTEP_CONSTANT, &
                            & zaxisCreate, zaxisDestroy, zaxisInqName, zaxisInqSize, &
                            & zaxisDefName, zaxisDefLevels, zaxisDefLbounds,  &
                            & zaxisDefUbounds, zaxisDefUnits, zaxisDefLongname 
  USE mo_cdi_constants, ONLY: GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL, ZA_SURFACE, ZA_DEPTH_BELOW_LAND
#else
  USE mo_linked_list,   ONLY: NETCDF2, NETCDF4, GRIB
  USE mo_netcdf,        ONLY: add_dim
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var, grib2_var
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
            ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND,                   &
            FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
            GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL, &
            ZA_SURFACE, ZA_DEPTH_BELOW_LAND, TSTEP_CONSTANT, &
            zaxisCreate, zaxisDestroy, zaxisInqName, zaxisInqSize, &
            zaxisDefName, zaxisDefLevels, zaxisDefLbounds,  &
            zaxisDefUbounds, zaxisDefUnits, zaxisDefLongname, &
            Create_zaxis, read_io_namelist


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

CONTAINS

#ifndef __ICON__

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

    SUBROUTINE Create_zaxis(cdiZaxisID, echamZaxisID, cdi_axis_type, length, name, &
    & levels, units, longname, lbounds, ubounds)

    INTEGER,           INTENT(inout) :: cdiZaxisID
    INTEGER,           INTENT(inout) :: echamZaxisID
    INTEGER,           INTENT(in)    :: cdi_axis_type
    INTEGER,           INTENT(in)    :: length
    CHARACTER(len=*),  INTENT(in)    :: name
    REAL(wp),          INTENT(in)    :: levels(:)
    CHARACTER(len=*),  INTENT(in)    :: units
    CHARACTER(len=*),  INTENT(in)    :: longname
    REAL(wp),          INTENT(in)    :: lbounds(:)
    REAL(wp),          INTENT(in)    :: ubounds(:)

    cdiZaxisID = zaxisCreate(cdi_axis_type, length)
    CALL zaxisDefName    (cdiZaxisID, TRIM(name))
    CALL zaxisDefLevels  (cdiZaxisID, levels(1:length))
    CALL zaxisDefLbounds (cdiZaxisID, lbounds(1:length))
    CALL zaxisDefUbounds (cdiZaxisID, ubounds(1:length))
    CALL zaxisDefUnits   (cdiZaxisID, TRIM(units))
    CALL zaxisDefLongname(cdiZaxisID, TRIM(longname))

#ifdef __ICON__
    echamZaxisID = -1
#else
    CALL add_dim(TRIM(name), length, &
      & longname=TRIM(longname), units=TRIM(units), VALUE=levels(1:length), &
      & indx=echamZaxisID)
#endif

  END SUBROUTINE Create_zaxis
  
  SUBROUTINE read_io_namelist(filename)

    USE mo_jsb_test, ONLY: read_interface_vars, write_interface_vars

#ifndef __ICON__
#ifdef STANDALONE
    USE mo_filename, ONLY: out_expname
#endif
#endif

    CHARACTER(LEN=*), INTENT(in) :: filename

    INTEGER :: nml_handler, nml_unit, istat

    NAMELIST /jsb_io_nml/ &
#ifdef STANDALONE
      out_expname,          &
      read_interface_vars,  &
#endif
      write_interface_vars

#ifdef __ICON__
    ! TODO: mo_jsb_test for now not inplemented with ICON
    read_interface_vars = .FALSE.
    write_interface_vars = .FALSE.
#else


    IF(my_process_is_stdio()) THEN
      ! default values
#ifdef STANDALONE
      out_expname = 'undefined_expname'
#endif
      read_interface_vars = .FALSE.
      write_interface_vars = .FALSE.

      nml_handler = open_nml(TRIM(filename))

      nml_unit =  position_nml('jsb_io_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_io_nml)

      CALL close_nml(nml_handler)

    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
#ifdef STANDALONE
      CALL p_bcast(out_expname, p_io)
#endif
      CALL p_bcast(read_interface_vars, p_io)
      CALL p_bcast(write_interface_vars, p_io)
    ENDIF
#endif
  END SUBROUTINE read_io_namelist


#endif
END MODULE mo_jsb_io_iface
