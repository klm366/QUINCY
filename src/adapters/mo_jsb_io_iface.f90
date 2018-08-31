!>
!! @brief Contains the interface to io
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                  by Reiner Schnur (2013-04-04)
!! IO namelist subroutine for jsb4 standalone     by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_io_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm
  USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, close_nml, position_nml
#ifdef __ICON__
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_grib2,         ONLY: t_grib2_var, grib2_var
  USE mo_io_config,     ONLY: lnetcdf_flt64_output
  USE mo_cdi,           ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24, &
                            & ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_LAND, &
                            & FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
                            & GRID_UNSTRUCTURED, &
                            & TSTEP_CONSTANT, TSTEP_INSTANT, &
                            & cdiDefMissval
  USE mo_cdi_constants, ONLY: GRID_CELL, GRID_UNSTRUCTURED_CELL
  USE mo_zaxis_type,    ONLY: zaxisTypeList, ZA_SURFACE
#else
  USE mo_linked_list,   ONLY: NETCDF2, NETCDF4, GRIB
  USE mo_netcdf,        ONLY: add_dim
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var, grib2_var
  PUBLIC :: Get_netcdf_precision
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
            ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_LAND,               &
            FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
            GRID_CELL, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, &
            TSTEP_CONSTANT, TSTEP_INSTANT, &
            Create_zaxis, cdiDefMissval, read_io_namelist
            ! Create_zaxis, Destroy_zaxis, cdiDefMissval, read_io_namelist


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
  INTEGER, PARAMETER :: GRID_CELL              =  1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL =  1
#endif

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

  SUBROUTINE Create_zaxis(ZaxisID, echamZaxisIdx, cdi_axis_type, length, name, &
    & levels, units, longname, lbounds, ubounds)

    INTEGER,           INTENT(out)   :: ZaxisID
    INTEGER,           INTENT(out)   :: echamZaxisIdx
    INTEGER,           INTENT(in)    :: cdi_axis_type
    INTEGER,           INTENT(in)    :: length
    CHARACTER(len=*),  INTENT(in)    :: name
    REAL(wp),          INTENT(in)    :: levels(:)
    CHARACTER(len=*),  INTENT(in)    :: units
    CHARACTER(len=*),  INTENT(in)    :: longname
    REAL(wp),          INTENT(in)    :: lbounds(:)
    REAL(wp),          INTENT(in)    :: ubounds(:)

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Create_zaxis'

#ifdef __ICON__
    echamZaxisIdx = -1

    IF (TRIM(name) == 'surface') THEN
      ZaxisID = ZA_SURFACE
      RETURN
    ELSE
      ZaxisID = zaxisTypeList%register(cdi_zaxis_type=cdi_axis_type, is_2D=(length==1))
    END IF
#else
    ZaxisID = zaxisCreate(cdi_axis_type, length)
    CALL zaxisDefName    (ZaxisID, TRIM(name))
    CALL zaxisDefLevels  (ZaxisID, levels(1:length))
    CALL zaxisDefLbounds (ZaxisID, lbounds(1:length))
    CALL zaxisDefUbounds (ZaxisID, ubounds(1:length))
    CALL zaxisDefUnits   (ZaxisID, TRIM(units))
    CALL zaxisDefLongname(ZaxisID, TRIM(longname))

    CALL add_dim(TRIM(name), length, &
      & longname=TRIM(longname), units=TRIM(units), VALUE=levels(1:length), &
      & indx=echamZaxisIdx)
#endif

  END SUBROUTINE Create_zaxis

!   SUBROUTINE Destroy_zaxis(ZaxisID)

!     INTEGER, INTENT(in) :: ZaxisID

!     ! In ICON, mo_vertical_axis_types:t_verticalAxisList_finalize takes care of destroying the vertical
!     ! axes and corresponding CDI zaxes. It is called in mo_name_list_output:close_name_list_output
! #ifndef __ICON__
!     CALL zaxisDestroy(ZaxisID)
! #endif

!   END SUBROUTINE Destroy_zaxis

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
      CALL p_bcast(out_expname, p_io, comm=mpi_comm)
#endif
      CALL p_bcast(read_interface_vars, p_io, comm=mpi_comm)
      CALL p_bcast(write_interface_vars, p_io, comm=mpi_comm)
    ENDIF
#endif

  END SUBROUTINE read_io_namelist

#endif
END MODULE mo_jsb_io_iface
