!>
!! @brief Contains routines for var_lists/streams
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-12-04)
!!
MODULE mo_jsb_varlist
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_varlist_iface, ONLY: t_var_list,                       &
                                  jsb_get_var_list => get_var_list, &
                                  jsb_new_var_list => new_var_list, &
                                  t_var_metadata, t_list_element,   &
                                  add_var_list_element_r2d, add_var_list_element_r3d
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, GRID_UNSTRUCTURED_CELL, GRID_REFERENCE, GRID_CELL
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var

  IMPLICIT NONE
  PRIVATE

  INTERFACE jsb_add_var
    MODULE PROCEDURE jsb_add_var_list_element_r2d
    MODULE PROCEDURE jsb_add_var_list_element_r3d
  END INTERFACE jsb_add_var

  PUBLIC :: t_var_list
  PUBLIC :: jsb_new_var_list, jsb_get_var_list, jsb_add_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist'

CONTAINS

  SUBROUTINE jsb_add_var_list_element_r2d(this_list, name, ptr, &
    hgrid, vgrid, cf, grib1, grib2, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, tlev_source, info,            &
    verbose, new_element, cdiTimeID                             &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid type
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)

    TYPE(t_cf_var) :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    INTEGER :: ldims(2), gdims(2)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_add_var_list_element_r2d'

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = t_grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_REFERENCE, GRID_CELL)

    ldims(1) = hgrid%nproma
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name), ptr, GRID_UNSTRUCTURED_CELL, vgrid%cdiZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      lmiss=lmiss, missval_r=missval_r)

  END SUBROUTINE jsb_add_var_list_element_r2d

  SUBROUTINE jsb_add_var_list_element_r3d(this_list, name, ptr, &
    hgrid, vgrid, cf, grib1, grib2, extra_dim_n, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,        &
    verbose, new_element, cdiTimeID                             &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:)          ! reference to field
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid type
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in), OPTIONAL :: extra_dim_n         ! Extra (third) dimension (not levels)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element
    INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)

    TYPE(t_cf_var) :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    INTEGER :: ldims(3), gdims(3)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':jsb_add_var_list_element_r3d'

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = t_grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_REFERENCE, GRID_CELL)

    ldims(1) = hgrid%nproma
    gdims(1) = hgrid%dims_g(1)
    IF (PRESENT(extra_dim_n)) THEN
      IF (PRESENT(loutput)) THEN
        IF (loutput) CALL finish(TRIM(routine), 'extra_dim_n not possible with loutput=TRUE')
      ELSE
        CALL finish(TRIM(routine), 'extra_dim_n requires loutput=FALSE')
      END IF
      IF (vgrid%n_levels > 1) CALL finish(TRIM(routine), 'nlevels>1 and extra_dim_n not allowed together')
      ldims(2) = hgrid%nblks
      ldims(3) = extra_dim_n
      gdims(2) = hgrid%dims_g(2)
      gdims(3) = extra_dim_n
    ELSE
      ldims(2) = vgrid%n_levels
      ldims(3) = hgrid%nblks
      gdims(2) = vgrid%n_levels
      gdims(3) = hgrid%dims_g(2)
    END IF

    CALL add_var_list_element_r3d( &
      this_list, TRIM(name), ptr, GRID_UNSTRUCTURED_CELL, vgrid%cdiZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      lmiss=lmiss, missval_r=missval_r, p5=p5)

  END SUBROUTINE jsb_add_var_list_element_r3d

#endif
END MODULE mo_jsb_varlist
