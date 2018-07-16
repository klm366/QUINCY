!>
!! @brief Contains the interface to varlists/streams
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_varlist_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp, dp
  USE mo_exception,          ONLY: finish
#ifdef __ICON__
  USE mo_var_list,           ONLY: new_var_list_icon => new_var_list, &
                                   default_var_list_settings,         &
                                   get_var_list,                      &
                                   add_var_icon => add_var,           &
                                   find_list_element
  USE mo_var_metadata,       ONLY: groups
  USE mo_var_metadata_types, ONLY: t_var_metadata, VARNAME_LEN
  USE mo_linked_list,        ONLY: t_var_list, t_list_element
#else
  USE mo_memory_base,        ONLY: new_stream,                       &
                                   get_var_list => get_stream,       &
                                   add_stream_element,               &
                                   default_stream_setting
  USE mo_linked_list,        ONLY: t_var_list => t_stream,           &
                                   t_var_metadata => memory_info,    &
                                   t_list_element => list_element,   &
                                   find_list_element
#endif

  USE mo_jsb_io_iface, ONLY: t_cf_var, t_grib2_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: VARNAME_LEN
  PUBLIC :: t_var_list, t_var_metadata, t_list_element, get_var_list
  PUBLIC :: new_var_list
  PUBLIC :: add_var_list_element_r2d, add_var_list_element_r3d

#ifndef __ICON__
  INTEGER, PARAMETER :: VARNAME_LEN = 64
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist_iface'

CONTAINS

  SUBROUTINE new_var_list (this_list, name, patch_id, output_type, restart_type, &
       &                       post_suf, rest_suf, init_suf, loutput, lrestart,  &
       &                       linitial, table)
    !
    TYPE(t_var_list), POINTER, INTENT(inout) :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)             :: name         ! name of output var_list
    INTEGER,          INTENT(in)             :: patch_id     ! patch ID
    INTEGER,          INTENT(in), OPTIONAL   :: output_type  ! 'GRIB1', 'GRIB2' or 'NetCDF[12]'
    INTEGER,          INTENT(in), OPTIONAL   :: restart_type ! 'NetCDF[12]'
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: post_suf     ! suffix of output file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: rest_suf     ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: init_suf     ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL   :: loutput      ! write to  output file
    LOGICAL,          INTENT(in), OPTIONAL   :: lrestart     ! write to restart file
    LOGICAL,          INTENT(in), OPTIONAL   :: linitial     ! read from initial file
    INTEGER,          INTENT(in), OPTIONAL   :: table

#ifdef __ICON__
    CALL new_var_list_icon(this_list, name,                                    &
                      output_type=output_type,                                 &
                      post_suf=post_suf, rest_suf=rest_suf, init_suf=init_suf, &
                      loutput=loutput, lrestart=lrestart, linitial=linitial,   &
                      patch_id=patch_id                                        &
                     )
#else
    INTEGER :: table_loc

    IF(PRESENT(table)) THEN
      table_loc = table
    ELSE
      table_loc = 180
    END IF
    CALL new_stream(this_list, name,                                         &
                    filetype=output_type,                                    &
                    post_suf=post_suf, rest_suf=rest_suf, init_suf=init_suf, &
                    lpost=loutput, lpout=loutput, lrerun=lrestart            &
                   )

    !@todo: hardcoded FALSE for now to reduce output size when running with ECHAM
    !       until the "in_groups" behavior is also implemented for ECHAM!
    ! - now: for required vars loutput needs to be set to true! (+ see namelist selector in jsb3 sitelevel)
    !@todo: implement the "in_groups" behavior for ECHAM!
    CALL default_stream_setting(this_list, table=table_loc, lpost=.FALSE., lrerun=lrestart)
#endif

  END SUBROUTINE new_var_list

  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(2)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

#ifndef __ICON__
    REAL(wp), POINTER :: p4(:,:,:,:)
#endif

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r2d'

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
#ifndef __ICON__
      p4 => p5(:,:,:,:,1)
#endif
    END IF

#ifdef __ICON__
    IF (PRESENT(in_groups)) THEN
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, in_group=groups(in_groups), verbose=verbose, new_element=new_element)
    ELSE
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, new_element=new_element)
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%ndims = 2
    element%field%info%used_dimensions(1:2) = ldims(1:2)
#else
    IF (PRESENT(p5)) THEN
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lrerun=lrestart, contnorest=lrestart_cont, &
        lmiss=lmiss, missval=missval_r, p4=p4, verbose=verbose &
        )
    ELSE
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lrerun=lrestart, contnorest=lrestart_cont, &
        lmiss=lmiss, missval=missval_r, verbose=verbose &
        )
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%zaxisID = vgrid

    IF (PRESENT(initval_r)) THEN
      ptr(:,:) = initval_r
    END IF
#endif

  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(3)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             POINTER,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

#ifndef __ICON__
    REAL(wp), POINTER :: p4(:,:,:,:)
#endif

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r3d'

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
#ifndef __ICON__
      p4 => p5(:,:,:,:,1)
#endif
    END IF

#ifdef __ICON__
    IF (PRESENT(in_groups)) THEN
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, in_group=groups(in_groups), new_element=new_element)
    ELSE
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, new_element=new_element)
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%ndims = 3
    element%field%info%used_dimensions(1:3) = ldims(1:3)
#else
    IF (PRESENT(p5)) THEN
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lmiss=lmiss, missval=missval_r, lrerun=lrestart, contnorest=lrestart_cont, &
        p4=p4 &
      !mem_info=info &
        )
    ELSE
      CALL add_stream_element(this_list, TRIM(name), ptr, &
        ldims=ldims, gdims=gdims, leveltype=levelindx, &
        longname=cf%long_name, units=cf%units, &
        code=code, table=table, &
        lpost=loutput, lmiss=lmiss, missval=missval_r, lrerun=lrestart, contnorest=lrestart_cont &
        )
    END IF

    element => find_list_element(this_list, TRIM(name))
    element%field%info%zaxisID = vgrid

    IF (PRESENT(initval_r)) THEN
      ptr(:,:,:) = initval_r
    END IF
#endif

  END SUBROUTINE add_var_list_element_r3d

#endif
END MODULE mo_jsb_varlist_iface
