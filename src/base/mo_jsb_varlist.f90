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

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d, REAL2D, REAL3D
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_varlist_iface, ONLY: t_jsb_varlist => t_var_list,      &
                                  jsb_get_var_list => get_var_list, &
                                  jsb_new_var_list => new_var_list, &
                                  add_var_list_element_r2d, add_var_list_element_r3d, &
                                  VARNAME_LEN
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_CELL, missval, tables
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var, grib2_var
  USE mo_jsb_process_class, ONLY: Get_process_name

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_varlist
  PUBLIC :: jsb_add_var, Get_varlist
  PUBLIC :: VARNAME_LEN
  PUBLIC :: BASIC, MEDIUM, FULL

  INTERFACE jsb_add_var
    MODULE PROCEDURE jsb_add_var_list_element_r2d
    MODULE PROCEDURE jsb_add_var_list_element_r3d
    MODULE PROCEDURE jsb_add_var_list_element_r2d_tile
  END INTERFACE jsb_add_var

  ENUM, BIND(C)
    ENUMERATOR :: BASIC=1, MEDIUM, FULL
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist'

CONTAINS

  FUNCTION Get_varlist(name, patch_id) RESULT(return_ptr)

    CHARACTER(LEN=*), INTENT(in) :: name
    INTEGER,          INTENT(in) :: patch_id
    TYPE(t_jsb_varlist), POINTER :: return_ptr

    CALL jsb_get_var_list(return_ptr, TRIM(name))

    IF (.NOT. ASSOCIATED(return_ptr)) THEN
      ALLOCATE(return_ptr)
      CALL jsb_new_var_list(return_ptr, TRIM(name), patch_id, &
        & loutput=.TRUE., &
        & lrestart=.TRUE., table=tables(1))
    END IF

  END FUNCTION Get_varlist

  SUBROUTINE jsb_add_var_list_element_r2d(mem, name, var,       &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose                                       &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 ! memory structure
!!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET, INTENT(inout) :: var                 ! reference to jsbach var
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid type
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              ! Prefix for variable name
    CHARACTER(len=*),     INTENT(in)           :: suffix              ! Suffix for variable name
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     ! aggregate over all child tiles?
    INTEGER,              INTENT(in), OPTIONAL :: output_level        ! minimum level of output groups this var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    CHARACTER(len=:), ALLOCATABLE :: output_level_name

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r2d'

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ldims(1) = hgrid%nproma
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    ALLOCATE(var%field(ldims(1),ldims(2),1,1,1))
    var%ptr => var%field(:,:,1,1,1)
    !$acc enter data create(var%field(:,:,1,1,1))

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3+MERGE(2, 0, PRESENT(output_level))))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        output_level_name = Get_var_group_name(output_level)
        groups(4) = 'jsb_all_'//output_level_name ! var group jsb_all_<level>
        groups(5) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_'//output_level_name ! var group jsb_<process>_<level>
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name_loc), var%ptr, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, p5=var%field, in_groups=groups, &
      verbose=verbose)

    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%type = REAL2D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all

    var%owner_proc_id = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

  END SUBROUTINE jsb_add_var_list_element_r2d

  SUBROUTINE jsb_add_var_list_element_r3d(mem, name, var,        &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,              &
    extra_dim_n, loutput, lcontainer,                            &
    lrestart, lrestart_cont, initval_r, isteptype,               &
    resetval_r, lmiss, missval_r, l_aggregate_all,               &
    output_level, verbose                                        &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    TYPE(t_jsb_var_real3d), INTENT(inout), TARGET :: var                 ! reference to jsbach var
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              ! Prefix for variable name
    CHARACTER(len=*),     INTENT(in)           :: suffix              ! Suffix for variable name
    INTEGER,              INTENT(in), OPTIONAL :: extra_dim_n         ! Extra (third) dimension (not levels)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     ! aggregate over all child tiles?
    INTEGER,              INTENT(in), OPTIONAL :: output_level        ! minimum level of output groups this var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(3), gdims(3)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    CHARACTER(len=:), ALLOCATABLE :: output_level_name

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r3d'

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

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

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    ALLOCATE(var%field(ldims(1),ldims(2),ldims(3),1,1))
    var%ptr => var%field(:,:,:,1,1)
    !$acc enter data create(var%field(:,:,:,1,1))

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3+MERGE(2, 0, PRESENT(output_level))))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        output_level_name = Get_var_group_name(output_level)
        groups(4) = 'jsb_all_'//output_level_name ! var group jsb_all_<level>
        groups(5) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_'//output_level_name ! var group jsb_<process>_<level>
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r3d( &
      this_list, TRIM(name_loc), var%ptr, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, p5=var%field, in_groups=groups, verbose=verbose)

    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%type = REAL3D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all

    var%owner_proc_id = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

  END SUBROUTINE jsb_add_var_list_element_r3d

  SUBROUTINE jsb_add_var_list_element_r2d_tile(varlist_name, name, ptr, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r,                               &
    verbose                                                     &
    )

    CHARACTER(len=*),     INTENT(in)           :: varlist_name
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              ! Prefix for variable name
    CHARACTER(len=*),     INTENT(in)           :: suffix              ! Suffix for variable name
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information

!!$    REAL(wp), POINTER :: ptr(:,:,:,:,:)
    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r2d'

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ldims(1) = hgrid%nproma
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    this_list => Get_varlist(TRIM(varlist_name), hgrid%host_patch_id)

    ALLOCATE(groups(3))
    groups(1) = 'ALL'
    groups(2) = 'jsb_all'
    groups(3) = 'jsb_tile_fractions'

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name_loc), ptr, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, verbose=verbose)

  END SUBROUTINE jsb_add_var_list_element_r2d_tile

  FUNCTION Get_var_group_name(id) RESULT(group_name)

    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), ALLOCATABLE :: group_name

    SELECT CASE (id)
    CASE(BASIC)
      group_name = 'basic'
    CASE(MEDIUM)
      group_name = 'medium'
    CASE(FULL)
      group_name = 'full'
    END SELECT

  END FUNCTION

#endif
END MODULE mo_jsb_varlist
