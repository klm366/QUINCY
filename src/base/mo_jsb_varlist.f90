!>
!! @brief Contains routines for var_lists/streams
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-12-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_varlist
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: message, finish
  USE mo_jsb_control,       ONLY: jsbach_runs_standalone
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_var_class,     ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d, REAL2D, REAL3D
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_varlist_iface, ONLY: t_jsb_varlist => t_var_list,      &
                                  jsb_get_var_list => get_var_list, &
                                  jsb_new_var_list => new_var_list, &
                                  add_var_list_element_r2d, add_var_list_element_r3d, &
                                  is_variable_in_output, &
                                  VARNAME_LEN
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_CELL, missval, tables
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var, grib2_var
  USE mo_jsb_process_class, ONLY: Get_process_name

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)  
#else
#define __acc_attach(ptr)
#endif

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

      IF (jsbach_runs_standalone()) return_ptr%p%model_type = 'jsbach'

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
    CHARACTER(len=50) :: name_loc, name_loc_t
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2), ldims_t(3), gdims_t(3)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups                           ! number of additional groups this variable belongs to
    INTEGER           :: out_level, igroup                   ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable
    CLASS(t_jsb_var), POINTER :: var_in_parent => NULL()

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r2d'

    IF (vgrid%n_levels /= 1) CALL finish(routine, 'Only one vertical level in vgrid allowed!')

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    name_loc_t = TRIM(name_loc)//'_t'
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

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

    IF (PRESENT(output_level)) THEN
      no_groups = 2 * (4 - output_level)
    ELSE
      no_groups = 0
    END IF
    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3 + no_groups))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        igroup = 4
        DO out_level=output_level,3
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    IF (mem%single_tile_output .OR. ASSOCIATED(mem%children) .OR. mem%pos_in_leaves == 0) THEN

      CALL add_var_list_element_r2d( &
        this_list, TRIM(name_loc), var%ptr2d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
        cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
        ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
        loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
        initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
        ! lmiss=lmiss_loc, missval_r=missval_r_loc)
        lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, &
        verbose=verbose)
      __acc_attach(var%ptr2d)
      !$acc update device(var%ptr2d)
      var%ptr => var%ptr2d
      __acc_attach(var%ptr)
      ! CALL message(routine, 'Added var '//TRIM(name_loc))

    END IF

    IF (.NOT. mem%single_tile_output) THEN
  
      IF (.NOT. ASSOCIATED(mem%parent) .AND. mem%no_of_leaves > 1) THEN
        ldims_t(1) = ldims(1)
        ldims_t(2) = mem%no_of_leaves
        ldims_t(3) = ldims(2)
        ! gdims_t(1) = gdims(1)
        ! gdims_t(2) = mem%no_of_leaves
        ! gdims_t(3) = gdims(2)
        CALL add_var_list_element_r3d( &
          this_list, TRIM(name_loc_t), var%ptr2d_t, GRID_UNSTRUCTURED_CELL, mem%vgrid_t%ZaxisID, &
          cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
          ldims=ldims_t, gdims=gdims_t, levelindx=mem%vgrid_t%echamZaxisIdx, &
          loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
          isteptype=isteptype, resetval_r=resetval_r, &
          ! lmiss=lmiss_loc, missval_r=missval_r_loc)
          initval_r=missval_r_loc, & ! initval_r=initval_r, 
          lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, &
          verbose=verbose)
        __acc_attach(var%ptr2d_t)
        !$acc update device(var%ptr2d_t)
        var%ptr_t => var%ptr2d_t
        __acc_attach(var%ptr_t)
        ! CALL message(routine, 'Added var '//TRIM(name_loc_t)//' for process leaves on box tile')
      ELSE IF (mem%no_of_leaves > 1) THEN
        var_in_parent => mem%parent%Get_var(TRIM(name))
        var%ptr2d_t => var_in_parent%ptr2d_t
        IF (mem%pos_in_leaves > 0) THEN
          var%pos_in_leaves = mem%pos_in_leaves
          var%ptr => var%ptr2d_t(:,mem%pos_in_leaves,:)
          IF (PRESENT(initval_r)) THEN
            var%ptr(:,:) = initval_r
          END IF
          __acc_attach(var%ptr)
          ! CALL message(routine, 'Linked var '//TRIM(name_loc)//' to process leaves')
        END IF
      END IF
  
    END IF
      
    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%type = REAL2D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all
    IF (loutput_loc) THEN
      IF (mem%in_var_groups) THEN
        var%is_in_output = is_variable_in_output(var%full_name, groups(:))
      ELSE
        var%is_in_output = is_variable_in_output(var%full_name)
      END IF
    END IF
    var%is_in_restart = .TRUE.  ! All jsbach var_lists have lrestart=.TRUE. by default
    IF (PRESENT(lrestart)) var%is_in_restart = lrestart
    ! var%is_in_restart = var%is_in_restart .OR. .TRUE.  ! TODO

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
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(3), gdims(3)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups                           ! number of additional groups this variable belongs to
    INTEGER           :: out_level, igroup                   ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsb_add_var_list_element_r3d'

    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

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

    IF (PRESENT(output_level)) THEN
      no_groups = 2 * (4 - output_level)
    ELSE
      no_groups = 0
    END IF
    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3 + no_groups))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        igroup = 4
        DO out_level=output_level,3
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//Get_process_name(mem%owner_proc_id)//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r3d( &
      this_list, TRIM(name_loc), var%ptr3d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, verbose=verbose)

    var%ptr => var%ptr3d
    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%type = REAL3D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all
    IF (loutput_loc) THEN
      IF (mem%in_var_groups) THEN
        var%is_in_output = is_variable_in_output(var%full_name, groups(:))
      ELSE
        var%is_in_output = is_variable_in_output(var%full_name)
      END IF
    END IF
    var%is_in_restart = .TRUE.  ! All jsbach var_lists have lrestart=.TRUE. by default
    IF (PRESENT(lrestart)) var%is_in_restart = lrestart
    ! var%is_in_restart = var%is_in_restart .OR. .TRUE.  ! TODO

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

    !$acc enter data copyin(var)
    !$acc update device(var%ptr3d)
    __acc_attach(var%ptr3d)
    __acc_attach(var%ptr)

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
