!>
!! @brief Contains the class and methods for the JSBACH memory structure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_varlist,       ONLY: t_jsb_varlist, Get_varlist,                         &
#ifdef __ICON__
    &                             add_var_list_element_r1d,                           &
#endif
    &                             add_var_list_element_r2d, add_var_list_element_r3d, &
    &                             Get_var_group_name, VARNAME_LEN
  USE mo_jsb_impl_constants,ONLY: SHORT_NAME_LEN
  USE mo_jsb_varlist_iface, ONLY: is_variable_in_output
  USE mo_jsb_var_class,     ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d, &
    &                             REAL1D, REAL2D, REAL3D
  USE mo_jsb_pool_class,    ONLY: t_jsb_pool
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, &
    &                             GRID_LONLAT,GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_CELL, missval, tables
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var, grib2_var
  USE mo_util,              ONLY: int2string
  USE mo_jsb_utils_iface,   ONLY: assign_if_present

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory_p, t_jsb_memory

  !> Abstract type for the memory state of each tile and process (to be extended by each process)
  TYPE, ABSTRACT :: t_jsb_memory
    INTEGER :: id = 0
    INTEGER :: grid_id = 0                                    !< grid ID of defined memory
    INTEGER :: owner_model_id = 0                             !< model ID of defined memory
    INTEGER :: owner_proc_id = 0                              !< process ID of defined memory
    CHARACTER(LEN=:), ALLOCATABLE  :: owner_proc_name         !< process name of defined memory
    INTEGER                        :: no_of_children = 0      !< counter for the number of children of the tile on which the mem lives
    TYPE(t_jsb_memory_p), POINTER  :: children(:) => NULL()   !< vector of pointers for children memory
    CLASS(t_jsb_memory) , POINTER  :: parent => NULL()        !< points for parent memory
    INTEGER, ALLOCATABLE           :: owner_tile_path(:)      !< tile location (in tile tree) of defined memory
    CHARACTER(len=SHORT_NAME_LEN)  :: owner_tile_name         !< tile name of defined memory
    !
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)                 !< collection of variables in a memory
    INTEGER                        :: max_no_of_vars = 0      !< max number of variables allowed in collection
    INTEGER                        :: no_of_vars = 0          !< counter for how many variables actually are allocated for this memory
    !
    CHARACTER(LEN=:), ALLOCATABLE  :: varlist_name            !< Name of varlist for the variables in memory instance
    LOGICAL                        :: in_var_groups = .FALSE. !< Whether variables of this memory are in any var group
    !
    LOGICAL                        :: has_bgc_materials = .FALSE.
    CLASS(t_jsb_pool), POINTER     :: pools         => NULL() !< TODO? This is old pool implementation and will be modified by new updates (Jan Engle)
    CLASS(t_jsb_pool), POINTER     :: bgc_material  => NULL() !< TODO t_jsb_memory rename 't_jsb_pool' -> 't_lnd_bgc_material'
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init                          !< Initial function to be implemented in each instance of tile memory
                                                                                  !< PASS(mem) -> pass memory object as the first argument 'mem'
    PROCEDURE                                    :: Get_var
    PROCEDURE                                    :: t_jsb_memory_add_var_r1d      !< to add 1-dimension variable in a speified memory
    PROCEDURE                                    :: t_jsb_memory_add_var_r2d      !< to add 2-dimension variable in a speified memory
    PROCEDURE                                    :: t_jsb_memory_add_var_r3d      !< to add 3-dimension variable in a speified memory
    PROCEDURE                                    :: t_jsb_memory_add_pool         !< Add a (hierarchical) pool structure to memory
                                                                                  !< (sub-pools and elements)
    PROCEDURE                                    :: t_jsb_memory_add_pool_var_r2d !< Add a 2d variable to pool (non-elements)
    GENERIC                                      :: Add_var => t_jsb_memory_add_var_r1d, t_jsb_memory_add_var_r2d, &
      &                                                        t_jsb_memory_add_var_r3d, &
      &                                                        t_jsb_memory_add_pool, t_jsb_memory_add_pool_var_r2d
    PROCEDURE                                    :: Get_var_position              !< Get position of a variable in list of variables `vars`
    PROCEDURE, PUBLIC                            :: Collect_var_real2d            !< Collect 2d var from its child tiles
    PROCEDURE, PUBLIC                            :: Collect_var_real3d            !< Collect 3d var from its child tiles
    PROCEDURE, PUBLIC                            :: HandDown_var_real2d           !< Hand back 2d var to child tiles
    PROCEDURE                                    :: Average_over_children_real2d  !< Average 2d var over the child tiles
  END TYPE t_jsb_memory

  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ! Note: argument "mem" needs to have the TARGET attribute so that we can point to individual variables
  !       (of class t_jsb_var) contained in the memory structure.
  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix, lct_ids, model_id)
      IMPORT :: t_jsb_memory
      CLASS(t_jsb_memory), INTENT(inout), TARGET :: mem
      CHARACTER(LEN=*),    INTENT(in)    :: prefix, suffix
      INTEGER,             INTENT(in)    :: lct_ids(:)     !< Primary lct (1) and lcts of descendant tiles
      INTEGER,             INTENT(in)    :: model_id
    END SUBROUTINE Init_memory
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

CONTAINS

  SUBROUTINE t_jsb_memory_add_var_r1d(mem, name, var,           &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose                                       &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< memory structure
    CHARACTER(len=*),     INTENT(in)           :: name                !< name of variable
    TYPE(t_jsb_var_real1d), TARGET, INTENT(inout) :: var              !< reference to jsbach var
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< vg:  needed for patch id
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< vertical grid type
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate Forecast) related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix for variable name (indicate process name for output variable)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix for variable name (indicate tile name for output variable)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< container flag (for ICON variables with tile dimension)
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           !< initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          !< reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           !< missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< aggregate over all tiles including tiles which the variable do not exist
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< minimum level of output groups this var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< print information

! 1d list elements are not supported with ECHAM streams
#ifdef __ICON__
    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(1), gdims(1)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups                           ! number of additional groups this variable belongs to
    INTEGER           :: out_level, igroup                   ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r1d'

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

    ldims = (/1/)
    gdims = (/1/)

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
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        igroup = 4
        DO out_level=output_level,3
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    var%full_name = TRIM(name_loc)
    IF (loutput_loc) THEN
      IF (mem%in_var_groups) THEN
        var%is_in_output = is_variable_in_output(var%full_name, groups(:))
      ELSE
        var%is_in_output = is_variable_in_output(var%full_name)
      END IF
    END IF
    var%is_in_restart = .TRUE.  ! All jsbach var_lists have lrestart=.TRUE. by default, even though
                                ! the monitoring variables are not needed in restart files
    IF (PRESENT(lrestart)) var%is_in_restart = lrestart

    ! As the 1d variables are only intended as diagnostic output variables (monitoring), the corresponding
    ! stream elements only need to be added, if the variables are actually defined as output variables.
    IF (.NOT. var%is_in_output .AND. .NOT. var%is_in_restart) RETURN

    CALL add_var_list_element_r1d(                                                            &
      this_list, TRIM(name_loc), var%ptr1d, GRID_LONLAT, vgrid%ZaxisID,                       &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table,                             &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx,                                &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r,                        &
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups,                             &
      verbose=verbose)

    var%ptr => var%ptr1d
    var%name = TRIM(name)
    var%unit = TRIM(cf%units)
    var%type = REAL1D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! if we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      mem%no_of_vars = pos_of_var
    END IF

    !$acc enter data copyin(var)
    !$acc update device(var%ptr1d)
    __acc_attach(var%ptr1d)
    __acc_attach(var%ptr)

#endif

  END SUBROUTINE t_jsb_memory_add_var_r1d

  SUBROUTINE t_jsb_memory_add_var_r2d(mem, name, var,       &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer, l_conserve_quan, cons_quan_type_id,    &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose                                       &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< memory structure
  !!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                !< name of variable
    TYPE(t_jsb_var_real2d), TARGET, INTENT(inout) :: var              !< reference to jsbach var
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< vertical grid type
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate Forecast) related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix for variable name (indicate process name for output variable)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix for variable name (indicate tile name for output variable)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< container flag (for ICON variables with tile dimension)
    LOGICAL,              INTENT(in), OPTIONAL :: l_conserve_quan     !< if var carries a to-be-conserved quantity
    INTEGER,              INTENT(in), OPTIONAL :: cons_quan_type_id   !< type id of the conserved quantity (mo_jsb_cqt_class)
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           !< initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          !< reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           !< missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< aggregate over all tiles including tiles which the variable do not exist
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< minimum level of output groups this var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< print information

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups                           ! number of additional groups this variable belongs to
    INTEGER           :: out_level, igroup                   ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r2d'

    !< set output variable name
    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    !< set loutput
    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

    !< set missing value
    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    !< set CF speicification
    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    !< ldims : local dimension (dimension of chunk block)
    !< gdims : global dimension for output file (netcdf, hdf, grip ... etc) grid format
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
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        igroup = 4
        DO out_level=output_level,3
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name_loc), var%ptr2d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, &
      verbose=verbose)

    var%ptr => var%ptr2d
    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%unit = TRIM(cf%units)
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

    IF (PRESENT(l_conserve_quan)) var%is_conserved_quan = l_conserve_quan
    IF (PRESENT(cons_quan_type_id)) var%cons_quan_type_id = cons_quan_type_id

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! if we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      mem%no_of_vars = pos_of_var
    END IF

    !$acc enter data copyin(var)
    !$acc update device(var%ptr2d)
    __acc_attach(var%ptr2d)
    __acc_attach(var%ptr)

  END SUBROUTINE t_jsb_memory_add_var_r2d

  SUBROUTINE t_jsb_memory_add_var_r3d(mem, name, var,        &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,              &
    extra_dim_n, loutput, lcontainer,                            &
    l_conserve_quan, cons_quan_type_id,                          &
    lrestart, lrestart_cont, initval_r, isteptype,               &
    resetval_r, lmiss, missval_r, l_aggregate_all,               &
    output_level, verbose                                        &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< memory structure
    CHARACTER(len=*),     INTENT(in)           :: name                !< name of variable
    TYPE(t_jsb_var_real3d), INTENT(inout), TARGET :: var              !< reference to jsbach var
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate Forecast) related metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 related metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 related metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix for variable name (indicate process name for output variable)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix for variable name (indicate tile name for output variable)
    INTEGER,              INTENT(in), OPTIONAL :: extra_dim_n         !< Extra (third) dimension (not levels)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< container flag (for ICON variables with tile dimension)
    LOGICAL,              INTENT(in), OPTIONAL :: l_conserve_quan     !< if var carries a to-be-conserved quantity
    INTEGER,              INTENT(in), OPTIONAL :: cons_quan_type_id   !< type id of the conserved quantity (mo_jsb_cqt_class)
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           !< initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          !< reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           !< missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< aggregate over all tiles including tiles which the variable do not exist
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< minimum level of output groups this var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< print information

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


    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r3d'

    !> set output variable name
    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    !> set loutput
    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

    !> set missing value
    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    !> set CF specification
    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    
    
    ldims(1) = hgrid%nproma           !< local dimension (dimension of chunk block)
    gdims(1) = hgrid%dims_g(1)        !< global dimension for output file (netcdf, hdf, grip ... etc) grid format
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
      ldims(2) = vgrid%n_levels    !< number of levels
      ldims(3) = hgrid%nblks
      gdims(2) = vgrid%n_levels
      gdims(3) = hgrid%dims_g(2)
    END IF

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    !> If output_level is set determine number of output groups for this variable; 
    !> global max output level is currently 3 [defined as enum in mo_jsb_varlist]
    !> TODO: Hardcoded "4" in this block -> defined "invalid" output level as variable?
    IF (PRESENT(output_level)) THEN
      no_groups = 2 * (4 - output_level)    
    ELSE
      no_groups = 0
    END IF
    IF (mem%in_var_groups) THEN                                  !< Check if this memory is part of any variable group
      ALLOCATE(groups(3 + no_groups))
      !> By default all variables are put into three 'all' groups 
      groups(1) = 'ALL'                                          !< HOST model
      groups(2) = 'jsb_all'                                      !< land model
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'            !< land process
      IF (PRESENT(output_level)) THEN                            !< Variable specific setting of output level
        igroup = 4
        !> Any variable that occures in 'basic' output level has to be also present in 'medium' and 'full'
        DO out_level=output_level,3
          out_level_name = Get_var_group_name(out_level)                      !< resolve level code through output group registry
          groups(igroup)   = 'jsb_all_'//out_level_name                       !< land model
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name !< land process
          igroup = igroup + 2
        END DO
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'                                          !< only present in HOST model
    END IF

    !> Add 3D variable as element to the list
    !> Use the generic routines of HOST model in mo_jsb_varlist 
    CALL add_var_list_element_r3d( &
      this_list, TRIM(name_loc), var%ptr3d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, verbose=verbose)

    !> JSB specific variables have to be set explicitly
    var%ptr => var%ptr3d
    var%name = TRIM(name)
    var%full_name = TRIM(name_loc)
    var%unit = TRIM(cf%units)
    var%type = REAL3D
    var%len_dim3 = vgrid%n_levels
    var%missval = missval_r_loc
    
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all !< Aggregation is optional
    !> Test if variable is requested in output (defined in output namelist) else it is not put out
    IF (loutput_loc) THEN
      IF (mem%in_var_groups) THEN
        var%is_in_output = is_variable_in_output(var%full_name, groups(:))
      ELSE
        var%is_in_output = is_variable_in_output(var%full_name)
      END IF
    END IF

    var%is_in_restart = .TRUE.                                  !< All jsbach var_lists are in restart by default
    IF (PRESENT(lrestart)) var%is_in_restart = lrestart         !< explicitly exclude variable list from restart
    ! var%is_in_restart = var%is_in_restart .OR. .TRUE.  ! TODO

    IF (PRESENT(l_conserve_quan)) var%is_conserved_quan = l_conserve_quan       !< Set if quantity is conserved
    IF (PRESENT(cons_quan_type_id)) var%cons_quan_type_id = cons_quan_type_id   !< Set type of conservation

    var%owner_model_id = mem%owner_model_id                     !< model that owns variable is identical owner of memory
    var%owner_proc_id  = mem%owner_proc_id                      !< process that owns variable is identical to owner of memory

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))    !< allocate tile_path (structure in memory to reach tile) 
    var%owner_tile_path(:) = mem%owner_tile_path(:)             !< tile_path of variable is identical to memory

    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! if we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      mem%no_of_vars = pos_of_var
    END IF

    !> Send the data to the CPU
    !$acc enter data copyin(var)
    !$acc update device(var%ptr3d)
    __acc_attach(var%ptr3d)
    __acc_attach(var%ptr)

  END SUBROUTINE t_jsb_memory_add_var_r3d

  RECURSIVE SUBROUTINE t_jsb_memory_add_pool(mem, pool,         &
    hgrid, vgrid_2d, vgrid_3d, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose                                       &
    )

    CLASS(t_jsb_memory),   INTENT(inout)        :: mem                 ! memory structure
  !!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CLASS(t_jsb_pool),    TARGET,  INTENT(inout) :: pool                ! reference to pool
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid_2d            ! vertical grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid_3d            ! vertical grid type
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
    ! TYPE(t_cf_var)    :: cf_loc
    TYPE(t_cf)        :: cf
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lrestart_loc
    ! REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    ! INTEGER           :: pos_of_var
    ! CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    ! CHARACTER(len=:), ALLOCATABLE :: output_level_name
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_pool'

    lrestart_loc = .FALSE.
    CALL assign_if_present(lrestart_loc, lrestart)
    lrestart_loc = lrestart_loc .AND. .NOT. ASSOCIATED(pool%pool_list)  ! Only elements on leaves are only put into restart file

    pool%owner_model_id = mem%owner_model_id

    IF (ALLOCATED(pool%element_list)) THEN
      DO i=1,SIZE(pool%element_list)
        cf = t_cf(prefix//'_'//pool%shortname//pool%element_list(i)%p%element_shortname, pool%element_unit, &
          &       pool%element_list(i)%p%element_name//' in pool '//pool%path//' on '//suffix)
        ! cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
        ! grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)
        SELECT TYPE (selected_element => pool%element_list(i)%p)
        TYPE IS (t_jsb_var_real2d)
          CALL mem%Add_var(selected_element%element_shortname, selected_element,       &
            & hgrid, vgrid_2d, cf, grib1, grib2, prefix//'_'//pool%shortname, suffix,  &
            & loutput=loutput, lcontainer=lcontainer,                                        &
            & lrestart=lrestart_loc, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,          &
            & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
            & output_level=output_level, verbose=verbose)
        TYPE IS (t_jsb_var_real3d)
          CALL mem%Add_var(selected_element%element_shortname, selected_element,       &
            & hgrid, vgrid_3d, cf, grib1, grib2, prefix//'_'//pool%shortname, suffix,  &
            & loutput=loutput, lcontainer=lcontainer,                                        &
            & lrestart=lrestart_loc, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,          &
            & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
            & output_level=output_level, verbose=verbose)
        END SELECT
      END DO
    END IF

    IF (ASSOCIATED(pool%pool_list)) THEN
      DO i=1,SIZE(pool%pool_list)
        CALL mem%Add_var(pool%pool_list(i)%p,       &
          & hgrid, vgrid_2d, vgrid_3d, grib1, grib2, prefix//'_'//pool%shortname, suffix,             &
          & loutput=loutput, lcontainer=lcontainer,                                        &
          & lrestart=lrestart, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,              &
          & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
          & output_level=output_level, verbose=verbose)
      END DO
    END IF

  END SUBROUTINE t_jsb_memory_add_pool

  SUBROUTINE t_jsb_memory_add_pool_var_r2d(mem, pool, name, var,       &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose                                       &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 ! memory structure
  !!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CLASS(t_jsb_pool),    INTENT(inout)           :: pool                ! reference to pool
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET, INTENT(inout) :: var              ! reference to jsbach var
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
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(dp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    CHARACTER(len=:), ALLOCATABLE :: output_level_name

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_pool_var_r2d'

    name_loc = pool%name // '_' // TRIM(name)
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
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3+MERGE(2, 0, PRESENT(output_level))))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (PRESENT(output_level)) THEN                            ! Variable specific setting of output level
        output_level_name = Get_var_group_name(output_level)
        groups(4) = 'jsb_all_'//output_level_name ! var group jsb_all_<level>
        groups(5) = 'jsb_'//mem%owner_proc_name//'_'//output_level_name ! var group jsb_<process>_<level>
      END IF
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

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

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! pos_of_var = mem%no_of_vars + 1
    ! IF (pos_of_var > mem%max_no_of_vars) THEN
    !   CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ! ELSE
    !   mem%vars(pos_of_var)%p => var
    !   mem%no_of_vars = pos_of_var
    ! END IF
    CALL pool%Add_var(var,dim=2)

  END SUBROUTINE t_jsb_memory_add_pool_var_r2d

  FUNCTION Get_var(this, name) RESULT(var_ptr)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name     !< variable name to search
    CLASS(t_jsb_var),        POINTER    :: var_ptr  !< return variable

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    var_ptr => NULL()
    DO i=1,this%no_of_vars ! for all variables
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN ! found?
        var_ptr => this%vars(i)%p ! assign to variable
        EXIT
      END IF
    END DO

    ! Note: Not finding the variable in this memory, i.e. returning NULL(), should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

  END FUNCTION Get_var

  ! ====================================================================================================== !
  !
  !> It returns the position from the vector of variables
  !
  FUNCTION Get_var_position(this, name) RESULT(pos)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name
    INTEGER                             :: pos

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var_position'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    pos = 0
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        pos = i
        EXIT
      END IF
    END DO

    ! Note: Not finding the variable in this memory, i.e. returning 0, should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

  END FUNCTION Get_var_position

  ! ====================================================================================================== !
  !
  !> Collect 'var' data from all the childrens and place them into an array
  !
  FUNCTION Collect_var_real2d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        ! mem_process of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children ! of the parent tile
    TYPE(t_jsb_var_real2d), INTENT(in)     :: var         ! of the parent tile: tile%mem(iproc)%p%var
    REAL(wp), POINTER                      :: var_collected(:,:)

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real2d'

    ALLOCATE(var_collected(ice-ics+1,no_children))

    DO i=1,no_children ! for each child
       var_collected(:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice,iblk) ! get the data
    END DO

  END FUNCTION Collect_var_real2d


  ! ====================================================================================================== !
  !
  !> Collect 'var' data from all the childrens and place them into an array
  !> The 3d includes a dimensions with layers.
  !
  FUNCTION Collect_var_real3d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        !< mem_process of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children !< of the parent tile
    TYPE(t_jsb_var_real3d), INTENT(in)     :: var         !< of the parent tile: tile%mem(iproc)%p%var
    REAL(wp), POINTER                      :: var_collected(:,:,:)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real3d'

    ALLOCATE(var_collected(ice - ics + 1, &
      &                    SIZE(this%children(1)%p%vars(var%child_idx(1))%p%ptr3d, DIM=2), & ! number of layers
      &                    no_children))

    DO i=1,no_children ! for each child
       var_collected(:,:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr3d(ics:ice,:, iblk) ! get the data
    END DO

  END FUNCTION Collect_var_real3d


  ! ====================================================================================================== !
  !
  !> The array 'var_handdown' places its data into the variable 'var' for each children.
  !> It is the opposite of Collect_var_read2d.
  !
  SUBROUTINE HandDown_var_real2d(this, ics, ice, iblk, no_children,  var, var_handdown)

    CLASS(t_jsb_memory), INTENT(inout)  :: this               !< mem_process of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children        !< of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var                !< of the parent tile only to find child memories
    REAL(wp),            INTENT(in)     :: var_handdown(:,:)  !< e.g. corresponds to (ics:ice,no_children)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':HandDown_var_real2d'

    DO i=1,no_children ! for each children
      this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice, iblk) = var_handdown(:,i)
    END DO

  END SUBROUTINE HandDown_var_real2d


  ! TODO: not used anywhere, ask to julia if it should be deleted
  FUNCTION Average_over_children_real2d(this, ics, ice, iblk, no_children,  var, fract_st) RESULT(var_averaged)

    CLASS(t_jsb_memory), INTENT(in)     :: this          ! mem_process of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children   ! of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var           ! of the parent tile: tile%mem(iproc)%p%var
                                                         ! R: kann man auf TYPE(t_jsb_var_real2d) setzen
    REAL(wp),            INTENT(in)     :: fract_st(:,:) ! fract_st(:,:) =tile%aggregators(1)%p%fractions(ics:ice,iblk,:)
    REAL(wp), allocatable               :: var_averaged(:)

    ! local variables
    INTEGER :: i
    REAL(wp), allocatable               :: var_collected(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':Average_over_children'

    ALLOCATE(var_collected(ice-ics+1,no_children))
    ! Collect corresponding vars of subtiles
    DO i=1,no_children
      IF ( .NOT. ASSOCIATED(this%children(i)%p%vars(var%child_idx(i))%p) ) THEN
        CALL finish(TRIM(routine), 'This var does not exist on the child tile '//TRIM(int2string(i))//' .')
      END IF
      var_collected(:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice, iblk)
    END DO
    var_averaged(:) = SUM(var_collected(:,:) * fract_st(:,:) , DIM=2) / SUM(fract_st(:,:), DIM=2)
    DEALLOCATE(var_collected)

  END FUNCTION Average_over_children_real2d

#endif
END MODULE mo_jsb_memory_class
