!> Contains the class and methods for the JSBACH memory structure
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: lp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_varlist,       ONLY: t_jsb_varlist, Get_varlist,                         &
    &                             add_var_list_element_r1d,                           &
    &                             add_var_list_element_r2d, add_var_list_element_r3d, &
    &                             Get_var_group_name, VARNAME_LEN
  USE mo_jsb_impl_constants,ONLY: SHORT_NAME_LEN
  USE mo_jsb_varlist_iface, ONLY: is_variable_in_output
  USE mo_jsb_var_class,     ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d, &
    &                             REAL1D, REAL2D, REAL3D
  USE mo_lnd_bgcm_class,    ONLY: t_lnd_bgcm
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2, &
    &                             GRID_LONLAT, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, GRID_CELL, missval, tables
  USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var, grib2_var
  USE mo_util,              ONLY: int2string
  USE mo_fortran_tools,     ONLY: assign_if_present

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory_p, t_jsb_memory

  ! ======================================================================================================= !
  !>
  !>Abstract type for the memory state of each tile and process (to be extended by each process)
  !>
  TYPE, ABSTRACT :: t_jsb_memory
    INTEGER :: id = 0                                         !< ID of the memory (= ID of the owning process)
    INTEGER :: grid_id = 0                                    !< ID of the horizontal grid of the memory
    INTEGER :: owner_model_id = 0                             !< ID of the model owning the memory
    INTEGER :: owner_proc_id = 0                              !< ID of the process owning the memory
    CHARACTER(LEN=:), ALLOCATABLE  :: owner_proc_name         !< Name of the process owning the memory
    INTEGER                        :: no_of_children = 0      !< Number of children of the tile the memory lives on
    TYPE(t_jsb_memory_p), POINTER  :: children(:) => NULL()   !< Memories of this process on the child tiles
    CLASS(t_jsb_memory) , POINTER  :: parent => NULL()        !< Memory of this process on the parent tile
    INTEGER, ALLOCATABLE           :: owner_tile_path(:)      !< Path of the owning tile in the tile tree
    CHARACTER(len=SHORT_NAME_LEN)  :: owner_tile_name         !< Name of the tile owning the memory
    !
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)                 !< collection of variables in a memory
    INTEGER                        :: max_no_of_vars = 0      !< max number of variables allowed in collection
    INTEGER                        :: no_of_vars = 0          !< Number of variables actually added to this memory
    !
    CHARACTER(LEN=:), ALLOCATABLE  :: varlist_name            !< Name of varlist for the variables in memory instance
    LOGICAL                        :: in_var_groups = .FALSE. !< Whether variables of this memory are in any var group
    LOGICAL                        :: lrestart_cont = .FALSE. !< TRUE: Continue although variables are missing in restart file
    !
    LOGICAL :: has_bgc_materials = .FALSE.                    !< Whether the memory contains bgc materials (QUINCY)
    CLASS(t_lnd_bgcm), POINTER :: bgc_material  => NULL()     !< bgc material
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init  !< Initialize the memory (implemented by each process memory)
    PROCEDURE         :: Get_var                       !< Get pointer to a variable of the memory by its name
    PROCEDURE         :: t_jsb_memory_add_var_r1d      !< Add a 1d (monitoring) variable to memory
    PROCEDURE         :: t_jsb_memory_add_var_r2d      !< Add a 2d variable to memory
    PROCEDURE         :: t_jsb_memory_add_var_r3d      !< Add a 3d variable to memory
    PROCEDURE         :: t_jsb_memory_add_bgcm         !< Add a (hierarchical) bgcm structure to memory (sub-bgcm and elements)
    GENERIC           :: Add_var => t_jsb_memory_add_var_r1d, t_jsb_memory_add_var_r2d, &
      &                             t_jsb_memory_add_var_r3d, &
      &                             t_jsb_memory_add_bgcm
    PROCEDURE         :: Get_var_position              !< Get position of a variable in the list of variables `vars`
    PROCEDURE, PUBLIC :: Collect_var_real2d            !< Collect 2d var from its child tiles
    PROCEDURE, PUBLIC :: Collect_var_real3d            !< Collect 3d var from its child tiles
    PROCEDURE, PUBLIC :: HandDown_var_real2d           !< Hand back 2d var to child tiles
    PROCEDURE         :: Average_over_children_real2d  !< Average 2d var over the child tiles
  END TYPE t_jsb_memory

  !> Pointer to a process memory (e.g. for arrays of memories)
  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ! Note: argument "mem" needs to have the TARGET attribute so that we can point to individual variables
  !       (of class t_jsb_var) contained in the memory structure.
  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
      IMPORT :: t_jsb_memory
      CLASS(t_jsb_memory), INTENT(inout), TARGET :: mem
      CHARACTER(LEN=*),    INTENT(in)    :: prefix, suffix
      INTEGER,             INTENT(in)    :: lct_ids(:)     !< Primary lct (1) and lcts of descendant tiles
      INTEGER,             INTENT(in)    :: lib_id         !< id of primary lct in lctlib
      INTEGER,             INTENT(in)    :: model_id
    END SUBROUTINE Init_memory
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Add a 1d variable to the memory
  !>
  !> 1d variables hold global values (e.g. for monitoring). Their var_list element in the host model is only
  !> created if the variable is in the output or in the restart file.
  !> Called through the generic `Add_var`.
  !
  SUBROUTINE t_jsb_memory_add_var_r1d(mem, name, var,           &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, output_groups, verbose                        &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< Memory structure
    CHARACTER(len=*),     INTENT(in)           :: name                !< Name of the variable (without prefix and suffix)
    TYPE(t_jsb_var_real1d), TARGET, INTENT(inout) :: var              !< The variable to be added
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< Horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< Vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate and Forecast) metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix of the output variable name (process name)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix of the output variable name (tile name)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< Output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< Container flag (host variable with tile dimension)
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< Restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< Continue if variable is missing in restart file
    REAL(lp),             INTENT(in), OPTIONAL :: initval_r           !< Initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< Type of statistical processing
    REAL(lp),             INTENT(in), OPTIONAL :: resetval_r          !< Reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< Missing value flag
    REAL(lp),             INTENT(in), OPTIONAL :: missval_r           !< Missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< Aggregation: count child tiles lacking the var
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< Output level; var is in groups of this and higher levels
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: output_groups(:)    !< Additional output groups the var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< Print information

! 1d list elements are not supported with ECHAM streams
    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(lp)          :: missval_r_loc
    LOGICAL           :: lrestart_cont_loc
    INTEGER           :: ldims(1), gdims(1)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups, no_add_groups            ! number of standard/additional groups this variable belongs to
    INTEGER           :: out_level, igroup, add_group        ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r1d'

    ! Output variable name: <prefix>_<name>_<suffix>
    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    ! Default for lrestart_cont from process configuration
    lrestart_cont_loc = mem%lrestart_cont
    IF (PRESENT(lrestart_cont)) lrestart_cont_loc = lrestart_cont

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ldims = (/1/)
    gdims = (/1/)

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    IF (PRESENT(output_level)) THEN
      IF (output_level > 0) THEN                   ! Output level is not NONE (i.e. variable is assigned to a group)
        ! Number of standard output groups: 2 groups per output_level:
        ! jsb_all_<output_level> and jsb_<process>_<output_level>
        ! e.g. 8 groups for MINIMUM(output_level=1): 2xMINIMUM + 2xBASIC + 2xMEDIUM + 2xFULL
        !      2 groups for FULL   (output_level=4): 2xFULL
        no_groups = 2 * (5 - output_level)
      ELSE
        no_groups = 0
      END IF
    ELSE
      no_groups = 0
    END IF

    IF (PRESENT(output_groups)) THEN
      no_add_groups = SIZE(output_groups)
    ELSE
      no_add_groups = 0
    END IF

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3 + no_groups + no_add_groups))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (no_groups > 0) THEN                                    ! Variable specific setting of output level
        igroup = 4
        ! Variables with lower output level are included in higher level output groups: Variables with an output
        ! level of e.g. BASIC are included in "group:jsb_all_medium" and "group:jsb_all_full" and in the corresponding
        ! process groups.
        DO out_level=output_level,4
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
      DO add_group = 1, no_add_groups
        igroup = 3 + no_groups + add_group
        groups(igroup)  = 'jsb_'//TRIM(output_groups(add_group))  ! e.g. group jsb_mon for monitoring output
      END DO
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
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont_loc, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r,                        &
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups,                             &
      verbose=verbose)

    var%ptr => var%ptr1d
    var%name = TRIM(name)
    var%standard_name = TRIM(cf%standard_name)
    var%unit = TRIM(cf%units)
    var%type = REAL1D
    var%missval = missval_r_loc
    IF (PRESENT(l_aggregate_all)) var%l_aggregate_all = l_aggregate_all

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! If we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      __acc_attach(mem%vars(pos_of_var)%p)
      mem%no_of_vars = pos_of_var
    END IF

    !$ACC ENTER DATA COPYIN(var)
    !$ACC UPDATE DEVICE(var%ptr1d)
    __acc_attach(var%ptr1d)
    __acc_attach(var%ptr)

  END SUBROUTINE t_jsb_memory_add_var_r1d

  ! ====================================================================================================== !
  !
  !> Add a 2d variable to the memory
  !>
  !> Adds the variable to the var_list of the memory in the host model (via `mo_jsb_varlist`), determines
  !> the output groups the variable belongs to, sets the ICON-Land specific variable attributes and adds the
  !> variable to the list of variables `vars` of the memory.
  !> Called through the generic `Add_var`.
  !
  SUBROUTINE t_jsb_memory_add_var_r2d(mem, name, var,           &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,             &
    loutput, lcontainer, cons_quan_type_id, transfer_all,       &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, output_groups, verbose                        &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< Memory structure
  !!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                !< Name of the variable (without prefix and suffix)
    TYPE(t_jsb_var_real2d), TARGET, INTENT(inout) :: var              !< The variable to be added
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< Horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< Vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate and Forecast) metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix of the output variable name (process name)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix of the output variable name (tile name)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< Output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< Container flag (host variable with tile dimension)
    INTEGER,              INTENT(in), OPTIONAL :: cons_quan_type_id   !< Conserved quantity type for lcc (`mo_jsb_cqt_class`)
    LOGICAL,              INTENT(in), OPTIONAL :: transfer_all        !< transfer all or only a proportion of matter upon lcc
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< Restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< Continue if variable is missing in restart file
    REAL(lp),             INTENT(in), OPTIONAL :: initval_r           !< Initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< Type of statistical processing
    REAL(lp),             INTENT(in), OPTIONAL :: resetval_r          !< Reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< Missing value flag
    REAL(lp),             INTENT(in), OPTIONAL :: missval_r           !< Missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< Aggregation: count child tiles lacking the var
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< Output level; var is in groups of this and higher levels
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: output_groups(:)    !< Additional output groups the var belongs to
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< Print information

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(lp)          :: missval_r_loc
    LOGICAL           :: lrestart_cont_loc
    INTEGER           :: ldims(2), gdims(2)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups, no_add_groups            ! number of standard/additional groups this variable belongs to
    INTEGER           :: out_level, igroup, add_group        ! loop variables
    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r2d'

    ! Output variable name: <prefix>_<name>_<suffix>
    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    ! Default for lrestart_cont from process configuration
    lrestart_cont_loc = mem%lrestart_cont
    IF (PRESENT(lrestart_cont)) lrestart_cont_loc = lrestart_cont

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ! ldims: local dimensions (on the domain of this process), gdims: global dimensions
    ldims(1) = hgrid%nproma
    ldims(2) = hgrid%nblks
    gdims(1:2) = hgrid%dims_g(1:2)

    this_list => Get_varlist(TRIM(mem%varlist_name), hgrid%host_patch_id)

    IF (PRESENT(output_level)) THEN
      IF (output_level > 0) THEN                   ! Output level is not NONE (i.e. variable is assigned to a group)
        ! Number of standard output groups: 2 groups per output_level:
        ! jsb_all_<output_level> and jsb_<process>_<output_level>
        ! e.g. 8 groups for MINIMUM(output_level=1): 2xMINIMUM + 2xBASIC + 2xMEDIUM + 2xFULL
        !      2 groups for FULL   (output_level=4): 2xFULL
        no_groups = 2 * (5 - output_level)
      ELSE
        no_groups = 0
      END IF
    ELSE
      no_groups = 0
    END IF

    IF (PRESENT(output_groups)) THEN
      no_add_groups = SIZE(output_groups)
    ELSE
      no_add_groups = 0
    END IF

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3 + no_groups + no_add_groups))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (no_groups > 0) THEN                                    ! Variable specific setting of output level
        igroup = 4
        ! Variables with lower output level are included in higher level output groups: Variables with an output
        ! level of e.g. BASIC are included in "group:jsb_all_medium" and "group:jsb_all_full" and in the corresponding
        ! process groups.
        DO out_level=output_level,4
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
      DO add_group = 1, no_add_groups
        igroup = 3 + no_groups + add_group
        groups(igroup)  = 'jsb_'//TRIM(output_groups(add_group))  ! var group jsb_<group>
      END DO
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r2d( &
      this_list, TRIM(name_loc), var%ptr2d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%parameter, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont_loc, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, &
      verbose=verbose)

    var%ptr => var%ptr2d
    var%name = TRIM(name)
    var%standard_name = TRIM(cf%standard_name)
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

    IF (PRESENT(cons_quan_type_id)) THEN
      var%is_conserved_quan = .TRUE.
      var%cons_quan_type_id = cons_quan_type_id
    END IF
    IF (PRESENT(transfer_all)) var%transfer_all = transfer_all

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! If we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      __acc_attach(mem%vars(pos_of_var)%p)
      mem%no_of_vars = pos_of_var
    END IF

    !$ACC ENTER DATA COPYIN(var)
    !$ACC UPDATE DEVICE(var%ptr2d)
    __acc_attach(var%ptr2d)
    __acc_attach(var%ptr)

  END SUBROUTINE t_jsb_memory_add_var_r2d

  ! ====================================================================================================== !
  !
  !> Add a 3d variable to the memory
  !>
  !> Same as [[t_jsb_memory_add_var_r2d]], but for variables with a vertical dimension or, if `extra_dim_n`
  !> is given, with an extra third dimension (not for output).
  !> Called through the generic `Add_var`.
  !
  ! TODO: Can we remove the "extra_dim_n" option -- it does not seem to be used and complicates the code
  SUBROUTINE t_jsb_memory_add_var_r3d(mem, name, var,            &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix,              &
    extra_dim_n, loutput, lcontainer,                            &
    cons_quan_type_id, transfer_all,                             &
    lrestart, lrestart_cont, initval_r, isteptype,               &
    resetval_r, lmiss, missval_r, l_aggregate_all,               &
    output_level, output_groups, verbose                         &
    )

    CLASS(t_jsb_memory),  INTENT(inout)        :: mem                 !< Memory structure
    CHARACTER(len=*),     INTENT(in)           :: name                !< Name of the variable (without prefix and suffix)
    TYPE(t_jsb_var_real3d), INTENT(inout), TARGET :: var              !< The variable to be added
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< Horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               !< Vertical grid
    TYPE(t_cf),           INTENT(in)           :: cf                  !< CF (Climate and Forecast) metadata
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix of the output variable name (process name)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix of the output variable name (tile name)
    INTEGER,              INTENT(in), OPTIONAL :: extra_dim_n         !< Extra (third) dimension (not levels)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< Output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< Container flag (host variable with tile dimension)
    INTEGER,              INTENT(in), OPTIONAL :: cons_quan_type_id   !< Conserved quantity type for lcc (`mo_jsb_cqt_class`)
    LOGICAL,              INTENT(in), OPTIONAL :: transfer_all        !< transfer all or only a proportion of matter upon lcc
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< Restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< Continue if variable is missing in restart file
    REAL(lp),             INTENT(in), OPTIONAL :: initval_r           !< Initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< Type of statistical processing
    REAL(lp),             INTENT(in), OPTIONAL :: resetval_r          !< Reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< Missing value flag
    REAL(lp),             INTENT(in), OPTIONAL :: missval_r           !< Missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< Aggregation: count child tiles lacking the var
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< Output level; var is in groups of this and higher levels
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: output_groups(:)    !< Additional output groups the var belongs to

    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< Print information

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    LOGICAL           :: loutput_loc
    TYPE(t_cf_var)    :: cf_loc
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lmiss_loc
    REAL(lp)          :: missval_r_loc
    LOGICAL           :: lrestart_cont_loc
    INTEGER           :: ldims(3), gdims(3)
    INTEGER           :: pos_of_var
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    INTEGER           :: no_groups, no_add_groups            ! number of standard/additional groups this variable belongs to
    INTEGER           :: out_level, igroup, add_group        ! loop variables

    CHARACTER(len=:), ALLOCATABLE :: out_level_name          ! loop variable


    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_var_r3d'

    ! Output variable name: <prefix>_<name>_<suffix>
    name_loc = TRIM(name)
    IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
    IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)

    loutput_loc = .TRUE.
    IF (PRESENT(loutput)) loutput_loc = loutput

    lmiss_loc = .TRUE.
    IF (PRESENT(lmiss)) lmiss_loc = lmiss
    missval_r_loc = missval
    IF (PRESENT(missval_r)) missval_r_loc = missval_r

    ! Default for lrestart_cont from process configuration
    lrestart_cont_loc = mem%lrestart_cont
    IF (PRESENT(lrestart_cont)) lrestart_cont_loc = lrestart_cont

    cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
    grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)

    ! ldims: local dimensions (on the domain of this process), gdims: global dimensions
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
      IF (output_level > 0) THEN                   ! Output level is not NONE (i.e. variable is assigned to a group)
        ! Number of standard output groups: 2 groups per output_level:
        ! jsb_all_<output_level> and jsb_<process>_<output_level>
        ! e.g. 8 groups for MINIMUM(output_level=1): 2xMINIMUM + 2xBASIC + 2xMEDIUM + 2xFULL
        !      2 groups for FULL   (output_level=4): 2xFULL
        no_groups = 2 * (5 - output_level)
      ELSE
        no_groups = 0
      END IF
    ELSE
      no_groups = 0
    END IF

    IF (PRESENT(output_groups)) THEN
      no_add_groups = SIZE(output_groups)
    ELSE
      no_add_groups = 0
    END IF

    IF (mem%in_var_groups) THEN                                  ! Should this tile memory be part of any var groups?
      ALLOCATE(groups(3 + no_groups + no_add_groups))
      ! All variables are put into these three catch-all groups
      groups(1) = 'ALL'
      groups(2) = 'jsb_all'
      ! var group jsb_<process>_all
      groups(3) = 'jsb_'//mem%owner_proc_name//'_all'
      IF (no_groups > 0) THEN                                    ! Variable specific setting of output level
        igroup = 4
        ! Variables with lower output level are included in higher level output groups: Variables with an output
        ! level of e.g. BASIC are included in "group:jsb_all_medium" and "group:jsb_all_full" and in the corresponding
        ! process groups.
        DO out_level=output_level,4
          out_level_name = Get_var_group_name(out_level)
          groups(igroup)   = 'jsb_all_'//out_level_name ! var group jsb_all_<level>
          groups(igroup+1) = 'jsb_'//mem%owner_proc_name//'_'//out_level_name ! var group jsb_<process>_<level>
          igroup = igroup + 2
        END DO
      END IF
      DO add_group = 1, no_add_groups
        igroup = 3 + no_groups + add_group
        groups(igroup)  = 'jsb_'//TRIM(output_groups(add_group))  ! var group jsb_<group>
      END DO
    ELSE
      ALLOCATE(groups(1))
      groups(1) = 'ALL'
    END IF

    CALL add_var_list_element_r3d( &
      this_list, TRIM(name_loc), var%ptr3d, GRID_UNSTRUCTURED_CELL, vgrid%ZaxisID, &
      cf_loc, grib2_loc, code=grib1%PARAMETER, table=grib1%table, &
      ldims=ldims, gdims=gdims, levelindx=vgrid%echamZaxisIdx, &
      loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont_loc, &
      initval_r=initval_r, isteptype=isteptype, resetval_r=resetval_r, &
      ! lmiss=lmiss_loc, missval_r=missval_r_loc)
      lmiss=lmiss_loc, missval_r=missval_r_loc, in_groups=groups, verbose=verbose)

    var%ptr => var%ptr3d
    var%name = TRIM(name)
    var%standard_name = TRIM(cf%standard_name)
    var%full_name = TRIM(name_loc)
    var%unit = TRIM(cf%units)
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

    IF (PRESENT(cons_quan_type_id)) THEN
      var%is_conserved_quan = .TRUE.
      var%cons_quan_type_id = cons_quan_type_id
    END IF
    IF (PRESENT(transfer_all)) var%transfer_all = transfer_all

    var%owner_model_id = mem%owner_model_id
    var%owner_proc_id  = mem%owner_proc_id

    ALLOCATE(var%owner_tile_path(SIZE(mem%owner_tile_path)))
    var%owner_tile_path(:) = mem%owner_tile_path(:)

    ! If we want to add a variable we need to store it in the next position in vars
    pos_of_var = mem%no_of_vars + 1
    ! Assert that we still have space for another variable
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var ! Collect the pointer to the newly added variable
      __acc_attach(mem%vars(pos_of_var)%p)
      mem%no_of_vars = pos_of_var
    END IF

    !$ACC ENTER DATA COPYIN(var)
    !$ACC UPDATE DEVICE(var%ptr3d)
    __acc_attach(var%ptr3d)
    __acc_attach(var%ptr)

  END SUBROUTINE t_jsb_memory_add_var_r3d

  ! ====================================================================================================== !
  !
  !> Add a (hierarchical) bgc material to the memory
  !>
  !> Recursively adds the element variables of the bgc material and of its sub-bgcms with `Add_var`.
  !> Only the elements on the leaves of the bgcm tree go into the restart file.
  !
  RECURSIVE SUBROUTINE t_jsb_memory_add_bgcm(mem, bgcm,         &
    hgrid, vgrid_2d, vgrid_3d, grib1, grib2, prefix, suffix,    &
    loutput, lcontainer,                                        &
    lrestart, lrestart_cont, initval_r, isteptype,              &
    resetval_r, lmiss, missval_r, l_aggregate_all,              &
    output_level, verbose,                                      &
    cons_quan_type_id, transfer_all                             &
    )

    CLASS(t_jsb_memory),   INTENT(inout)        :: mem                 !< Memory structure
  !!$    TYPE(t_jsb_varlist),  INTENT(inout)        :: this_list           ! list
    CLASS(t_lnd_bgcm),    TARGET,  INTENT(inout) :: bgcm                !< The bgc material to be added
    TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               !< Horizontal grid
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid_2d            !< Vertical grid for 2d elements
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid_3d            !< Vertical grid for 3d elements
    TYPE(t_grib1),        INTENT(in)           :: grib1               !< GRIB1 metadata
    TYPE(t_grib2),        INTENT(in)           :: grib2               !< GRIB2 metadata
    CHARACTER(len=*),     INTENT(in)           :: prefix              !< Prefix of the output variable name (process name)
    CHARACTER(len=*),     INTENT(in)           :: suffix              !< Suffix of the output variable name (tile name)
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             !< Output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          !< Container flag (host variable with tile dimension)
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            !< Restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       !< Continue if variable is missing in restart file
    REAL(lp),             INTENT(in), OPTIONAL :: initval_r           !< Initial value
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           !< Type of statistical processing
    REAL(lp),             INTENT(in), OPTIONAL :: resetval_r          !< Reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               !< Missing value flag
    REAL(lp),             INTENT(in), OPTIONAL :: missval_r           !< Missing value
    LOGICAL,              INTENT(in), OPTIONAL :: l_aggregate_all     !< Aggregation: count child tiles lacking the var
    INTEGER,              INTENT(in), OPTIONAL :: output_level        !< Output level; var is in groups of this and higher levels
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             !< Print information
    INTEGER,              INTENT(in), OPTIONAL :: cons_quan_type_id   !< Conserved quantity type for lcc (`mo_jsb_cqt_class`)
    LOGICAL,              INTENT(in), OPTIONAL :: transfer_all        !< transfer all or only a proportion of matter upon lcc

    TYPE(t_jsb_varlist), POINTER :: this_list           ! list
    CHARACTER(len=50) :: name_loc
    ! TYPE(t_cf_var)    :: cf_loc
    TYPE(t_cf)        :: cf
    TYPE(t_grib2_var) :: grib2_loc
    LOGICAL           :: lrestart_loc
    ! REAL(lp)          :: missval_r_loc
    INTEGER           :: ldims(2), gdims(2)
    ! INTEGER           :: pos_of_var
    ! CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: groups(:)     ! groups to which this variable belongs to
    ! CHARACTER(len=:), ALLOCATABLE :: output_level_name
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_memory_add_bgcm'

    lrestart_loc = .FALSE.
    CALL assign_if_present(lrestart_loc, lrestart)
    lrestart_loc = lrestart_loc .AND. .NOT. ASSOCIATED(bgcm%bgcm_list)  ! Only elements on leaves are only put into restart file

    bgcm%owner_model_id = mem%owner_model_id

    IF (ALLOCATED(bgcm%element_list)) THEN
      DO i=1,SIZE(bgcm%element_list)
        cf = t_cf(prefix//'_'//bgcm%shortname//bgcm%element_list(i)%p%element_shortname, bgcm%element_unit, &
          &       bgcm%element_list(i)%p%element_name//' in bgcm '//bgcm%path//' on '//suffix)
        ! cf_loc = t_cf_var(cf%standard_name, cf%units, cf%long_name, cf%datatype)
        ! grib2_loc = grib2_var(grib2%discipline, grib2%category, grib2%number, grib2%bits, GRID_UNSTRUCTURED, GRID_CELL)
        SELECT TYPE (selected_element => bgcm%element_list(i)%p)
        TYPE IS (t_jsb_var_real2d)
          CALL mem%Add_var(selected_element%element_shortname, selected_element,       &
            & hgrid, vgrid_2d, cf, grib1, grib2, prefix//'_'//bgcm%shortname, suffix,  &
            & loutput=loutput, lcontainer=lcontainer,                                        &
            & lrestart=lrestart_loc, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,          &
            & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
            & output_level=output_level, verbose=verbose, cons_quan_type_id=cons_quan_type_id, transfer_all=transfer_all )
        TYPE IS (t_jsb_var_real3d)
          CALL mem%Add_var(selected_element%element_shortname, selected_element,       &
            & hgrid, vgrid_3d, cf, grib1, grib2, prefix//'_'//bgcm%shortname, suffix,  &
            & loutput=loutput, lcontainer=lcontainer,                                        &
            & lrestart=lrestart_loc, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,          &
            & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
            & output_level=output_level, verbose=verbose, cons_quan_type_id=cons_quan_type_id, transfer_all=transfer_all)
        END SELECT
      END DO
    END IF

    IF (ASSOCIATED(bgcm%bgcm_list)) THEN
      DO i=1,SIZE(bgcm%bgcm_list)
        CALL mem%Add_var(bgcm%bgcm_list(i)%p,       &
          & hgrid, vgrid_2d, vgrid_3d, grib1, grib2, prefix//'_'//bgcm%shortname, suffix,             &
          & loutput=loutput, lcontainer=lcontainer,                                        &
          & lrestart=lrestart, lrestart_cont=lrestart_cont, initval_r=initval_r, isteptype=isteptype,              &
          & resetval_r=resetval_r, lmiss=lmiss, missval_r=missval_r, l_aggregate_all=l_aggregate_all,              &
          & output_level=output_level, verbose=verbose, cons_quan_type_id=cons_quan_type_id, transfer_all=transfer_all)
      END DO
    END IF

  END SUBROUTINE t_jsb_memory_add_bgcm


  ! ====================================================================================================== !
  !
  !> Get a pointer to the variable `name` of the memory
  !>
  !> Returns NULL() if the variable is not found.
  !
  FUNCTION Get_var(this, name) RESULT(var_ptr)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name     !< Name of the variable
    CLASS(t_jsb_var),        POINTER    :: var_ptr  !< Pointer to the variable

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var'

    IF (.NOT. ASSOCIATED(this%vars(1)%p)) CALL finish(TRIM(routine), 'Invalid vector of vars')

    var_ptr => NULL()
    DO i=1,this%no_of_vars
      IF (TRIM(this%vars(i)%p%name) == TRIM(name)) THEN
        var_ptr => this%vars(i)%p
        EXIT
      END IF
    END DO

    ! Note: Not finding the variable in this memory, i.e. returning NULL(), should only happen because
    !       the variable belongs to an LCT that's not handled by this tile.

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var

  ! ====================================================================================================== !
  !
  !> Get the position of the variable `name` in the list of variables `vars` of the memory
  !>
  !> Returns 0 if the variable is not found.
  !
  FUNCTION Get_var_position(this, name) RESULT(pos)

    CLASS(t_jsb_memory),     INTENT(in) :: this
    CHARACTER(len=*),        INTENT(in) :: name     !< Name of the variable
    INTEGER                             :: pos      !< Position of the variable in `vars`

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_var_position'

    ! Return in case this tile's memory is empty
    IF (this%no_of_vars == 0) RETURN

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

    !IF (.NOT. ASSOCIATED(var_ptr)) &
    !  & CALL finish(TRIM(routine), 'Variable '//TRIM(name)//' on tile '//TRIM(this%owner_tile_name)//' not found')

  END FUNCTION Get_var_position

  ! ====================================================================================================== !
  !
  !> Collect the values of 2d variable `var` from all child tiles into an array (cells, children)
  !
  FUNCTION Collect_var_real2d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        !< Process memory of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children !< Number of children of the parent tile
    TYPE(t_jsb_var_real2d), INTENT(in)     :: var         !< Variable on the parent tile
    REAL(lp), POINTER                      :: var_collected(:,:)

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real2d'

#ifdef _OPENACC
    CALL finish(routine, 'GPU port not implemented yet')
#endif

    ALLOCATE(var_collected(ice-ics+1,no_children))

    DO i=1,no_children
       var_collected(:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice,iblk)
    END DO

  END FUNCTION Collect_var_real2d


  ! ====================================================================================================== !
  !
  !> Collect the values of 3d variable `var` from all child tiles into an array (cells, levels, children)
  !
  FUNCTION Collect_var_real3d(this, ics, ice, iblk, no_children,  var) RESULT(var_collected)

    CLASS(t_jsb_memory),    INTENT(in)     :: this        !< Process memory of the parent tile
    INTEGER,                INTENT(in)     :: ics, ice, iblk
    INTEGER,                INTENT(in)     :: no_children !< Number of children of the parent tile
    TYPE(t_jsb_var_real3d), INTENT(in)     :: var         !< Variable on the parent tile
    REAL(lp), POINTER                      :: var_collected(:,:,:)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_var_real3d'

#ifdef _OPENACC
    CALL finish(routine, 'GPU port not implemented yet')
#endif

    ALLOCATE(var_collected(ice - ics + 1, &
      &                    SIZE(this%children(1)%p%vars(var%child_idx(1))%p%ptr3d, DIM=2), &
      &                    no_children))

    DO i=1,no_children
       var_collected(:,:,i) = this%children(i)%p%vars(var%child_idx(i))%p%ptr3d(ics:ice,:, iblk)
    END DO

  END FUNCTION Collect_var_real3d


  ! ====================================================================================================== !
  !
  !> Hand down the values in `var_handdown` (cells, children) to 2d variable `var` on the child tiles
  !>
  !> Opposite of [[Collect_var_real2d]].
  !
  SUBROUTINE HandDown_var_real2d(this, ics, ice, iblk, no_children,  var, var_handdown)

    CLASS(t_jsb_memory), INTENT(inout)  :: this        !< Process memory of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children !< Number of children of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var         !< Variable on the parent tile (to find the child variables)
    REAL(lp),            INTENT(in)     :: var_handdown(:,:) !< Values for the child tiles (ics:ice,no_children)

    INTEGER :: i !, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':HandDown_var_real2d'

#ifdef _OPENACC
    CALL finish(routine, 'GPU port not implemented yet')
#endif

    DO i=1,no_children
      this%children(i)%p%vars(var%child_idx(i))%p%ptr2d(ics:ice, iblk) = var_handdown(:,i)
    END DO

  END SUBROUTINE HandDown_var_real2d


  ! ====================================================================================================== !
  !
  !> Average 2d variable `var` over the child tiles, weighted with the fractions `fract_st`
  !>
  !> Note: currently not used
  !
  FUNCTION Average_over_children_real2d(this, ics, ice, iblk, no_children,  var, fract_st) RESULT(var_averaged)

    CLASS(t_jsb_memory), INTENT(in)     :: this          !< Process memory of the parent tile: tile%mem(iproc)%p
    INTEGER,             INTENT(in)     :: ics, ice, iblk
    INTEGER,             INTENT(in)     :: no_children   !< Number of children of the parent tile
    CLASS(t_jsb_var),    INTENT(in)     :: var           !< Variable on the parent tile
    REAL(lp),            INTENT(in)     :: fract_st(:,:) !< Fractions of the child tiles (ics:ice,no_children)
    REAL(lp), allocatable               :: var_averaged(:)

    ! local variables
    INTEGER :: i
    REAL(lp), allocatable               :: var_collected(:,:)
    CHARACTER(len=*), PARAMETER :: routine = modname//':Average_over_children'

#ifdef _OPENACC
    CALL finish(routine, 'GPU port not implemented yet')
#endif

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
