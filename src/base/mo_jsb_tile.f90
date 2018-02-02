!>
!! @brief Contains concrete tile class that contains the surface structure and memory state.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! This had to be separated from mo_jsb_tile_class in order to avoid circular dependencies.
!!
MODULE mo_jsb_tile
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator_p, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process_p, &
    &                               SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  USE mo_jsb_process_factory, ONLY: max_no_of_processes
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task !, t_jsb_task_msg
  USE mo_jsb_config_class,    ONLY: t_jsb_config, t_jsb_config_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile
    LOGICAL :: must_create_memory_parents(max_no_of_processes)    !< helper variable
    TYPE(t_jsb_model),         POINTER :: model => NULL()
    LOGICAL :: in_var_groups
  CONTAINS
    PROCEDURE :: Init      => Init_tile
    PROCEDURE :: Init_vars => Init_vars_on_tile
    PROCEDURE :: Print     => Print_tile
    PROCEDURE :: Handler   => Process_task_on_tile
  END TYPE t_jsb_tile

  TYPE t_jsb_tile_p
    TYPE(t_jsb_tile), POINTER :: p
  END TYPE t_jsb_tile_p

  ! Derived-type constructor
  INTERFACE t_jsb_tile
    PROCEDURE Construct_tile
  END INTERFACE t_jsb_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile'

CONTAINS

  FUNCTION Construct_tile(name, description, processes, lct, &
    & process_actions, parent, model_id, grid_id, fract_filename, fract_varname) RESULT(return_ptr)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_util,                ONLY: one_of
    USE mo_jsb_lct_class,       ONLY: t_jsb_lct, LAND_TYPE, LAKE_TYPE, Contains_lct
    USE mo_jsb_process_class,   ONLY: SEB_, A2L_, L2A_, Get_process_name

    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in), OPTIONAL :: processes(:)
    TYPE(t_jsb_lct),   INTENT(in), OPTIONAL :: lct
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id       !< Only needed if parent is not present
    INTEGER,           INTENT(in), OPTIONAL :: grid_id        !< Only needed if parent is not present
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_filename !< Name of file with tile fractions
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_varname  !< Name of fraction variable in fract_filename
    CLASS(t_jsb_tile), POINTER              :: return_ptr

    TYPE(t_jsb_model),            POINTER :: model
    TYPE(t_jsb_config_p),         POINTER :: configs(:)
    CHARACTER(LEN=15)                     :: process_name
    INTEGER                               :: i, j, iproc
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
    CLASS(t_State),               POINTER :: parent_tmp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Construct_tile'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting construction of tile '//TRIM(name))

    ! Allocate and initialize base structure (=state) for tile
    ALLOCATE(t_jsb_tile::return_ptr)
    IF (PRESENT(parent)) THEN
      parent_tmp => parent
      CALL return_ptr%Init_state(name, description, parent_tmp)
    ELSE
      CALL return_ptr%Init_state(name, description)
    END IF

    ! Initialize further components of t_jsb_tile_abstract extended type

    ALLOCATE(return_ptr%process_action(max_no_of_processes))

    IF (PRESENT(processes) .AND. PRESENT(process_actions)) THEN
      IF (SIZE(processes) /= SIZE(process_actions)) &
        & CALL finish(TRIM(routine), 'Size of processes and process_actions must be equal!')
    END IF

    IF (PRESENT(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile)
        return_ptr%parent_tile => parent
        return_ptr%grid_id = parent%grid_id
        return_ptr%owner_model_id = parent%owner_model_id
      END SELECT
    ELSE
      IF (PRESENT(model_id)) THEN
        return_ptr%owner_model_id = model_id
      ELSE
        CALL finish(TRIM(routine), 'Need model_id')
      END IF
      IF (PRESENT(grid_id)) THEN
        return_ptr%grid_id = grid_id
      ELSE
        CALL finish(TRIM(routine), 'Need grid_id')
      END IF
    END IF

    model => Get_model(return_ptr%owner_model_id)
    return_ptr%model => model

    IF (ASSOCIATED(return_ptr%parent_tile)) THEN
      IF (.NOT. ASSOCIATED(return_ptr%parent_tile%first_child_tile)) THEN
          return_ptr%parent_tile%first_child_tile => return_ptr
      END IF
    END IF
    IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
      SELECT TYPE (prev_sibling => return_ptr%prev_sibling)
      CLASS IS (t_jsb_tile)
        prev_sibling%next_sibling_tile => return_ptr
      END SELECT
    END IF

    IF (PRESENT(lct)) THEN
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        IF (ASSOCIATED(return_ptr%parent_tile%parent_tile) .AND. (lct%id == LAND_TYPE .OR. lct%id == LAKE_TYPE)) &
          & CALL finish(TRIM(routine), 'Tiles with lct LAND_TYPE or LAKE_TYPE must be direct children of root tile.')
      END IF
      CALL return_ptr%Add_lct(lct)
      CALL return_ptr%Add_lct_to_parents(lct)
    ELSE
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        CALL return_ptr%Add_lct(return_ptr%parent_tile%lcts(1))
      END IF
    END IF

    IF (PRESENT(fract_filename)) THEN
      return_ptr%fract_filename = TRIM(fract_filename)
    ELSE
      return_ptr%fract_filename = TRIM(model%config%fract_filename)
    END IF

    IF (PRESENT(fract_varname)) THEN
      return_ptr%fract_varname = TRIM(fract_varname)
    ELSE
      return_ptr%fract_varname = ''
    END IF

    ! Loop over all processes and set actions for this tile
    return_ptr%must_create_memory_parents(:) = .FALSE.
    DO i=1,max_no_of_processes
      return_ptr%process_action(i) = SKIP_
      iproc = -1
      j = 0
      IF (PRESENT(processes)) THEN
        IF (SIZE(processes) > 0) THEN
          j = one_of(i, processes) ! i counts through the number of all processes; j is the number of the i-process within
                                   ! the explicitly defined input processes for this tile (e.g. defined in init_usecase)
        END IF
      END IF
      IF (j > 0) THEN ! the i-process was explicitly defined for this tile
        iproc = processes(j)

        IF (iproc == SEB_ .AND. ALLOCATED(return_ptr%lcts)) THEN
          IF (return_ptr%lcts(1)%id /= LAND_TYPE .AND. return_ptr%lcts(1)%id /= LAKE_TYPE) THEN
            CALL finish(TRIM(routine), 'SEB_ process only allowed on tiles with lct LAND_TYPE or LAKE_TYPE.')
          END IF
        END IF

        IF (PRESENT(process_actions)) THEN
          CALL return_ptr%Set_process_action(iproc, process_actions(j), action_to_parents=AGGREGATE_)
        ELSE          ! use defaults if process_actions are not given
          IF (iproc == A2L_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE IF (iproc == L2A_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE
            CALL return_ptr%Set_process_action(iproc, ON_LEAFS_, action_to_parents=AGGREGATE_)
          END IF
        END IF
!!$      ELSE IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
!!$        iproc = i
!!$        SELECT TYPE (sibling => return_ptr%prev_sibling)
!!$        CLASS IS (t_jsb_tile_abstract)
!!$          return_ptr%process_action(iproc) = sibling%process_action(iproc)
!!$        END SELECT
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN  ! the i-process was NOT explicitly defined for this tile
        !IF (return_ptr%parent_tile%process_action(i) > AGGREGATE_) THEN
          iproc = i
          SELECT CASE (return_ptr%parent_tile%process_action(iproc))
          CASE (SKIP_, INHERIT_, ON_SUBTREE_, ON_LEAFS_)
            CALL return_ptr%Set_process_action(iproc, return_ptr%parent_tile%process_action(iproc))
            !return_ptr%parent_tile%process_action(iproc) = AGGREGATE_
          CASE (ON_TILE_)
            IF (iproc == L2A_) THEN
              CALL return_ptr%Set_process_action(iproc, SKIP_)
            ELSE
              CALL return_ptr%Set_process_action(iproc, INHERIT_)
            END IF
          END SELECT
        !END IF
      END IF

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'Finished construction of tile '//TRIM(name))

  END FUNCTION Construct_tile

  ! Note: It is important that Init_tile is only called for the tiles after the hierarchical tile tree has
  !       been completely constructed. Therefore, Init_tile is called from init_model_common.
  !
  SUBROUTINE Init_tile(this, varlist_name, prefix, suffix, l_fixed_fractions, in_var_groups)

    USE mo_jsb_process_factory, ONLY: Create_process_memory
    USE mo_jsb_process_class,   ONLY: Get_process_name
    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_tile_class,      ONLY: t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_varlist,         ONLY: t_jsb_varlist, jsb_add_var, VARNAME_LEN
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_INSTANT, TSTEP_CONSTANT
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

    CLASS(t_jsb_tile), INTENT(inout) :: this
    CHARACTER(len=*),  INTENT(in)    :: varlist_name
    CHARACTER(len=*),  INTENT(in)    :: prefix         ! Prefix for variable names
    CHARACTER(len=*),  INTENT(in)    :: suffix
    LOGICAL,           INTENT(in)    :: l_fixed_fractions  ! Tile fractions can change or not
    LOGICAL,           INTENT(in)    :: in_var_groups

    INTEGER                               :: iproc
    CLASS(t_jsb_memory), POINTER          :: tmp_ptr
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
    CLASS(t_State), POINTER               :: tmp_child
    CHARACTER(LEN=15)                     :: process_name

    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator_weighted_by_fract
    CLASS(t_jsb_aggregator),                  POINTER :: aggregator
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    CHARACTER(len=10)          :: parent_tile_name
    INTEGER                    :: table
    TYPE(t_input_file)         :: input_file
    REAL(wp),          POINTER :: return_pointer(:,:), tmp_fract(:,:)
    CHARACTER(len=50)          :: name_loc
    CHARACTER(len=VARNAME_LEN) :: name
    INTEGER                    :: ilct, i, no_of_children, child_no

    INTEGER :: isteptype

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_tile'

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile '//TRIM(this%name))

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(this%mem(max_no_of_processes))

    ! Loop over all processes and initialize memory for this tile
    DO iproc=1,max_no_of_processes

      this%mem(iproc)%p => NULL()

      IF (this%process_action(iproc) == SKIP_) CYCLE

      IF (this%process_action(iproc) == ON_LEAFS_) THEN
        IF (ASSOCIATED(this%parent_tile)) this%parent_tile%process_action(iproc) = AGGREGATE_
        IF (.NOT. this%Has_children()) this%process_action(iproc) = ON_TILE_
      END IF

      process_name = Get_process_name(iproc)

      IF (this%process_action(iproc) == INHERIT_) THEN
        this%mem(iproc)%p => this%parent_tile%mem(iproc)%p
        this%mem(iproc)%p%in_var_groups = in_var_groups
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Linked memory for process '//TRIM(process_name)//&
          &                             ' on tile '//TRIM(this%name)//' to parent')
      ELSE
        no_of_children = this%Get_no_of_children()
        child_no = this%path(this%level)
        ! this%param(iproc)%p => t_jsb_param(iproc)
        tmp_ptr => Create_process_memory(iproc)
        IF (ASSOCIATED(this%parent_tile)) THEN
          IF (.NOT. ASSOCIATED(this%parent_tile%mem(iproc)%p%children) .AND. this%parent_tile%Get_no_of_children() > 0) &
            & ALLOCATE(this%parent_tile%mem(iproc)%p%children(this%parent_tile%Get_no_of_children()))
          this%parent_tile%mem(iproc)%p%children(child_no)%p => tmp_ptr
        END IF
        this%mem(iproc)%p => tmp_ptr
        this%mem(iproc)%p%id = iproc
        this%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
        this%mem(iproc)%p%grid_id = this%grid_id
        this%mem(iproc)%p%owner_model_id = this%owner_model_id
        this%mem(iproc)%p%owner_proc_id = iproc
        ALLOCATE(this%mem(iproc)%p%owner_tile_path(SIZE(this%path)))
        this%mem(iproc)%p%owner_tile_path(:) = this%path(:)
        this%mem(iproc)%p%owner_tile_name = TRIM(this%name)
        this%mem(iproc)%p%no_of_children = no_of_children
        this%mem(iproc)%p%in_var_groups = in_var_groups    ! Must be before CALL to %Init !
        IF (prefix /= '') THEN
          CALL this%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
        ELSE
          CALL this%mem(iproc)%p%Init(TRIM(process_name), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
        END IF
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name))

      END IF

      IF (this%must_create_memory_parents(iproc)) THEN
        child_no = this%path(this%level)
        tmp_tile => this%parent_tile
        DO WHILE (ASSOCIATED(tmp_tile))
          no_of_children = tmp_tile%Get_no_of_children()
          IF (.NOT. ASSOCIATED(tmp_tile%mem(iproc)%p)) THEN
            tmp_tile%mem(iproc)%p => Create_process_memory(iproc)
            IF (.NOT. ASSOCIATED(tmp_tile%mem(iproc)%p%children) .AND. no_of_children > 0) &
              & ALLOCATE(tmp_tile%mem(iproc)%p%children(no_of_children))
            tmp_tile%mem(iproc)%p%children(child_no)%p => tmp_ptr
            tmp_ptr => tmp_tile%mem(iproc)%p
            tmp_tile%mem(iproc)%p%id = iproc
            tmp_tile%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
            tmp_tile%mem(iproc)%p%grid_id = this%grid_id
            tmp_tile%mem(iproc)%p%owner_model_id = this%owner_model_id
            ALLOCATE(tmp_tile%mem(iproc)%p%owner_tile_path(SIZE(tmp_tile%path)))
            tmp_tile%mem(iproc)%p%owner_tile_path(:) = tmp_tile%path(:)
            tmp_tile%mem(iproc)%p%owner_tile_name = TRIM(tmp_tile%name)
            tmp_tile%mem(iproc)%p%no_of_children = no_of_children
            IF (debug_on()) &
              & CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(tmp_tile%name))
          END IF
          ! Note: Init for mem(iproc) might be called several times for the same tile and must make sure that only
          ! those parts of the memory (for certain lct%id) are created that haven't been initialized before.
          IF (prefix /= '') THEN
            CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(tmp_tile%name), &
              &                             tmp_tile%lcts(:)%id, this%owner_model_id)
          ELSE
            CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name), '', tmp_tile%lcts(:)%id, this%owner_model_id)
          END IF

          child_no = tmp_tile%path(tmp_tile%level)
          tmp_tile => tmp_tile%parent_tile
        END DO
        NULLIFY(tmp_tile)
      END IF

    END DO

    hgrid   => Get_grid(this%grid_id)
    surface => Get_vgrid('surface')

    table = tables(1)

    parent_tile_name = this%Get_parent_name()
    IF (parent_tile_name == '') parent_tile_name = 'grid box  '

    IF (l_fixed_fractions) THEN
      isteptype = TSTEP_CONSTANT
    ELSE
      isteptype = TSTEP_INSTANT
    END IF

    ! Add fraction(s) of this tile to stream
    CALL jsb_add_var(TRIM(varlist_name), 'fract', tmp_fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction', '',                                                                          &
      &      'fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                     &
      & initval_r=0.0_wp, isteptype=isteptype)
    CALL this%Set_fract_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_max', tmp_fract,                         &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_max', '',                                                                          &
      &      'max fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                      &
      & initval_r=0.0_wp, isteptype=isteptype)
    CALL this%Set_fract_max_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_old', tmp_fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_old', '',                                                                          &
      &      'old fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                      &
      & initval_r=0.0_wp, isteptype=isteptype)
    CALL this%Set_fract_old_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_max_old', tmp_fract,                         &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_max_old', '',                                                                          &
      &      'old max fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                      &
      & initval_r=0.0_wp, isteptype=isteptype)
    CALL this%Set_fract_max_old_ptr(tmp_fract)
    NULLIFY(tmp_fract)

    ALLOCATE(this%l_fract_children_changed(hgrid%nproma, hgrid%nblks))
    this%l_fract_children_changed(:,:) = .TRUE.

    ! Initialize fraction(s) of this tile from input file
    ALLOCATE(tmp_fract(hgrid%nproma, hgrid%nblks))
    IF (this%fract_varname == 'equal') THEN
      IF (ASSOCIATED(this%parent)) THEN
        tmp_fract(:,:) = 1._wp / REAL(this%parent%Get_no_of_children(), wp)
      ELSE
        tmp_fract(:,:) = 1._wp
      END IF
      CALL this%Set_fraction(fract=tmp_fract(:,:))
      tmp_fract(:,:) = 1._wp
      CALL this%Set_fraction(fract_max=tmp_fract(:,:))
    ELSEIF (this%fract_varname == 'allFirst') THEN !TODO-JN-COHORTS
      IF(ASSOCIATED(this%parent)) THEN
        IF (this%Is_first_child()) THEN
          tmp_fract(:,:) = 1._wp !- 0.005_wp * (REAL(this%parent%Get_no_of_children(), wp) - 1._wp)
        ELSE
          tmp_fract(:,:) = 0._wp !0.005_wp
        END IF
      ELSE
        tmp_fract(:,:) = 1._wp
      END IF
      CALL this%Set_fraction(fract=tmp_fract(:,:))
      tmp_fract(:,:) = 1._wp
      CALL this%Set_fraction(fract_max=tmp_fract(:,:))
    ELSE
      input_file = jsb_netcdf_open_input(TRIM(this%fract_filename), this%owner_model_id)

      IF (this%fract_varname == '') THEN
        name_loc = 'fract'
        IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
        IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)
      ELSE
        name_loc = TRIM(this%fract_varname)
      END IF
      IF (input_file%Has_var(TRIM(name_loc))) THEN
        return_pointer => input_file%Read_2d(variable_name=TRIM(name_loc))
        CALL this%Set_fraction(fract=return_pointer(:,:))
      ELSE
        CALL finish(TRIM(routine), 'Fraction '//TRIM(name_loc)//' for tile '//TRIM(this%name)// &
          &                        ' not found in '//TRIM(this%fract_filename))
      END IF

      IF (this%fract_varname == '') THEN
        name_loc = 'fract_max'
        IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
        IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)
      ELSE
        name_loc = TRIM(this%fract_varname)//'_max'
      END IF
      IF (input_file%Has_var(TRIM(name_loc))) THEN
        return_pointer => input_file%Read_2d(variable_name=TRIM(name_loc))
        CALL this%Set_fraction(fract_max=return_pointer(:,:))
      ELSE
        tmp_fract(:,:) = 1._wp
        CALL this%Set_fraction(fract_max=tmp_fract(:,:))
      END IF

      CALL input_file%Close()

    END IF

    ! If tile has children, create and register aggregator
    IF (this%Has_children()) THEN
      aggregator_weighted_by_fract => t_jsb_aggregator_weighted_by_fract(this%grid_id, this%Get_no_of_children())
      CALL this%Register_aggregator(aggregator_weighted_by_fract)
    END IF

    ! If tile has parent, put fractions of current tile into corresponding aggregator of parent
    IF (ASSOCIATED(this%parent)) THEN
      i = this%Get_pos_of_child()
      SELECT TYPE (parent => this%parent)
      CLASS IS (t_jsb_tile)
        aggregator => parent%Get_aggregator('weighted_by_fract')
        IF (ASSOCIATED(aggregator)) THEN
          SELECT TYPE (aggregator)
          TYPE IS (t_jsb_aggregator_weighted_by_fract)
            CALL this%Get_fraction(aggregator%fractions(:,:,i))
          END SELECT
        END IF
      END SELECT
    END IF

    ! Loop over lcts of this tile and add variables for fractions for each lct
    DO ilct=1,SIZE(this%lcts)

      CALL jsb_add_var(TRIM(varlist_name), 'fract_'//TRIM(this%lcts(ilct)%name), this%lcts(ilct)%fract,           &
        & hgrid, surface,                                                                                         &
        & t_cf('lct_fraction', '',                                                                                &
        &      'fraction of '//TRIM(this%name)//' covered by LCT '//TRIM(this%lcts(ilct)%name)),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                       &
        & initval_r=0.0_wp, isteptype=isteptype)

    END DO

    DEALLOCATE(tmp_fract)

  END SUBROUTINE Init_tile

  SUBROUTINE Init_vars_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_process_class,   ONLY: Get_process_name
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN

    CLASS(t_jsb_tile), INTENT(inout) :: this

    CLASS(t_jsb_memory), POINTER  :: mem
    CHARACTER(len=VARNAME_LEN)    :: name
    CHARACTER(LEN=15)             :: process_name
    INTEGER                       :: iproc, i, no_of_children, child_no

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_vars_on_tile'

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile structure on vars '//TRIM(this%name))

    ! Loop over all processes and initialize memory for this tile
    DO iproc=1,max_no_of_processes

      process_name = Get_process_name(iproc)
      ! if the current process is not on leafs of current tile exit DO loop
      SELECT CASE (this%process_action(iproc))
      CASE (SKIP_, ON_TILE_, INHERIT_)
        CYCLE
      END SELECT

      no_of_children = this%Get_no_of_children()
      mem => this%mem(iproc)%p
      DO i=1,mem%no_of_vars
        name = mem%vars(i)%p%name
        IF (no_of_children > 0) THEN
          IF (.NOT. ALLOCATED(mem%vars(i)%p%child_idx)) &
            & ALLOCATE(mem%vars(i)%p%child_idx(no_of_children))
          DO child_no=1,no_of_children
            IF (ASSOCIATED(mem%children(child_no)%p)) THEN
              mem%vars(i)%p%child_idx(child_no) = mem%children(child_no)%p%Get_var_position(TRIM(name))
            ELSE
              mem%vars(i)%p%child_idx(child_no) = 0
            END IF
          END DO
        END IF
      END DO

    END DO

  END SUBROUTINE Init_vars_on_tile

  ! Task handler function for tile (HSM message handler)
  FUNCTION Process_task_on_tile(this, msg_in) RESULT(return_ptr)

    USE mo_util, ONLY: one_of, logical2string, int2string
    USE mo_jsb_parallel, ONLY: Get_omp_thread, Get_omp_no_of_threads
    USE mo_jsb_control,  ONLY: l_timer, timer_process_task
    USE mo_timer,        ONLY: timer_start, timer_stop

    CLASS(t_jsb_tile),  INTENT(inout) :: this
    CLASS(t_Message),   INTENT(in)    :: msg_in
    CLASS(t_Message),   POINTER       :: return_ptr

    TYPE(t_jsb_model), POINTER :: model
    !CLASS(t_jsb_task_msg), POINTER :: msg
    TYPE(t_Message), POINTER :: msg
    INTEGER :: iproc
    LOGICAL :: l_debug

    CLASS(t_jsb_process_task), POINTER :: task

    INTEGER :: iblk, no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_task_on_tile'

    no_omp_thread = Get_omp_thread()

    IF (l_timer .AND. debug_on()) CALL timer_start(timer_process_task(this%owner_model_id))

    ! Gfortran crashes on the next line
    ! ALLOCATE(msg, source=msg_in)
    ALLOCATE(msg)
    msg%name = msg_in%name
    msg%action = msg_in%action

    model => this%model

    iblk = model%options(no_omp_thread)%iblk
    l_debug = debug_on('hsm') .AND. iblk == 1

    task => model%current_task(no_omp_thread)%p

    IF (l_debug) THEN
      !IF (one_of(msg_in%action%name, (/'ENTER','START','EXIT '/)) < 0) THEN
        CALL message(TRIM(routine), 'Processing task '//TRIM(task%name)//'/action '//&
          &                         TRIM(msg%action%name)//' on tile '//TRIM(this%name)//&
          &                         '(visited='//logical2string(this%visited(no_omp_thread))//')')
!!$        PRINT*,'AAA ', TRIM(routine), 'Processing task '//TRIM(task%name)//'/action '//&
!!$          &                         TRIM(msg%action%name)//' on tile '//TRIM(this%name)
!!$        CALL this%Print()
      !END IF
    ELSE IF (debug_on('detail') .AND. iblk == 1) THEN
      IF (one_of(msg_in%action%name, (/'ENTER','START','EXIT '/)) < 0) THEN
        CALL message(TRIM(routine), 'Processing task '//TRIM(task%name)//'/action '//&
          &                         TRIM(msg%action%name)//' on tile '//TRIM(this%name))
      END IF
    END IF

    SELECT CASE (TRIM(msg%action%name))
    CASE ('START')
      IF (.NOT. this%visited(no_omp_thread)) THEN
        IF (this%process_action(task%process_id) >= AGGREGATE_) THEN
          IF (this%Has_children() .AND. &
            & this%process_action(task%process_id) /= ON_SUBTREE_ .AND. &
            & this%process_action(task%process_id) /= ON_TILE_) THEN
            CALL model%Start_state(this%first_child_tile)
          END IF
        ELSE IF (ASSOCIATED(this%next_sibling_tile)) THEN
          CALL model%Start_state(this%next_sibling_tile)
        END IF
      END IF
    CASE ('INTEGRATE')
      IF (this%process_action(task%process_id) >= AGGREGATE_) THEN
        IF (this%Has_children() .AND. &
          & this%process_action(task%process_id) /= ON_SUBTREE_ .AND. &
          & this%process_action(task%process_id) /= ON_TILE_) THEN
          CALL model%Goto(this%first_child_tile, l_debug)
        ELSE
          CALL task%Do_it(msg, this, model%options(no_omp_thread))
          IF (ASSOCIATED(this%next_sibling_tile)) THEN
            !IF (this%next_sibling%Has_children()) THEN
            !  CALL model%Goto(this%next_sibling%first_child, l_debug)
            !ELSE
              CALL model%Goto(this%next_sibling_tile, l_debug)
            !END IF
          ELSE IF (ASSOCIATED(this%parent_tile)) THEN
            !IF (ASSOCIATED(msg%action)) DEALLOCATE(msg%action)
            msg%action = model%Get_action('AGGREGATE')
            CALL model%Goto(this%parent_tile, l_debug)
          ELSE
            this%visited(no_omp_thread) = .FALSE.
            IF (ASSOCIATED(msg)) DEALLOCATE(msg)
            msg => NULL()
          END IF
        END IF
      ELSE ! this%process_action(task%process_id) < AGGREGATE_ ; e.g. SKIP_ or INHERIT_
        IF (ASSOCIATED(this%next_sibling_tile)) THEN
          ! IF (this%next_sibling_tile%Has_children()) THEN
            ! CALL model%Goto(this%next_sibling_tile%first_child_tile, l_debug)
          ! ELSE
            CALL model%Goto(this%next_sibling_tile, l_debug)
          ! END IF
        ELSE IF (ASSOCIATED(this%parent_tile)) THEN
          !IF (ASSOCIATED(msg%action)) DEALLOCATE(msg%action)
          msg%action = model%Get_action('AGGREGATE')
          CALL model%Goto(this%parent_tile, l_debug)
        ELSE
          this%visited(no_omp_thread) = .FALSE.
          IF (ASSOCIATED(msg)) DEALLOCATE(msg)
          msg => NULL()
        END IF
      END IF
    CASE ('AGGREGATE')
      IF (this%process_action(task%process_id) >= AGGREGATE_ .AND. this%visited(no_omp_thread)) THEN
        CALL task%Do_it(msg, this, model%options(no_omp_thread))
      END IF
      IF (ASSOCIATED(this%next_sibling_tile)) THEN
        msg%action = model%Get_action('INTEGRATE')
        CALL model%Goto(this%next_sibling_tile, l_debug)
      ELSE IF (this == model%top) THEN
        this%visited(no_omp_thread) = .FALSE.
        IF (ASSOCIATED(msg)) DEALLOCATE(msg)
        msg => NULL()
      ELSE
        msg%action = model%Get_action('AGGREGATE')
        CALL model%Goto(this%parent_tile, l_debug)
      END IF
    END SELECT

    !return_ptr => msg
    IF (ASSOCIATED(msg)) THEN
      ALLOCATE(return_ptr, source=msg)
!      return_ptr%action = msg%action
    ELSE
      return_ptr => NULL()
    END IF

    IF (l_debug) THEN
      IF (ASSOCIATED(msg)) THEN
        IF (Get_omp_no_of_threads() > 1) THEN
          CALL message(TRIM(routine),'Finished (new message action is '//TRIM(msg%action%name)//')'//&
            &                                               ', thread '//TRIM(int2string(no_omp_thread)))
        ELSE
          CALL message(TRIM(routine),'Finished (new message action is '//TRIM(msg%action%name)//')')
        END IF
      ELSE
        IF (Get_omp_no_of_threads() > 1) THEN
          CALL message(TRIM(routine),'Finished (new message is NULL)'//', thread '//TRIM(int2string(no_omp_thread)))
        ELSE
          CALL message(TRIM(routine),'Finished (new message is NULL)')
        END IF
      END IF
    ELSE IF (debug_on('detail') .AND. iblk == 1) THEN
      IF (one_of(msg_in%action%name, (/'ENTER','START','EXIT '/)) < 0) THEN
        IF (Get_omp_no_of_threads() > 1) THEN
          CALL message(TRIM(routine),'Finished, thread '//TRIM(int2string(no_omp_thread)))
        ELSE
          CALL message(TRIM(routine),'Finished')
        END IF
      END IF
    END IF

    IF (ASSOCIATED(msg)) DEALLOCATE(msg)

    IF (l_timer .AND. debug_on()) CALL timer_stop(timer_process_task(this%owner_model_id))

  END FUNCTION Process_task_on_tile

  SUBROUTINE Print_tile(this)

    USE mo_jsb_process_class,  ONLY: Get_process_name, Get_action_name

    CLASS(t_jsb_tile), INTENT(in) :: this

    INTEGER :: ilct, iproc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Print_tile'

    ! Print for basic type
    CALL this%Print_state()

    IF (this%fract_filename /= '') THEN
      CALL message('     fract_filename', TRIM(this%fract_filename))
    END IF

    ! Additional prints for tile type
    IF (ALLOCATED(this%lcts)) THEN
      DO ilct=1,SIZE(this%lcts)
        IF (this%lcts(ilct)%lib_id > 0) THEN
          WRITE(message_text,*) TRIM(this%lcts(ilct)%Get_name()),' ... ',TRIM(this%lcts(ilct)%Get_longname()), &
            &                   ' (',this%lcts(ilct)%lib_id,')'
        ELSE
          WRITE(message_text,*) TRIM(this%lcts(ilct)%Get_name())//' ... '//TRIM(this%lcts(ilct)%Get_longname())
        END IF
        CALL message('     LCT', message_text)
      END DO
    END IF

    DO iproc=1,max_no_of_processes
      IF (this%process_action(iproc) > 0) THEN
        CALL message('     Process '//Get_process_name(iproc), 'Action '//Get_action_name(this%process_action(iproc)))
      END IF
    END DO

  END SUBROUTINE Print_tile

#endif
END MODULE mo_jsb_tile
