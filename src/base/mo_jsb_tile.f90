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
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_msg
  USE mo_jsb_config_class,    ONLY: t_jsb_config, t_jsb_config_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile
    LOGICAL :: must_create_memory_parents(max_no_of_processes)    !< helper variable
  CONTAINS
    PROCEDURE :: Init    => Init_tile
    PROCEDURE :: Print   => Print_tile
    PROCEDURE :: Handler => Process_task_on_tile
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
    USE mo_jsb_process_class,   ONLY: SEB, ATM2LAND, LAND2ATM, Get_process_name

    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in)           :: processes(:)
    TYPE(t_jsb_lct),   INTENT(in), OPTIONAL :: lct
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id       !< Only needed if parent is not present
    INTEGER,           INTENT(in), OPTIONAL :: grid_id        !< Only needed if parent is not present
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_filename !< Name of file with tile fractions
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_varname  !< Name of fraction variable in fract_filename
    CLASS(t_jsb_tile), POINTER              :: return_ptr

    TYPE(t_jsb_config_p),         POINTER :: configs(:)
    CHARACTER(LEN=15)                     :: process_name
    INTEGER                               :: i, j, iproc
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
    CLASS(t_State),               POINTER :: parent_tmp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Construct_tile'

    CALL message(TRIM(routine), 'Starting construction of tile '//TRIM(name))

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

    IF (PRESENT(lct)) THEN

      IF (PRESENT(parent)) THEN
        IF (ASSOCIATED(parent%parent) .AND. (lct%id == LAND_TYPE .OR. lct%id == LAKE_TYPE)) &
          & CALL finish(TRIM(routine), 'Tiles with lct LAND_TYPE or LAKE_TYPE must be direct children of root tile.')
      END IF

      CALL return_ptr%Add_lct(lct)
      CALL return_ptr%Add_lct_to_parents(lct)

    ELSE

      IF (PRESENT(parent)) THEN
        SELECT TYPE (parent)
        CLASS IS (t_jsb_tile)
          CALL return_ptr%Add_lct(parent%lcts(1))
        END SELECT
      END IF

    END IF

    IF (PRESENT(process_actions)) THEN
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

    IF (PRESENT(fract_filename)) THEN
      return_ptr%fract_filename = TRIM(fract_filename)
    ELSE
      return_ptr%fract_filename = 'bc_land_frac.nc'
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
      IF (SIZE(processes) > 0) THEN
        j = one_of(i, processes) ! i counts through the number of all processes, j is the number of the i-process
                                 ! within the processes vector (input) 
      ELSE
        j = 0
      END IF
      IF (j > 0) THEN
        iproc = processes(j)

        IF (iproc == SEB .AND. ALLOCATED(return_ptr%lcts)) THEN
          IF (.NOT. Contains_lct(return_ptr%lcts, LAND_TYPE) .AND. .NOT. Contains_lct(return_ptr%lcts, LAKE_TYPE)) THEN
            CALL finish(TRIM(routine), 'SEB process only allowed on tiles with lct LAND_TYPE or LAKE_TYPE.')
          END IF
        END IF

        IF (PRESENT(process_actions)) THEN
          return_ptr%process_action(iproc) = process_actions(j)
        ELSE
          IF (iproc == ATM2LAND) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE IF (iproc == LAND2ATM) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE
            return_ptr%process_action(iproc) = ON_LEAFS_
          END IF
        END IF
        ! All parents of this tile aggregate memory for this new process
        CALL return_ptr%Set_process_action_parents(iproc, AGGREGATE_)
        return_ptr%must_create_memory_parents(i) = .TRUE.
!!$      ELSE IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
!!$        iproc = i
!!$        SELECT TYPE (sibling => return_ptr%prev_sibling)
!!$        CLASS IS (t_jsb_tile_abstract)
!!$          return_ptr%process_action(iproc) = sibling%process_action(iproc)
!!$        END SELECT
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        IF (return_ptr%parent_tile%process_action(i) > AGGREGATE_) THEN
          iproc = i
          SELECT CASE (return_ptr%parent_tile%process_action(iproc))
          CASE (SKIP_)
            return_ptr%process_action(iproc) = SKIP_
          CASE (INHERIT_)
            return_ptr%process_action(iproc) = INHERIT_
          CASE (ON_SUBTREE_)
            return_ptr%process_action(iproc) = ON_SUBTREE_
          CASE (ON_LEAFS_)
            return_ptr%process_action(iproc) = ON_LEAFS_
            !return_ptr%parent_tile%process_action(iproc) = AGGREGATE_
          CASE (ON_TILE_)
            IF (iproc == LAND2ATM) THEN
              return_ptr%process_action(iproc) = SKIP_
            ELSE
              return_ptr%process_action(iproc) = INHERIT_
            END IF
          END SELECT
        END IF
      END IF

!!$      IF (iproc < 1) CYCLE
    END DO

    CALL message(TRIM(routine), 'Finished construction of tile '//TRIM(name))

  END FUNCTION Construct_tile

  ! Note: It is important that Init_tile is only called for the tiles after the hierarchical tile tree has
  !       been completely constructed. Therefore, Init_tile is called from model_init_common.
  !
  SUBROUTINE Init_tile(this, varlist_name, prefix, suffix)

    USE mo_jsb_process_factory, ONLY: Create_process_memory
    USE mo_jsb_process_class,   ONLY: SEB, ATM2LAND, LAND2ATM, Get_process_name
    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_tile_class,      ONLY: t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_varlist,         ONLY: t_jsb_varlist, Get_varlist, jsb_add_var
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input, ncfile_has_var

    CLASS(t_jsb_tile), INTENT(inout) :: this
    CHARACTER(len=*),  INTENT(in)    :: varlist_name
    CHARACTER(len=*),  INTENT(in)    :: prefix         ! Prefix for variable names
    CHARACTER(len=*),  INTENT(in)    :: suffix

    INTEGER                               :: iproc
    CLASS(t_jsb_memory), POINTER          :: tmp_ptr
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
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
    INTEGER                    :: ilct, i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_tile'

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(this%mem(max_no_of_processes))

    ! Loop over all processes and initialize memory for this tile
    DO iproc=1,max_no_of_processes

      this%mem(iproc)%p => NULL()

      IF (this%process_action(iproc) == SKIP_) CYCLE

      process_name = Get_process_name(iproc)

      IF (this%process_action(iproc) == INHERIT_) THEN
        this%mem(iproc)%p => this%parent_tile%mem(iproc)%p
        CALL message(TRIM(routine), 'Linked memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name)//' to parent')
      ELSE
        !this%param(iproc)%p => t_jsb_param(iproc)
        tmp_ptr => Create_process_memory(iproc)
        SELECT TYPE (tmp_ptr)
        CLASS IS (t_jsb_memory)
          this%mem(iproc)%p => tmp_ptr
          this%mem(iproc)%p%id = iproc
          this%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
          this%mem(iproc)%p%grid_id = this%grid_id
          this%mem(iproc)%p%owner_model_id = this%owner_model_id
          this%mem(iproc)%p%owner_proc_id = iproc
          ALLOCATE(this%mem(iproc)%p%owner_tile_path(SIZE(this%path)))
          this%mem(iproc)%p%owner_tile_path(:) = this%path(:)
          this%mem(iproc)%p%owner_tile_name = TRIM(this%name)
          IF (prefix /= '') THEN
            CALL this%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
          ELSE
            CALL this%mem(iproc)%p%Init(TRIM(process_name), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
          END IF
          CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name))
        END SELECT
        NULLIFY(tmp_ptr)
      END IF

      IF (this%must_create_memory_parents(iproc)) THEN
        tmp_tile => this%parent_tile
        DO WHILE (ASSOCIATED(tmp_tile))
          IF (.NOT. ASSOCIATED(tmp_tile%mem(iproc)%p)) THEN
            tmp_ptr => Create_process_memory(iproc)
            SELECT TYPE (tmp_ptr)
            CLASS IS (t_jsb_memory)
              tmp_tile%mem(iproc)%p => tmp_ptr
              tmp_tile%mem(iproc)%p%id = iproc
              tmp_tile%mem(iproc)%p%varlist_name = TRIM(varlist_name)//'_'//TRIM(process_name)
              tmp_tile%mem(iproc)%p%grid_id = this%grid_id
              tmp_tile%mem(iproc)%p%owner_model_id = this%owner_model_id
              ALLOCATE(tmp_tile%mem(iproc)%p%owner_tile_path(SIZE(tmp_tile%path)))
              tmp_tile%mem(iproc)%p%owner_tile_path(:) = tmp_tile%path(:)
              tmp_tile%mem(iproc)%p%owner_tile_name = TRIM(tmp_tile%name)
              CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(tmp_tile%name))
            END SELECT
            NULLIFY(tmp_ptr)
          END IF
          ! Note: Init for mem(iproc) might be called several times for the same tile and must make sure that only
          ! those parts of the memory (for certain lct%id) are created that haven't been initialized before.
          IF (prefix /= '') THEN
            CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(tmp_tile%name), &
              &                             tmp_tile%lcts(:)%id, this%owner_model_id)
          ELSE
            CALL tmp_tile%mem(iproc)%p%Init(TRIM(process_name), '', tmp_tile%lcts(:)%id, this%owner_model_id)
          END IF
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

    ! Add fraction(s) of this tile to stream
    CALL jsb_add_var(TRIM(varlist_name), 'fract', tmp_fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction', '',                                                                          &
      &      'fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)
    CALL this%Set_fract_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_max', tmp_fract,                         &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_max', '',                                                                          &
      &      'max fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)
    CALL this%Set_fract_max_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_old', tmp_fract,                                 &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_old', '',                                                                          &
      &      'old fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)
    CALL this%Set_fract_old_ptr(tmp_fract)

    CALL jsb_add_var(TRIM(varlist_name), 'fract_max_old', tmp_fract,                         &
      & hgrid, surface,                                                                                         &
      & t_cf('tile_fraction_max_old', '',                                                                          &
      &      'old max fraction of '//TRIM(parent_tile_name)//' covered by tile '//TRIM(this%name)),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp)
    CALL this%Set_fract_max_old_ptr(tmp_fract)
    NULLIFY(tmp_fract)

    ALLOCATE(this%l_fract_children_changed(hgrid%nproma, hgrid%nblks))
    this%l_fract_children_changed(:,:) = .TRUE.

    ! Initialize fraction(s) of this tile from input file
    ALLOCATE(tmp_fract(hgrid%nproma, hgrid%nblks))
    IF (this%fract_filename == '') THEN
      IF (ASSOCIATED(this%parent)) THEN
        tmp_fract(:,:) = 1._wp / REAL(this%parent%Get_no_of_children(), wp)
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
      IF (ncfile_has_var(TRIM(this%fract_filename), TRIM(name_loc))) THEN
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
      IF (ncfile_has_var(TRIM(this%fract_filename), TRIM(name_loc))) THEN
        return_pointer => input_file%Read_2d(variable_name=TRIM(name_loc))
        CALL this%Set_fraction(fract_max=return_pointer(:,:))
      ELSE
        tmp_fract(:,:) = 1._wp
        CALL this%Set_fraction(fract_max=tmp_fract(:,:))
      END IF

      CALL input_file%Close()

    END IF
    DEALLOCATE(tmp_fract)

    ! If tile has children, create and register aggregator
    IF (this%Has_children()) THEN
      aggregator_weighted_by_fract => t_jsb_aggregator_weighted_by_fract(this)
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

    ! Loop over lcts of this tile and add fractions for each lct
    DO ilct=1,SIZE(this%lcts)
      CALL jsb_add_var(TRIM(varlist_name), 'fract_'//TRIM(this%lcts(ilct)%name), this%lcts(ilct)%fract,           &
        & hgrid, surface,                                                                                         &
        & t_cf('lct_fraction', '',                                                                                &
        &      'fraction of '//TRIM(this%name)//' covered by LCT '//TRIM(this%lcts(ilct)%name)),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & initval_r=0.0_wp)
    END DO

  END SUBROUTINE Init_tile

  ! Task handler function for tile (HSM message handler)
  FUNCTION Process_task_on_tile(this, hsm, msg_in) RESULT(return_ptr)

    USE mo_util, ONLY: one_of

    CLASS(t_jsb_tile),  INTENT(inout) :: this
    CLASS(t_Hsm),       INTENT(inout) :: hsm
    CLASS(t_Message),   INTENT(in)    :: msg_in
    CLASS(t_Message),   POINTER       :: return_ptr

    CLASS(t_jsb_task_msg), POINTER :: msg
    INTEGER :: iproc

    CLASS(t_jsb_process_task), POINTER :: task

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_task_on_tile'

    ALLOCATE(msg, source=msg_in)
    !return_ptr => NULL()

    SELECT TYPE(model=>hsm)
    CLASS IS (t_jsb_model)

      task => model%current_task

!      IF (debug_on() .AND. model%options%iblk == 1) THEN
!        IF (one_of(msg_in%action%name, (/'ENTER','START','EXIT '/)) < 0) THEN
!          CALL message(TRIM(routine), 'Processing task '//TRIM(task%name)//'/action '//&
!            &                         TRIM(msg%action%name)//' on tile '//TRIM(this%name))
!!$        PRINT*,'AAA ', TRIM(routine), 'Processing task '//TRIM(task%name)//'/action '//&
!!$          &                         TRIM(msg%action%name)//' on tile '//TRIM(this%name)
!!$        CALL this%Print()
!        END IF
!      END IF

      SELECT CASE (TRIM(msg%action%name))
      CASE ('START')
        IF (.NOT. this%visited) THEN
          IF (this%process_action(task%process_id) >= AGGREGATE_) THEN
            IF (this%Has_children() .AND. &
              & this%process_action(task%process_id) /= ON_SUBTREE_ .AND. &
              & this%process_action(task%process_id) /= ON_TILE_) THEN
              CALL model%Start_state(this%first_child)
            END IF
          ELSE IF (ASSOCIATED(this%next_sibling)) THEN
            CALL model%Start_state(this%next_sibling)
          END IF
        END IF
      CASE ('INTEGRATE')
        IF (this%process_action(task%process_id) >= AGGREGATE_) THEN
          IF (this%Has_children() .AND. &
            & this%process_action(task%process_id) /= ON_SUBTREE_ .AND. &
            & this%process_action(task%process_id) /= ON_TILE_) THEN
            CALL model%GOTO(this%first_child)
          ELSE
            task%msg => msg
            CALL task%Do_it(this, model%options)
            IF (ASSOCIATED(this%next_sibling)) THEN
              IF (this%next_sibling%Has_children()) THEN
                CALL model%Goto(this%next_sibling%first_child)
              ELSE
                CALL model%Goto(this%next_sibling)
              END IF
            ELSE IF (ASSOCIATED(this%parent)) THEN
              msg%action => model%Get_action('AGGREGATE')
              CALL model%Goto(this%parent)
            ELSE
              this%visited = .FALSE.
              msg => NULL()
            END IF
          END IF
        ELSE
          IF (ASSOCIATED(this%next_sibling)) THEN
            IF (this%next_sibling%Has_children()) THEN
              CALL model%Goto(this%next_sibling%first_child)
            ELSE
              CALL model%Goto(this%next_sibling)
            END IF
          ELSE IF (ASSOCIATED(this%parent)) THEN
            msg%action => model%Get_action('AGGREGATE')
            CALL model%Goto(this%parent)
          ELSE
            this%visited = .FALSE.
            msg => NULL()
          END IF
        END IF
      CASE ('AGGREGATE')
        IF (this%process_action(task%process_id) >= AGGREGATE_ .AND. this%visited) THEN
          task%msg => msg
          CALL task%Do_it(this, model%options)
        END IF
        IF (ASSOCIATED(this%next_sibling)) THEN
          msg%action => model%Get_action('INTEGRATE')
          CALL model%Goto(this%next_sibling)
        ELSE IF (this == model%top) THEN
          this%visited = .FALSE.
          msg => NULL()
        ELSE
          msg%action => model%Get_action('AGGREGATE')
          CALL model%Goto(this%parent)
        END IF
      END SELECT
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Error running task')
    END SELECT

    return_ptr => msg

    !SELECT TYPE(model=>hsm)
    !CLASS IS (t_jsb_model)
    !  IF (debug_on() .AND. model%options%iblk == 1) THEN
    !    IF (one_of(msg_in%action%name, (/'ENTER','START','EXIT '/)) < 0) THEN
    !      CALL message(TRIM(routine),'Finished')
    !    END IF
    !  END IF
    !END SELECT

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
        CALL message('     LCT', TRIM(this%lcts(ilct)%Get_name())//' ... '//TRIM(this%lcts(ilct)%Get_longname()))
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
