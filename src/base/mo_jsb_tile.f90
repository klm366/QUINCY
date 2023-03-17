!>
!> @brief Contains concrete tile class that contains the surface structure and memory state.
!>
!> @author
!>  Reiner Schnur, MPI-M Hamburg
!>
!> @par Revision History
!> First version by Reiner Schnur (2013-04-04)
!>
!> This had to be separated from mo_jsb_tile_class in order to avoid circular dependencies.
!>
!> @par Copyright
!> This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!> Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!>

MODULE mo_jsb_tile
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_control,         ONLY: debug_on, get_no_of_models

  USE mo_hsm_class,           ONLY: t_State, t_Message
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,   ONLY: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  USE mo_jsb_process_factory, ONLY: max_no_of_processes
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task !, t_jsb_task_msg
  USE mo_jsb_cqt_class,       ONLY: t_jsb_consQuan

#ifdef _OPENACC
  USE openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)  
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile           !< provides the components defined in tile_class abstract object
    TYPE(t_jsb_model),         POINTER :: model => NULL()    !< connects tile to the model object
    LOGICAL :: in_var_groups                                 !< defines whether tile is part of a output group
  CONTAINS
    !> These deferred routines must match the interfaces provided in jsb_tile_class
    PROCEDURE :: Init           => Init_tile
    PROCEDURE :: Init_fractions => Init_fractions_on_tile
    PROCEDURE :: Init_vars      => Init_vars_on_tile
    PROCEDURE :: Print          => Print_tile
    PROCEDURE :: Handler        => Process_task_on_tile
    PROCEDURE :: Is_process_calculated       => Is_process_calculated_on_tile
    PROCEDURE :: Is_process_active           => Is_process_active_on_tile
    PROCEDURE :: Has_process_memory          => Has_process_memory_on_tile
    PROCEDURE :: Count_conserved_quantities  => Count_conserved_quantities_on_tile
    PROCEDURE :: Collect_conserved_variables => Collect_conserved_variables_on_tile
    PROCEDURE :: Count_and_classify_bgc_materials => Count_and_classify_bgc_materials_on_tile
    PROCEDURE :: Collect_bgc_materials            => Collect_bgc_materials_on_tile
  END TYPE t_jsb_tile

  !> Used to create a vector of tiles
  ! (Define a type only containing a pointer to enable referring to objects of t_jsb_tile in vectors)
  TYPE t_jsb_tile_p
    TYPE(t_jsb_tile), POINTER :: p
  END TYPE t_jsb_tile_p

  !> Derived-type constructor used to create the tile instance
  INTERFACE t_jsb_tile
    PROCEDURE Construct_tile
  END INTERFACE t_jsb_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile'

CONTAINS

  ! ======================================================================================================= !
  !> 
  !> Constructs the tile object, as part of the tile tree, and returns the pointer to the tile
  !>
  FUNCTION Construct_tile(name, description, processes, lct, &
    & process_actions, parent, model_id, fract_filename, fract_varname) RESULT(return_ptr)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_util,                ONLY: one_of
    USE mo_jsb_cqt_class,       ONLY: Get_number_of_types
    USE mo_jsb_lct_class,       ONLY: t_jsb_lct, LAND_TYPE, LAKE_TYPE
    USE mo_jsb_process_class,   ONLY: SEB_, A2L_, L2A_, ANTHRO2L_     ! < Rackhun Son >

    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in), OPTIONAL :: processes(:)
    TYPE(t_jsb_lct),   INTENT(in), OPTIONAL :: lct            !< primary lct of this tile
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id       !< Only needed if parent is not present
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_filename !< Name of file with tile fractions
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_varname  !< Name of fraction variable in fract_filename
    CLASS(t_jsb_tile), POINTER              :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),            POINTER :: model
    INTEGER                               :: i, j, iproc
    CLASS(t_State),               POINTER :: parent_tmp       !< Temporary pointer to the parent tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Construct_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting construction of tile '//TRIM(name))

    ! Allocate and initialize base structure (=state) for tile
    ALLOCATE(t_jsb_tile::return_ptr)

    ! Add the parent tile to the actual tile (for all but the box tile)
    ! init the parent with information (properties) from the state machine only
    IF (PRESENT(parent)) THEN
      parent_tmp => parent
      CALL return_ptr%Init_state(name, description, parent_tmp)
    ELSE
      CALL return_ptr%Init_state(name, description)
    END IF

    ! Initialize further components of t_jsb_tile_abstract extended type
    return_ptr%self_tile => return_ptr

    ! allocate parent tile
    ALLOCATE(return_ptr%process_action(max_no_of_processes))
    ! as defined in usecase, check if each process has one action (e.g., on_tile etc.)
    IF (PRESENT(processes) .AND. PRESENT(process_actions)) THEN
      IF (SIZE(processes) /= SIZE(process_actions)) &
        & CALL finish(TRIM(routine), 'Size of processes and process_actions must be equal!')
    END IF
    ! give owner_model_id to child tiles from parent tile, hence, only box tile needs to get a owner_model_id
    IF (PRESENT(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile)
        return_ptr%parent_tile => parent
        return_ptr%owner_model_id = parent%owner_model_id
      END SELECT
    ELSE
      ! box tile already has a owner_model_id
      IF (PRESENT(model_id)) THEN
        return_ptr%owner_model_id = model_id
      ELSE
        CALL finish(TRIM(routine), 'Need model_id')
      END IF
    END IF

    ! get mode from model_id
    model => Get_model(return_ptr%owner_model_id)
    return_ptr%model => model
    ! ..
    IF (ASSOCIATED(return_ptr%parent_tile)) THEN
      IF (.NOT. ASSOCIATED(return_ptr%parent_tile%first_child_tile)) THEN
          return_ptr%parent_tile%first_child_tile => return_ptr
      END IF
    END IF
    IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
      ! to avoid further SELECT TYPE, just define the CLASS/TYPE here
      SELECT TYPE (prev_sibling => return_ptr%prev_sibling)
      CLASS IS (t_jsb_tile)
        prev_sibling%next_sibling_tile => return_ptr
      END SELECT
    END IF
    ! define lct - land cover type
    ! each tile can have multiple lct, teh 1st one (1st element of the below vector) is the important one
    !  other lct are "inherited" from leave tiles
    !  if an lct is not provided for a tile, this tile gets the same lct as the parent
    IF (PRESENT(lct)) THEN
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        IF (ASSOCIATED(return_ptr%parent_tile%parent_tile) .AND. (lct%id == LAND_TYPE .OR. lct%id == LAKE_TYPE)) &
          & CALL finish(TRIM(routine), 'Tiles with lct LAND_TYPE or LAKE_TYPE must be direct children of root tile.')
      END IF
      ! add lct
      CALL return_ptr%Add_lct(lct)
      ! add the lct of this tile to all parent tiles -> needed for aggregation of variables
      CALL return_ptr%Add_lct_to_parents(lct)
    ELSE
      ! If no lct is passed upon tile construction the primary lct of the parent tile is used
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        CALL return_ptr%Add_lct(return_ptr%parent_tile%lcts(1))
      END IF
    END IF
    ! filename of the fract file
    IF (PRESENT(fract_filename)) THEN
      ! use this filename
      return_ptr%fract_filename = TRIM(fract_filename)
    ELSE
      ! use filename from namelist
      return_ptr%fract_filename = TRIM(model%config%fract_filename)
    END IF
    ! what variable name to read the fraction from
    ! 'fract_varname' is define in the usecase
    ! if not defined in the usecase, it is set to "empty" here, but changed later
    IF (PRESENT(fract_varname)) THEN
      return_ptr%fract_varname = TRIM(fract_varname)
    ELSE
      return_ptr%fract_varname = ''
    END IF

    ! - Init the list of conserved quantities (common for all LCC processes working with the tile)
    ALLOCATE(return_ptr%conserved_quantities(Get_number_of_types()))

    ! Loop over all processes and set actions for this tile
    ! (init is going from box to bottom tile)
    ! process action of each process, as defined in the usecase
    DO i=1,max_no_of_processes
      ! default is SKIP, i.e., update_*() routine may not do anything
      return_ptr%process_action(i) = SKIP_
      ! if process and config is associated
      IF (ASSOCIATED(model%processes(i)%p)) THEN
        IF (ASSOCIATED(model%processes(i)%p%config)) THEN
          ! go on only if process is set ACTIVE, the default is FALSE; needs namelist setting to activate
          IF (.NOT. model%processes(i)%p%config%active) CYCLE
        END IF
      END IF
      iproc = -1
      j = 0
      IF (PRESENT(processes)) THEN
        IF (SIZE(processes) > 0) THEN
          ! j is the process that was explicitly defined in the usecase fir this tile
          j = one_of(i, processes) ! i counts through the number of all processes; j is the number of the i-process within
                                   ! the explicitly defined input processes for this tile (e.g. defined in init_usecase)
        END IF
      END IF
      ! the i-process was explicitly defined for this tile
      IF (j > 0) THEN
        iproc = processes(j)
        ! TODO hard coded checks of process dependencies - should not be handled here
        IF (iproc == SEB_ .AND. ALLOCATED(return_ptr%lcts)) THEN
          IF (return_ptr%lcts(1)%id /= LAND_TYPE .AND. return_ptr%lcts(1)%id /= LAKE_TYPE) THEN
            CALL finish(TRIM(routine), 'SEB_ process only allowed on tiles with lct LAND_TYPE or LAKE_TYPE.')
          END IF
        END IF
        ! if there is a process action at this tile for the process
        !  (would probably be not if process action is simply inherited from parent tile)
        IF (PRESENT(process_actions)) THEN
          CALL return_ptr%Set_process_action(iproc, process_actions(j), action_to_parents=AGGREGATE_)
        ELSE          ! use defaults if process_actions are not given
          IF (iproc == A2L_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE IF (iproc == L2A_) THEN
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)
          ELSE IF (iproc == ANTHRO2L_) THEN     ! < Rackhun Son >
            CALL return_ptr%Set_process_action(iproc, ON_TILE_, action_to_parents=AGGREGATE_)

          ! if there is no process action defined, then we assume that the parent tile has defined it as on_leaves
          !   otherwise the process would 
          ! after finish contructing the tiles we know how the tile tree is structured -> ON_LEAFS_ is not needed anymore then
          !   then, after tile tree construction (with init tile), 
          !      the AGGREGATE ON_TILE INHERIT SKIP are the process actions used for the tiles
          !   i.e., the SKIP and ON_LEAFS_ process actions are transformed in one of the above for four process actions
          
          ELSE
            CALL return_ptr%Set_process_action(iproc, ON_LEAFS_, action_to_parents=AGGREGATE_)
          END IF
        END IF
      ! the i-process was NOT explicitly defined for this tile
      ! "ASSOCIATED(return_ptr%parent_tile)" asks whether we are at the box tile or not
      ! then look at the process action of the parent tile, and by default INHERIT 
      !   the tile would use the memory of the parent tile
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN  
        !IF (return_ptr%parent_tile%process_action(i) > AGGREGATE_) THEN
          iproc = i
          SELECT CASE (return_ptr%parent_tile%process_action(iproc))
          CASE (SKIP_, INHERIT_, ON_SUBTREE_, ON_LEAFS_)
            CALL return_ptr%Set_process_action(iproc, return_ptr%parent_tile%process_action(iproc))
            !return_ptr%parent_tile%process_action(iproc) = AGGREGATE_
          CASE (ON_TILE_)
            ! especially for L2A_ we want to do nothing on the child tiles of the box tile
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

  ! ======================================================================================================= !
  !> 
  !> Init tile - in particular initialise memory for processes on the tile
  !>
  ! Note: It is important that Init_tile is only called for the tiles after the hierarchical tile tree has
  !       been completely constructed. Init_tile is called for each tile from jsbach_init (CALL tile%Init,
  !       in mo_jsb_model_init.f90). The corresponding Init_interface is defined in mo_jsb_tile_class.
  SUBROUTINE Init_tile(this, varlist_name, prefix, suffix, grid_id, in_var_groups)

    USE mo_jsb_process_factory, ONLY: Create_process_memory
    USE mo_jsb_process_class,   ONLY: Get_process_name
    USE mo_jsb_memory_class,    ONLY: t_jsb_memory

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this              !< Current tile
    CHARACTER(len=*),  INTENT(in)    :: varlist_name      !< Model shortname
    CHARACTER(len=*),  INTENT(in)    :: prefix            !< Prefix for variable names (generally empty)
    CHARACTER(len=*),  INTENT(in)    :: suffix            !< Suffix for variable names
    INTEGER,           INTENT(in)    :: grid_id           !< Grid id
    LOGICAL,           INTENT(in)    :: in_var_groups     !< Should output be written for this tile?
    ! -------------------------------------------------------------------------------------------------- !

    INTEGER                               :: iproc
    TYPE(t_jsb_model), POINTER            :: model
    CLASS(t_jsb_memory), POINTER          :: tmp_ptr
    CLASS(t_jsb_tile_abstract), POINTER   :: parent_tile => NULL()
    CHARACTER(LEN=:), ALLOCATABLE         :: process_name

    CHARACTER(len=:), ALLOCATABLE :: varlist_name_loc     ! Local varlist name
    CHARACTER(len=3)              :: model_id_suffix      ! Suffix for domain if more than one model
    INTEGER                       :: no_of_children, child_no

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile '//TRIM(this%name))

    WRITE(model_id_suffix, '(A1,I2.2)') 'D', this%owner_model_id

    model => this%model
    this%grid_id = grid_id                ! Defined for convenience

    ! Allocate pointer array to collect lcc structures for lcc procs
    ALLOCATE(this%lcc_processes(max_no_of_processes))

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(this%mem(max_no_of_processes))
    !$acc enter data create(this%mem)

    ! Loop over all processes and initialize memory for this tile
    parent_tile => this%parent_tile       ! Defined for convenience
    DO iproc=1,max_no_of_processes

      this%mem(iproc)%p => NULL()

      IF (this%process_action(iproc) == SKIP_) CYCLE      ! Nothing needs to be done

      ! If current tile is a leaf, we replace ON_LEAF by ON_TILE and set the parent tile to AGGREGATE
      IF (this%process_action(iproc) == ON_LEAFS_) THEN
        IF (ASSOCIATED(parent_tile)) parent_tile%process_action(iproc) = AGGREGATE_
        IF (.NOT. this%Has_children()) this%process_action(iproc) = ON_TILE_
      END IF

      ! Not all processes need memory (compare process_factory):
      ! Only create / link memory for processes that have a memory
      IF (.NOT. model%processes(iproc)%p%has_memory) CYCLE

      process_name = Get_process_name(iproc)

      IF (this%process_action(iproc) == INHERIT_) THEN
        ! We do not need to generate memory, just use the memory of the parent.
        this%mem(iproc)%p => parent_tile%mem(iproc)%p
        ! No output is needed for the tile, as the parent output is identical (???).
        this%mem(iproc)%p%in_var_groups = in_var_groups
        __acc_attach(this%mem(iproc)%p)
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Linked memory for process '//TRIM(process_name)//&
          &                             ' on tile '//TRIM(this%name)//' to parent')
      ELSE
        ! this%param(iproc)%p => t_jsb_param(iproc)
        tmp_ptr => Create_process_memory(iproc) ! Create memory for the specific process
        __acc_attach(tmp_ptr)
        this%mem(iproc)%p => tmp_ptr     !TODO: Still needed? Same as this%mem(iproc)%p%owner_proc_id  (see below)
        __acc_attach(this%mem(iproc)%p)  !      (Needs to be tested on gpu)
        this%mem(iproc)%p%id = iproc

        ! Preparations needed for the Init process memory function called below
        ! Generate variable list name. Each process has its own variable list; used in add_var
        IF (get_no_of_models() > 1) THEN
          varlist_name_loc = TRIM(varlist_name)//'_'//TRIM(process_name)//'_'//model_id_suffix
        ELSE
          varlist_name_loc = TRIM(varlist_name)//'_'//TRIM(process_name)
        END IF
        this%mem(iproc)%p%varlist_name = TRIM(varlist_name_loc)
  
        ! Also the memory of the parent tile is available on the current tile (if it has a parent).
        IF (ASSOCIATED(parent_tile)) this%mem(iproc)%p%parent => parent_tile%mem(iproc)%p

        this%mem(iproc)%p%grid_id = grid_id
        this%mem(iproc)%p%owner_model_id = this%owner_model_id
        this%mem(iproc)%p%owner_proc_id = iproc
        this%mem(iproc)%p%owner_proc_name = process_name
        ALLOCATE(this%mem(iproc)%p%owner_tile_path(SIZE(this%path)))
        this%mem(iproc)%p%owner_tile_path(:) = this%path(:)
        this%mem(iproc)%p%owner_tile_name = TRIM(this%name)
        this%mem(iproc)%p%in_var_groups = in_var_groups    ! Must be before CALL to %Init !
        !$acc update device(this%mem(iproc)%p)
        !$acc enter data create(this%mem(iproc)%p%vars)
        IF (prefix /= '') THEN
          !TODO: check if prefix is still needed ...
          CALL this%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
        ELSE
          CALL this%mem(iproc)%p%Init(TRIM(process_name), TRIM(suffix), this%lcts(:)%id, this%owner_model_id)
        END IF
        !$acc update device(this%mem(iproc)%p%vars)
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name))

        ! Make the current tile's memory available from the parent tile
        IF (ASSOCIATED(parent_tile)) THEN                               ! Tile has a parent (i.e. this is not the box tile)
          no_of_children = parent_tile%Get_no_of_children()             ! Total number of siblings
          IF (.NOT. ASSOCIATED(parent_tile%mem(iproc)%p%children) .AND. no_of_children > 0) THEN
            parent_tile%mem(iproc)%p%no_of_children = no_of_children    !TODO: this line is probably not needed
            ALLOCATE(parent_tile%mem(iproc)%p%children(no_of_children)) ! Allocate memory for all siblings
            !$acc enter data create(parent_tile%mem(iproc)%p%children)
          END IF
          child_no = this%path(this%level)                ! Index of the child (i.e. the current tile; comes from state machine)
          parent_tile%mem(iproc)%p%children(child_no)%p => this%mem(iproc)%p
          !$acc update device(parent_tile%mem(iproc)%p%children(child_no)%p)
          __acc_attach(parent_tile%mem(iproc)%p%children(child_no)%p)
        END IF

      END IF

    END DO
    !$acc update device(this%mem)

  END SUBROUTINE Init_tile

  ! ======================================================================================================= !
  !> 
  !> Add and initialise fraction variables for this tile - either from a file or according to specification 
  !>
  SUBROUTINE Init_fractions_on_tile(this, varlist_name, prefix, suffix, l_fixed_fractions, l_rel_fractions)

    USE mo_jsb_tile_class,      ONLY: t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_varlist,         ONLY: jsb_add_var, VARNAME_LEN
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_INSTANT, TSTEP_CONSTANT
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this               !< Current tile
    CHARACTER(len=*),  INTENT(in)    :: varlist_name       !< Model shortname
    CHARACTER(len=*),  INTENT(in)    :: prefix             !< Prefix for variable names (generally empty)
    CHARACTER(len=*),  INTENT(in)    :: suffix             !< Suffix for variable names
    LOGICAL,           INTENT(in)    :: l_fixed_fractions
      !< If tile fractions can change (active lcc process). I.a. determines if fracts go to restart file
    LOGICAL,           INTENT(in)    :: l_rel_fractions
      !< Tile fractions in input file are usually relative with respect to parent (historical reasons)
      !< while the fractions in the code are absolute (with respect to the box tile)
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator_weighted_by_fract
    CLASS(t_jsb_aggregator),                  POINTER :: aggregator
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER                    :: table
    TYPE(t_input_file)         :: input_file
    REAL(wp),          POINTER :: return_pointer(:,:), tmp_fract(:,:)
    CHARACTER(len=50)          :: name_loc
    CHARACTER(len=:), ALLOCATABLE :: varlist_name_loc
    CHARACTER(len=3)           :: model_id_suffix                    ! Suffix for domain if more than one model

    INTEGER :: isteptype, i, ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fractions_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing fractions on tile '//TRIM(this%name))

    hgrid   => Get_grid(this%grid_id)
    surface => Get_vgrid('surface')

    table = tables(1)

    ! -------------------------------------------------------------------------------------------------- !
    ! Add fraction(s) of this tile to varlist

    ! isteptype us used in the add_var below 
    ! - if constant cdi knows that its not per timestep / i.e. not instantaneous
    IF (l_fixed_fractions) THEN
      isteptype = TSTEP_CONSTANT
    ELSE
      isteptype = TSTEP_INSTANT
    END IF

    WRITE(model_id_suffix, '(A1,I2.2)') 'D', this%owner_model_id
    IF (get_no_of_models() > 1) THEN
      varlist_name_loc = TRIM(varlist_name)//'_'//model_id_suffix
    ELSE
      varlist_name_loc = TRIM(varlist_name)
    END IF

    !$acc enter data copyin(this)

    ! Note: uses a much simpler add_var routine than process variables do
    ! Further, the fractions are direct components of the tile class and do not 'live' on a process memory
    ! Therefore, as opposed to the process variables a temporal pointer (tmp_fract) is used here,
    ! where variables on the process memory already use the corresponding pointer on the memory
    ! In addition, there is also no process prefix on this variable in the output
    CALL jsb_add_var(varlist_name_loc, 'fract', tmp_fract,                                                    &
      & hgrid, surface,                                                                                       &
      & t_cf('tile_fraction', '',                                                                             &
      &      'fraction of grid box covered by tile '//TRIM(this%name)),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
      & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                     &
      & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])

    ! Now we need to set the tile component fract pointer to the just allocated memory
    CALL this%Set_fract_ptr(tmp_fract)
    NULLIFY(tmp_fract)

    ! And reuse tmp_fract
    CALL jsb_add_var(varlist_name_loc, 'fract_old', tmp_fract,                                                 &
      & hgrid, surface,                                                                                        &
      & t_cf('tile_fraction_old', '',                                                                          &
      &      'old fraction of grid box covered by tile '//TRIM(this%name)),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                     &
      & prefix, suffix, loutput=.NOT. l_fixed_fractions, lrestart=.NOT. l_fixed_fractions,                     &
      & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])
    CALL this%Set_fract_old_ptr(tmp_fract)
    NULLIFY(tmp_fract)


    ! -------------------------------------------------------------------------------------------------- !
    ! Initialize fraction(s) of this tile from input file
  
    ! **** here we continue next time ****
    ALLOCATE(this%l_fract_children_changed(hgrid%nproma, hgrid%nblks))
    this%l_fract_children_changed(:,:) = .TRUE.

    ALLOCATE(tmp_fract(hgrid%nproma, hgrid%nblks))
    IF (TRIM(this%fract_varname) == 'equal') THEN
      IF (ASSOCIATED(this%parent_tile)) THEN
        CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        tmp_fract(:,:) = tmp_fract(:,:) / REAL(this%parent_tile%Get_no_of_children(), wp)
      ELSE
        tmp_fract(:,:) = 1._wp
      END IF
      CALL this%Set_fraction(fract=tmp_fract(:,:))
    ELSEIF (this%fract_varname == 'allFirst') THEN
        IF(ASSOCIATED(this%parent)) THEN
          CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
          IF (.NOT. this%Is_first_child()) THEN
            tmp_fract(:,:) = 0._wp
          END IF
        ELSE
          tmp_fract(:,:) = 1._wp
        END IF
        CALL this%Set_fraction(fract=tmp_fract(:,:))
    ELSE
      input_file = jsb_netcdf_open_input(TRIM(this%fract_filename), this%grid_id)

      IF (TRIM(this%fract_varname) == '') THEN
        name_loc = 'fract'
        IF (prefix /= '') name_loc = TRIM(prefix)//'_'//TRIM(name_loc)
        IF (suffix /= '') name_loc = TRIM(name_loc)//'_'//TRIM(suffix)
      ELSE
        name_loc = TRIM(this%fract_varname)
      END IF
      IF (input_file%Has_var(TRIM(name_loc))) THEN
        return_pointer => input_file%Read_2d(variable_name=TRIM(name_loc))

        !
        ! TODO: For security set all values of the box fractions to 1 that are close to 1, and all box fractions to 0
        ! that are close to 0 (similar to mo_echam_phy_init) in order to account for precision errors
        ! IF (TRIM(this%name) == 'box') THEN
        !   return_pointer(:,:) = MERGE(1._wp, return_pointer(:,:), return_pointer(:,:) > 1._wp - 10._wp*EPSILON(1._wp))
        !   return_pointer(:,:) = MERGE(0._wp, return_pointer(:,:), return_pointer(:,:) <         10._wp*EPSILON(1._wp))
        ! END IF
        !
        ! TODO: For security, as a workaround, set lake fractions in coastal points that are 
        ! partially land and partially ocean to zero, and corresponding land fraction to one.
        ! IF (TRIM(this%name) == 'lake') THEN
        !   CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        !   return_pointer(:,:) = MERGE(return_pointer(:,:), 0._wp, tmp_fract(:,:) >= 1._wp)
        ! END IF
        ! IF (TRIM(this%name) == 'land') THEN
        !   CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        !   return_pointer(:,:) = MERGE(return_pointer(:,:), 1._wp, tmp_fract(:,:) >= 1._wp)
        ! END IF
      ELSE
        CALL finish(TRIM(routine), 'Fraction '//TRIM(name_loc)//' for tile '//TRIM(this%name)// &
          &                        ' not found in '//TRIM(this%fract_filename))
      END IF
      IF (ASSOCIATED(this%parent_tile) .AND. l_rel_fractions) THEN
        CALL this%parent_tile%Get_fraction(fract=tmp_fract(:,:))
        CALL this%Set_fraction(fract=MAX(0._wp, return_pointer(:,:)*tmp_fract(:,:)))
      ELSE
        CALL this%Set_fraction(fract=return_pointer(:,:))
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
    !$acc enter data create(this%lcts)
    DO ilct=1,SIZE(this%lcts)

      !! !$acc enter data create(this%lcts(ilct)%fract)
      CALL jsb_add_var(varlist_name_loc, 'fract_'//TRIM(this%lcts(ilct)%name), this%lcts(ilct)%fract,             &
        & hgrid, surface,                                                                                         &
        & t_cf('lct_fraction', '',                                                                                &
        &      'fraction of '//TRIM(this%name)//' covered by LCT '//TRIM(this%lcts(ilct)%name)),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix, lrestart=.NOT. l_fixed_fractions,                                                       &
        & initval_r=0.0_wp, isteptype=isteptype, groups=[character(len=VARNAME_LEN) :: 'jsb_tile_fractions'])
      __acc_attach(this%lcts(ilct)%fract)

    END DO

    DEALLOCATE(tmp_fract)

  END SUBROUTINE Init_fractions_on_tile

  ! ======================================================================================================= !
  !> 
  !> Init variabls on nodes (not-leave tiles); called during initialization from model_init
  !>
  SUBROUTINE Init_vars_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER  :: mem
    CHARACTER(len=VARNAME_LEN)    :: name
    INTEGER                       :: iproc, i, no_of_children, child_no

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_vars_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile structure on vars '//TRIM(this%name))

    no_of_children = this%Get_no_of_children()
    IF (no_of_children == 0) CALL finish(routine, 'This should not happen')

    ! Loop over all processes and initialize memory for this tile
    DO iproc=1,max_no_of_processes

      ! if the current process is not on leaves of current tile exit DO loop
      SELECT CASE (this%process_action(iproc))
      CASE (SKIP_, ON_TILE_, INHERIT_)
        CYCLE
      END SELECT

      IF (.NOT. ASSOCIATED(this%mem(iproc)%p)) CYCLE

      mem => this%mem(iproc)%p
      DO i=1,mem%no_of_vars
        name = mem%vars(i)%p%name
        IF (.NOT. ALLOCATED(mem%vars(i)%p%child_idx)) &
          & ALLOCATE(mem%vars(i)%p%child_idx(no_of_children))
        DO child_no=1,no_of_children
          IF (ASSOCIATED(mem%children(child_no)%p)) THEN
            mem%vars(i)%p%child_idx(child_no) = mem%children(child_no)%p%Get_var_position(TRIM(name))
          ELSE
            mem%vars(i)%p%child_idx(child_no) = 0
          END IF
        END DO
        !$acc enter data copyin(mem%vars(i)%p%child_idx)
        ! __acc_attach(mem%vars(i)%p%child_idx)
      END DO

    END DO

  END SUBROUTINE Init_vars_on_tile

  ! ====================================================================================================== !
  !> 
  !> Check if a given process is to be integrated on the tile: 
  !>  i.e. if the process action on this tile is one of ON_TILE_ or ON_SUBTREE_
  !>
  LOGICAL FUNCTION Is_process_calculated_on_tile(this, iproc)

    USE mo_util,                ONLY: one_of
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile on which the process activity is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which the activity is to be checked
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Is_process_calculated_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    ! Note: It's not necessary to test for ON_LEAFS_ because, after the tile tree has been completely
    !       set up, ON_LEAFS_ on a tile has been replaced by AGGREGATE_ and ON_TILE_ on the sub-tree.
    Is_process_calculated_on_tile = (one_of(this%process_action(iproc), (/ON_TILE_, ON_SUBTREE_/)) > 0)
    
  END FUNCTION Is_process_calculated_on_tile

  ! ====================================================================================================== !
  !>
  !> Check if a given process is to be active (integrated or aggregated) on the tile: 
  !>  i.e. if the process action on this tile is not SKIP_ or INHERIT_
  !>
  LOGICAL FUNCTION Is_process_active_on_tile(this, iproc)

    USE mo_util,                ONLY: one_of
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile on which the process activity is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which the activity is to be checked
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Is_process_active_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    Is_process_active_on_tile = (one_of(this%process_action(iproc), (/SKIP_, INHERIT_/)) < 0)
    
  END FUNCTION Is_process_active_on_tile

  ! ====================================================================================================== !
  !
  !> Check if the tile has memory for the process: i.e. if the process action on this tile is not SKIP_
  !
  LOGICAL FUNCTION Has_process_memory_on_tile(this, iproc)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this  !< tile for which the process memory is to be checked
    INTEGER,           INTENT(in) :: iproc !< process for which it is checked if memory is available
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Has_process_memory_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    Has_process_memory_on_tile = (.NOT. this%process_action(iproc) == SKIP_) 
  END FUNCTION Has_process_memory_on_tile

  ! ======================================================================================================= !
  !>
  !> Identifies the number and types of bgc materials on this tile
  !> NOTE: Currently we are not really generic here but only expect 1l or 2l bgcms! 
  !>       If something more generic is required in the future it'll need to be implemented then
  !>       or alternatively if the bgcm stores as required for quincy are not required but "deeper" bgcm
  !>       structures are, the collection of the pointers in the bgcm store might be surpressed if using
  !>       bgcms but not Quincy...
  !>
  SUBROUTINE Count_and_classify_bgc_materials_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_pool_class,      ONLY: t_jsb_pool
    USE mo_lnd_bgcm_store,      ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_class,      ONLY: assert_no_compartment_bgc_materials, get_dimension_of_element_variables
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< tile for which the bgc materials are checked
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER     :: mem
    CLASS(t_jsb_pool), POINTER       :: bgcm
    CLASS(t_lnd_bgcm_store), POINTER :: bgcm_store
    INTEGER                          :: i_proc, i_bgcm, i_compartment
    INTEGER                          :: nr_of_bgcms_this_proc, nr_of_compartments
    INTEGER                          :: this_dim, common_dim

    CHARACTER(len=*), PARAMETER :: routine = modname//':Count_and_classify_bgc_materials_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on())     CALL message(TRIM(routine), TRIM(this%name))
    
    ! Loop over all processes
    DO i_proc = 1, max_no_of_processes

      ! if the current process is not on the current tile exit DO loop
      SELECT CASE (this%process_action(i_proc))
      CASE (SKIP_, ON_LEAFS_, INHERIT_, AGGREGATE_, ON_SUBTREE_)
        CYCLE
      END SELECT

      ! also exit the Do loop if this process has no memory
      IF (.NOT. ASSOCIATED(this%mem(i_proc)%p)) CYCLE

      ! Only processes are relevant here that have any bgc material in their process memory
      mem => this%mem(i_proc)%p
      IF(mem%has_bgc_materials) THEN

        ! If this is the first memory in which we find a bgc material the store needs to be initialised
        IF (.NOT. ASSOCIATED(this%bgcm_store)) THEN
          ALLOCATE(bgcm_store)
          this%bgcm_store => bgcm_store
        ENDIF

        nr_of_bgcms_this_proc = SIZE(mem%bgc_material%pool_list)
        this%bgcm_store%nr_of_bgc_materials = this%bgcm_store%nr_of_bgc_materials + nr_of_bgcms_this_proc
        DO i_bgcm = 1, nr_of_bgcms_this_proc
          bgcm => mem%bgc_material%pool_list(i_bgcm)%p

          IF(bgcm%contains_elements) THEN
            ! i.e. directly elements and no compartments
            CALL assert_no_compartment_bgc_materials(bgcm, this%name, routine)
            this_dim = get_dimension_of_element_variables(bgcm, this%name, routine)
            IF(this_dim == 2) THEN
              this%bgcm_store%nr_of_1l_2d_bgc_materials = this%bgcm_store%nr_of_1l_2d_bgc_materials + 1
            ELSEIF(this_dim == 3) THEN
              this%bgcm_store%nr_of_1l_3d_bgc_materials = this%bgcm_store%nr_of_1l_3d_bgc_materials + 1
            ENDIF
          ELSE
            ! i.e. no elements but one layer of compartments which are then expected to have elements
            ! but no further compartments! 
            nr_of_compartments = SIZE(bgcm%pool_list)
            common_dim = -1
            DO i_compartment = 1, nr_of_compartments
              IF(.NOT. bgcm%pool_list(i_compartment)%p%contains_elements) THEN
                CALL finish(TRIM(routine), 'Found bgcm '//TRIM(bgcm%name)//' in tile '//TRIM(this%name)// &
                  & ' that contains a compartment without elements '//                                    &
                  &  TRIM(bgcm%pool_list(i_compartment)%p%name)//                                         &
                  & ' but no behaviour is implemented for this so far, please check!')
              ENDIF
              CALL assert_no_compartment_bgc_materials(bgcm%pool_list(i_compartment)%p, this%name, routine)
              this_dim = get_dimension_of_element_variables(bgcm%pool_list(i_compartment)%p, this%name, routine)
              IF(i_compartment == 1) THEN
                common_dim = this_dim
              ELSEIF(common_dim /= this_dim) THEN
                CALL finish(TRIM(routine), &
                  & 'Found compartments with elements with different dimensions for the bgcm '// &
                  & TRIM(bgcm%name)//' in tile '//TRIM(this%name)//                              &
                  & ', but no behaviour is implemented for this so far, please check!')
              ENDIF
            ENDDO ! For each compartment
            IF(common_dim == 2) THEN
              this%bgcm_store%nr_of_2l_2d_bgc_materials = this%bgcm_store%nr_of_2l_2d_bgc_materials + 1
            ELSEIF(common_dim == 3) THEN
              this%bgcm_store%nr_of_2l_3d_bgc_materials = this%bgcm_store%nr_of_2l_3d_bgc_materials + 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE Count_and_classify_bgc_materials_on_tile

  ! ======================================================================================================= !
  !>
  !> Collects the pointers to the variables that contain bgc elements
  !>
  SUBROUTINE Collect_bgc_materials_on_tile(this)
    USE mo_jsb_memory_class,     ONLY: t_jsb_memory
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< tile for which the conserved quantities are collected
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER  :: mem
    INTEGER                       :: i_proc
    INTEGER                       :: ind_1l_2d_bgcm, ind_2l_2d_bgcm, ind_1l_3d_bgcm, ind_2l_3d_bgcm
    INTEGER                       :: ind_bgcm

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_bgc_materials_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on())     CALL message(TRIM(routine), TRIM(this%name))

    ! If there is any process on this tile that has bgc materials, an instance of the bgcm store
    ! would have been created on this tile in the collection routine
    IF(ASSOCIATED(this%bgcm_store)) THEN
      CALL this%bgcm_store%init()

      ! and init indices for bookkeeping
      ind_bgcm = 0
      ind_1l_2d_bgcm = 0
      ind_2l_2d_bgcm = 0
      ind_1l_3d_bgcm = 0
      ind_2l_3d_bgcm = 0
    ELSE
      ! If there is no bgcm store, then there are no active processes with bgcm material on this tile
      RETURN
    ENDIF

    ! Loop over all processes
    DO i_proc=1,max_no_of_processes
      ! if the current process is not on the current tile exit DO loop
      SELECT CASE (this%process_action(i_proc))
      CASE (SKIP_, ON_LEAFS_, INHERIT_, AGGREGATE_, ON_SUBTREE_)
        CYCLE
      END SELECT

      ! also exit the DO loop if this process has no memory
      IF (.NOT. ASSOCIATED(this%mem(i_proc)%p)) CYCLE

      ! Only processes are relevant here that have any bgc material in their process memory
      mem => this%mem(i_proc)%p
      IF(mem%has_bgc_materials) THEN
        CALL this%bgcm_store%append_bgcms(mem%bgc_material, this%name, &
          & ind_bgcm, ind_1l_2d_bgcm, ind_2l_2d_bgcm, ind_1l_3d_bgcm, ind_2l_3d_bgcm)
      ENDIF
    ENDDO

  END SUBROUTINE Collect_bgc_materials_on_tile

  ! ======================================================================================================= !
  !>
  !> Identifies the conserved quantity types (CQTs) in use on this tile and counts conserved quantities (CQs) of each CQT
  !>
  SUBROUTINE Count_conserved_quantities_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_cqt_class,       ONLY: Get_cqt_name, max_cqt_name_length
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< tile for which the conserved quantities are to be counted
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER   :: mem
    INTEGER                        :: iproc, i, id, ind
    LOGICAL                        :: consQuanTypeInUse
    TYPE(t_jsb_consQuan), POINTER  :: consQuanType
    CHARACTER(len=max_cqt_name_length) :: thisCQTName

    CHARACTER(len=*), PARAMETER :: routine = modname//':Count_conserved_quantities_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on())     CALL message(TRIM(routine), TRIM(this%name))
    
    ! Loop over all processes
    DO iproc=1,max_no_of_processes

      ! if the current process is not on the current tile exit DO loop
      SELECT CASE (this%process_action(iproc))
      CASE (SKIP_, ON_LEAFS_, INHERIT_, AGGREGATE_, ON_SUBTREE_)
        CYCLE
      END SELECT

      IF (.NOT. ASSOCIATED(this%mem(iproc)%p)) CYCLE

      ! loop over all variables in the memory of the process and count number of vars for each used type
      mem => this%mem(iproc)%p
      DO i=1,mem%no_of_vars

        ! check those variables that have a conservative quantity type
        IF ( mem%vars(i)%p%is_conserved_quan ) THEN
          id = mem%vars(i)%p%cons_quan_type_id

          ! check if this type is already accounted for
          consQuanTypeInUse = .FALSE.
          DO ind=1,this%nr_of_cqts
            IF (this%conserved_quantities(ind)%p%type_id .EQ. id) THEN
              consQuanTypeInUse = .TRUE.
              this%conserved_quantities(ind)%p%no_of_vars = this%conserved_quantities(ind)%p%no_of_vars + 1             
              EXIT
            END IF
          END DO

          IF (.NOT. consQuanTypeInUse) THEN
            ! needs to be one of those in mo_jsb_cqt_class
            thisCQTName = Get_cqt_name(id)
            IF (debug_on())  CALL message(TRIM(routine), '... vars of type '//TRIM(thisCQTName))

            this%nr_of_cqts = this%nr_of_cqts + 1

            ALLOCATE(consQuanType)
            this%conserved_quantities(this%nr_of_cqts)%p => consQuanType
            this%conserved_quantities(this%nr_of_cqts)%p%type_id = id
            this%conserved_quantities(ind)%p%no_of_vars = this%conserved_quantities(ind)%p%no_of_vars + 1
          ENDIF

        END IF
      END DO
    END DO

  END SUBROUTINE Count_conserved_quantities_on_tile




  ! ======================================================================================================= !
  !>
  !> Collects the pointers to the variables that are conserved quantities on the given tile
  !>
  SUBROUTINE Collect_conserved_variables_on_tile(this)

    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< tile for which the conserved quantities are collected
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER  :: mem
    INTEGER                       :: iproc, i, id, ind, last_index, nrOfConsQuans
    LOGICAL                       :: consQuanTypeInUse

    CHARACTER(len=*), PARAMETER :: routine = modname//':Collect_conserved_variables_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    IF (debug_on())     CALL message(TRIM(routine), TRIM(this%name))
    
    ! Loop over all processes
    DO iproc=1,max_no_of_processes

      ! if the current process is not on the current tile exit DO loop
      SELECT CASE (this%process_action(iproc))
      CASE (SKIP_, ON_LEAFS_, INHERIT_, AGGREGATE_, ON_SUBTREE_)
        CYCLE
      END SELECT

      IF (.NOT. ASSOCIATED(this%mem(iproc)%p)) CYCLE

      ! loop over all variables in the memory of the process
      mem => this%mem(iproc)%p
      DO i=1,mem%no_of_vars

        ! only address those variables that have a conservative quantity type
        IF ( mem%vars(i)%p%is_conserved_quan ) THEN
          id = mem%vars(i)%p%cons_quan_type_id
          consQuanTypeInUse = .FALSE.

          ! find container for this type 
          DO ind=1,this%nr_of_cqts
            IF (this%conserved_quantities(ind)%p%type_id .EQ. id) THEN
              consQuanTypeInUse = .TRUE.
              nrOfConsQuans = this%conserved_quantities(ind)%p%no_of_vars

              ! If this is the first collected variable: allocate collector for var pointers 
              IF (.NOT. (ALLOCATED(this%conserved_quantities(ind)%p%cq_vars_2D))) THEN
                ALLOCATE(this%conserved_quantities(ind)%p%associated_process(nrOfConsQuans))
                ALLOCATE(this%conserved_quantities(ind)%p%cq_vars_2D(nrOfConsQuans))
              ENDIF

              this%conserved_quantities(ind)%p%last_index_used = this%conserved_quantities(ind)%p%last_index_used + 1
              last_index = this%conserved_quantities(ind)%p%last_index_used
              ! Assert: last_index > no_of_vars should not be possible
              IF ( nrOfConsQuans < last_index) &
                & CALL finish(TRIM(routine), 'Assigned conserved quantities exceed expected number.')

              ! Add current variable
              this%conserved_quantities(ind)%p%cq_vars_2D(last_index)%p => mem%vars(i)%p
              this%conserved_quantities(ind)%p%associated_process(last_index) = mem%owner_proc_id

              ! Since the type has been found: exit do loop now
              EXIT
            END IF
          END DO

          ! Assert: type should have been found
          IF (.NOT. consQuanTypeInUse) &
            & CALL finish(TRIM(routine), 'Conserved quantity type not found.')

        END IF
      END DO
    END DO

  END SUBROUTINE Collect_conserved_variables_on_tile

  ! ======================================================================================================= !
  !>
  !> Task handler function for tile (HSM message handler)
  !>
  FUNCTION Process_task_on_tile(this, msg_in) RESULT(return_ptr)

    USE mo_util, ONLY: one_of, logical2string, int2string
    USE mo_jsb_parallel, ONLY: Get_omp_thread, Get_omp_no_of_threads
    USE mo_jsb_control,  ONLY: timer_on, timer_process_task
    USE mo_timer,        ONLY: timer_start, timer_stop
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile),  INTENT(inout) :: this
    CLASS(t_Message),   INTENT(in)    :: msg_in
    CLASS(t_Message),   POINTER       :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    !CLASS(t_jsb_task_msg), POINTER :: msg
    TYPE(t_Message), POINTER :: msg
    LOGICAL :: l_debug

    CLASS(t_jsb_process_task), POINTER :: task

    INTEGER :: iblk, no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = modname//':Process_task_on_tile'
    ! -------------------------------------------------------------------------------------------------- !
    no_omp_thread = Get_omp_thread()

    IF (timer_on('detail')) CALL timer_start(timer_process_task(this%owner_model_id))

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
      return_ptr => msg
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

    IF (timer_on('detail')) CALL timer_stop(timer_process_task(this%owner_model_id))

  END FUNCTION Process_task_on_tile

  ! ======================================================================================================= !
  !>
  !> Print routine for the tile
  !>
  SUBROUTINE Print_tile(this)

    USE mo_jsb_process_class,  ONLY: Get_process_name, Get_action_name
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(in) :: this !< tile to be printed
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: ilct, iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_tile'
    ! -------------------------------------------------------------------------------------------------- !
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
