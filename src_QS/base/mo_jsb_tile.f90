!> Contains concrete tile class that contains the surface structure and memory state.
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
!>### Definition of concrete class and methods for tiles
!> This module implements 't_jsb_tile' as an extension of [[mo_jsb_tile_class:t_jsb_tile_abstract]] and
!>   implements, among others, the methods
!>
!> - Init => [[Init_tile]]
!> - Init_fractions => [[Init_fractions_on_tile]]
!> - Init_vars => [[Init_vars_on_tile]]
!> - Handler => [[Process_task_on_tile]]
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_tile
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator_p, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process_p, &
                                    SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  USE mo_jsb_process_factory, ONLY: max_no_of_processes

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile, t_jsb_tile_p

  TYPE, EXTENDS(t_jsb_tile_abstract) :: t_jsb_tile
    LOGICAL :: must_create_memory_parents(max_no_of_processes)    !< helper variable
    TYPE(t_jsb_model),         POINTER :: model => NULL()
    LOGICAL                            :: in_var_groups           !< defines whether tile is part of a output group
  CONTAINS
    PROCEDURE :: Init                             => Init_tile
    PROCEDURE :: Is_process_calculated            => Is_process_calculated_on_tile
    PROCEDURE :: Is_process_active                => Is_process_active_on_tile
    PROCEDURE :: Has_process_memory               => Has_process_memory_on_tile
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

  FUNCTION Construct_tile(name, description, processes, lct, &
    & process_actions, parent, model_id, grid_id, fract_filename, fract_varname) RESULT(return_ptr)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_util,                ONLY: one_of
    USE mo_jsb_lct_class,       ONLY: t_jsb_lct, LAND_TYPE, LAKE_TYPE
    USE mo_jsb_process_class,   ONLY: SEB_, A2L_, L2A_, Get_process_name

    CHARACTER(len=*),  INTENT(in)           :: name
    CHARACTER(len=*),  INTENT(in)           :: description
    INTEGER,           INTENT(in), OPTIONAL :: processes(:)
    TYPE(t_jsb_lct),   INTENT(in), OPTIONAL :: lct            !< primary lct of this tile
    INTEGER,           INTENT(in), OPTIONAL :: process_actions(:)
    CLASS(t_jsb_tile), POINTER,    OPTIONAL :: parent
    INTEGER,           INTENT(in), OPTIONAL :: model_id       !< Only needed if parent is not present
    INTEGER,           INTENT(in), OPTIONAL :: grid_id        !< Only needed if parent is not present
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_filename !< Name of file with tile fractions
    CHARACTER(len=*),  INTENT(in), OPTIONAL :: fract_varname  !< Name of fraction variable in fract_filename
    CLASS(t_jsb_tile), POINTER              :: return_ptr

    TYPE(t_jsb_model),            POINTER :: model
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
      ! ! add the lct of this tile to all parent tiles -> needed for aggregation of variables
      ! CALL return_ptr%Add_lct_to_parents(lct)
    ELSE
      ! If no lct is passed upon tile construction the primary lct of the parent tile is used
      IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        CALL return_ptr%Add_lct(return_ptr%parent_tile%lcts(1))
      END IF
    END IF

      return_ptr%fract_filename = TRIM(fract_filename)

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
          j = one_of(i, processes) ! i counts through the number of all processes, j is the number of the i-process
                                   ! within the processes vector (input)
        END IF
      END IF
      IF (j > 0) THEN
        iproc = processes(j)

        IF (PRESENT(process_actions)) THEN
          return_ptr%process_action(iproc) = process_actions(j)
        ELSE
          IF (iproc == A2L_) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE IF (iproc == L2A_) THEN
            return_ptr%process_action(iproc) = ON_TILE_
          ELSE
            return_ptr%process_action(iproc) = ON_LEAFS_
          END IF
        END IF
        ! All parents of this tile aggregate memory for this new process
        return_ptr%must_create_memory_parents(i) = .TRUE.
!!$      ELSE IF (ASSOCIATED(return_ptr%prev_sibling)) THEN
!!$        iproc = i
!!$        SELECT TYPE (sibling => return_ptr%prev_sibling)
!!$        CLASS IS (t_jsb_tile_abstract)
!!$          return_ptr%process_action(iproc) = sibling%process_action(iproc)
!!$        END SELECT
      ELSE IF (ASSOCIATED(return_ptr%parent_tile)) THEN
        !IF (return_ptr%parent_tile%process_action(i) > AGGREGATE_) THEN
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
            IF (iproc == L2A_) THEN
              return_ptr%process_action(iproc) = SKIP_
            ELSE
              return_ptr%process_action(iproc) = INHERIT_
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
  SUBROUTINE Init_tile(this, varlist_name, prefix, suffix, grid_id, in_var_groups)

    USE mo_jsb_process_factory, ONLY: Create_process_memory
    USE mo_jsb_process_class,   ONLY: Get_process_name
    USE mo_jsb_memory_class,    ONLY: t_jsb_memory
    USE mo_jsb_tile_class,      ONLY: t_jsb_aggregator, t_jsb_aggregator_weighted_by_fract
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables

    CLASS(t_jsb_tile), INTENT(inout) :: this
    CHARACTER(len=*),  INTENT(in)    :: varlist_name
    CHARACTER(len=*),  INTENT(in)    :: prefix         ! Prefix for variable names
    CHARACTER(len=*),  INTENT(in)    :: suffix
    INTEGER,           INTENT(in)    :: grid_id           !< Grid id
    LOGICAL,           INTENT(in)    :: in_var_groups     !< Should output be written for this tile?

    INTEGER                               :: iproc
    CLASS(t_jsb_memory), POINTER          :: tmp_ptr
    CLASS(t_jsb_tile_abstract),   POINTER :: tmp_tile
    CLASS(t_State), POINTER               :: tmp_child
    CLASS(t_jsb_tile_abstract), POINTER   :: parent_tile => NULL()
    CHARACTER(LEN=15)                     :: process_name

    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator_weighted_by_fract
    CLASS(t_jsb_aggregator),                  POINTER :: aggregator
    REAL(wp),          POINTER :: return_pointer(:,:)
    CHARACTER(len=50)          :: name_loc
    CHARACTER(len=VARNAME_LEN) :: name
    INTEGER                    :: ilct, i, no_of_children, child_no

    INTEGER :: isteptype

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_tile'

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing tile '//TRIM(this%name))

    ! Allocate all memory components for new tile and initialize their data part with NULL()
    ALLOCATE(this%mem(max_no_of_processes))

    ! Loop over all processes and initialize memory for this tile
    parent_tile => this%parent_tile       ! Defined for convenience
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
        ! We do not need to generate memory, just use the memory of the parent.
        this%mem(iproc)%p => parent_tile%mem(iproc)%p
        ! No output is needed for the tile, as the parent output is identical (???).
        this%mem(iproc)%p%in_var_groups = in_var_groups
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Linked memory for process '//TRIM(process_name)//&
          &                             ' on tile '//TRIM(this%name)//' to parent')
      ELSE
        no_of_children = this%Get_no_of_children()
        child_no = this%path(this%level)
        ! this%param(iproc)%p => t_jsb_param(iproc)
        tmp_ptr => Create_process_memory(iproc)
        IF (ASSOCIATED(parent_tile)) THEN
          IF (.NOT. ASSOCIATED(parent_tile%mem(iproc)%p%children) .AND. parent_tile%Get_no_of_children() > 0) &
            & ALLOCATE(parent_tile%mem(iproc)%p%children(parent_tile%Get_no_of_children()))
          parent_tile%mem(iproc)%p%children(child_no)%p => tmp_ptr
        END IF
        this%mem(iproc)%p => tmp_ptr
        this%mem(iproc)%p%id = iproc
        this%mem(iproc)%p%grid_id = this%grid_id
        this%mem(iproc)%p%owner_model_id = this%owner_model_id
        this%mem(iproc)%p%owner_proc_id = iproc
        this%mem(iproc)%p%no_of_children = no_of_children
        this%mem(iproc)%p%in_var_groups = in_var_groups    ! Must be before CALL to %Init !
        ! IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
        !   IF (this%Has_process_memory(iproc)) THEN
        !     this%mem(iproc)%p%lrestart_cont = model%processes(iproc)%p%config%lrestart_cont
        !   END IF
        ! END IF
        IF (prefix /= '') THEN
          CALL this%mem(iproc)%p%Init(TRIM(process_name)//'_'//TRIM(prefix), &
            &     TRIM(suffix), this%lcts(:)%id, this%lcts(1)%lib_id, this%owner_model_id)
        ELSE
          CALL this%mem(iproc)%p%Init(TRIM(process_name), &
            &     TRIM(suffix), this%lcts(:)%id, this%lcts(1)%lib_id, this%owner_model_id)
        END IF
        IF (debug_on()) &
          & CALL message(TRIM(routine), 'Created memory for process '//TRIM(process_name)//' on tile '//TRIM(this%name))
      END IF

    END DO

    ! If tile has children, create and register aggregator
!    IF (this%Has_children()) THEN
      aggregator_weighted_by_fract => t_jsb_aggregator_weighted_by_fract(this)
      CALL this%Register_aggregator(aggregator_weighted_by_fract)
!    END IF

  END SUBROUTINE Init_tile

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
    USE mo_jsb_grid_class,       ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,             ONLY: Get_grid, Get_vgrid
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile), INTENT(inout) :: this !< tile for which the conserved quantities are collected

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_memory), POINTER :: mem
    TYPE(t_jsb_grid),    POINTER :: hgrid
    TYPE(t_jsb_vgrid),   POINTER :: vgrid
    INTEGER                      :: i_proc
    INTEGER                      :: ind_1l_2d_bgcm, ind_2l_2d_bgcm, ind_1l_3d_bgcm, ind_2l_3d_bgcm
    INTEGER                      :: ind_bgcm
    INTEGER                      :: nblks, nproma, nsoil

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

    hgrid  => Get_grid(this%grid_id)
    nproma = hgrid%nproma
    nblks  = hgrid%nblks
    vgrid  => Get_vgrid('soil_depth_water')
    nsoil  =  vgrid%n_levels

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
        CALL this%bgcm_store%append_bgcms(mem%bgc_material, nblks, nproma, nsoil, this%name, &
          & ind_bgcm, ind_1l_2d_bgcm, ind_2l_2d_bgcm, ind_1l_3d_bgcm, ind_2l_3d_bgcm)
      ENDIF
    ENDDO

  END SUBROUTINE Collect_bgc_materials_on_tile

#endif
END MODULE mo_jsb_tile
