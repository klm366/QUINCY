!>
!! @brief Contains abstract tile class that contains the surface structure and memory state
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
                                    LAND_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,       ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_abstract
  PUBLIC :: t_jsb_aggregator, t_jsb_aggregator_p, t_jsb_aggregator_weighted_by_fract


!!$  INTEGER, PARAMETER :: max_depth_of_tile_tree = 3

  ! Declarations for aggregators
  TYPE, ABSTRACT :: t_jsb_aggregator
    CHARACTER(len=25) :: name
  CONTAINS
    PROCEDURE(aggregate_interface_2d), DEFERRED             :: Aggregate_2d
    PROCEDURE(aggregate_interface_3d), DEFERRED             :: Aggregate_3d
    GENERIC                                                 :: Aggregate => Aggregate_2d, Aggregate_3d
  END TYPE t_jsb_aggregator

  TYPE t_jsb_aggregator_p
    CLASS(t_jsb_aggregator), POINTER :: p => NULL()
  END TYPE t_jsb_aggregator_p

  ! Declarations for tile
  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct)                        :: lcts                 !< List of land cover types for tile
                                                                   !! The first entry is always the primary LCT for this tile;
                                                                   !! the remaining lcts are those inherited for aggregation
    LOGICAL                                :: is_vegetation = .FALSE.
    LOGICAL                                :: is_bare = .FALSE.
    LOGICAL                                :: is_land = .FALSE.
    LOGICAL                                :: is_glacier = .FALSE.
    LOGICAL                                :: is_lake = .FALSE.
    LOGICAL                                :: contains_vegetation = .FALSE.
    LOGICAL                                :: contains_bare = .FALSE.
    LOGICAL                                :: contains_land = .FALSE.
    LOGICAL                                :: contains_glacier = .FALSE.
    LOGICAL                                :: contains_lake = .FALSE.
    LOGICAL                                :: contains_soil = .FALSE.
    INTEGER,                ALLOCATABLE    :: process_action(:)
    REAL(wp), PRIVATE,          POINTER    :: fract     (:,:)      !< Fraction of this tile (relative to parent tile)
    REAL(wp),                   POINTER    :: fract_max (:,:)      !< Maximum potential fraction for this tile (rel. to parent)
    REAL(wp), PRIVATE,          POINTER    :: fract_old     (:,:)  !< Old fract
    REAL(wp),                   POINTER    :: fract_max_old (:,:)  !< Old fract_max
    CHARACTER(LEN=filename_max)            :: fract_filename       !< Name of file to read tile fractions from. The variable names
                                                                   !! must be either "fract[_max]_"//tile%name, or specified by 
                                                                   !! fract_varname.
    CHARACTER(LEN=20)                      :: fract_varname        !< Name of variable containing fractions for this tile
    LOGICAL,                ALLOCATABLE    :: l_fract_children_changed(:,:)
    TYPE(t_jsb_aggregator_p), ALLOCATABLE  :: aggregators(:)
    !CLASS(t_jsb_param_p),       POINTER    :: param(:)
    CLASS(t_jsb_memory_p),      POINTER    :: mem(:)
    ! The next three pointers basically repeat what's already in the base type, but they save time-consuming "SELECT TYPE"s
    CLASS(t_jsb_tile_abstract), POINTER    :: parent_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: first_child_tile => NULL()
    CLASS(t_jsb_tile_abstract), POINTER    :: next_sibling_tile => NULL()
    !
    INTEGER                                :: grid_id = -1
    INTEGER                                :: owner_model_id = 0
  CONTAINS
    PROCEDURE                              :: Is_process_active
    PROCEDURE                              :: Register_aggregator
    PROCEDURE                              :: Get_aggregator
    PROCEDURE                              :: Get_fraction_1d
    PROCEDURE                              :: Get_fraction_2d
    GENERIC                                :: Get_fraction => Get_fraction_1d, Get_fraction_2d
    PROCEDURE(Init_interface),      DEFERRED :: Init
  END type t_jsb_tile_abstract

  ABSTRACT INTERFACE
    SUBROUTINE Init_interface(this, varlist_name, prefix, suffix, l_fixed_fractions)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
      LOGICAL,                    INTENT(in)    :: l_fixed_fractions
    END SUBROUTINE Init_interface
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE aggregate_interface_2d(this, tile, var, ics, ice, iblk, caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real2d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real2d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
      CHARACTER(len=*),           INTENT(in)    :: caller
    END SUBROUTINE aggregate_interface_2d
    SUBROUTINE aggregate_interface_3d(this, tile, var, ics, ice, iblk,caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real3d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real3d),     INTENT(inout) :: var
      INTEGER,                    INTENT(in)    :: ics, ice, iblk
      CHARACTER(len=*),           INTENT(in)    :: caller
    END SUBROUTINE aggregate_interface_3d
  END INTERFACE

  TYPE, EXTENDS(t_jsb_aggregator) :: t_jsb_aggregator_weighted_by_fract
    REAL(wp), ALLOCATABLE :: fractions(:,:,:)
                                                    !  variable under consideration
  CONTAINS
    PROCEDURE :: Set_fractions => Set_fractions_aggregator
    PROCEDURE :: Aggregate_2d                  => Aggregate_weighted_by_fract_2d
    PROCEDURE :: Aggregate_3d                  => Aggregate_weighted_by_fract_3d
  END TYPE t_jsb_aggregator_weighted_by_fract

  INTERFACE t_jsb_aggregator_weighted_by_fract
    PROCEDURE Create_aggregator_weighted_by_fract
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  FUNCTION Is_process_active(this, iproc) RESULT(return_value)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    INTEGER,                    INTENT(in) :: iproc
    LOGICAL                                :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Is_process_active'

    return_value = .FALSE.

    SELECT TYPE (this)
    CLASS IS (t_jsb_tile_abstract)
      return_value = ASSOCIATED(this%mem(iproc)%p)
    END SELECT

  END FUNCTION Is_process_active

  SUBROUTINE Get_fraction_1d(this, ics, ice, iblk, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:),        &
      & fract_max(:),    &
      & fract_old(:),    &
      & fract_max_old(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_1d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:) = this%fract(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max(:) = this%fract_max(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:) = this%fract_old(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max_old(:) = this%fract_max_old(ics:ice, iblk)
    END IF

  END SUBROUTINE Get_fraction_1d

  SUBROUTINE Get_fraction_2d(this, fract, fract_max, fract_old, fract_max_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:,:),        &
      & fract_max(:,:),    &
      & fract_old(:,:),    &
      & fract_max_old(:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:,:) = this%fract(:,:)
    END IF

    IF (PRESENT(fract_max)) THEN
      IF (SIZE(fract_max,1) /= SIZE(this%fract_max,1) .OR. SIZE(fract_max,2) /= SIZE(this%fract_max,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max(:,:) = this%fract_max(:,:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:,:) = this%fract_old(:,:)
    END IF

    IF (PRESENT(fract_max_old)) THEN
      IF (SIZE(fract_max_old,1) /= SIZE(this%fract_max_old,1) .OR. SIZE(fract_max_old,2) /= SIZE(this%fract_max_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_max_old(:,:) = this%fract_max_old(:,:)
    END IF

  END SUBROUTINE Get_fraction_2d

  SUBROUTINE Register_aggregator(this, aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)        :: this
    CLASS(t_jsb_aggregator),    INTENT(in),   TARGET :: aggregator

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_aggregator'


    IF (.NOT. ALLOCATED(this%aggregators)) THEN
      this%aggregators = (/ t_jsb_aggregator_p(aggregator) /)
    ELSE
      this%aggregators = (/ this%aggregators, t_jsb_aggregator_p(aggregator) /)
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Registered aggregator '//TRIM(aggregator%name)//' for tile '//TRIM(this%name))

  END SUBROUTINE Register_aggregator

  FUNCTION Get_aggregator(this, name, copy) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    CHARACTER(len=*),           INTENT(in) :: name
    LOGICAL, OPTIONAL,          INTENT(in) :: copy
    CLASS(t_jsb_aggregator),    POINTER    :: aggregator

    INTEGER :: i
    LOGICAL :: l_copy

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_aggregator'

    l_copy = .FALSE.
    IF (PRESENT(copy)) l_copy = copy

    aggregator => NULL()
    IF (ALLOCATED(this%aggregators)) THEN
      DO i=1,SIZE(this%aggregators)
        IF (TRIM(this%aggregators(i)%p%name) == TRIM(name)) THEN
          IF (l_copy) THEN
            ALLOCATE(aggregator, source=this%aggregators(i)%p)
          ELSE
            aggregator => this%aggregators(i)%p
          END IF
          EXIT
        END IF
      END DO
    END IF

    IF (.NOT.(ASSOCIATED(aggregator))) &
      & CALL finish(TRIM(routine), 'Aggregator "'//TRIM(name)//'" not found on tile '//TRIM(this%name))

  END FUNCTION Get_aggregator

  FUNCTION Create_aggregator_weighted_by_fract(tile) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)          :: tile
    CLASS(t_jsb_aggregator_weighted_by_fract), POINTER :: aggregator

    INTEGER                   :: no_children

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Create_aggregator_weighted_by_fract'

    ALLOCATE(t_jsb_aggregator_weighted_by_fract::aggregator)
    aggregator%name = 'weighted_by_fract'

    IF (debug_on()) CALL message(TRIM(routine), 'Created aggregator '//TRIM(aggregator%name)//' on tile '//TRIM(tile%name))

  END FUNCTION Create_aggregator_weighted_by_fract

  SUBROUTINE Set_fractions_aggregator(this, tile, ics, ice, iblk)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk

    INTEGER :: no_children, i
    CLASS(t_State), POINTER :: current
    REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_fractions_aggregator'

    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(ics:ice,iblk))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child 
      DO i=1,no_children
        SELECT TYPE (current)
        CLASS IS (t_jsb_tile_abstract)
          CALL current%Get_fraction(ics, ice, iblk, fract=this%fractions(ics:ice,iblk,i))
        END SELECT
        current => current%next_sibling
      END DO

      ALLOCATE(fract(ice-ics+1))
      ALLOCATE(fract_sum(ice-ics+1))
      fract(:) = tile%fract(ics:ice,iblk)
      fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      !     &          all_print=.TRUE.)
      !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! END IF
      DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(ics:ice,iblk) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator

  SUBROUTINE Aggregate_weighted_by_fract_2d(this, tile, var, ics, ice, iblk, caller)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real2d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk
    CHARACTER(len=*),                          INTENT(in)    :: caller

    INTEGER :: no_children, i, j, iproc, nc
    ! CLASS(t_jsb_tile_abstract), POINTER :: current
    ! TYPE(t_jsb_var_real2d), POINTER :: var_tmp
    CLASS(t_jsb_var), POINTER            :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp)              :: fract_sum(ice-ics+1)
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'
    END IF

    IF (.NOT. ALLOCATED(var%ptr)) RETURN

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()
    if (no_children == 0) return

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ALLOCATE(data     (nc,no_children))
    ALLOCATE(fractions(nc,no_children))

    j = 0
    fractions(:,:) = 0._wp
    fract_sum(:)   = 0._wp
    data     (:,:) = 0._wp

    ! current => tile%first_child_tile
    DO i=1,no_children
      ! IF (ASSOCIATED(current%mem(iproc)%p)) THEN  ! Process is active on child tile
      !   SELECT TYPE (v => current%mem(iproc)%p%Get_var_real2d(TRIM(var%name)))
      !   TYPE IS (t_jsb_var_real2d)
      !     var_tmp => v
      !   END SELECT
      ! IF (var%child_idx(i) > 0) THEN
      !   var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
      ! ELSE
      !   var_tmp => NULL()
      ! END IF
      ! IF (ASSOCIATED(var_tmp)) THEN
      IF (var%child_idx(i) > 0 .AND. ASSOCIATED(tile%mem(iproc)%p)) THEN
        IF (.NOT. ASSOCIATED(tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p)) THEN
          CALL finish(TRIM(routine), 'ERROR - this should not happen!')
        ELSE
          var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
        END IF
        fractions(:,i) = this%fractions(ics:ice,iblk,i)
        IF (.NOT. var%l_aggregate_all) fract_sum(:)   = fract_sum(:) + fractions(:,i)
        data(:,i) = var_tmp%field(ics:ice,iblk,1,1,1)
        ! data(:,i) = var_tmp%ptr(ics:ice,iblk)
        j = j + 1
      END IF
      ! END IF
      ! current => current%next_sibling_tile
    END DO
    j=no_children
    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   ! print*, fract_sum
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    IF (j == 1 .AND. .NOT. var%l_aggregate_all) THEN
      ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction 
      var%ptr(ics:ice,iblk) = SUM(data(:,:), DIM=2)
    ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
      ! For tiles where variable is not defined data(:,tile) is 0
      var%ptr(ics:ice,iblk) = SUM(fractions(:,:) * data(:,:), DIM=2)
    ELSE
      WHERE (fract_sum(:) > 0._wp)
        var%ptr(ics:ice,iblk) = SUM(fractions(:,:) * data(:,:), DIM=2) / fract_sum(:)
      !ELSE WHERE
      !  var%ptr(ics:ice,iblk) = missval
      END WHERE
    END IF

    DEALLOCATE(data, fractions)

    IF (var%l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_2d

  SUBROUTINE Aggregate_weighted_by_fract_3d(this, tile, var, ics, ice, iblk, caller)

    USE mo_hsm_class,      ONLY: t_State

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile
    TYPE(t_jsb_var_real3d),                    INTENT(inout) :: var
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk
    CHARACTER(len=*),                          INTENT(in)    :: caller

    INTEGER :: no_children, i, j, iproc, nc, ilev
    ! CLASS(t_jsb_tile_abstract), POINTER :: current
    ! TYPE(t_jsb_var_real3d), POINTER :: var_tmp
    CLASS(t_jsb_var), POINTER :: var_tmp
    REAL(wp), ALLOCATABLE :: data(:,:,:), fractions(:,:)    ! Collected data and fractions from children
    REAL(wp)              :: fract_sum(ice-ics+1)
    REAL(dp) :: missval

    CHARACTER(len=256) :: routine

    IF (debug_on()) THEN
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = TRIM(modname)//':Aggregate_weighted_by_fract_2d'
    END IF

    IF (.NOT. ALLOCATED(var%ptr)) RETURN

    CALL this%Set_fractions(tile, ics, ice, iblk)

    nc = ice - ics + 1

    missval = var%missval

    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    no_children = tile%Get_no_of_children()

    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ALLOCATE(data     (nc, SIZE(var%ptr,2), no_children))
    ALLOCATE(fractions(nc,                  no_children))

    j = 0
    fractions(:,:)   = 0._wp
    fract_sum(:  )   = 0._wp
    data     (:,:,:) = 0._wp
    ! current => tile%first_child_tile
    DO i=1,no_children
      ! IF (ASSOCIATED(current%mem(iproc)%p)) THEN   ! Process is active on child tile
      !   SELECT TYPE (v => current%mem(iproc)%p%Get_var_real2d(TRIM(var%name)))
      !   TYPE IS (t_jsb_var_real3d)
      !     var_tmp => v
      !   END SELECT
        ! var_tmp => current%mem(iproc)%p%Get_var_real3d(TRIM(var%name))
        ! IF (ASSOCIATED(var_tmp)) THEN
      IF (var%child_idx(i) > 0 .AND. ASSOCIATED(tile%mem(iproc)%p)) THEN
        IF (.NOT. ASSOCIATED(tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p)) THEN
          CALL finish(TRIM(routine), 'ERROR - this should not happen!')
        ELSE
          var_tmp => tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p
        END IF
        fractions(:,i) = this%fractions(ics:ice,iblk,i)
        IF (.NOT. var%l_aggregate_all) fract_sum(:)   = fract_sum(:) + fractions(:,i)
        data(:,:,i) = var_tmp%field(ics:ice,:,iblk,1,1)
        ! data(:,:,i) = var_tmp%ptr(ics:ice,:,iblk)
        j = j + 1
      END IF
      ! END IF
      ! current => current%next_sibling_tile
    END DO
    j=no_children
    IF (j == 0) CALL finish(TRIM(routine), 'Variable '//TRIM(var%name)//' not found in any of children of tile '//&
      &                                    TRIM(tile%name)//' - this should not happen!')

    ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
    !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
    !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
    !     &          all_print=.TRUE.)
    !   CALL finish(TRIM(routine),  '    ... this should not happen')
    ! END IF

    DO ilev=1,SIZE(data,2)
      IF (j == 1 .AND. .NOT. var%l_aggregate_all) THEN
        ! Only one tile to aggregate, avoid multiplying and dividing by a potentially small fraction 
        var%ptr(ics:ice,ilev,iblk) = SUM(data(:,ilev,:), DIM=2)
      ELSE IF (var%l_aggregate_all) THEN              ! Sum of fractions is 1
        ! For tiles where variable is not defined data(:,ilev,tile) is 0
        var%ptr(ics:ice,ilev,iblk) = SUM(fractions(:,:) * data(:,ilev,:), DIM=2)
      ELSE
        WHERE (fract_sum(:) > 0._wp)
          var%ptr(ics:ice,ilev,iblk) = SUM(fractions(:,:) * data(:,ilev,:), DIM=2) / fract_sum(:)
        ! ELSE WHERE
          ! var%ptr(ics:ice,ilev,iblk) = missval
        END WHERE
      END IF
    END DO

    DEALLOCATE(data, fractions)

    IF (var%l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_3d

#endif
END MODULE mo_jsb_tile_class
