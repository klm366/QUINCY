!> Summary: Methods for lcc processes
!>
!>### Contains methods required to deal with actively or passively relocated conserved quantities upon lcc
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-04-12 - present)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_jsb_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: finish, message
  USE mo_util,                ONLY: int2string
  USE mo_kind,                ONLY: wp

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN

  USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: end_lcc_reloc, start_lcc_reloc, init_lcc_reloc, init_lcc, transfer_from_active_to_passive_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lcc'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initializes/sets up the lcc instance for and on the given process
  !
  SUBROUTINE init_lcc(process_name, tile, activeCQTs, passiveCQTs, involved_child_tiles)

    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
    USE mo_jsb_process_class, ONLY: Get_process_id
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_jsb_cqt_class,     ONLY: Get_cqt_name
    USE mo_jsb_lcc_class,     ONLY: t_jsb_lcc_var

    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*),  INTENT(in) :: process_name     !< name of the encompassing lcc process  
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< current tile on which the lcc process runs
    INTEGER,           INTENT(in) :: activeCQTs(:)    !< cqts that are actively relocated by this lcc process
    INTEGER,           INTENT(in) :: passiveCQTs(:)   !< cqts that are passively relocated by this lcc process
    CHARACTER(len=*),  INTENT(in) :: involved_child_tiles(:) !< child tiles involved in this lcc process
    ! JN-TODO: probably not only direct child tiles, i.e. descendants?
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc'

    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid),  POINTER          :: hgrid  
    TYPE(t_jsb_lcc_proc), POINTER       :: lcc_relocations
    CLASS(t_jsb_tile_abstract), POINTER :: child
    INTEGER :: process_id, nr_of_children, nr_of_active_vars, nr_of_passive_vars, nr_of_tiles, nr_of_vars
    INTEGER :: ind, i_child, i_cqt, i_var, i_active, i_active_this_tile, i_passive, i_passive_this_tile, i_tile
    INTEGER :: nproma, nblks
    LOGICAL :: is_involved_tile, has_searched_var
    LOGICAL, ALLOCATABLE :: is_child_tile(:)
    TYPE(t_jsb_lcc_var), POINTER  :: a_lcc_var
    CHARACTER(len=VARNAME_LEN) :: thisName

    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    hgrid => Get_grid(model%grid_id)
    process_id = Get_process_id(process_name)
    lcc_relocations => model%processes(process_id)%p%lcc_relocations

    IF (debug_on()) CALL message(TRIM(routine), 'For proc '//TRIM(process_name)//' ...')

    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    nr_of_tiles = SIZE(involved_child_tiles(:))
    !>
    !> Assertion: number of involved child tiles should be larger than one
    !>
    IF (.NOT. nr_of_tiles > 1) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: number of tiles involved in a lcc process need to be larger than one! '&
        & //'But was: '// int2string(nr_of_tiles) //'. Please check.')
    ENDIF

    lcc_relocations%nr_of_tiles = nr_of_tiles
    ALLOCATE(is_child_tile(nr_of_tiles))
    is_child_tile = .FALSE.

    lcc_relocations%name = TRIM(process_name) // "_lcc"
    lcc_relocations%tile_names = involved_child_tiles(:)

    !>
    !> Assertion: active and passive cqts are not allowed to overlap
    !>
    DO i_cqt = 1,SIZE(activeCQTs(:))
      IF (ANY(activeCQTs(i_cqt) == passiveCQTs(:) )) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Conceptually each CQT type can only be active or passive ' &
            & //'for one lcc proc! If further distinctions are required new sub-types need to be implemented. ' &
            & //'CQ type which was active and passive: '// Get_cqt_name(activeCQTs(i_cqt)) //'. Please check.')
      ENDIF
    ENDDO

    nr_of_active_vars = 0
    nr_of_passive_vars = 0

    ! Iterate over the children of the tile to collect the CQ variables of 
    ! the actively and passively relocated CQTs
    nr_of_children = tile%Get_no_of_children()
    child => tile%Get_first_child_tile()
    i_tile = 0
    DO i_child = 1,nr_of_children

      ! only child tiles listed as involved are used
      is_involved_tile = .FALSE.
      DO ind = 1,nr_of_tiles
        IF (TRIM(child%name) .EQ. TRIM(involved_child_tiles(ind))) THEN
          i_tile = i_tile + 1
          is_involved_tile = .TRUE.
          is_child_tile(ind) = .TRUE.
          EXIT
        ENDIF
      ENDDO

      ! The lcc structure is created according to active and passive cqs on the first involved tile
      IF (is_involved_tile .AND. (i_tile == 1)) THEN

        ! In the first step only collect the number of variables
        DO i_cqt = 1, child%numberOfConservedQuanTypes
          IF (ANY(child%conserved_quantities(i_cqt)%p%type_id == activeCQTs(:) )) THEN
            nr_of_active_vars = nr_of_active_vars &
              & + child%conserved_quantities(i_cqt)%p%no_of_vars
          ELSEIF (ANY(child%conserved_quantities(i_cqt)%p%type_id == passiveCQTs(:) )) THEN
            nr_of_passive_vars = nr_of_passive_vars &
              & + child%conserved_quantities(i_cqt)%p%no_of_vars
          ENDIF
        ENDDO

        lcc_relocations%nr_of_active_vars = nr_of_active_vars
        lcc_relocations%nr_of_passive_vars = nr_of_passive_vars

        ! allocate
        ALLOCATE(lcc_relocations%active_vars(nr_of_active_vars))
        ALLOCATE(lcc_relocations%active_vars_cqt(nr_of_active_vars))
        ALLOCATE(lcc_relocations%passive_vars(nr_of_passive_vars))
        ALLOCATE(lcc_relocations%passive_vars_cqt(nr_of_passive_vars))

        ! and then collect vars and pointers for the first involved child
        i_active = 0
        i_passive = 0
        DO i_cqt = 1, child%numberOfConservedQuanTypes
          ! for the active cqts
          IF (ANY(child%conserved_quantities(i_cqt)%p%type_id == activeCQTs(:) )) THEN
            nr_of_vars = child%conserved_quantities(i_cqt)%p%no_of_vars
            DO i_var = 1,nr_of_vars
              i_active = i_active + 1

              lcc_relocations%active_vars_cqt(i_active) &
                & = child%conserved_quantities(i_cqt)%p%type_id

              ALLOCATE(a_lcc_var)
              lcc_relocations%active_vars(i_active)%p => a_lcc_var
              lcc_relocations%active_vars(i_active)%p%name &
                & = child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p%name

              ALLOCATE(lcc_relocations%active_vars(i_active)%p%tile_pointer(nr_of_tiles))
              lcc_relocations%active_vars(i_active)%p%tile_pointer(i_tile)%p &
                & => child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p

              ALLOCATE(lcc_relocations%active_vars(i_var)%p%relocate_this(nr_of_tiles, nproma, nblks))
              lcc_relocations%active_vars(i_var)%p%relocate_this(:,:,:) = 0.0_wp

            ENDDO 
          ! and separately for the passive cqts
          ELSEIF (ANY(child%conserved_quantities(i_cqt)%p%type_id == passiveCQTs(:) )) THEN
            nr_of_vars = child%conserved_quantities(i_cqt)%p%no_of_vars
            DO i_var = 1,nr_of_vars
              i_passive = i_passive + 1

              lcc_relocations%passive_vars_cqt(i_passive) &
                & = child%conserved_quantities(i_cqt)%p%type_id

              ALLOCATE(a_lcc_var)
              lcc_relocations%passive_vars(i_passive)%p => a_lcc_var
              lcc_relocations%passive_vars(i_passive)%p%name &
                & = child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p%name

              ALLOCATE(lcc_relocations%passive_vars(i_passive)%p%tile_pointer(nr_of_tiles))
              lcc_relocations%passive_vars(i_passive)%p%tile_pointer(i_tile)%p &
                & => child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p

              ALLOCATE(lcc_relocations%passive_vars(i_passive)%p%relocate_this(nr_of_tiles, nproma, nblks))
              lcc_relocations%passive_vars(i_passive)%p%relocate_this(:,:,:) = 0.0_wp

            ENDDO 
          ENDIF
        ENDDO
      ENDIF

      ! For all other involved tiles (i.e. other than the first) only the pointers need to be collected
      ! + assert that the demanded active and passive variables are part of all involved tiles
      IF (is_involved_tile .AND. (i_tile > 1)) THEN

        i_active_this_tile = 0
        i_passive_this_tile = 0

        DO i_cqt = 1, child%numberOfConservedQuanTypes
          IF (ANY(child%conserved_quantities(i_cqt)%p%type_id == activeCQTs(:) )) THEN
            nr_of_vars = child%conserved_quantities(i_cqt)%p%no_of_vars
            DO i_var = 1,nr_of_vars
              i_active_this_tile = i_active_this_tile + 1

              thisName = child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p%name
              has_searched_var = .FALSE.
              DO i_active = 1,nr_of_active_vars
                IF(TRIM(lcc_relocations%active_vars(i_active)%p%name) .EQ. TRIM(thisName)) THEN
                  has_searched_var = .TRUE.
                  lcc_relocations%active_vars(i_active)%p%tile_pointer(i_tile)%p &
                    & => child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p
                  EXIT
                ENDIF
              ENDDO
              !>
              !> Assertion: all passive and active variables need to exist on all tiles
              !>
              IF (.NOT. has_searched_var) THEN
                CALL finish(TRIM(routine), &
                  & 'Violation of assertion: Conceptually each tile needs to contain the same active and passive vars. ' &
                  & //'Passive var ' // TRIM(thisName) // ' does not exist on tile ' &
                  & //TRIM(involved_child_tiles(1))//'. Please check.' )
              ENDIF
            ENDDO
          ELSEIF(ANY(child%conserved_quantities(i_cqt)%p%type_id == passiveCQTs(:) )) THEN
            nr_of_vars = child%conserved_quantities(i_cqt)%p%no_of_vars
            DO i_var = 1,nr_of_vars
              i_passive_this_tile = i_passive_this_tile + 1

              thisName = child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p%name
              has_searched_var = .FALSE.
              DO i_passive = 1,nr_of_passive_vars
                IF(TRIM(lcc_relocations%passive_vars(i_passive)%p%name) .EQ. TRIM(thisName)) THEN
                  has_searched_var = .TRUE.
                  lcc_relocations%passive_vars(i_passive)%p%tile_pointer(i_tile)%p &
                    & => child%conserved_quantities(i_cqt)%p%consQuanVars2D(i_var)%p
                  EXIT
                ENDIF
              ENDDO

              ! Assertion: all passive and active variables need to exist on all tiles
              IF (.NOT. has_searched_var) THEN
                CALL finish(TRIM(routine), &
                  & 'Violation of assertion: Conceptually each tile needs to contain the same active and passive vars. ' &
                  & //'Active var ' // TRIM(thisName) // ' does not exist on tile ' &
                  & //TRIM(involved_child_tiles(1))//'. Please check.' )
              ENDIF
            ENDDO
          ENDIF
        ENDDO 

        ! Assertion: all variables found on the first tile need to also be cqs for all other tiles
        IF (.NOT. i_active_this_tile == nr_of_active_vars ) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Not all active variables found on first tile '//TRIM(involved_child_tiles(1)) &
            & //' were also active cqs of this tile '//TRIM(involved_child_tiles(i_tile))//'. Found ' &
            & //int2string(nr_of_active_vars)//' on first and '//int2string(i_active_this_tile)//' on this. Please check.')
        ENDIF
        IF (.NOT. i_passive_this_tile == nr_of_passive_vars ) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Not all passive variables found on first tile '//TRIM(involved_child_tiles(1)) &
            & //' were also passive cqs of this tile '//TRIM(involved_child_tiles(i_tile))//'. Found ' &
            & //int2string(nr_of_passive_vars)//' on first and '//int2string(i_passive_this_tile)//' on this. Please check.')
        ENDIF
      ENDIF

      child => child%Get_next_sibling_tile()
    ENDDO

    !>
    !> Assertion: all involved tiles are actually child tiles
    !>
    ! JN-TODO: also allow descendant tiles
    IF (.NOT. ALL(is_child_tile)) THEN
      DO i_tile = 1,nr_of_tiles
        IF (.NOT. is_child_tile(i_tile)) THEN
          CALL finish(TRIM(routine), &
            & 'Violation of assertion: Conceptually each involved tile needs to be a child tile of the tile on '&
            & // 'which the lcc proc is running -- did not found ' &
            & // TRIM(involved_child_tiles(i_tile))//'. Please check.' )
        ENDIF
      ENDDO
    ENDIF

    DEALLOCATE(is_child_tile)

  END SUBROUTINE init_lcc

  ! ====================================================================================================== !
  !
  !> Init the relocation upon lcc process task update (currently: collect initial areas and reset relocation arrays)
  !
  SUBROUTINE init_lcc_reloc(lcc_relocations, options, tile, initial_area)

    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations  !< lcc structure of the calling lcc process
    TYPE(t_jsb_task_options),   INTENT(in) :: options             !< run-time parameters
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile               
        !< tile on which the lcc process is running for which the lcc structure is initialized here
    REAL(wp),                   INTENT(inout) :: initial_area(:,:)  !< current area of involved tiles 
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_lcc_reloc'

    CLASS(t_jsb_tile_abstract), POINTER :: current_tile
    REAL(wp), DIMENSION(options%nc) :: current_fract
    INTEGER :: i_tile, i_var, iblk, ics, ice, nc, nr_of_tiles
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles

    !>
    !> Assertion: initial area size and number of tiles need to suit 
    !>
    IF ((SIZE(initial_area,DIM=1) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=2) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area array has unexpected size. Please check!')
    ENDIF

    ! Collect initial areas of all involved tiles
    current_tile => tile%first_child_tile
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      IF (ANY(current_tile%name .EQ. lcc_relocations%tile_names)) THEN !JN-TODO: does this work or TRIM and DO loop?
        i_tile = i_tile + 1
        
        CALL current_tile%Get_fraction(ics, ice, iblk, fract=current_fract(:))

        initial_area(i_tile, :) = current_fract(:)
      ENDIF

      current_tile => current_tile%next_sibling_tile
    ENDDO

    ! reset relocation arrays
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      lcc_relocations%active_vars(i_var)%p%relocate_this(:,1:nc,iblk) = 0.0_wp
    ENDDO
    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      lcc_relocations%passive_vars(i_var)%p%relocate_this(:,1:nc,iblk) = 0.0_wp
    ENDDO

  END SUBROUTINE init_lcc_reloc

  ! ====================================================================================================== !
  !
  !> Start the relocation upon lcc (collect to be transferred matter from all active and passive vars)
  !
  SUBROUTINE start_lcc_reloc(lcc_relocations, options, lost_area, initial_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for collecting to be transferred matter 
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< run-time parameters
    REAL(wp),                   INTENT(in)    :: lost_area(:,:)    !< area moved from involved tiles
    REAL(wp),                   INTENT(in)    :: initial_area(:,:) !< initial area of involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':start_lcc_reloc'

    INTEGER :: i_var, i_tile, iblk, ics, ice, nc, nr_of_tiles
    REAL(wp) :: to_be_transferred(options%nc)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles
    !>
    !> Assertion: dimensions of lost_area and initial_area need to match
    !>
    IF ((SIZE(lost_area,DIM=1) /= nr_of_tiles) .OR. (SIZE(initial_area,DIM=1) /= nr_of_tiles) &
      & .OR. (SIZE(lost_area,DIM=2) /= nc) .OR. (SIZE(initial_area,DIM=2) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area arrays have unexpected size. Please check!')
    ENDIF
    !>
    !> Assertion: lost area is >= 0.0_wp and <= 1.0_wp
    !>
    IF (ANY(lost_area < 0.0_wp) .OR. ANY(lost_area > 1.0_wp) &
      & .OR. ANY(initial_area < 0.0_wp) .OR. ANY(initial_area > 1.0_wp)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area arrays have values out of range ([0.0, 1.0]). Please check!')
    ENDIF
    !>
    !> Assertion: lost area <= initial_area 
    !>
    ! -> JN-TODO: we might want/need to discuss the concept of the e-10 areas here
    IF (ANY(lost_area > initial_area)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Lost area exceeds available area! Please check!')
    ENDIF
    ! Iterate over all active and passive pools to identify the matter that needs to be transferred
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      DO i_tile = 1,nr_of_tiles
        to_be_transferred(:) = 0.0_wp
        WHERE(lost_area(i_tile, :) > 0.0_wp)
          ! collect to be transferred matter
          to_be_transferred(:) = lcc_relocations%active_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            &   * (lost_area(i_tile, :) / initial_area(i_tile,:))
          lcc_relocations%active_vars(i_var)%p%relocate_this(i_tile,1:nc,iblk) = to_be_transferred(:)
          ! and subtract it 
          lcc_relocations%active_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            & = lcc_relocations%active_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            &   - to_be_transferred(:)
        END WHERE
      ENDDO
    ENDDO
    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      DO i_tile = 1,nr_of_tiles
        to_be_transferred(:) = 0.0_wp
        WHERE(lost_area(i_tile, :) > 0.0_wp)
          ! collect to be transferred matter
          to_be_transferred(:) = lcc_relocations%passive_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
          &   * (lost_area(i_tile, :) / initial_area(i_tile,:))
          lcc_relocations%passive_vars(i_var)%p%relocate_this(i_tile,1:nc,iblk) = to_be_transferred(:)
          ! and subtract it 
          lcc_relocations%passive_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            & = lcc_relocations%passive_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            &   - to_be_transferred(:)
        END WHERE
      ENDDO
    ENDDO
    
  END SUBROUTINE start_lcc_reloc

  ! ====================================================================================================== !
  !
  !> End the relocation upon lcc (distribute all to be transferred matter - passive relocation only!)
  !>
  !> Note: the concept here is that the to be distributed area is collected of ALL shrinking tiles
  !>       i.e. if there are several shrinking tiles in the call
  !>            then there is no directed movement from tile_x -> tile_y
  !>            but sum (shrinking tiles) -> per extending tiles prop. to sum (extending tiles)
  !
  SUBROUTINE end_lcc_reloc(lcc_relocations, options, gained_area)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing to be relocated matter 
    TYPE(t_jsb_task_options),   INTENT(in)    :: options           !< run-time parameters
    REAL(wp),                   INTENT(in)    :: gained_area(:,:)  !< area moved to involved tiles
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':end_lcc_reloc'

    INTEGER :: i_var, i_tile, iblk, ics, ice, nc, nr_of_tiles
    REAL(wp) :: passive_var_sum(options%nc), gained_area_sum(options%nc)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    nr_of_tiles = lcc_relocations%nr_of_tiles
    gained_area_sum = SUM(gained_area(:, :), DIM=1)

    !>
    !> Assertion: dimensions of gained_area need to be tiles x nc
    !>
    IF ((SIZE(gained_area,DIM=1) /= nr_of_tiles) .OR. (SIZE(gained_area,DIM=2) /= nc)) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Area array has unexpected size. Please check!')
    ENDIF
    !>
    !> Assertion: reloc matter of active variables expected to be zero (i.e. relocated before routine is called)
    !>
    DO i_var = 1,lcc_relocations%nr_of_active_vars
      DO i_tile = 1,nr_of_tiles
        IF(ANY(lcc_relocations%active_vars(i_var)%p%relocate_this(i_tile,1:nc,iblk) /= 0.0_wp)) THEN
          CALL finish(TRIM(routine), 'Violation of assertion: All active variables need to already have been ' &
            & // 'relocated upon call of this routine. However, found != 0.0 values for active var '&
            & // TRIM(lcc_relocations%active_vars(i_var)%p%name) // ' on tile ' &
            & // TRIM(lcc_relocations%tile_names(i_tile)) // '. Please check!')
        ENDIF
      ENDDO
    ENDDO

    ! Passive variables are proportionally distributed
    DO i_var = 1,lcc_relocations%nr_of_passive_vars
      passive_var_sum = SUM(lcc_relocations%passive_vars(i_var)%p%relocate_this(:,1:nc,iblk), DIM=1)
      DO i_tile = 1,nr_of_tiles
        WHERE(gained_area(i_tile, :) > 0.0_wp)
          lcc_relocations%passive_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            & = lcc_relocations%passive_vars(i_var)%p%tile_pointer(i_tile)%p%ptr2d(ics:ice,iblk) &
            & + passive_var_sum(:) * (gained_area(i_tile,:) / gained_area_sum(:))
        END WHERE
      ENDDO
      lcc_relocations%passive_vars(i_var)%p%relocate_this(:,1:nc,iblk) = 0.0_wp
    ENDDO

  END SUBROUTINE end_lcc_reloc

  ! ====================================================================================================== !
  !
  !> Relocates content of a given active to a given passive var
  !
  SUBROUTINE transfer_from_active_to_passive_var(lcc_relocations, options, active_var_name, passive_var_name)

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(inout) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var 
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< run-time parameters
    CHARACTER(len=*), INTENT(in)    :: active_var_name   !< (active) source var
    CHARACTER(len=*), INTENT(in)    :: passive_var_name  !< (passive) sink var
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':transfer_from_active_to_passive_var'

    INTEGER :: i_active, i_passive, i_tile, nc, ics, ice, iblk
    LOGICAL :: is_active, is_passive
    ! -------------------------------------------------------------------------------------------------- !
    !JN-TODO: is it less expensive to only do it for tiles with actual area changes?
    nc = options%nc
    ics = options%ics
    ice = options%ice
    iblk = options%iblk

    is_active = .FALSE.
    is_passive = .FALSE.
    DO i_active = 1,lcc_relocations%nr_of_active_vars
      IF(lcc_relocations%active_vars(i_active)%p%name .EQ. active_var_name) THEN
        is_active = .TRUE.
        DO i_passive = 1,lcc_relocations%nr_of_passive_vars
          IF(lcc_relocations%passive_vars(i_passive)%p%name .EQ. passive_var_name) THEN
            is_passive = .TRUE.
              DO i_tile = 1,lcc_relocations%nr_of_tiles
                ! Add active content to passive var 
                lcc_relocations%passive_vars(i_passive)%p%relocate_this(i_tile,1:nc,iblk) &
                  & = lcc_relocations%passive_vars(i_passive)%p%relocate_this(i_tile,1:nc,iblk) &
                  & + lcc_relocations%active_vars(i_active)%p%relocate_this(i_tile,1:nc,iblk)
                ! And set active to zero
                lcc_relocations%active_vars(i_active)%p%relocate_this(i_tile,1:nc,iblk) = 0.0_wp
              ENDDO

            EXIT
          ENDIF
        ENDDO

        EXIT
      ENDIF
    ENDDO

    !>
    !> Assertion source variable needs to be active variable for this lcc process
    !>
    IF (.NOT. is_active) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed active var ' // TRIM(active_var_name) &
        & // ' is not an active var in lcc structure for ' // lcc_relocations%name )
    ENDIF
    !>
    !> Assertion sink variable needs to be passive variable for this lcc process
    !>
    IF (.NOT. is_passive) THEN
      CALL finish(TRIM(routine), 'Violation of assertion: Passed passive var ' // TRIM(passive_var_name) &
        & // ' is not a passive var in lcc structure for ' // lcc_relocations%name )
    ENDIF

  END SUBROUTINE transfer_from_active_to_passive_var

#endif
END MODULE mo_jsb_lcc
