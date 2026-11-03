!> Contains class definitions for biogeochemical material (bgcm)
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
!>### define bgcm class, indices (enumerator) and helper functions for bgc material and their components
!>
MODULE mo_lnd_bgcm_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: lp
  USE mo_util,                ONLY: int2string
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_jsb_var_class,       ONLY: t_jsb_var, t_jsb_var_p, t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lnd_bgcm, t_lnd_bgcm_p
  PUBLIC :: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, FIRST_ELEM_ID, LAST_ELEM_ID

  PUBLIC :: get_name_for_element_id, get_max_number_of_elements, get_element_index, get_short_name_for_element_id
  PUBLIC :: assert_no_compartment_bgc_materials, get_dimension_of_element_variables

  ! ======================================================================================================= !
  !>
  !>TYPE for bgcm - biogeochemical materials
  !>
  TYPE, ABSTRACT :: t_lnd_bgcm
    CHARACTER(LEN=:),  ALLOCATABLE :: shortname                                ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: name                                     ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: path                                     !
    CHARACTER(LEN=:),  ALLOCATABLE :: element_unit                             ! unit of bgcm element variables (used for output stream)
    INTEGER                        :: id                            = -1
      !< a bgcm id that - if set /= -1 - should be unique (from mo_lnd_bgcm_store_class)
    INTEGER                        :: owner_model_id                = 0        ! ID of the jsbach model this bgcm is running at
    LOGICAL                        :: contains_elements             = .FALSE.  ! contains element variables (e.g., Carbon, Nitrogen)
    LOGICAL                        :: contains_bgcm                 = .FALSE.  ! contains bgc_material (t_lnd_bgcm)
    LOGICAL                        :: contains_elem_1d_ptr_to_chunk = .FALSE.  ! t_lnd_bgcm_elem2d -> set TRUE in the bgcm init routines in the PROC_memory_class
    LOGICAL                        :: contains_elem_2d_ptr_to_chunk = .FALSE.  ! t_lnd_bgcm_elem3d -> set TRUE in the bgcm init routines in the PROC_memory_class
    CLASS(t_lnd_bgcm_p),   POINTER :: bgcm_list(:) => NULL()                   ! children
    CLASS(t_lnd_bgcm_p),   POINTER :: leaf_list(:) => NULL()                   ! leaves
    TYPE(t_jsb_var_p), ALLOCATABLE :: element_list(:)                          ! elements of this bgcm (elements used in PROC are defined in 'all_elements(:)' in PROC_memory_class)
    TYPE(t_jsb_var_p), ALLOCATABLE :: var_list(:)                              ! other associated variables - currently not used
  CONTAINS
    PROCEDURE(init_iface), DEFERRED :: Init
    PROCEDURE(get_el_name_iface_id), DEFERRED :: Get_element_name_by_id
    GENERIC                         :: Get_element_name => Get_element_name_by_id
    PROCEDURE                       :: Add_element           => t_lnd_bgcm_add_element
    PROCEDURE                       :: Add_var               => t_lnd_bgcm_add_var
    PROCEDURE                       :: Has_elements          => t_lnd_bgcm_has_elements
    PROCEDURE                       :: Add_bgcm              => t_lnd_bgcm_add_bgcm
    PROCEDURE                       ::                          t_lnd_bgcm_find_element_position_by_id
    PROCEDURE                       ::                          t_lnd_bgcm_find_element_position_by_name
    PROCEDURE                       :: Set_paths             => t_lnd_bgcm_set_paths
    GENERIC                         :: Find_element_position => t_lnd_bgcm_find_element_position_by_id,   &
      &                                                         t_lnd_bgcm_find_element_position_by_name
  END TYPE t_lnd_bgcm

  ! ======================================================================================================= !
  !>
  !>TYPE for bgcm pointer
  !>
  TYPE :: t_lnd_bgcm_p
    CLASS(t_lnd_bgcm), POINTER :: p
  END TYPE t_lnd_bgcm_p

  ! ======================================================================================================= !
  !>
  !>Interface for init routine of the bgcm TYPE
  !>
  ABSTRACT INTERFACE
    SUBROUTINE init_iface(this, name, shortname, elements_index_map, l_elements, element_list, element_unit)
      IMPORT t_lnd_bgcm
      CLASS(t_lnd_bgcm),          INTENT(inout) :: this
      CHARACTER(LEN=*),           INTENT(in)    :: name
      CHARACTER(LEN=*),           INTENT(in)    :: shortname
      INTEGER,                    INTENT(in)    :: elements_index_map(:)
      LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements      ! Add elements yes or no?
      INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:) ! Select elements to add
      CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit
    END SUBROUTINE init_iface
    FUNCTION get_el_name_iface_id(this, id) RESULT(name)
      IMPORT t_lnd_bgcm
      CLASS(t_lnd_bgcm),          INTENT(in) :: this
      INTEGER,           INTENT(in) :: id
      CHARACTER(LEN=:), ALLOCATABLE :: name
    END FUNCTION get_el_name_iface_id
  END INTERFACE

  ! ======================================================================================================= !
  !>
  !>IDs of elements, Carbon, Nitrogen, Phosphorous and their isotopes
  !>
  !>  (since elements are dynamic these IDs can only be used as idx in arrays with all possible elements)
  !>
  ENUM, BIND(C)
    ENUMERATOR ::           &
      & ELEM_C_ID = 1,      &
      & ELEM_N_ID,          &
      & ELEM_P_ID,          &
      & ELEM_C13_ID,        &
      & ELEM_C14_ID,        &
      & ELEM_N15_ID,        &
      & ONEPLUS_LAST_ELEM_ID  ! needs to be the last -- it is only used to determine the max number of elements
  END ENUM

  !> ID of first and last bgcm elements
  INTEGER, PARAMETER :: FIRST_ELEM_ID = ELEM_C_ID
  INTEGER, PARAMETER :: LAST_ELEM_ID  = ONEPLUS_LAST_ELEM_ID - 1

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_class'

CONTAINS

  ! ====================================================================================================== !
  !>
  !> Returns the number of elements defined in the enumerator
  !>
  FUNCTION get_max_number_of_elements() RESULT(last_type_id)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: last_type_id
    ! -------------------------------------------------------------------------------------------------- !
    ! The elem id enum  starts with 1 to be able to use the ids as indices in arrays with all possible elements
    last_type_id = LAST_ELEM_ID

  END FUNCTION get_max_number_of_elements


  ! ======================================================================================================= !
  !>
  !> Helper: assert that bgc material with elements does not contain further bgc materials
  !>
  SUBROUTINE assert_no_compartment_bgc_materials(this, tile_name, routine)
    USE mo_exception,           ONLY: finish
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm), POINTER, INTENT(IN) :: this      !< bgc material to test assertion for
    CHARACTER(len=*),           INTENT(IN) :: tile_name !< name of the tile with this bgc material
    CHARACTER(len=*),           INTENT(IN) :: routine   !< name of the calling routine
    ! -------------------------------------------------------------------------------------------------- !
    IF (ASSOCIATED(this%bgcm_list)) THEN
      IF (.NOT. SIZE(this%bgcm_list) .EQ. 0) THEN
        CALL finish(TRIM(routine), 'Found bgcm '//TRIM(this%name)//' in tile '//TRIM(tile_name)// &
          & ' with elements and compartments, but no behaviour is implemented for this so far, please check!')
      END IF
    END IF

  END SUBROUTINE assert_no_compartment_bgc_materials

  ! ======================================================================================================= !
  !>
  !> returns the dimension of the element variables contained in the given bgc material
  !> asserts that it is a dimension which can be handled by the current implementation
  !> and asserts that all elements have the same dimension
  !>
  FUNCTION get_dimension_of_element_variables(this, tile_name, routine) RESULT(this_dim)
    USE mo_exception,           ONLY: finish
    USE mo_jsb_var_class,       ONLY: REAL2D, REAL3D
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm), POINTER, INTENT(IN) :: this      !< bgc material to test assertion for
    CHARACTER(len=*),           INTENT(IN) :: tile_name !< name of the tile with this bgc material
    CHARACTER(len=*),           INTENT(IN) :: routine   !< name of the calling routine
    INTEGER                                :: this_dim  !< dimension
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: nr_of_elements, i_element
    ! -------------------------------------------------------------------------------------------------- !
    nr_of_elements = SIZE(this%element_list)
    this_dim = -1

    ! assert that all variables have the same dimensionality
    DO i_element = 1, nr_of_elements
      IF (this_dim == -1) THEN
        this_dim = this%element_list(i_element)%p%type
      ELSEIF (this_dim /= this%element_list(i_element)%p%type) THEN
        CALL finish(TRIM(routine), 'Found elements with different dimensions for bgcm '//TRIM(this%name)// &
          & ' in tile '//TRIM(tile_name)//', but no behaviour is implemented for this so far, please check!')
      ENDIF
    ENDDO

    IF (this_dim /= REAL2D .AND. this_dim /= REAL3D) THEN
      CALL finish(TRIM(routine), 'Found elements with other than 2D or 3D type for bgcm '//TRIM(this%name)// &
        & ' in tile '//TRIM(tile_name)//' for this so far no behaviour is implemented, please check!')
    ENDIF

  END FUNCTION get_dimension_of_element_variables


  ! ======================================================================================================= !
  !>
  !> get the short name for this element ID
  !>
  FUNCTION get_short_name_for_element_id(id_element) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)           :: id_element
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id_element)
    CASE (ELEM_C_ID)
      return_value = 'C'
    CASE (ELEM_N_ID)
      return_value = 'N'
    CASE (ELEM_P_ID)
      return_value = 'P'
    CASE (ELEM_C13_ID)
      return_value = 'C13'
    CASE (ELEM_C14_ID)
      return_value = 'C14'
    CASE (ELEM_N15_ID)
      return_value = 'N15'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION get_short_name_for_element_id


  ! ======================================================================================================= !
  !>
  !> get the name for this element ID
  !>
  FUNCTION get_name_for_element_id(id_element) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: id_element
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id_element)
    CASE (ELEM_C_ID)
      return_value = 'carbon'
    CASE (ELEM_N_ID)
      return_value = 'nitrogen'
    CASE (ELEM_P_ID)
      return_value = 'phosphorus'
    CASE (ELEM_C13_ID)
      return_value = 'carbon13'
    CASE (ELEM_C14_ID)
      return_value = 'carbon14'
    CASE (ELEM_N15_ID)
      return_value = 'nitrogen15'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION get_name_for_element_id


  ! ======================================================================================================= !
  !>
  !> Returns (and checks!) the element index in the index map for the queried id
  !>
  FUNCTION get_element_index(element_id, elements_index_map, caller) RESULT(element_idx)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,          INTENT(in) :: element_id            !< id of the queried element
    INTEGER,          INTENT(in) :: elements_index_map(:) !< id-index map
    CHARACTER(len=*), INTENT(in) :: caller                !< calling routine
    INTEGER                      :: element_idx           !< index of the queried element
    ! -------------------------------------------------------------------------------------------------- !

    element_idx = elements_index_map(element_id)

    ! Assert that the model is running with this element
    IF (element_idx == 0) THEN
      CALL finish(TRIM(caller), 'Unknown or unused element id: ' // int2string(element_id) &
        & // ' - element name: "' // get_name_for_element_id(element_id) // '". Please check!' )
    ENDIF
  END FUNCTION get_element_index

  ! ======================================================================================================= !
  !>
  !>bgcm add element
  !>
  SUBROUTINE t_lnd_bgcm_add_element(this, name, shortname, id, dim)

    CLASS(t_lnd_bgcm), INTENT(inout) :: this
    CHARACTER(LEN=*),  INTENT(in)    :: name
    CHARACTER(LEN=*),  INTENT(in)    :: shortname
    INTEGER,           INTENT(in)    :: id
    INTEGER,           INTENT(in)    :: dim

    CLASS(t_jsb_var),  POINTER     :: element
    TYPE(t_jsb_var_p), ALLOCATABLE :: element_list(:)
    INTEGER :: n_elements, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_lnd_bgcm_add_element'

    IF (dim==2) THEN
      ALLOCATE(t_jsb_var_real2d :: element)
    ELSE IF (dim==3) THEN
      ALLOCATE(t_jsb_var_real3d :: element)
    END IF
    element%element_name      = name
    element%element_shortname = shortname
    element%element_id        = id

    IF (ALLOCATED(this%element_list)) THEN
      n_elements = SIZE(this%element_list)
    ELSE
      n_elements = 0
    END IF

    ALLOCATE(element_list(n_elements+1))
    DO i=1,n_elements
      element_list(i)%p => this%element_list(i)%p
    END DO
    element_list(n_elements+1)%p => element

    CALL MOVE_ALLOC(element_list, this%element_list)

  END SUBROUTINE t_lnd_bgcm_add_element

  ! ======================================================================================================= !
  !>
  !>bgcm has elements
  !>
  FUNCTION t_lnd_bgcm_has_elements(this) RESULT(has_elements)

    CLASS(t_lnd_bgcm), INTENT(in) :: this
    LOGICAL                       :: has_elements

    has_elements = ALLOCATED(this%element_list)

  END FUNCTION t_lnd_bgcm_has_elements

  ! ======================================================================================================= !
  !>
  !>bgcm get element position
  !>
  FUNCTION t_lnd_bgcm_find_element_position_by_id(this, element_id) RESULT(pos)

    CLASS(t_lnd_bgcm), INTENT(in) :: this
    INTEGER,           INTENT(in) :: element_id
    INTEGER                       :: pos

    INTEGER :: i

    pos = 0

    IF (.NOT. this%Has_elements()) RETURN

    i = 0
    DO
      i = i + 1
      IF (.NOT. ASSOCIATED(this%element_list(i)%p)) THEN   ! End of list
        EXIT
      ELSE IF (this%element_list(i)%p%element_id == element_id) THEN
        pos = i
        EXIT
      END IF
    END DO

  END FUNCTION t_lnd_bgcm_find_element_position_by_id

  ! ======================================================================================================= !
  !>
  !>bgcm get element position by element name
  !>
  FUNCTION t_lnd_bgcm_find_element_position_by_name(this, element_name) RESULT(pos)

    USE mo_util,       ONLY: tolower

    CLASS(t_lnd_bgcm), INTENT(in) :: this
    CHARACTER(LEN=*),  INTENT(in) :: element_name
    INTEGER                       :: pos

    INTEGER :: i

    pos = 0

    IF (.NOT. ALLOCATED(this%element_list)) RETURN

    i = 0
    DO
      i = i + 1
      IF (tolower(this%element_list(i)%p%element_name) == tolower(element_name)) THEN
        pos = i
        EXIT
      END IF
    END DO

  END FUNCTION t_lnd_bgcm_find_element_position_by_name

  ! ======================================================================================================= !
  !>
  !>bgcm add variable
  !>
  SUBROUTINE t_lnd_bgcm_add_var(this, var, dim)

    CLASS(t_lnd_bgcm),        INTENT(inout) :: this
    CLASS(t_jsb_var), TARGET, INTENT(in)    :: var
    INTEGER,                  INTENT(in)    :: dim

    TYPE(t_jsb_var_p), ALLOCATABLE :: var_list(:)
    INTEGER :: n_vars, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_lnd_bgcm_add_var'

    IF (ALLOCATED(this%var_list)) THEN
      n_vars = SIZE(this%var_list)
    ELSE
      n_vars = 0
    END IF

    ALLOCATE(var_list(n_vars+1))
    DO i=1,n_vars
      var_list(i)%p => this%var_list(i)%p
    END DO
    var_list(n_vars+1)%p => var

    CALL MOVE_ALLOC(var_list, this%var_list)

  END SUBROUTINE t_lnd_bgcm_add_var

  ! ======================================================================================================= !
  !>
  !>bgcm add bgcm
  !>
  SUBROUTINE t_lnd_bgcm_add_bgcm(this, bgcm, id)

    CLASS(t_lnd_bgcm),         INTENT(inout) :: this
    CLASS(t_lnd_bgcm), TARGET, INTENT(inout) :: bgcm
    INTEGER, OPTIONAL,         INTENT(in)    :: id  !< unique id of this bgcm

    TYPE(t_lnd_bgcm_p), POINTER :: bgcm_list(:), leaf_list(:)
    INTEGER :: n_bgcm_main, n_leaves, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_lnd_bgcm_add_bgcm'

    ! if bgc materials should be used within the bgcm store, they need to have a unique id for quick association
    IF (PRESENT(id)) bgcm%id = id

    IF (ASSOCIATED(this%bgcm_list)) THEN
      n_bgcm_main = SIZE(this%bgcm_list)
    ELSE
      n_bgcm_main = 0
    END IF

    ALLOCATE(bgcm_list(n_bgcm_main + 1))
    DO i = 1, n_bgcm_main
      bgcm_list(i)%p => this%bgcm_list(i)%p
    END DO
    bgcm_list(n_bgcm_main + 1)%p => bgcm

    IF (ASSOCIATED(this%bgcm_list)) DEALLOCATE(this%bgcm_list)
    this%bgcm_list => bgcm_list

    IF (ASSOCIATED(this%leaf_list)) THEN
      n_leaves = SIZE(this%leaf_list)
    ELSE
      n_leaves = 0
    END IF

    IF (ASSOCIATED(bgcm%leaf_list)) THEN
      ALLOCATE(leaf_list(n_leaves + SIZE(bgcm%leaf_list)))
    ELSE
      ALLOCATE(leaf_list(n_leaves + 1))
    END IF

    DO i = 1, n_leaves
      leaf_list(i)%p => this%leaf_list(i)%p
    END DO

    IF (ASSOCIATED(bgcm%leaf_list)) THEN
      DO i = 1, SIZE(bgcm%leaf_list)
        leaf_list(n_leaves+i)%p => bgcm%leaf_list(i)%p
      END DO
    ELSE
      leaf_list(n_leaves + 1)%p => bgcm
    END IF

    IF (ASSOCIATED(this%leaf_list)) DEALLOCATE(this%leaf_list)
    this%leaf_list => leaf_list

  END SUBROUTINE t_lnd_bgcm_add_bgcm

  ! ======================================================================================================= !
  !>
  !>bgcm set paths (memory management)
  !>
  RECURSIVE SUBROUTINE t_lnd_bgcm_set_paths(this, path)

    CLASS(t_lnd_bgcm), INTENT(inout) :: this
    CHARACTER(LEN=*),     OPTIONAL   :: path

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_lnd_bgcm_set_paths'

    IF (PRESENT(path)) THEN
      IF (path /= '') THEN
        this%path = path // ':' // this%shortname
      ELSE
        this%path = this%shortname
      END IF
    ELSE
      this%path = this%shortname
    END IF
    ! print*,'DDD ', this%shortname, path, this%path

    IF (ASSOCIATED(this%bgcm_list)) THEN
      DO i=1,SIZE(this%bgcm_list)
        CALL this%bgcm_list(i)%p%Set_paths(this%path)
      END DO
    END IF

  END SUBROUTINE t_lnd_bgcm_set_paths

#endif
END MODULE mo_lnd_bgcm_class
