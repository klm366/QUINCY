!> Contains class definitions for pools
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
MODULE mo_jsb_pool_class
#ifndef __NO_JSBACH__

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message
  USE mo_util,          ONLY: int2string
  USE mo_jsb_var_class, ONLY: t_jsb_var, t_jsb_var_p, t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_lnd_bgcm, t_lnd_bgcm_p
  PUBLIC :: ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15  ! element_id (element enumerator) for pools

  TYPE, ABSTRACT :: t_lnd_bgcm
    CHARACTER(LEN=:),  ALLOCATABLE :: shortname                                ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: name                                     ! No spaces!
    CHARACTER(LEN=:),  ALLOCATABLE :: path                                     !
    CHARACTER(LEN=:),  ALLOCATABLE :: element_unit                             ! unit of pool element variables (used for output stream)
    INTEGER                        :: id                            = -1
      !< a pool id that - if set /= -1 - should be unique (from mo_lnd_bgcm_store_class)
    INTEGER                        :: owner_model_id                = 0        ! ID of the jsbach model this pool is running at
    LOGICAL                        :: contains_elements             = .FALSE.  ! contains element variables
    LOGICAL                        :: contains_bgcm                 = .FALSE.  ! contains bgc_material (t_lnd_bgcm)
    LOGICAL                        :: contains_elem_1d_ptr_to_chunk = .FALSE.  ! t_lnd_bgcm_elem2d -> set TRUE in the pool init routines in the PROC_memory_class
    LOGICAL                        :: contains_elem_2d_ptr_to_chunk = .FALSE.  ! t_lnd_bgcm_elem3d -> set TRUE in the pool init routines in the PROC_memory_class
    CLASS(t_lnd_bgcm_p),   POINTER :: pool_list(:) => NULL()                   ! children
    CLASS(t_lnd_bgcm_p),   POINTER :: leaf_list(:) => NULL()                   ! leaves
    TYPE(t_jsb_var_p), ALLOCATABLE :: element_list(:)                          ! elements of this pool (elements used in PROC are defined in 'all_elements(:)' in PROC_memory_class)
    TYPE(t_jsb_var_p), ALLOCATABLE :: var_list(:)                              ! other associated variables - currently not used
  CONTAINS
    PROCEDURE(init_iface), DEFERRED :: Init
    PROCEDURE(get_el_name_iface_id), DEFERRED :: Get_element_name_by_id
    GENERIC                         :: Get_element_name => Get_element_name_by_id
    PROCEDURE                       :: Add_element           => t_lnd_bgcm_add_element
    PROCEDURE                       :: Add_var               => t_lnd_bgcm_add_var
    PROCEDURE                       :: Has_elements          => t_lnd_bgcm_has_elements
    PROCEDURE                       :: Add_pool              => t_lnd_bgcm_add_pool
    PROCEDURE                       ::                          t_lnd_bgcm_find_element_position_by_id
    PROCEDURE                       ::                          t_lnd_bgcm_find_element_position_by_name
    PROCEDURE                       :: Set_paths             => t_lnd_bgcm_set_paths
    GENERIC                         :: Find_element_position => t_lnd_bgcm_find_element_position_by_id,   &
      &                                                         t_lnd_bgcm_find_element_position_by_name
  END TYPE t_lnd_bgcm

  TYPE :: t_lnd_bgcm_p
    CLASS(t_lnd_bgcm), POINTER :: p
  END TYPE t_lnd_bgcm_p

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

  !>element_id (element enumerator) for bgc_material
  !>
  !! used in the PROC_memory_class for initialization of pools
  !!
  !! names of element variables used in JSBACH are:
  !! (sequence accordingly to the below enumerator)
  !!    carbon, nitrogen, phosphorus, carbon13, carbon14, nitrogen15
  !!
  !! but see element ID ENUM in mo_lnd_bgcm_class
  !!
  ENUM, BIND(C)
    ENUMERATOR :: ELEM_C = 1, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_pool_class'

CONTAINS

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

  FUNCTION t_lnd_bgcm_has_elements(this) RESULT(has_elements)

    CLASS(t_lnd_bgcm), INTENT(in) :: this
    LOGICAL                       :: has_elements

    has_elements = ALLOCATED(this%element_list)

  END FUNCTION t_lnd_bgcm_has_elements

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

  SUBROUTINE t_lnd_bgcm_add_pool(this, pool, id)

    CLASS(t_lnd_bgcm),         INTENT(inout) :: this
    CLASS(t_lnd_bgcm), TARGET, INTENT(inout) :: pool
    INTEGER, OPTIONAL,         INTENT(in)    :: id  !< unique id of this bgcm

    TYPE(t_lnd_bgcm_p), POINTER :: pool_list(:), leaf_list(:)
    INTEGER :: n_pools, n_leaves, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_lnd_bgcm_add_pool'

    ! if bgc materials should be used within the bgcm store, they need to have a unique id for quick association
    IF( PRESENT(id)) pool%id = id

    IF (ASSOCIATED(this%pool_list)) THEN
      n_pools = SIZE(this%pool_list)
    ELSE
      n_pools = 0
    END IF

    ALLOCATE(pool_list(n_pools+1))
    DO i=1,n_pools
      pool_list(i)%p => this%pool_list(i)%p
    END DO
    pool_list(n_pools+1)%p => pool

    IF (ASSOCIATED(this%pool_list)) DEALLOCATE(this%pool_list)
    this%pool_list => pool_list

    IF (ASSOCIATED(this%leaf_list)) THEN
      n_leaves = SIZE(this%leaf_list)
    ELSE
      n_leaves = 0
    END IF

    IF (ASSOCIATED(pool%leaf_list)) THEN
      ALLOCATE(leaf_list(n_leaves + SIZE(pool%leaf_list)))
    ELSE
      ALLOCATE(leaf_list(n_leaves + 1))
    END IF

    DO i=1,n_leaves
      leaf_list(i)%p => this%leaf_list(i)%p
    END DO

    IF (ASSOCIATED(pool%leaf_list)) THEN
      DO i=1,SIZE(pool%leaf_list)
        leaf_list(n_leaves+i)%p => pool%leaf_list(i)%p
      END DO
    ELSE
      leaf_list(n_leaves+1)%p => pool
    END IF

    IF (ASSOCIATED(this%leaf_list)) DEALLOCATE(this%leaf_list)
    this%leaf_list => leaf_list

  END SUBROUTINE t_lnd_bgcm_add_pool

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

    IF (ASSOCIATED(this%pool_list)) THEN
      DO i=1,SIZE(this%pool_list)
        CALL this%pool_list(i)%p%Set_paths(this%path)
      END DO
    END IF

  END SUBROUTINE t_lnd_bgcm_set_paths

#endif
END MODULE mo_jsb_pool_class
