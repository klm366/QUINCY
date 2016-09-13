!>
!! @brief 
!!  In SLM: contains a small version of mo_jsb_varlist
!!  In JSBACH4: Contains routines for var_lists/streams
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20160816)
!!
MODULE mo_slm_varlist

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish

  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                  t_jsb_var_real1d, &
                                  t_jsb_struct, t_jsb_n_struct
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: slm_add_var

  INTERFACE slm_add_var
    MODULE PROCEDURE slm_add_struct               ! *slmjsb4  a TYPE with only variables included
    MODULE PROCEDURE slm_add_n_struct             ! *slmjsb4  a TYPE with nested TYPEs
    MODULE PROCEDURE slm_add_l2_struct            ! *slmjsb4  a TYPE with only variables included within a TYPE
    MODULE PROCEDURE slm_add_struct_element_r1d   ! *slmjsb4 *slmtemp | works with struct & l2_struct !!
    MODULE PROCEDURE slm_add_var_list_element_r1d
    MODULE PROCEDURE slm_add_var_list_element_r2d
    MODULE PROCEDURE slm_add_var_list_element_r3d
!    MODULE PROCEDURE slm_add_var_list_element_r2d_tile
  END INTERFACE slm_add_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_slm_varlist'

CONTAINS

! *slmtemp -> 1D to 2D
! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> add 1D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real1d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r1d( &
    struct, name, var, &
    ics, ice, iblk, &
    initval_r) 

    ! 0.1 InOut
    CLASS(t_jsb_struct),  INTENT(inout)           :: struct              ! struct & l2_struct
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    TYPE(t_jsb_var_real1d), TARGET, INTENT(inout) :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    INTEGER,             INTENT(in)               :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: iblk  ! model%options -> for SLM compatibility
    REAL(wp),             INTENT(in)              :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER           :: pos_of_var

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_real'

    var%name = TRIM(name)
    !var%owner_proc_id = struct%owner_proc_id

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr)
    var%ptr = initval_r

  END SUBROUTINE slm_add_struct_element_r1d


! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> add structure to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_struct:struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct( &
    mem, name, struct) 

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)    :: struct              ! struct
    ! 0.2 Local
    INTEGER           :: pos_of_struct

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct'
    struct%name = TRIM(name)
    pos_of_struct = mem%no_of_structs + 1
    IF (pos_of_struct > mem%max_no_of_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of structs exceeded: '//TRIM(mem%structlist_name)//', '//TRIM(name))
    ELSE
      mem%structs(pos_of_struct)%p => struct
      mem%no_of_structs = pos_of_struct
    END IF

    ALLOCATE(struct%vars(struct%max_no_of_vars))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_struct  ! *slmtemp

  END SUBROUTINE slm_add_struct


! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> add nested structure to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_n_struct:n_struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_n_struct( &
    mem, name, n_struct) 

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_n_struct), TARGET, INTENT(inout)  :: n_struct            ! n_struct
    ! 0.2 Local
    INTEGER           :: pos_of_n_struct

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_n_struct'
    n_struct%name = TRIM(name)
    pos_of_n_struct = mem%no_of_n_structs + 1
    IF (pos_of_n_struct > mem%max_no_of_n_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of n_structs exceeded: '//TRIM(mem%n_structlist_name)//', '//TRIM(name))
    ELSE
      mem%n_structs(pos_of_n_struct)%p => n_struct
      mem%no_of_n_structs = pos_of_n_struct
    END IF

    ALLOCATE(n_struct%l2_structs(n_struct%max_no_of_l2_structs))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_n_struct  ! *slmtemp

  END SUBROUTINE slm_add_n_struct


! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> add structure to a nested structure that is located directly in the process memory
  !!
  !!  t_jsb_n_struct:n_struct <- t_jsb_struct:l2_struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_l2_struct( &
    n_struct, name, l2_struct) 

    ! 0.1 InOut
    CLASS(t_jsb_n_struct),  INTENT(inout)       :: n_struct            ! memory structure
    CHARACTER(len=*),     INTENT(in)            :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)  :: l2_struct           ! l2_struct
    ! 0.2 Local
    INTEGER           :: pos_of_l2_struct

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_l2_struct'
    l2_struct%name = TRIM(name)
    pos_of_l2_struct = n_struct%no_of_l2_structs + 1
    IF (pos_of_l2_struct > n_struct%max_no_of_l2_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of l2_structs exceeded: '//TRIM(n_struct%l2_structlist_name)//', '//TRIM(name))
    ELSE
      n_struct%l2_structs(pos_of_l2_struct)%p => l2_struct
      n_struct%no_of_l2_structs = pos_of_l2_struct
    END IF

    ALLOCATE(l2_struct%vars(l2_struct%max_no_of_vars))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_struct  ! *slmtemp

  END SUBROUTINE slm_add_l2_struct


! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> add 1D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_real1d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r1d( &
    mem, name, var, &
    ics, ice, iblk, &
    initval_r) 

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    TYPE(t_jsb_var_real1d), TARGET, INTENT(inout) :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    INTEGER,             INTENT(in)               :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: iblk  ! model%options -> for SLM compatibility
    REAL(wp),             INTENT(in)              :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER           :: pos_of_var

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_real'

    var%name = TRIM(name)
    !var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)
    ALLOCATE(var%ptr)
    var%ptr = initval_r

  END SUBROUTINE slm_add_var_list_element_r1d

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> add 2D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_real2d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r2d( &
    mem, name, var, &
    ics, ice, iblk, &
!    hgrid, vgrid, &
!    cf, &
!    grib1, grib2, &
!    prefix, suffix, &
!    !loutput, lcontainer, &
!    !lrestart, lrestart_cont, initval_r, isteptype,              &
!    lrestart, initval_r
!    !resetval_r, lmiss, missval_r, tlev_source, info,            &
!    !verbose, new_element, cdiTimeID                             &
    initval_r) 

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET, INTENT(inout) :: var                 ! reference to jsbach var
    INTEGER,             INTENT(in)               :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: iblk  ! model%options -> for SLM compatibility
    REAL(wp),             INTENT(in)              :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER           :: pos_of_var

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r2d'

    var%name = TRIM(name)
    !var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)
    ALLOCATE(var%ptr(ics:ice,iblk))
    var%ptr(ics:ice,iblk) = initval_r


  END SUBROUTINE slm_add_var_list_element_r2d


! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> add 3D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_real3d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r3d( &
    mem, name, var, &
    ics, ice, iblk, &
!    hgrid, vgrid, &
!    cf, &
!    grib1, grib2, &
!    prefix, suffix, &
!    !loutput, lcontainer, &
!    !lrestart, lrestart_cont, initval_r, isteptype,              &
!    lrestart, initval_r
!    !resetval_r, lmiss, missval_r, tlev_source, info,            &
!    !verbose, new_element, cdiTimeID                             &
    initval_r) 

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    TYPE(t_jsb_var_real3d), TARGET, INTENT(inout) :: var                 ! reference to jsbach var
    INTEGER,             INTENT(in)               :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)               :: iblk  ! model%options -> for SLM compatibility
    REAL(wp),             INTENT(in)              :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER           :: pos_of_var

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r3d'

    var%name = TRIM(name)
    !var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)
    ALLOCATE(var%ptr(ics:ice,1:1,iblk))     ! *slmtemp added 1:1 to make allocation work
    var%ptr(ics:ice,1:1,iblk) = initval_r   ! *slmtemp added 1:1 to make allocation work


  END SUBROUTINE slm_add_var_list_element_r3d



END MODULE mo_slm_varlist
