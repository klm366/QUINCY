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

  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: slm_add_var

  INTERFACE slm_add_var
    MODULE PROCEDURE slm_add_var_list_element_r1d
    MODULE PROCEDURE slm_add_var_list_element_r2d
!    MODULE PROCEDURE slm_add_var_list_element_r3d
!    MODULE PROCEDURE slm_add_var_list_element_r2d_tile
  END INTERFACE slm_add_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_slm_varlist'

CONTAINS

  SUBROUTINE slm_add_var_list_element_r1d( &
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


END MODULE mo_slm_varlist
