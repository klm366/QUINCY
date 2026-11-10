!> Contains basic definitions for variable class
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
!> Variables of ICON-Land are objects of the types defined here, one for each rank (1d, 2d, 3d). They are part
!> of the process memories (see `mo_jsb_memory_class`) and are created with `mem%Add_var(...)` in the memory
!> initialization of each process, which also adds them to the var_list of the host model and optionally
!> to output and restart files (1d variables are only added if requested for output or restart).
!>
MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: lp, dp
  USE mo_exception,         ONLY: finish, message
  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p, REAL1D, REAL2D, REAL3D

  !> Type (rank) of a variable, see [[t_jsb_var:type]]
  ENUM, BIND(C)
    ENUMERATOR :: REAL1D=1, REAL2D, REAL3D
  END ENUM

  !================================================================================================================================
  !>
  !> Abstract base type of ICON-Land variables
  !>
  !> Holds the meta data of a variable and pointers to its data. The data arrays themselves are allocated in
  !> the var_list of the host model; depending on the rank, `ptr1d`, `ptr2d` or `ptr3d` points to them, the
  !> other two are not associated, respectively. The extended types add a rank-specific pointer `ptr` to the
  !> same data. Most components are set in the `Add_var` routines of `mo_jsb_memory_class`.
  !>
  TYPE, ABSTRACT :: t_jsb_var
    REAL(lp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr1d(:) => NULL()          !< Data of 1d variables (one global value)
    REAL(lp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr2d(:,:) => NULL()        !< Data of 2d variables (nproma, nblks)
    REAL(lp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr3d(:,:,:) => NULL()      !< Data of 3d variables (used with: nproma, nlev, nblks)
    CHARACTER(len=VARNAME_LEN)    :: name               !< Name of variable (without prefix and suffix)
    CHARACTER(len=VARNAME_LEN)    :: standard_name      !< Standard name of the variable (CF)
    CHARACTER(len=VARNAME_LEN)    :: full_name          !< Full name of variable in varlist, including prefix
                                                        !< and suffix (e.g. `<process>_<name>_<tile>`), as in
                                                        !< output files
    CHARACTER(len=:), ALLOCATABLE :: unit               !< Unit of the variable (CF)
    INTEGER                       :: type               !< REAL1D / REAL2D / REAL3D
    INTEGER                       :: subset_type        !< Currently not used (never set)
    CHARACTER(len=:), ALLOCATABLE :: element_name       !< Only for elements of bgc materials (QUINCY):
                                                        !< name of the chemical element, e.g. 'carbon'
    CHARACTER(len=:), ALLOCATABLE :: element_shortname  !< Only for elements of bgc materials: short name of
                                                        !< the element, e.g. 'C' (part of the variable name)
    INTEGER                       :: element_id         !< Only for element variables of bgc materials: ID of
                                                        !< the element (`ELEM_*_ID`, see `mo_lnd_bgcm_class`)
    INTEGER                       :: owner_model_id     !< ID of the model owning the variable
    INTEGER                       :: owner_proc_id = -1 !< ID of the process owning the variable
    INTEGER, ALLOCATABLE          :: owner_tile_path(:) !< Path of the owning tile in the tile tree
    INTEGER, ALLOCATABLE          :: child_idx(:)       !< Position-numbers of the same var (as this var) on
                                                        !< the children of the current tile (i.e. within
                                                        !< mem%vars(:) of that child)
    REAL(dp)                      :: missval            !< Missing value
    LOGICAL :: l_aggregate_all = .FALSE.                !< When aggregating, consider fractions of all child
                                                        !< tiles even of those where variable doesn't exist
                                                        !< (e.g. for fluxes per grid box); see
                                                        !< `Aggregate_weighted_by_fract_2d/3d`
    LOGICAL :: is_conserved_quan = .FALSE.              !< Whether variable carries a to-be-conserved quantity
    INTEGER :: cons_quan_type_id = -1                   !< ID identifying the type of the conserved quantity
    LOGICAL :: transfer_all = .FALSE.                   !< Flag indicating if all matter should be collected
                                                        !< (true) or if matter should only be collected
                                                        !< proportionally to the lost area
    LOGICAL :: is_in_output  = .FALSE.                  !< Whether variable is requested to be written to
                                                        !< output (by name or group in the output namelists)
    LOGICAL :: is_in_restart = .FALSE.                  !< Whether variable is written to restart file
                                                        !< (`lrestart` of `Add_var`, default: true)
  END TYPE t_jsb_var

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"

  !================================================================================================================================
  !>
  !> 1d variable (e.g. global values for monitoring)
  !>
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real1d
    REAL(lp), POINTER &
      :: ptr(:) => NULL()          !< Associated with ptr1d(:)
  END TYPE t_jsb_var_real1d

  !================================================================================================================================
  !>
  !> 2d variable (horizontal grid)
  !>
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(lp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
      :: ptr(:,:) => NULL()    !< Associated with ptr2d(:,:)
  END TYPE t_jsb_var_real2d

  !================================================================================================================================
  !>
  !> 3d variable (horizontal grid and vertical)
  !>
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(lp), POINTER :: ptr(:,:,:) => NULL()     !< Associated with ptr3d(:,:,:)
  END TYPE t_jsb_var_real3d

!#pragma GCC diagnostic pop

  !================================================================================================================================
  !>
  !> Pointer to a variable (e.g. for the list of variables `vars` of a memory)
  !>
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p => NULL()    !< Pointer to the variable
  END TYPE t_jsb_var_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'  !< Module name

#endif
END MODULE mo_jsb_var_class
