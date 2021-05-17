!> Summary: Types for lcc processes
!>
!>### Contains types required for dealing with actively or passively relocated conserved quantities upon lcc
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
MODULE mo_jsb_lcc_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
    
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_lcc_proc, t_jsb_lcc_var

  !> Type used to enclose pointers and values for each conserved quantity with active or passive relocation 
  TYPE :: t_jsb_lcc_var
    CHARACTER(len=VARNAME_LEN)    :: name               !< Name of the variable
!    INTEGER :: type_id = -1     !< one of the CQ_TYPEs ( mo_jsb_cqt_class ) -- JN-TODO: required -> currently in encompassing type...
    TYPE(t_jsb_var_p), ALLOCATABLE :: tile_pointer(:)   !< One pointer per tile involved in the encompassing lcc process
    REAL(wp), ALLOCATABLE :: relocate_this(:,:,:)       !< One array (on domain) per tile
  END TYPE t_jsb_lcc_var
  TYPE :: t_jsb_lcc_var_p
    TYPE(t_jsb_lcc_var), POINTER :: p => NULL()
  END TYPE t_jsb_lcc_var_p

  !> Type used on lcc processes to collect conserved quantities with active and passive relocations
  TYPE :: t_jsb_lcc_proc
    CHARACTER(len=SHORT_NAME_LEN) :: name !< Name of the structure: process + _lcc 
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: tile_names(:) !< Names of the tiles involved with this lcc process
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: active_vars_names(:) ! JN-TODO maybe not required?
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: passive_vars_names(:) ! JN-TODO maybe not required?
    INTEGER, ALLOCATABLE :: active_vars_cqt(:)  ! JN-TODO: required? 
    INTEGER, ALLOCATABLE :: passive_vars_cqt(:) ! JN-TODO: required? 
    INTEGER :: nr_of_tiles = 0          !< Number of tiles involved with this lcc process
    INTEGER :: nr_of_active_vars = 0    !< Number of variables with active relocation in the encompassing lcc process
    INTEGER :: nr_of_passive_vars = 0   !< Number of variables with passive relocation in the encompassing lcc process
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: active_vars(:)
        !< Collection of lcc var types that are actively relocated for the encompassing lcc process
    TYPE(t_jsb_lcc_var_p), ALLOCATABLE :: passive_vars(:) 
        !< Collection of lcc var types that are passively relocated for the encompassing lcc process
  END TYPE t_jsb_lcc_proc

CONTAINS

#endif
END MODULE mo_jsb_lcc_class
