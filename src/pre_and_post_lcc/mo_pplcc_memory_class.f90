!> Summary: pplcc (pre- and post lcc) memory class (currently of no use beyond being required by the infrastructure)
!>
!>### Could contain memory definitions for pplcc, but currently of no use beyond infrastructural requirements
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-02-02)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_pplcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pplcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_pplcc_memory
    
  CONTAINS
    PROCEDURE :: Init => Init_pplcc_memory
  END TYPE t_pplcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_memory_class'

CONTAINS

  SUBROUTINE Init_pplcc_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_pplcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pplcc_memory'
    ! -------------------------------------------------------------------------------------------------- !
    
    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    !This process currently has no own variables 

  END SUBROUTINE Init_pplcc_memory

#endif
END MODULE mo_pplcc_memory_class
