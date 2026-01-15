!> Summary: alcc (anthropogenic lcc) memory class
!>
!>### Contains memory definitions for alcc
!>
!>
!> Copyright (c) 2018 Rainer Schneck
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> - First version   by Rainer Schneck (2018-12-03)
!> - condensed       by Julia Nabel    (2021-06-23)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_alcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: FULL

  ! Use of processes in this module
  ! dsl4jsb_Use_processes ALCC_

  ! Use process configurations
  ! dsl4jsb_Use_config(ALCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_alcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  TYPE, EXTENDS(t_jsb_memory) :: t_alcc_memory

    TYPE(t_jsb_var_real3d) :: cf_current_year
    !< this years target cfs derived from annually read land use data file. [-]

  CONTAINS
    PROCEDURE :: Init => Init_alcc_memory
  END TYPE t_alcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_memory_class'

CONTAINS

  SUBROUTINE Init_alcc_memory(mem, prefix, suffix, lct_ids, model_id)

    ! Declarations
    USE mo_jsb_io,             ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    !USE mo_jsb_model_usecases, ONLY: npft ! R: Note, this leads to circular dependency
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_alcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)             :: prefix
    CHARACTER(len=*),      INTENT(in)             :: suffix
    INTEGER,               INTENT(in)             :: lct_ids(:)
    INTEGER,               INTENT(in)             :: model_id
    ! -------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Def_config(ALCC_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid               ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: pft_vgrid           ! Vertical grid
    INTEGER                     :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_memory'
    ! -------------------------------------------------------------------------------------------------- !

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid   => Get_grid(mem%grid_id)
    pft_vgrid    => Get_vgrid('pfts')

!    dsl4jsb_Get_config(ALCC_)

    !JN-TODO: if changes are to be done on a daily basis maybe other field more interesting?
    !           e.g. target_cf ? or better remember something else e.g. "required_cf_change"?
    !         then also change out_level to BASIC and lrestart to TRUE
    CALL mem%Add_var('cf_current_year', mem%cf_current_year,                              &
      & hgrid, pft_vgrid,                                                                 &
      & t_cf('This years target cfs derived from an annually read land use data file.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=FULL,                                                                &
      & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)

  END SUBROUTINE Init_alcc_memory

#endif
END MODULE mo_alcc_memory_class