!> alcc (anthropogenic lcc) memory class
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains memory definitions for alcc
!>
MODULE mo_alcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: dp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d
  USE mo_jsb_varlist,            ONLY: FULL

  ! Use of processes in this module
  ! dsl4jsb_Use_processes ALCC_

  ! Use process configurations
  ! dsl4jsb_Use_config(ALCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_alcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  ! TO-DO, associate the Transmatrix with the girds, so iblk (<= nblk) and ic (<= nc) in the jsb_interface
  ! comments from Julia N: The update_alcc routine is already only called for one chunk, 
  ! thus you do not have to care about the block number (iblk) but only about the cells on this chunk (ic <= nc). 
  ! So You would loop over nc which is already available in the update_all routine (nc  = options%nc).
  
  TYPE, EXTENDS(t_jsb_memory) :: t_alcc_memory
  

    TYPE(t_jsb_var_real3d) :: &
      & cf_current_year,       & !< This years cover fractions derived from annually read land use data [-].
      & cf_day_delta,          & !< Daily change in cover fractions calculated from cf_current_year.
      & cf_transition,         & ! Mingyue Zhang:define transition gain v_grid = 3*3
      & cf_transition_pfts,    & ! Mingyue Zhang:define transition loss v_grid = npft*npft
      & cf_transition_gain,    & ! Mingyue Zhang:define transition matrix v_grid = 11
      & cf_transition_loss       ! Mingyue Zhang:define transition matrix v_grid = 11

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
    CLASS(t_alcc_memory), INTENT(inout), TARGET   :: mem

    CHARACTER(len=*),      INTENT(in)             :: prefix
    CHARACTER(len=*),      INTENT(in)             :: suffix
    INTEGER,               INTENT(in)             :: lct_ids(:)
    INTEGER,               INTENT(in)             :: model_id
    ! -------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Def_config(ALCC_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                       ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: pft_vgrid                   ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER  :: trans_vgrid                 ! Mingyue: Vertical grid for C,N,P transitions
    TYPE(t_jsb_vgrid), POINTER  :: trans_pft_vgrid             ! Mingyue: Vertical grid for 11 pfts transitions
    INTEGER                     :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_memory'
    ! -------------------------------------------------------------------------------------------------- !

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid   => Get_grid(mem%grid_id)
    pft_vgrid    => Get_vgrid('pfts')
    trans_vgrid    => Get_vgrid('transition')
    trans_pft_vgrid    => Get_vgrid('transition_pfts')

!    dsl4jsb_Get_config(ALCC_)

    CALL mem%Add_var('cf_current_year', mem%cf_current_year,                              &
      & hgrid, pft_vgrid,                                                                 &
      & t_cf('This years cover fractions read from forcing.'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.FALSE., initval_r=-0.1_wp)

    CALL mem%Add_var('cf_day_delta', mem%cf_day_delta,                                    &
      & hgrid, pft_vgrid,                                                                 &
      & t_cf('Daily change in cover fractions.'),                                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=FULL,                                                                &
      & lrestart=.TRUE., lrestart_cont=.TRUE., initval_r=0.0_wp)

     ! Mingyue Zhang: add variables for transitions type start 
    CALL mem%Add_var('cf_transition', mem%cf_transition,                                   &
      & hgrid, trans_vgrid,                                                                &
      & t_cf('transitions between crop, natural and pasture'),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                                    &
      & output_level=FULL,                                                                 &
      & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
     ! Mingyue Zhang: add variables for transitions type finish

     ! Mingyue Zhang: add variables for transitions type start 
    CALL mem%Add_var('cf_transition_pfts', mem%cf_transition_pfts,                          &
      & hgrid, trans_pft_vgrid,                                                             &
      & t_cf('transitions between pfts'),                                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
      & prefix, suffix,                                                                     &
      & output_level=FULL,                                                                  &
      & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
     ! Mingyue Zhang: add variables for transitions type finish
 
     ! Mingyue Zhang: add variables for transitions gain start 
     ! TODO: the vgrid(trans_vgrid) should be modified
    CALL mem%Add_var('cf_transition_gain', mem%cf_transition_gain,                         &
      & hgrid, pft_vgrid,                                                                  &
      & t_cf('transitions gain_pfts'),                                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                                    &
      & output_level=FULL,                                                                 &
      & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
     ! Mingyue Zhang: add variables for transitions gain finish

     ! Mingyue Zhang: add variables for transitions loss start 
     ! TODO: the vgrid(trans_vgrid) should be modified
    CALL mem%Add_var('cf_transition_loss', mem%cf_transition_loss,                          &
      & hgrid, pft_vgrid,                                                                   &
      & t_cf('transitions loss_pfts'),                                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
      & prefix, suffix,                                                                     &
      & output_level=FULL,                                                                  &
      & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
     ! Mingyue Zhang: add variables for transitions loss finish

  END SUBROUTINE Init_alcc_memory

#endif
END MODULE mo_alcc_memory_class
