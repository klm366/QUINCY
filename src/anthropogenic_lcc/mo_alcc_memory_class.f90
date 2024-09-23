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

   ! TO-DO, associate the Transmatrix with the girds, so iblk (<= nblk) and ic (<= nc) in the jsb_interface
  ! comments from Julia N: The update_alcc routine is already only called for one chunk, 
  ! thus you do not have to care about the block number (iblk) but only about the cells on this chunk (ic <= nc). 
  ! So You would loop over nc which is already available in the update_all routine (nc  = options%nc).
  ! Mingyue Zhang: define TransMtrx start
  TYPE MatrixPtr
    real(dp), pointer :: ptr(:,:)
  end TYPE MatrixPtr
  
  real(dp), pointer :: TransMtrx(:)
  
  ! Mingyue Zhang: Set the number of matrices dynamically
  INTEGER, PARAMETER :: numMatrices = 9 ! Adjust as needed according to the size of transition matrix
  
  allocate(TransMtrx(numMatrices))
  
  do i = 1, numMatrices   ! Adjust 'rows' and 'cols' as needed, should be same as the simulation domian lat/lon box
    allocate(TransMtrx(i)%ptr(rows, cols))
  end do

  ! Accessing elements of the matrices 
  ! Use a do loop for the read in LUH2 transition data: 
  ! crop2pasture pasture2crop natural2crop natural2pasture crop2natural pasture2natural
  TransMtrx(1)%ptr(i, j) = 1.0_dp ! Example assignment (now go to init.f90 to read in data and then come back to here)
  ! Mingyue Zhang: define TransMtrx end 

  ! Mingyue Zhang: defeine transition types start  CHARACTER(len=*) :: transitions(9)

    ! Initialize the array with the character strings
    transitions = (/ 'crop2crop','crop2pasture', 'crop2natural',&
                     'pasture2crop','pasture2pasture', 'pasture2natural',&
                     'natural2crop', 'natural2pasture', 'natural2natural'/)
    ! Output to verify the array
    print *, transitions
  ! Mingyue Zhang: defeine transition types end


  TYPE, EXTENDS(t_jsb_memory) :: t_alcc_memory

    TYPE(t_jsb_var_real3d) :: &
      & cf_current_year,      & !< This years cover fractions derived from annually read land use data [-].
      & cf_day_delta            !< Daily change in cover fractions calculated from cf_current_year.

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

    CALL mem%Add_var('cf_current_year', mem%cf_current_year,                              &
      & hgrid, pft_vgrid,                                                                 &
      & t_cf('This years cover fractions read from forcing.'),  &
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

    ! Mingyue Zhang: add variables for transitions matrix start 
    CALL mem%Add_var('TransMtrx', mem%TransMtrx,                 &
    & hgrid, pft_vgrid,                                                                            &
    & prefix, suffix,                                                                              &
    & output_level=FULL,                                                                          &
    & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
    ! Mingyue Zhang: add variables for transitions matrix finish
    ! Mingyue Zhang: add variables for transitions type start 
    CALL mem%Add_var('transitions', mem%transitions,                 &
    & hgrid, pft_vgrid,                                                                            &
    & prefix, suffix,                                                                              &
    & output_level=FULL,                                                                          &
    & loutput=.FALSE., lrestart=.FALSE., initval_r=-0.1_wp)
    ! Mingyue Zhang: add variables for transitions type finish


  END SUBROUTINE Init_alcc_memory

#endif
END MODULE mo_alcc_memory_class
