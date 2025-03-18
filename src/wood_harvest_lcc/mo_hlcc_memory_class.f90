! Notes: 
! 27.1. changed vgrid to forest age classes 
!
! TO-DO:
! when this is "generalised" for non-AC cases, need to have a variable *vgrid*
!
!
!> hlcc (harvesting lcc) memory class
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
!>#### Contains memory definitions for hlcc
!>
MODULE mo_hlcc_memory_class
    #ifndef __NO_JSBACH__
    
      USE mo_kind, ONLY: wp
    
      USE mo_jsb_model_class,        ONLY: t_jsb_model
      USE mo_jsb_class,              ONLY: Get_model
      USE mo_jsb_memory_class,       ONLY: t_jsb_memory
      USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d
      USE mo_jsb_varlist,            ONLY: FULL
    
      ! Use of processes in this module
      ! dsl4jsb_Use_processes HLCC_
    
      ! Use process configurations
      ! dsl4jsb_Use_config(HLCC_)
    
      IMPLICIT NONE
      PRIVATE
    
      PUBLIC :: t_hlcc_memory, max_no_of_vars
    
      INTEGER, PARAMETER :: max_no_of_vars = 5
    
      TYPE, EXTENDS(t_jsb_memory) :: t_hlcc_memory
    
        TYPE(t_jsb_var_real3d) :: &
          & harvest_fract_year!,      & !< This years cover fractions derived from annually read land use data [-].
          !& harvest_fract_day_delta            !< Daily change in cover fractions calculated from harvest_fract_year. !MR-TODO: doesn't exist yet
    
      CONTAINS
        PROCEDURE :: Init => Init_hlcc_memory
      END TYPE t_hlcc_memory
    
      CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_memory_class'
    
    CONTAINS
    
      SUBROUTINE Init_hlcc_memory(mem, prefix, suffix, lct_ids, model_id)
    
        ! Declarations
        USE mo_jsb_io,             ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
        USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
        USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
        !USE mo_jsb_model_usecases, ONLY: npft ! R: Note, this leads to circular dependency
        ! -------------------------------------------------------------------------------------------------- !
        CLASS(t_hlcc_memory), INTENT(inout), TARGET :: mem
        CHARACTER(len=*),      INTENT(in)             :: prefix
        CHARACTER(len=*),      INTENT(in)             :: suffix
        INTEGER,               INTENT(in)             :: lct_ids(:)
        INTEGER,               INTENT(in)             :: model_id
        ! -------------------------------------------------------------------------------------------------- !
    !    dsl4jsb_Def_config(HLCC_)
    
        TYPE(t_jsb_model), POINTER  :: model
        TYPE(t_jsb_grid),  POINTER  :: hgrid               ! Horizontal grid
        TYPE(t_jsb_vgrid), POINTER  :: ac_grid           ! Vertical grid - MR: called ac-Grid in mo_fage_interface, so kept the same name here instead of e.g. ac_vgrid
        INTEGER                     :: table
    
        CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hlcc_memory'
        ! -------------------------------------------------------------------------------------------------- !
    
        model => Get_model(model_id)
    
        table = tables(1)
    
        ! Set pointers to variables in the "grid"
        hgrid   => Get_grid(mem%grid_id)
        ac_grid    => Get_vgrid('age_classes_grid') ! MR-TODO: is this found? - vgrid is the forest age class tiles  MR: (Get_grid is a public function in mo_jsb_grid.f90...
        !> I think as long as mo_fage_config_class is run before this (registering this grid) this should work...)
    
    !    dsl4jsb_Get_config(HLCC_) !MR-TODO: maybe this is needed if a new vgrid declared in config?
    
        CALL mem%Add_var('harvest_fract_year', mem%harvest_fract_year,                              &
          & hgrid, ac_grid,                                                                 &
          & t_cf('This years harvested fractions read from forcing.'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                & ! MR-TODO: disturb has output stuff... but probs don't want that here
          & prefix, suffix,                                                                   &
          & lrestart=.FALSE., initval_r=-0.1_wp) !---------------------MR-TODO: for alcc (cf_current_year) initval=-1.0_wp ... If the code checks for negatives, this should be a'right
    
!        CALL mem%Add_var('harvest_fract_day_delta', mem%harvest_fract_day_delta,                     & ! MR-TODO: not in use yet!
!          & hgrid, ac_grid,                                                                 &
!          & t_cf('Daily rate of harvested cover fractions.'),                                        &
!          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
!          & prefix, suffix,                                                                   &
!          & output_level=FULL,                                                                &
!          & lrestart=.TRUE., lrestart_cont=.TRUE., initval_r=0.0_wp)   !-----------MR-TODO: note that initval is 0 here - daily: only "filled" on needs basis?

    
      END SUBROUTINE Init_hlcc_memory
    
    #endif
    END MODULE mo_hlcc_memory_class