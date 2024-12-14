! Meri's to-be-deleted notes: this is now modeled after flcc, but see also alcc memory
!
!
!
!---------------------------------------------
!> Contains the memory class for the hlcc process.
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
MODULE mo_hlcc_memory_class
    #ifndef __NO_JSBACH__
    
      USE mo_kind, ONLY: wp
    
      USE mo_jsb_model_class,        ONLY: t_jsb_model
      USE mo_jsb_class,              ONLY: Get_model
      USE mo_jsb_memory_class,       ONLY: t_jsb_memory
      USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
    
      dsl4jsb_Use_processes HLCC_
    
      IMPLICIT NONE
      PRIVATE
    
      PUBLIC :: t_hlcc_memory, max_no_of_vars
    
      INTEGER, PARAMETER :: max_no_of_vars = 5 ! MR-TODO : significance?
    
      !> Type definition for memory
      TYPE, EXTENDS(t_jsb_memory) :: t_hlcc_memory
    
      TYPE(t_jsb_var_real2d) ::  &
        harvested_fract,         &              ! harvested fraction of tile
        cconservation_harvest                   ! MR: I think this is just for checking for the conservation (used in connection w CALL check_carbon_conservation...) 
    
      CONTAINS
        PROCEDURE :: Init => Init_hlcc_memory
      END TYPE t_hlcc_memory
    
    
      CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_memory_class'
    
    CONTAINS
    
      SUBROUTINE Init_hlcc_memory(mem, prefix, suffix, lct_ids, model_id)
    
        USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
        USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
        USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
        USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    
        CLASS(t_hlcc_memory), INTENT(inout), TARGET :: mem
        CHARACTER(len=*),    INTENT(in)    :: prefix
        CHARACTER(len=*),    INTENT(in)    :: suffix
        INTEGER,             INTENT(in)    :: lct_ids(:)
        INTEGER,             INTENT(in)    :: model_id
    
        TYPE(t_jsb_model), POINTER :: model
        TYPE(t_jsb_grid),  POINTER :: hgrid           ! Horizontal grid
        TYPE(t_jsb_vgrid), POINTER :: surface         ! Vertical grids
        INTEGER :: table
    
        CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hlcc_memory'
    
        IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
        IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used
    
        model => Get_model(model_id)
        table = tables(1)
    
        hgrid   => Get_grid(mem%grid_id)
        surface => Get_vgrid('surface')
    
        CALL mem%Add_var( 'harvested_fract', mem%harvested_fract,                                 &
          & hgrid, surface,                                                                 &
          & t_cf('harvested_fract', 'frac', 'Fraction of tile that is harvested.'),               &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
          & prefix, suffix,                                                                 &
          & output_level=BASIC,                                                             &
          & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3
    
        IF (.NOT. model%processes(HLCC_)%p%config%active) THEN ! MR-TODO : maybe a good check? (copied as well from disturb. mem)
          CALL mem%Add_var( 'cconservation_harvest', mem%cconservation_harvest,         &
            & hgrid, surface,                                                     &
            & t_cf('cconservation_harvest', 'mol(C) m-2(canopy)',                    &
            &  'Test for c conservation on relocation for harvested area.'),         &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
            & prefix, suffix,                                                     &
            & loutput = .TRUE., lrestart=.FALSE., initval_r=0.0_wp )
        END IF
    
      END SUBROUTINE Init_hlcc_memory
    
    #endif
    END MODULE mo_hlcc_memory_class