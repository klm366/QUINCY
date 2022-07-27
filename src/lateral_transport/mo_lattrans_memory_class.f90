!>
!! @brief Contains the memory class for the lateral transport process.
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-07-22)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lattrans_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL
  ! USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lattrans_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 10

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_lattrans_memory
    
    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      & wtr_transport,        & !< Lateral water transport flux    [kg m-2 s-1]
      & wtr_transstor           !< Lateral transport water storage [m]

  CONTAINS
    PROCEDURE :: Init => Init_lattrans_memory
  END TYPE t_lattrans_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lattrans_memory_class'

CONTAINS

  SUBROUTINE Init_lattrans_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_lattrans_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lattrans_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    IF (TRIM(mem%owner_tile_name) == 'box') THEN
      grib2_desc = t_grib2(1,0,202, grib_bits)
    ELSE
      grib2_desc = t_grib2(255, 255, 255, grib_bits)
    END IF


    ! Common variables

    CALL mem%Add_var('wtr_transport', mem%wtr_transport,                             &
      & hgrid, surface,                                                              &
      & t_cf('wtr_transport', 'kg m-2 s-1', 'water flow between tiles'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var('wtr_transstor', mem%wtr_transstor,                             &
      & hgrid, surface,                                                              &
      & t_cf('wtr_transstor', 'm (water equivalent)', 'intermediate lateral flow storage'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_lattrans_memory

#endif
END MODULE mo_lattrans_memory_class
