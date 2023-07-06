! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format e.g. "2016-04-30"                              ===
! === ... lentic with your process name e.g. "seb" for surface energy balance                              ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementations points (Implementation:) that you    ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================


!>
!! @brief Contains the memory class for the lentic process.
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (22-05-2023)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lentic_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL
  !USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt
  
  ! Use of prcesses in this module
  dsl4jsb_Use_processes LENTIC_

  ! Use process configurations
  dsl4jsb_Use_config(LENTIC_)
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lentic_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_lentic_memory

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & T_sfc, &  !< lake surface temperature wtr
      & h_ice, &  !< lake ice height
      & h_ML,  &  !< lake mixed-layer height
      & h_snow    !< lake snow-layer height
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_lentic_memory
  END TYPE t_lentic_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lentic_memory_class'

CONTAINS

  SUBROUTINE Init_lentic_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_lentic_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lentic_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0) THEN

      ! CALL mem%Add_var( 'depth_w', mem%depth_w,                                            &
      !   & hgrid, surface,                                                              &
      !   & t_cf('lake_depth', 'm', 'lake mean wtr depth'),                  &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      !   & prefix, suffix,                                                              &
      !   & output_level=BASIC, initval_r=280.0_wp )

      CALL mem%Add_var( 'h_ice', mem%h_ice,                                          &
        & hgrid, surface,                                                              &
        & t_cf('lice_height', 'm', 'lake ice height'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=280.0_wp )

      CALL mem%Add_var( 'h_ML', mem%h_ML,                                          &
        & hgrid, surface,                                                              &
        & t_cf('lake_mixedLayer', 'm', 'lake mixed-layer height'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=280.0_wp )
        
      CALL mem%Add_var( 'h_snow', mem%h_snow,                                          &
        & hgrid, surface,                                                              &
        & t_cf('lake_snowthickness', 'm', 'lake snow-layer height'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=280.0_wp )

      CALL mem%Add_var( 'T_sfc', mem%T_sfc,                                          &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=280.0_wp )
    END IF

  END SUBROUTINE Init_lentic_memory

#endif
END MODULE mo_lentic_memory_class