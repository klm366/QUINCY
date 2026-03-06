!> Contains the memory class for the land_ice process.
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
!>#### Contains memory definitions for the land-ice coupling process
!>
MODULE mo_land_ice_memory_class
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_LAND_ICE__)

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_land_ice_memory, max_no_of_land_ice_vars

  INTEGER, PARAMETER :: max_no_of_land_ice_vars = 20

  !> Type definition for land_ice memory
  TYPE, EXTENDS(t_jsb_memory) :: t_land_ice_memory

    !-- Variables sent to the ice-sheet model (accumulated/averaged over coupling interval)
    TYPE(t_jsb_var_real2d) :: &
      & surface_temp,             & !< Surface temperature on glacier tile                    [K]
      & surface_mass_balance,     & !< Surface mass balance (snowfall - sublimation - melt)   [kg m-2 s-1]
      & ground_heat_flux,         & !< Ground heat flux into glacier ice                      [W m-2]
      & rainfall_rate,            & !< Rainfall rate on glacier                               [kg m-2 s-1]
      & snowfall_rate               !< Snowfall rate on glacier                               [kg m-2 s-1]

    !-- Variables received from the ice-sheet model
    TYPE(t_jsb_var_real2d) :: &
      & ice_thickness,            & !< Ice-sheet thickness                                    [m]
      & ice_surface_elevation,    & !< Ice-sheet surface elevation                            [m]
      & glacier_fraction,         & !< Updated glacier extent from ice dynamics               []
      & basal_melt_rate,          & !< Basal melting rate                                     [kg m-2 s-1]
      & calving_flux                !< Calving flux at ice margins                            [kg m-2 s-1]

  CONTAINS
    PROCEDURE :: Init => Init_land_ice_memory
  END TYPE t_land_ice_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_land_ice_memory_class'

CONTAINS

  SUBROUTINE Init_land_ice_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_land_ice_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_land_ice_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    ! ------------------------------------------------------------------
    ! Variables sent to the ice-sheet model
    ! ------------------------------------------------------------------

    IF (One_of(GLACIER_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var( 'lice_surface_temp', mem%surface_temp,                             &
        & hgrid, surface,                                                                  &
        & t_cf('lice_surface_temp', 'K', 'Surface temperature on glacier (for ice-sheet coupling)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=273.15_wp )

      CALL mem%Add_var( 'lice_smb', mem%surface_mass_balance,                              &
        & hgrid, surface,                                                                  &
        & t_cf('lice_smb', 'kg m-2 s-1', 'Surface mass balance on glacier'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_ground_hflx', mem%ground_heat_flux,                          &
        & hgrid, surface,                                                                  &
        & t_cf('lice_ground_hflx', 'W m-2', 'Ground heat flux into glacier ice'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_rainfall', mem%rainfall_rate,                                &
        & hgrid, surface,                                                                  &
        & t_cf('lice_rainfall', 'kg m-2 s-1', 'Rainfall rate on glacier'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_snowfall', mem%snowfall_rate,                                &
        & hgrid, surface,                                                                  &
        & t_cf('lice_snowfall', 'kg m-2 s-1', 'Snowfall rate on glacier'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      ! ------------------------------------------------------------------
      ! Variables received from the ice-sheet model
      ! ------------------------------------------------------------------

      CALL mem%Add_var( 'lice_ice_thickness', mem%ice_thickness,                           &
        & hgrid, surface,                                                                  &
        & t_cf('lice_ice_thickness', 'm', 'Ice-sheet thickness'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_ice_elevation', mem%ice_surface_elevation,                   &
        & hgrid, surface,                                                                  &
        & t_cf('lice_ice_elevation', 'm', 'Ice-sheet surface elevation'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_glacier_fract', mem%glacier_fraction,                        &
        & hgrid, surface,                                                                  &
        & t_cf('lice_glacier_fract', '', 'Glacier fraction from ice-sheet model'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_basal_melt', mem%basal_melt_rate,                            &
        & hgrid, surface,                                                                  &
        & t_cf('lice_basal_melt', 'kg m-2 s-1', 'Basal melting rate from ice-sheet model'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

      CALL mem%Add_var( 'lice_calving', mem%calving_flux,                                  &
        & hgrid, surface,                                                                  &
        & t_cf('lice_calving', 'kg m-2 s-1', 'Calving flux at ice margins'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                  &
        & initval_r=0.0_wp )

    END IF

  END SUBROUTINE Init_land_ice_memory

#endif
END MODULE mo_land_ice_memory_class
