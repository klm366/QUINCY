!>
!! @brief Contains the memory class for the flake process.
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (22-05-2023)
!!cd
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_memory_class
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
  dsl4jsb_Use_processes FLAKE_

  ! Use process configurations
  dsl4jsb_Use_config(FLAKE_)
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_flake_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_flake_memory

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & c_t,          &  !< lake tc shape factor
      & T_lice,       &  !< lake ice temperature 
      ! & T_sfc,        &  !< lake surface temperature wtr
      & t_bot,        &  !< lake bottom temperature at water/sediment interface
      & t_wml,        &  !< lake mixed-layer temperaturer wtr
      & t_mnw,        &  !< lake mean wtr temperature
      & h_ice,        &  !< lake ice height
      & h_ML,         &  !< lake mixed-layer height
      & h_snow,       &  !< lake snow-layer height
      & fetch,        &  !< lake fetch
      & w_star,       &  !< convective velocity scale 
      & z0t_sf, &
      & z0u_sf, &
      & Q_bot!, & !< ! Heat flux through the water-bottom sediment interface [W m^{-2}]
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_flake_memory
  END TYPE t_flake_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_memory_class'

CONTAINS

  SUBROUTINE Init_flake_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_flake_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: lib_id
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_flake_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var( 'c_t', mem%c_t,                                                &
        & hgrid, surface,                                                              &
        & t_cf('l_tc_shapefactor', 'm', 'shape factor (thermocline)'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'h_ice', mem%h_ice,                                            &
        & hgrid, surface,                                                              &
        & t_cf('lice_height', 'm', 'lake ice height'),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'h_ML', mem%h_ML,                                              &
        & hgrid, surface,                                                              &
        & t_cf('lake_mixedLayer', 'm', 'lake mixed-layer height'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )
        
      CALL mem%Add_var( 'h_snow', mem%h_snow,                                          &
        & hgrid, surface,                                                              &
        & t_cf('lake_snowthickness', 'm', 'lake snow-layer height'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )
        
      ! CALL mem%Add_var( 'T_sfc', mem%T_sfc,                                            &
      !   & hgrid, surface,                                                              &
      !   & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr'),                   &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      !   & prefix, suffix,                                                              &
      !   & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'T_lice', mem%T_lice,                                          &
        & hgrid, surface,                                                              &
        & t_cf('lice_temp', 'K', 'lake ice temperature'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'T_bot', mem%T_bot,                                            &
        & hgrid, surface,                                                              &
        & t_cf('lake_temp_bot', 'K', 'lake bottom temperature wtr/sediment'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )      
        
      CALL mem%Add_var( 'T_wML', mem%T_wML,                                            &
        & hgrid, surface,                                                              &
        & t_cf('lake_temp_wtr_ML', 'K', 'lake mixed-layer temperature wtr'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )  

      CALL mem%Add_var( 'T_mnw', mem%T_mnw,                                            &
        & hgrid, surface,                                                              &
        & t_cf('lake_temp_wtr_mean', 'K', 'lake mean temperature wtr'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp ) 

      CALL mem%Add_var( 'fetch', mem%fetch,                                            &
        & hgrid, surface,                                                              &
        & t_cf('fetch', 'm', 'lake fetch'),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'Q_bot', mem%Q_bot,                                            &
        & hgrid, surface,                                                              &
        & t_cf('Q_bot', 'W m-2', 'Heat flux through water-sediment interface (lakes)'),&
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

      CALL mem%Add_var( 'w_star', mem%w_star,                                          &
        & hgrid, surface,                                                              &
        & t_cf('w_star', 'm s-1', 'lake convective velocity scale'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=BASIC, initval_r=0.0_wp )

        CALL mem%Add_var( 'z0t_sf', mem%z0t_sf,                                        &
        & hgrid, surface,                                                              &
        & t_cf('z0t_sf', 'm', 'lake roughness length (heat)'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL, initval_r=0.0_wp )

        CALL mem%Add_var( 'z0u_sf', mem%z0u_sf,                                        &
        & hgrid, surface,                                                              &
        & t_cf('z0u_sf', 'm', 'lake roughness length (momentum)'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL, initval_r=0.0_wp )
    END IF

  END SUBROUTINE Init_flake_memory

#endif
END MODULE mo_flake_memory_class