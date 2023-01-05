!>
!! @brief Contains types for the anthropogenic forcing interface memory.
!!
!! @author
!!  Rackhun Son, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rackhun Son (2022-12-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_anthro2l_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_anthro2l_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  ! ======================================================================================================= !
  !>Type definition for anthropogenic2land memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_anthro2l_memory

    !------------------------------------------------------------------------------------------------------ !
    !> for DL-fire model in disturbance

    TYPE(t_jsb_var_real2d) ::     &
      & DEBUG_VAR,                &  !< 
      & a_temperature,            & 
      & a_rel_humidity,           &
      & population,               &
      & gdp,                      &
      & hdi,                      &
      & land_use_primf,           &
      & land_use_c4per,           &
      & land_use_urban,           &
      & land_use_secdf,           &
      & land_use_secma,           &
      & land_use_primn,           &
      & land_use_secdn,           &
      & land_use_c4ann,           &
      & land_use_pastr,           &
      & land_use_range,           &
      & land_use_c3ann,           &
      & land_use_c3nfx,           &
      & land_use_c3per,           &
      & land_use_secmb  

  CONTAINS
    PROCEDURE :: Init => Init_anthro2l_memory
  END TYPE t_anthro2l_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_anthro2l_memory_class'


CONTAINS


  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: anthropogenic2land
  !>
  SUBROUTINE Init_anthro2l_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
    USE mo_jsb_class,         ONLY: Get_model

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_anthro2l_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),                INTENT(in)            :: prefix
    CHARACTER(len=*),                INTENT(in)            :: suffix
    INTEGER,                         INTENT(in)            :: lct_ids(:)
    INTEGER,                         INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),   POINTER :: hgrid
    TYPE(t_jsb_vgrid),  POINTER :: surface
    TYPE(t_jsb_model),  POINTER :: model
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_anthro2l_memory'


    ! ----------------------------------------------------------------------------------------------------- !
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used
    ! ----------------------------------------------------------------------------------------------------- !
    model        => Get_model(model_id)
    table        = tables(1)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')

    ! ----------------------------------------------------------------------------------------------------- !
    ! DL-fire model
    !IF (model%config%model_scheme == MODEL_JSBACH_DLfire) THEN

    CALL mem%Add_var( 'DEBUG_VAR', mem%DEBUG_VAR,                            &
      & hgrid, surface,                                                      &
      & t_cf('DEBUG_VAR', '', ''),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC, &
      & initval_r=777.7_wp )

    CALL mem%Add_var( 'a_temperature', mem%a_temperature,                    &
      & hgrid, surface,                                                      &
      & t_cf('near surface temperature (anthropogenic)', '', ''),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'a_rel_humidity', mem%a_rel_humidity,                  &
      & hgrid, surface,                                                      &
      & t_cf('relative humidity (anthropogenic)', '', ''),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'population', mem%population,                          &
      & hgrid, surface,                                                      &
      & t_cf('population density (anthropogenic)', '', ''),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'gdp', mem%gdp,                                        &
      & hgrid, surface,                                                      &
      & t_cf('GDP (anthropogenic)', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'hdi', mem%hdi,                                        &
      & hgrid, surface,                                                      &
      & t_cf('HDI (anthropogenic)', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_primf', mem%land_use_primf,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_primf (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_c4per', mem%land_use_c4per,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_c4per (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_urban', mem%land_use_urban,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_urban (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_secdf', mem%land_use_secdf,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_secdf (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_secma', mem%land_use_secma,                      &
      & hgrid, surface,                                                      &
      & t_cf('land_use_secma (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_primn', mem%land_use_primn,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_primn (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_secdn', mem%land_use_secdn,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_secdn (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_c4ann', mem%land_use_c4ann,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_c4ann (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_pastr', mem%land_use_pastr,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_pastr (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_range', mem%land_use_range,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_range (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_c3ann', mem%land_use_c3ann,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_c3ann  (anthropogenic)', '', ''),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )
    
    CALL mem%Add_var( 'land_use_c3nfx', mem%land_use_c3nfx,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_c3nfx (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_c3per', mem%land_use_c3per,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_c3per (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'land_use_secmb', mem%land_use_secmb,                  &
      & hgrid, surface,                                                      &
      & t_cf('land_use_secmb (anthropogenic)', '', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    !END IF ! IF (model%config%model_scheme == MODEL_JSBACH_DLfire) THEN

  END SUBROUTINE Init_anthro2l_memory

#endif
END MODULE mo_anthro2l_memory_class
