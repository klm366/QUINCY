! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format e.g. "2016-04-30"                              ===
! === ... <PROCESS_NAME_LOWER_CASE> with your process name e.g. "seb" for surface energy balance                              ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementations points (Implementation:) that you    ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================


!>
!! @brief Contains the memory class for the dynveg process.
!!
!! @author
!!  Thomas Raddatz, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Thomas Raddatz (24.03.2020)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_dynveg_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL
  USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt

  ! Use of prcesses in this module
  dsl4jsb_Use_processes DYNVEG_

  ! Use process configurations
  dsl4jsb_Use_config(DYNVEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dynveg_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_dynveg_memory
    
    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      & seconds_day, &
      & seconds_month, &
      & temp_sum_day, &
      & temp_sum_month, &
      & min_mmtemp_of_yr, &
      & max_mmtemp_of_yr, &
      & min_mmtemp20, &
      & max_mmtemp20, &
      & gdd_sum_year, &
      & gdd_prev_year, &
      & bare_fpc, &
      & sum_green_bio_memory, &
      & desert_fpc
    TYPE(t_jsb_var_real3d) :: &
      & act_fpc, &
      & pot_fpc, &
      & cover_fract_pot, &
      & cover_fract, &
      & bio_exist
    ! Additional variables for LAND lct_type
!!$  TR  TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & s_star,               & !< Surface dry static energy (s^star, see manual)    [m2 s-2]

      
    ! Additional variables for LAKE lct_type
!!$ TR    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_dynveg_memory
  END TYPE t_dynveg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dynveg_memory_class'

CONTAINS

  SUBROUTINE Init_dynveg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_dynveg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, pft_vgrid           ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_dynveg_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid     => Get_grid(mem%grid_id)
    surface   => Get_vgrid('surface')
    pft_vgrid => Get_vgrid('pfts')

    CALL mem%Add_var('seconds_day', mem%seconds_day,                          &
      & hgrid, surface,                                                       &
      & t_cf('seconds_day', 'sum of seconds in day', 'time sum day'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('seconds_month', mem%seconds_month,                      &
      & hgrid, surface,                                                       &
      & t_cf('seconds_month', 'sum of seconds in month', 'time sum month'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('temp_sum_day', mem%temp_sum_day,                        &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_day', 'sum of air temperature in day', 'temperature sum day'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('temp_sum_month', mem%temp_sum_month,                    &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_month', 'sum of air temperature in month', 'temperature sum month'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('min_mmtemp_of_yr', mem%min_mmtemp_of_yr,                &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp_of_yr', 'min temperature of year', 'temperature min'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp_of_yr', mem%max_mmtemp_of_yr,                      &
      & hgrid, surface,                                                             &
      & t_cf('max_mmtemp_of_yr', 'max temperature of year', 'temperature max'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('min_mmtemp20', mem%min_mmtemp20,                  &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp20', 'min temperature climatology', 'temp min climate'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp20', mem%max_mmtemp20,                  &
      & hgrid, surface,                                                       &
      & t_cf('max_mmtemp20', 'max temperature climatology', 'temp max climate'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('gdd_sum_year', mem%gdd_sum_year,                  &
      & hgrid, surface,                                                       &
      & t_cf('gdd_sum_year', 'gdd current year', 'gdd sum year'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('gdd_prev_year', mem%gdd_prev_year,                &
      & hgrid, surface,                                                       &
      & t_cf('gdd_prev_year', 'gdd previous year', 'gdd prev year'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('act_fpc', mem%act_fpc,                            &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('act_fpc', 'actual foliage projective cover', 'act fpc'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bare_fpc', mem%bare_fpc,                          &
      & hgrid, surface,                                                       &
      & t_cf('bare_fpc', 'bare foliage projective cover', 'bare fpc'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('pot_fpc', mem%pot_fpc,                            &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('pot_fpc', 'potential foliage projective cover', 'pot fpc'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('cover_fract_pot', mem%cover_fract_pot,            &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract_pot', 'potential cover fractions', 'cover fract pot'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bio_exist', mem%bio_exist,                        &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('bio_exist', 'bio-climatic limits', 'bio exist'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('sum_green_bio_memory', mem%sum_green_bio_memory,  &
      & hgrid, surface,                                                       &
      & t_cf('sum_green_bio_memory', 'desert memory', 'sum_green_bio_memory'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('desert_fpc', mem%desert_fpc,                      &
      & hgrid, surface,                                                       &
      & t_cf('desert_fpc', 'desert foliage projective cover', 'desert fpc'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

! diagnostic output variables
    CALL mem%Add_var('cover_fract', mem%cover_fract,                    &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract', 'PFT cover fractions', 'cover fract'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    !X add your common variables here e.g:
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

    !X  CALL mem%Add_var('t_old', mem%t_old,                                             &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

!!$    END IF

    ! Additional variables for LAND lct
    !X add your LAND-lct variables here e.g:
    !X IF (One_of(LAND_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%s_star%ptr)) THEN

    !X  CALL mem%Add_var('s_sfc', mem%s_star,                                                  &
    !X    & hgrid, surface,                                                                    &
    !X    & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
    !X    & prefix, suffix,                                                                    &
    !X    & initval_r=2.9E5_wp )

    !X  CALL mem%Add_var('qsat_sfc', mem%qsat_star,                                            &
    !X    & hgrid, surface,                                                                    &
    !X    & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
    !X    & prefix, suffix,                                                                    &
    !X    & initval_r=0.0075_wp )

    !X END IF

    ! Additional variables for lakes
    !X add your LAKE-lct variables here e.g:
    !XIF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_wtr%ptr)) THEN

    !X  CALL mem%Add_var('t_wtr', mem%t_wtr,                                             &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

    !X  CALL mem%Add_var('t_ice', mem%t_ice,                                             &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=tmelt )

    !X END IF

  END SUBROUTINE Init_dynveg_memory

#endif
END MODULE mo_dynveg_memory_class
