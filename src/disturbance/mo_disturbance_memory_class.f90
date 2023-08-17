!>
!! @brief Contains the memory class for the disturb process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-6)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d    ! <Rackhun Son>

  dsl4jsb_Use_processes WLCC_, FLCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_disturb_memory, max_no_of_vars

  !INTEGER, PARAMETER :: max_no_of_vars = 40
  INTEGER, PARAMETER :: max_no_of_vars = 90   ! < Rackhun Son >

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_disturb_memory

  TYPE(t_jsb_var_real2d) ::   &
    cconservation_fire,       &
    cconservation_wind,       &
    burned_fract,             & ! burned fraction of tile
    damaged_fract,            & ! windbreak fraction of tile
    q_rel_air_climbuf,        & ! Relative humidity of lowest air layer smoothed in time.
    q_rel_air_climbuf_yday,   & ! Relative humidity of lowest air layer smoothed in time at the end of the last day.
                                ! This is only needed for cbalone forcing!
    prev_day_max_wind_10m,    & ! for windbreak
    max_wind_10m,             & ! for windbreak
    max_wind_10m_act,         & ! for windbreak

    koppen,                   & ! for Fire <Rackhun Son>
    latitude,                 & ! for Fire <Rackhun Son>
    burned_fract_infer,       & ! for Fire <Rackhun Son>
    burned_fract_test,        & ! for Fire <Rackhun Son>
    t_air_ave,                & ! for Fire <Rackhun Son>
    q_ave,                    & ! for Fire <Rackhun Son>
    q_rel_ave,                & ! for Fire <Rackhun Son>
    ws10_ave,                 & ! for Fire <Rackhun Son>
    rain_ave,                 & ! for Fire <Rackhun Son>
    swl1_ave,                 & ! for Fire <Rackhun Son>
    swl2_ave,                 & ! for Fire <Rackhun Son>
    swl3_ave,                 & ! for Fire <Rackhun Son>
    swl4_ave,                 & ! for Fire <Rackhun Son>
    swl5_ave,                 & ! for Fire <Rackhun Son>
    lai_ave                     ! for Fire <Rackhun Son>

  TYPE(t_jsb_var_real3d) ::   &
    t_air_annual,             & ! for Fire <Rackhun Son>
    rain_annual,              & ! for Fire <Rackhun Son>
    t_air_clim,               & ! for Fire <Rackhun Son>
    q_clim,                   & ! for Fire <Rackhun Son>
    q_rel_clim,               & ! for Fire <Rackhun Son>
    ws10_clim,                & ! for Fire <Rackhun Son>
    swl1_clim,                & ! for Fire <Rackhun Son>
    swl2_clim,                & ! for Fire <Rackhun Son>
    swl3_clim,                & ! for Fire <Rackhun Son>
    swl4_clim,                & ! for Fire <Rackhun Son>
    swl5_clim,                & ! for Fire <Rackhun Son>
    lai_clim,                 & ! for Fire <Rackhun Son>

    V1_memory,                & ! for Fire <Rackhun Son>
    V2_memory,                & ! for Fire <Rackhun Son>
    V3_memory,                & ! for Fire <Rackhun Son>
    V4_memory,                & ! for Fire <Rackhun Son>
    V5_memory,                & ! for Fire <Rackhun Son>
    V6_memory,                & ! for Fire <Rackhun Son>
    V7_memory,                & ! for Fire <Rackhun Son>
    V8_memory,                & ! for Fire <Rackhun Son>
    V9_memory,                & ! for Fire <Rackhun Son>
    V10_memory,               & ! for Fire <Rackhun Son>
    V11_memory,               & ! for Fire <Rackhun Son>
    V12_memory,               & ! for Fire <Rackhun Son>
    V13_memory,               & ! for Fire <Rackhun Son>
    V14_memory,               & ! for Fire <Rackhun Son>
    V15_memory,               & ! for Fire <Rackhun Son>
    V16_memory,               & ! for Fire <Rackhun Son>
    V17_memory,               & ! for Fire <Rackhun Son>
    V18_memory,               & ! for Fire <Rackhun Son>
    V19_memory,               & ! for Fire <Rackhun Son>
    V20_memory,               & ! for Fire <Rackhun Son>
    V21_memory,               & ! for Fire <Rackhun Son>
    V22_memory,               & ! for Fire <Rackhun Son>
    V23_memory,               & ! for Fire <Rackhun Son>
    V24_memory,               & ! for Fire <Rackhun Son>
    V25_memory,               & ! for Fire <Rackhun Son>
    V26_memory,               & ! for Fire <Rackhun Son>
    V27_memory,               & ! for Fire <Rackhun Son>
    V28_memory,               & ! for Fire <Rackhun Son>
    V29_memory,               & ! for Fire <Rackhun Son>
    V30_memory,               & ! for Fire <Rackhun Son>
    V31_memory,               & ! for Fire <Rackhun Son>
    V32_memory,               & ! for Fire <Rackhun Son>
    V33_memory,               & ! for Fire <Rackhun Son>
    V34_memory,               & ! for Fire <Rackhun Son>
    V35_memory,               & ! for Fire <Rackhun Son>
    V36_memory,               & ! for Fire <Rackhun Son>
    V37_memory,               & ! for Fire <Rackhun Son>
    V38_memory,               & ! for Fire <Rackhun Son>
    V39_memory,               & ! for Fire <Rackhun Son>
    Y_memory,                 & ! for Fire <Rackhun Son>
    RHC_memory,               & ! for Fire <Rackhun Son>
    R_memory,                 & ! for Fire <Rackhun Son>
    Koppen_memory               ! for Fire <Rackhun Son>


  REAL(wp), POINTER :: lightning_climatology(:,:,:)   => NULL()   ! lightning climatology (nproma,nblks,366) <Rackhun Son>
  
  REAL(wp), POINTER :: slope(:,:)                     => NULL()   ! slope                     (nproma,nblks) <Rackhun Son>
  REAL(wp), POINTER :: roughness(:,:)                 => NULL()   ! roughness                 (nproma,nblks) <Rackhun Son>
  REAL(wp), POINTER :: road_density(:,:)              => NULL()   ! road_density              (nproma,nblks) <Rackhun Son>
  REAL(wp), POINTER :: gfed_region(:,:)               => NULL()   ! GFED4 reference regions   (nproma,nblks) <Rackhun Son>
 
  CONTAINS
    PROCEDURE :: Init => Init_disturb_memory
  END TYPE t_disturb_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_memory_class'

CONTAINS

  SUBROUTINE Init_disturb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_disturb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid           ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface         ! Vertical grids
    TYPE(t_jsb_vgrid), POINTER :: daily           ! daily time grids <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: memory          ! memory           <Rackhun Son>
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)
    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    daily   => Get_vgrid('daily')      ! 366 days <Rackhun Son>
    memory  => Get_vgrid('memory')     ! training memory <Rackhun Son>

    ALLOCATE(mem%lightning_climatology(hgrid%nproma, hgrid%nblks, 1:366))   ! lightning climatology <Rackhun Son>

    ALLOCATE(mem%slope                     (hgrid%nproma, hgrid%nblks))   ! slope <Rackhun Son>
    ALLOCATE(mem%roughness                 (hgrid%nproma, hgrid%nblks))   ! roughness <Rackhun Son>
    ALLOCATE(mem%road_density              (hgrid%nproma, hgrid%nblks))   ! road_density <Rackhun Son>
    ALLOCATE(mem%gfed_region               (hgrid%nproma, hgrid%nblks))   ! GFED4 reference regions <Rackhun Son>


    CALL mem%Add_var( 'burned_fract', mem%burned_fract,                                 &
      & hgrid, surface,                                                                 &
      & t_cf('burned_fract', 'frac', 'Fraction of tile that is burned.'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'damaged_fract', mem%damaged_fract,                                &
      & hgrid, surface,                                                                  &
      & t_cf('damaged_fract', 'frac', 'Fraction of tile that is damaged by windbreak.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'q_air_clim', mem%q_rel_air_climbuf,                                            &
      & hgrid, surface,                                                                               &
      & t_cf('q_rel_air_climbuf', '-', 'Relative humiditiy of lowest air layer smoothed in time.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
      & prefix, suffix,                                                                               &
      & lrestart=.TRUE., initval_r=50.0_wp ) ! initval taken from JS3; lrestart=true only necessary if we restart within a day

    CALL mem%Add_var( 'q_air_clim_yday', mem%q_rel_air_climbuf_yday,                                    &
      & hgrid, surface,                                                                                 &
      & t_cf('q_rel_air_climbuf_yday', '-', 'Relative humiditiy of lowest air layer smoothed in time    &
                                               & at the end of the last day.'),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=BASIC,                                                                             &
      & lrestart=.FALSE., initval_r=50.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'prev_day_max_wind_10m', mem%prev_day_max_wind_10m,         &
      & hgrid, surface,                                                           &
      & t_cf('prev_day_max_wind_10m', '-', 'Previous day maximum wind speed.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m', mem%max_wind_10m,                                                &
      & hgrid, surface,                                                                                &
      & t_cf('max_wind_10m', '-', 'Maximum wind speed smoothed in time at the end of the last day..'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m_act', mem%max_wind_10m_act,             &
      & hgrid, surface,                                                     &
      & t_cf('max_wind_10m_act', '-', 'This days maximum wind speed.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp )

    IF (.NOT. model%processes(WLCC_)%p%config%active) THEN
      CALL mem%Add_var( 'cconservation_wind', mem%cconservation_wind,         &
        & hgrid, surface,                                                     &
        & t_cf('cconservation_wind', 'mol(C) m-2(canopy)',                    &
        &   'Test for c conservation on relocation for wind damaged area.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput = .TRUE., lrestart=.FALSE., initval_r=0.0_wp )
    ENDIF

    IF (.NOT. model%processes(FLCC_)%p%config%active) THEN
      CALL mem%Add_var( 'cconservation_fire', mem%cconservation_fire,         &
        & hgrid, surface,                                                     &
        & t_cf('cconservation_fire', 'mol(C) m-2(canopy)',                    &
        &  'Test for c conservation on relocation for burned area.'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput = .TRUE., lrestart=.FALSE., initval_r=0.0_wp )
    ENDIF

    ! additional variable for fire <Rackhun Son>
    CALL mem%Add_var('t_air_ave', mem%t_air_ave,                            &
      & hgrid, surface,                                                     &
      & t_cf('t_air_ave', 'K',                                              &
      &  'temperature daily average'),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)
    
    CALL mem%Add_var('t_air_annual', mem%t_air_annual,                      &
      & hgrid, daily,                                                       &
      & t_cf('t_air_annual', 'K',                                           &
      &  'temperature annual data'),                                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 293.15_wp)

    CALL mem%Add_var('q_ave', mem%q_ave,                                    &
      & hgrid, surface,                                                     &
      & t_cf('q_ave', 'g/g',                                                &
      &  'specific humidity daily average'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('q_rel_ave', mem%q_rel_ave,                            &
      & hgrid, surface,                                                     &
      & t_cf('q_rel_ave', '%',                                              &
      &  'relative humidity daily average'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('ws10_ave', mem%ws10_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('ws10_ave', 'm/s',                                             &
      &  'wind speed daily average'),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)
    
    CALL mem%Add_var('rain_ave', mem%rain_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('rain_ave', '',                                                &
      &  'precipitation daily average'),                                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('rain_annual', mem%rain_annual,                        &
      & hgrid, daily,                                                       &
      & t_cf('rain_annual', '',                                             &
      &  'precipitation annual data'),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl1_ave', mem%swl1_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('swl1_ave', 'm3/m3',                                           &
      &  'daily average of volumetric soil water content in the 1st layer'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl2_ave', mem%swl2_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('swl2_ave', 'm3/m3',                                           &
      &  'daily average of volumetric soil water content in the 2nd layer'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp) 

    CALL mem%Add_var('swl3_ave', mem%swl3_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('swl3_ave', 'm3/m3',                                           &
      &  'daily average of volumetric soil water content in the 3rd layer'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp) 

    CALL mem%Add_var('swl4_ave', mem%swl4_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('swl4_ave', 'm3/m3',                                           &
      &  'daily average of volumetric soil water content in the 4th layer'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl5_ave', mem%swl5_ave,                              &
      & hgrid, surface,                                                     &
      & t_cf('swl5_ave', 'm3/m3',                                           &
      &  'daily average of volumetric soil water content in the 5th layer'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('lai_ave', mem%lai_ave,                                &
      & hgrid, surface,                                                     &
      & t_cf('lai_ave', '-',                                                &
      &  'LAI daily average'),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp) 

    CALL mem%Add_var('t_air_clim', mem%t_air_clim,                          &
      & hgrid, daily,                                                       &
      & t_cf('t_air_clim', 'K',                                             &
      &  'temperature daily climatology'),                                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('q_clim', mem%q_clim,                                  &
      & hgrid, daily,                                                       &
      & t_cf('q_clim', 'g/g',                                               &
      &  'specific humidity daily climatology'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('q_rel_clim', mem%q_rel_clim,                          &
      & hgrid, daily,                                                       &
      & t_cf('q_rel_clim', '%',                                             &
      &  'relative humidity daily climatology'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('ws10_clim', mem%ws10_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('ws10_clim', 'm/s',                                            &
      &  'wind speed daily climatology'),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl1_clim', mem%swl1_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('swl1_clim', 'm3/m3',                                          &
      &  'daily climatology of soil water content in the 1st layer'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl2_clim', mem%swl2_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('swl2_clim', 'm3/m3',                                          &
      &  'daily climatology of soil water content in the 2nd layer'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl3_clim', mem%swl3_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('swl3_clim', 'm3/m3',                                          &
      &  'daily climatology of soil water content in the 3rd layer'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl4_clim', mem%swl4_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('swl4_clim', 'm3/m3',                                          &
      &  'daily climatology of soil water content in the 4th layer'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('swl5_clim', mem%swl5_clim,                            &
      & hgrid, daily,                                                       &
      & t_cf('swl5_clim', 'm3/m3',                                          &
      &  'daily climatology of soil water content in the 5th layer'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('lai_clim', mem%lai_clim,                              &
      & hgrid, daily,                                                       &
      & t_cf('lai_clim', '-',                                               &
      &  'LAI daily climatology'),                                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V1_memory', mem%V1_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V1_memory', '-',                                              &
      &  'V1 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V2_memory', mem%V2_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V2_memory', '-',                                              &
      &  'V2 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V3_memory', mem%V3_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V3_memory', '-',                                              &
      &  'V3 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V4_memory', mem%V4_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V4_memory', '-',                                              &
      &  'V4 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V5_memory', mem%V5_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V5_memory', '-',                                              &
      &  'V5 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V6_memory', mem%V6_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V6_memory', '-',                                              &
      &  'V6 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V7_memory', mem%V7_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V7_memory', '-',                                              &
      &  'V7 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V8_memory', mem%V8_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V8_memory', '-',                                              &
      &  'V8 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V9_memory', mem%V9_memory,                            &
      & hgrid, memory,                                                      &
      & t_cf('V9_memory', '-',                                              &
      &  'V9 memory'),                                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V10_memory', mem%V10_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V10_memory', '-',                                             &
      &  'V10 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V11_memory', mem%V11_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V11_memory', '-',                                             &
      &  'V11 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V12_memory', mem%V12_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V12_memory', '-',                                             &
      &  'V12 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V13_memory', mem%V13_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V13_memory', '-',                                             &
      &  'V13 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V14_memory', mem%V14_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V14_memory', '-',                                             &
      &  'V14 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V15_memory', mem%V15_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V15_memory', '-',                                             &
      &  'V15 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V16_memory', mem%V16_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V16_memory', '-',                                             &
      &  'V16 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V17_memory', mem%V17_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V17_memory', '-',                                             &
      &  'V17 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V18_memory', mem%V18_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V18_memory', '-',                                             &
      &  'V18 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V19_memory', mem%V19_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V19_memory', '-',                                             &
      &  'V19 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V20_memory', mem%V20_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V20_memory', '-',                                             &
      &  'V20 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V21_memory', mem%V21_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V21_memory', '-',                                             &
      &  'V21 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V22_memory', mem%V22_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V22_memory', '-',                                             &
      &  'V22 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V23_memory', mem%V23_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V23_memory', '-',                                             &
      &  'V23 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V24_memory', mem%V24_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V24_memory', '-',                                             &
      &  'V24 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V25_memory', mem%V25_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V25_memory', '-',                                             &
      &  'V25 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V26_memory', mem%V26_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V26_memory', '-',                                             &
      &  'V26 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V27_memory', mem%V27_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V27_memory', '-',                                             &
      &  'V27 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V28_memory', mem%V28_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V28_memory', '-',                                             &
      &  'V28 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V29_memory', mem%V29_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V29_memory', '-',                                             &
      &  'V29 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V30_memory', mem%V30_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V30_memory', '-',                                             &
      &  'V30 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V31_memory', mem%V31_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V31_memory', '-',                                             &
      &  'V31 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V32_memory', mem%V32_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V32_memory', '-',                                             &
      &  'V32 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V33_memory', mem%V33_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V33_memory', '-',                                             &
      &  'V33 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V34_memory', mem%V34_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V34_memory', '-',                                             &
      &  'V34 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V35_memory', mem%V35_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V35_memory', '-',                                             &
      &  'V35 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V36_memory', mem%V36_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V36_memory', '-',                                             &
      &  'V36 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V37_memory', mem%V37_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V37_memory', '-',                                             &
      &  'V37 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V38_memory', mem%V38_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V38_memory', '-',                                             &
      &  'V38 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('V39_memory', mem%V39_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('V39_memory', '-',                                             &
      &  'V39 memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('RHC_memory', mem%RHC_memory,                          &
      & hgrid, memory,                                                      &
      & t_cf('RHC_memory', '-',                                             &
      &  'RHC memory'),                                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('Y_memory', mem%Y_memory,                              &
      & hgrid, memory,                                                      &
      & t_cf('Y_memory', '-',                                               &
      &  'Y memory'),                                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('R_memory', mem%R_memory,                              &
      & hgrid, memory,                                                      &
      & t_cf('R_memory', '-',                                               &
      &  'R memory'),                                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = -999.0_wp)

    CALL mem%Add_var('Koppen_memory', mem%Koppen_memory,                    &
      & hgrid, memory,                                                      &
      & t_cf('Koppen_memory', '-',                                          &
      &  'Koppen memory'),                                                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var( 'koppen', mem%koppen,                                             &
      & hgrid, surface,                                                                 &
      & t_cf('koppen', '-', 'Koppen-Geiger classification'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var( 'latitude', mem%latitude,                                         &
      & hgrid, surface,                                                                 &
      & t_cf('latitude', '-', 'latitude'),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var( 'burned_fract_infer', mem%burned_fract_infer,                     &
      & hgrid, surface,                                                                 &
      & t_cf('burned_fract_infer', 'frac', 'Fraction of tile that is burned.'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) 

    CALL mem%Add_var( 'burned_fract_test', mem%burned_fract_test,                       &
      & hgrid, surface,                                                                 &
      & t_cf('burned_fract_test', 'frac', 'Fraction of tile that is burned.'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) 

  END SUBROUTINE Init_disturb_memory

#endif
END MODULE mo_disturb_memory_class
