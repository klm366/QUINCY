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

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d    ! <Rackhun Son>

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_disturb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

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

    t_air_ave,                & ! for Fire <Rackhun Son>
    q_ave,                    & ! for Fire <Rackhun Son>
    q_rel_ave,                & ! for Fire <Rackhun Son>
    ws10_ave,                 & ! for Fire <Rackhun Son>
    rain_ave,                 & ! for Fire <Rackhun Son>
    swl1_ave,                 & ! for Fire <Rackhun Son>
    swl2_ave,                 & ! for Fire <Rackhun Son>
    swl3_ave,                 & ! for Fire <Rackhun Son>
    swl4_ave,                 & ! for Fire <Rackhun Son>
    lai_ave                     ! for Fire <Rackhun Son>
    !fuel_yearly                 ! for Fire <Rackhun Son>

  TYPE(t_jsb_var_real3d) ::   &
    t_air_clim,               & ! for Fire <Rackhun Son>
    q_clim,                   & ! for Fire <Rackhun Son>
    q_rel_clim,               & ! for Fire <Rackhun Son>
    ws10_clim,                & ! for Fire <Rackhun Son>
    swl1_clim,                & ! for Fire <Rackhun Son>
    swl2_clim,                & ! for Fire <Rackhun Son>
    swl3_clim,                & ! for Fire <Rackhun Son>
    swl4_clim,                & ! for Fire <Rackhun Son>
    lai_clim,                 & ! for Fire <Rackhun Son>
    weather_lstm_hidden,      & ! for Fire <Rackhun Son>
    weather_lstm_cell,        & ! for Fire <Rackhun Son>
    land_lstm_hidden,         & ! for Fire <Rackhun Son>
    land_lstm_cell              ! for Fire <Rackhun Son>

  REAL(wp), POINTER :: lightning_climatology(:,:,:)   => NULL()   ! lightning climatology (nproma,nblks,366) <Rackhun Son>
  
  REAL(wp), POINTER :: slope(:,:)                     => NULL()   ! slope                     (nproma,nblks) <Rackhun Son>
  REAL(wp), POINTER :: roughness(:,:)                 => NULL()   ! roughness                 (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: topographic_position_idx(:,:)  => NULL()   ! topographic_position_idx  (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: terrain_ruggedness_idx(:,:)    => NULL()   ! terrain_ruggedness_idx    (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: vector_ruggedness_measure(:,:) => NULL()   ! vector_ruggedness_measure (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: profile_curvature(:,:)         => NULL()   ! profile_curvature         (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: tangential_curvature(:,:)      => NULL()   ! tangential_curvature      (nproma,nblks) <Rackhun Son>

  !REAL(wp), POINTER :: population_density(:,:)        => NULL()   ! population_density        (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: gdp(:,:)                       => NULL()   ! gdp                       (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: hdi(:,:)                       => NULL()   ! hdi                       (nproma,nblks) <Rackhun Son>
  REAL(wp), POINTER :: road_density(:,:)              => NULL()   ! road_density              (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_c3ann(:,:)      => NULL()   ! land_use_state_c3ann      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_c3nfx(:,:)      => NULL()   ! land_use_state_c3nfx      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_c3per(:,:)      => NULL()   ! land_use_state_c3per      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_c4ann(:,:)      => NULL()   ! land_use_state_c4ann      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_c4per(:,:)      => NULL()   ! land_use_state_c4per      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_pastr(:,:)      => NULL()   ! land_use_state_pastr      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_primf(:,:)      => NULL()   ! land_use_state_primf      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_primn(:,:)      => NULL()   ! land_use_state_primn      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_range(:,:)      => NULL()   ! land_use_state_range      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_secdf(:,:)      => NULL()   ! land_use_state_secdf      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_secdn(:,:)      => NULL()   ! land_use_state_secdn      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_secma(:,:)      => NULL()   ! land_use_state_secma      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_secmb(:,:)      => NULL()   ! land_use_state_secmb      (nproma,nblks) <Rackhun Son>
  !REAL(wp), POINTER :: land_use_state_urban(:,:)      => NULL()   ! land_use_state_urban      (nproma,nblks) <Rackhun Son>

  REAL(wp), POINTER :: gfed_region(:,:)               => NULL()   ! GFED4 reference regions   (nproma,nblks) <Rackhun Son>

  !REAL(wp), POINTER :: population_density(:,:,:)    => NULL()   ! population density        (nproma,nblks,366) <Rackhun Son>

  REAL(wp), POINTER :: t_air_local_clim(:,:,:)        => NULL()   ! temperature           (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: q_rel_local_clim(:,:,:)        => NULL()   ! relative humidity     (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: q_local_clim(:,:,:)            => NULL()   ! specific humidity     (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: ws10_local_clim(:,:,:)         => NULL()   ! wind speed            (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: swvl1_local_clim(:,:,:)        => NULL()   ! soil water content    (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: swvl2_local_clim(:,:,:)        => NULL()   ! soil water content    (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: swvl3_local_clim(:,:,:)        => NULL()   ! soil water content    (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: swvl4_local_clim(:,:,:)        => NULL()   ! soil water content    (nproma,nblks,366)   <Rackhun Son>
  REAL(wp), POINTER :: lai_local_clim(:,:,:)          => NULL()   ! lai                   (nproma,nblks,366)   <Rackhun Son>
 
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

    TYPE(t_jsb_grid),  POINTER :: hgrid             ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface           ! Vertical grids
    TYPE(t_jsb_vgrid), POINTER :: daily             ! daily time grids <Rackhun Son>
    !TYPE(t_jsb_vgrid), POINTER :: monthly           ! monthlly time grids <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: d_weather_lstm    ! LSTM (weather) hidden grids <Rackhun Son>
    TYPE(t_jsb_vgrid), POINTER :: d_land_lstm       ! LSTM (land) hidden grids <Rackhun Son>
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    daily   => Get_vgrid('daily')      ! 366 days <Rackhun Son>
    !monthly  => Get_vgrid('monthly')    ! 12 months <Rackhun Son>
    d_weather_lstm => Get_vgrid('d_weather_lstm')  ! LSTM (weather) hidden dimension <Rackhun Son>
    d_land_lstm    => Get_vgrid('d_land_lstm')     ! LSTM (land) hidden dimension <Rackhun Son>

    ALLOCATE(mem%lightning_climatology(hgrid%nproma, hgrid%nblks, 1:366))   ! lightning climatology <Rackhun Son>

    ALLOCATE(mem%slope                     (hgrid%nproma, hgrid%nblks))   ! slope <Rackhun Son>
    ALLOCATE(mem%roughness                 (hgrid%nproma, hgrid%nblks))   ! roughness <Rackhun Son>
    !ALLOCATE(mem%topographic_position_idx  (hgrid%nproma, hgrid%nblks))   ! topographic position idx <Rackhun Son>
    !ALLOCATE(mem%terrain_ruggedness_idx    (hgrid%nproma, hgrid%nblks))   ! terrain ruggedness idx <Rackhun Son>
    !ALLOCATE(mem%vector_ruggedness_measure (hgrid%nproma, hgrid%nblks))   ! vector ruggedness measure <Rackhun Son>
    !ALLOCATE(mem%profile_curvature         (hgrid%nproma, hgrid%nblks))   ! profile curvature <Rackhun Son>
    !ALLOCATE(mem%tangential_curvature      (hgrid%nproma, hgrid%nblks))   ! tangential curvature <Rackhun Son>

    !ALLOCATE(mem%population_density        (hgrid%nproma, hgrid%nblks))   ! population density <Rackhun Son>
    !ALLOCATE(mem%gdp                       (hgrid%nproma, hgrid%nblks))   ! gdp <Rackhun Son>
    !ALLOCATE(mem%hdi                       (hgrid%nproma, hgrid%nblks))   ! hdi <Rackhun Son>
    ALLOCATE(mem%road_density              (hgrid%nproma, hgrid%nblks))   ! road_density <Rackhun Son>
    !ALLOCATE(mem%land_use_state_c3ann      (hgrid%nproma, hgrid%nblks))   ! land_use_state_c3ann <Rackhun Son>
    !ALLOCATE(mem%land_use_state_c3nfx      (hgrid%nproma, hgrid%nblks))   ! land_use_state_c3nfx <Rackhun Son>
    !ALLOCATE(mem%land_use_state_c3per      (hgrid%nproma, hgrid%nblks))   ! land_use_state_c3per <Rackhun Son>
    !ALLOCATE(mem%land_use_state_c4ann      (hgrid%nproma, hgrid%nblks))   ! land_use_state_c4ann <Rackhun Son>
    !ALLOCATE(mem%land_use_state_c4per      (hgrid%nproma, hgrid%nblks))   ! land_use_state_c4per <Rackhun Son>
    !ALLOCATE(mem%land_use_state_pastr      (hgrid%nproma, hgrid%nblks))   ! land_use_state_pastr <Rackhun Son>
    !ALLOCATE(mem%land_use_state_primf      (hgrid%nproma, hgrid%nblks))   ! land_use_state_primf <Rackhun Son>
    !ALLOCATE(mem%land_use_state_primn      (hgrid%nproma, hgrid%nblks))   ! land_use_state_primn <Rackhun Son>
    !ALLOCATE(mem%land_use_state_range      (hgrid%nproma, hgrid%nblks))   ! land_use_state_range <Rackhun Son>
    !ALLOCATE(mem%land_use_state_secdf      (hgrid%nproma, hgrid%nblks))   ! land_use_state_secdf <Rackhun Son>
    !ALLOCATE(mem%land_use_state_secdn      (hgrid%nproma, hgrid%nblks))   ! land_use_state_secdn <Rackhun Son>
    !ALLOCATE(mem%land_use_state_secma      (hgrid%nproma, hgrid%nblks))   ! land_use_state_secma <Rackhun Son>
    !ALLOCATE(mem%land_use_state_secmb      (hgrid%nproma, hgrid%nblks))   ! land_use_state_secmb <Rackhun Son>
    !ALLOCATE(mem%land_use_state_urban      (hgrid%nproma, hgrid%nblks))   ! land_use_state_urban <Rackhun Son>

    ALLOCATE(mem%gfed_region               (hgrid%nproma, hgrid%nblks))   ! GFED4 reference regions <Rackhun Son>

    !ALLOCATE(mem%population_density   (hgrid%nproma, hgrid%nblks, 1:366))   ! population density    <Rackhun Son>

    ALLOCATE(mem%t_air_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! temperature           <Rackhun Son>
    ALLOCATE(mem%q_rel_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! relative humidity     <Rackhun Son>
    ALLOCATE(mem%q_local_clim         (hgrid%nproma, hgrid%nblks, 1:366))   ! specific humidity     <Rackhun Son>
    ALLOCATE(mem%ws10_local_clim      (hgrid%nproma, hgrid%nblks, 1:366))   ! wind speed            <Rackhun Son>
    ALLOCATE(mem%swvl1_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! soil water content    <Rackhun Son>
    ALLOCATE(mem%swvl2_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! soil water content    <Rackhun Son>
    ALLOCATE(mem%swvl3_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! soil water content    <Rackhun Son>
    ALLOCATE(mem%swvl4_local_clim     (hgrid%nproma, hgrid%nblks, 1:366))   ! soil water content    <Rackhun Son>
    ALLOCATE(mem%lai_local_clim       (hgrid%nproma, hgrid%nblks, 1:366))   ! lai                   <Rackhun Son>
    !ALLOCATE(mem%t_air_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! temperature           <Rackhun Son>
    !ALLOCATE(mem%q_rel_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! relative humidity     <Rackhun Son>
    !ALLOCATE(mem%q_local_clim         (hgrid%nproma, hgrid%nblks, 1:12))   ! specific humidity     <Rackhun Son>
    !ALLOCATE(mem%ws10_local_clim      (hgrid%nproma, hgrid%nblks, 1:12))   ! wind speed            <Rackhun Son>
    !ALLOCATE(mem%swvl1_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! soil water content    <Rackhun Son>
    !ALLOCATE(mem%swvl2_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! soil water content    <Rackhun Son>
    !ALLOCATE(mem%swvl3_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! soil water content    <Rackhun Son>
    !ALLOCATE(mem%swvl4_local_clim     (hgrid%nproma, hgrid%nblks, 1:12))   ! soil water content    <Rackhun Son>
    !ALLOCATE(mem%lai_local_clim       (hgrid%nproma, hgrid%nblks, 1:12))   ! lai                   <Rackhun Son>


    CALL mem%Add_var( 'burned_fract', mem%burned_fract,                                 &
      & hgrid, surface,                                                                 &
      & t_cf('burned_fract', 'frac', 'Fraction of tile that is burned.'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'damaged_fract', mem%damaged_fract,                                &
      & hgrid, surface,                                                                  &
      & t_cf('damaged_fract', 'frac', 'Fraction of tile that is damaged by windbreak.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

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

    CALL mem%Add_var( 'cconservation_wind', mem%cconservation_wind,         &
      & hgrid, surface,                                                     &
      & t_cf('cconservation_wind', 'mol(C) m-2(canopy)',                    &
      &   'Test for c conservation on relocation for wind damaged area.'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'cconservation_fire', mem%cconservation_fire,         &
      & hgrid, surface,                                                     &
      & t_cf('cconservation_fire', 'mol(C) m-2(canopy)',                    &
      &  'Test for c conservation on relocation for burned area.'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp )

    ! additional variable for fire <Rackhun Son>
    CALL mem%Add_var('t_air_ave', mem%t_air_ave,                            &
      & hgrid, surface,                                                     &
      & t_cf('t_air_ave', 'K',                                              &
      &  'temperature daily average'),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

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

    CALL mem%Add_var('lai_clim', mem%lai_clim,                              &
      & hgrid, daily,                                                       &
      & t_cf('lai_clim', '-',                                               &
      &  'LAI daily climatology'),                                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    !CALL mem%Add_var('fuel_yearly', mem%fuel_yearly,                        &
    !  & hgrid, surface,                                                     &
    !  & t_cf('fuel_yearly', '-',                                            &
    !  &  'Fuel yearly mean'),                                               &
    !  & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !  & prefix, suffix,                                                     &
    !  & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('weather_lstm_hidden', mem%weather_lstm_hidden,        &
      & hgrid, d_weather_lstm,                                              &
      & t_cf('weather_lstm_hidden', '-',                                    &
      &  'weather lstm hidden'),                                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('weather_lstm_cell', mem%weather_lstm_cell,            &
      & hgrid, d_weather_lstm,                                              &
      & t_cf('weather_lstm_cell', '-',                                      &
      &  'weather lstm cell'),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('land_lstm_hidden', mem%land_lstm_hidden,              &
      & hgrid, d_land_lstm,                                                 &
      & t_cf('land_lstm_hidden', '-',                                       &
      &  'land lstm hidden'),                                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

    CALL mem%Add_var('land_lstm_cell', mem%land_lstm_cell,                  &
      & hgrid, d_land_lstm,                                                 &
      & t_cf('land_lstm_cell', '-',                                         &
      &  'land lstm cell'),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart = .TRUE., initval_r = 0.0_wp)

  END SUBROUTINE Init_disturb_memory

#endif
END MODULE mo_disturb_memory_class
