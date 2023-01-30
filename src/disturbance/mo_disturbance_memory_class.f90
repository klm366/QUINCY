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

  TYPE(t_jsb_var_real3d) ::   &
    t_air_clim,               & ! for Fire <Rackhun Son>
    q_clim,                   & ! for Fire <Rackhun Son>
    q_rel_clim,               & ! for Fire <Rackhun Son>
    ws10_clim,                & ! for Fire <Rackhun Son>
    swl1_clim,                & ! for Fire <Rackhun Son>
    swl2_clim,                & ! for Fire <Rackhun Son>
    swl3_clim,                & ! for Fire <Rackhun Son>
    swl4_clim,                & ! for Fire <Rackhun Son>
    lai_clim                    ! for Fire <Rackhun Son>
    !lstm_hidden,             & ! for Fire <Rackhun Son>
    !lstm_cell,               & ! for Fire <Rackhun Son>

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
    TYPE(t_jsb_vgrid), POINTER :: daily             ! daily time grids <Rackhun Son>
    !TYPE(t_jsb_vgrid), POINTER :: d_lstm            ! LSTM hidden grids <Rackhun Son>
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)
    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    daily   => Get_vgrid('daily')      ! 366 days <Rackhun Son>
    !d_lstm  => Get_vgrid('d_lstm')     ! LSTM hidden dimension <Rackhun Son>

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

    !CALL mem%Add_var('lstm_hidden', mem%lstm_hidden,                        &
    !  & hgrid, d_lstm,                                                      &
    !  & t_cf('lstm_hidden', '-',                                            &
    !  &  'lstm hidden'),                                                    &
    !  & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !  & prefix, suffix,                                                     &
    !  & lrestart = .TRUE., initval_r = 0.0_wp)

    !CALL mem%Add_var('lstm_cell', mem%lstm_cell,                            &
    !  & hgrid, d_lstm,                                                      &
    !  & t_cf('lstm_cell', '-',                                              &
    !  &  'lstm cell'),                                                      &
    !  & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
    !  & prefix, suffix,                                                     &
    !  & lrestart = .TRUE., initval_r = 0.0_wp)

  END SUBROUTINE Init_disturb_memory

#endif
END MODULE mo_disturb_memory_class
