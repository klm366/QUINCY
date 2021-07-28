!>
!! @brief Contains the memory class for the surface energy balance process.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version in mo_srf_memory_class.f90 by Reiner Schnur (2013-04-04)
!! Separated out the surface energy balance part by Reiner Schnur (2016-03-21)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_seb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
  USE mo_jsb_physical_constants, ONLY: tmelt

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_seb_memory
    !
    ! Common variables
    TYPE(t_jsb_var_real2d) ::   &
      & t,                      & !< Surface temperature                                             [K]
      & t_old,                  & !< Surface temperature at previous timestep                        [K]
      & t_unfilt,               & !< Surface temperature, unfiltered (before Asselin)                [K]
      & heat_cap,               & !< Heat capacity of surface layer                                  [J m-2 K-1]
      & heat_cap_old,           & !< Heat capacity of surface layer at previous timestep             [J m-2 K-1]
      & latent_hflx,            & !< Latent heat flux                                                [W m-2]
      & sensible_hflx,          & !< Sensible heat flux                                              [W m-2]
      & forc_hflx                 !< Additional heat flux at surface                                 [W m-2]
                                  !! (additional to radiative fluxes and latent and sensible
                                  !!  heat fluxes, e.g. ground heat flux)
    !
    ! Additional variables for LAND lct_type
    TYPE(t_jsb_var_real2d) ::   &
      & s_star,                 & !< Surface dry static energy (s^star, see manual)                  [m2 s-2]
                                  !! This is the s that is used in the vertical diffusion
                                  !! scheme in the atmosphere
      & qsat_star,              & !< Saturation specific humidity at surface (qsat^star, see manual) [kg kg-1]
                                  !! This is the qsat_srf that is used in the vertical diffusion
                                  !! scheme in the atmosphere
      & t_unfilt_old,           & !< Surface temperature at previous timestep, unfiltered            [K]
                                  !! (before Asselin)
      & t_filt,                 & !< Surface temperature, filtered (after Asselin)                   [K]
      & t_eff,                  & !< Effective surface temperature (see manual)                      [K]
                                  !! This is the surface temperature used in the atmospheric
                                  !! radiative heating
      & latent_hflx_lnd,        & !< Latent heat flux                                                [W m-2]
      & sensible_hflx_lnd,      & !< Sensible heat flux                                              [W m-2]
      & previous_day_temp_mean, & ! mean day air-temperature [Celsius] of previous day,
                                  ! because not known for ..
                                  ! current day; updated at the beginning of each new day (stream)
      & previous_day_temp_min,  & ! minimum air-temperature [Celsius] of previous day, because not known for  ..
                                  ! current day; updated at the beginning of each new day (stream)
      & previous_day_temp_max,  & ! maximum air-temperature [Celsius] of previous day, because not known for   ..
                                  ! current day; updated at the beginning of each new day (stream)
      & pseudo_soil_temp,       & ! Pseudo-soil-temperature [Celsius]: a weighted running mean of the land ...
                                  ! surface temperature, to simulate soil temperature. Updated every time step ..
                                  ! by the subroutine "pseudo_soil_temp()". (stream)
      & day_temp_sum,           & ! sum of air temperatures since midnight (needed to compute mean air
                                  ! temperature) [Celsius * 1] (stream)
      & day_temp_min,           & ! minimum of air temperatures since midnight (needed to compute growing
                                  ! degree days GDD) [Celsius * 1] (stream)
      & day_temp_max,           & ! maximum of air temperatures since midnight (needed to compute growing
                                  ! degree days GDD) [Celsius * 1] (stream)
      & N_pseudo_soil_temp,     & ! Normalization for computing the pseudo soil temperature
      & F_pseudo_soil_temp        ! exponential factor used for updating the pseudo soil temperature
    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & t_lwtr,               & !< Lake surface temperature (water)                                [K]
      & t_lice,               & !< Lake surface temperature (ice)                                  [K]
!!$      & t_lwtr_old,            & !< Lake surface temperature (water) at previous timestep           [K]
!!$      & t_lice_old,            & !< Lake surface temperature (ice) at previous timestep             [K]
      & s_lwtr,               & !< Surface dry static energy over lake water                       [m2 s-2]
      & qsat_lwtr,            & !< Saturation specific humidity at surface over lake water        [kg kg-1]
      & s_lice,               & !< Surface dry static energy over lake ice                         [m2 s-2]
      & qsat_lice,            & !< Saturation specific humidity at surface over lake ice           [kg kg-1]
      & latent_hflx_wtr,      &
      & sensible_hflx_wtr,    &
      & latent_hflx_ice,      &
      & sensible_hflx_ice,    &
      & depth_lice,           & !< Thickness of ice on lake                                        [m]
      & fract_lice,           & !< Fraction of lake ice                                            []
      & hflx_form_ice,        & !< Heat flux for ice formation                                     [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux for ice on lakes
      & hflx_melt_lice,       & !< Heat flux from melting ice on lakes
      & hflx_melt_snow          !< Heat flux from melting snow on lake ice

  CONTAINS
    PROCEDURE :: Init => Init_seb_memory
  END TYPE t_seb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_memory_class'

CONTAINS

  SUBROUTINE Init_seb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var, BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables !, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_seb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_memory'
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    !$acc enter data pcreate(mem)
    dsl4jsb_Get_config(SEB_)


    !$acc enter data pcreate(mem%t,mem%t_old,mem%t_unfilt,mem%heat_cap,     &
    !$acc mem%heat_cap_old,mem%latent_hflx,mem%sensible_hflx,mem%forc_hflx)
    CALL jsb_add_var( mem, 't', mem%t,                                               &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 't_old', mem%t_old,                                       &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 't_unfilt', mem%t_unfilt,                                 &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_unfilt', 'K', 'unfiltered surface temperature'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 'heat_cap', mem%heat_cap,                                 &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap', 'J m-2 K-1)', 'surface layer heat capacity'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'heat_cap_old', mem%heat_cap_old,                         &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap_old', 'J m-2 K-1)', 'srf layer heat cap. prev timestep'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'latent_hflx', mem%latent_hflx,                           &
      & hgrid, surface,                                                              &
      & t_cf('latent_hflx', 'W m-2', 'latent heat flux density at surface'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL jsb_add_var( mem, 'sensible_hflx', mem%sensible_hflx,                       &
      & hgrid, surface,                                                              &
      & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL jsb_add_var( mem, 'forc_hflx', mem%forc_hflx,                               &
      & hgrid, surface,                                                              &
      & t_cf('forc_hflx', 'W m-2)', 'additional heat flux at surface'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    ! Additional variables for LAND lct
    IF (One_of(LAND_TYPE, lct_ids(:)) > 0 ) THEN

      !$acc enter data pcreate(mem%s_star,mem%qsat_star,mem%t_unfilt_old, &
      !$acc mem%t_filt,mem%t_eff,mem%latent_hflx_lnd,mem%sensible_hflx_lnd, &
      !$acc mem%previous_day_temp_mean,mem%previous_day_temp_min, &
      !$acc mem%previous_day_temp_max,mem%pseudo_soil_temp,mem%day_temp_sum, &
      !$acc mem%day_temp_min,mem%day_temp_max,mem%N_pseudo_soil_temp,mem%F_pseudo_soil_temp)

      CALL jsb_add_var( mem, 's_star', mem%s_star,                                           &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & output_level=BASIC, &
        & initval_r=2.9E5_wp )

      CALL jsb_add_var( mem, 'qsat_star', mem%qsat_star,                                     &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & output_level=BASIC, &
        & initval_r=0.0075_wp )

      CALL jsb_add_var( mem, 't_unfilt_old', mem%t_unfilt_old,                               &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_unfilt_old', 'K', 'unf. surface temperature at previous timestep'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_filt', mem%t_filt,                                           &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_filt', 'K', 'filtered surface temperature'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_eff', mem%t_eff,                                             &
        & hgrid, surface,                                                                    &
        & t_cf('surface_temperature_eff', 'K', 'effective surface temperature'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 'latent_hflx_lnd', mem%latent_hflx_lnd,                         &
        & hgrid, surface,                                                                    &
        & t_cf('latent_hflx_lnd', 'W m-2', 'latent heat flux density at surface (lnd)'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'sensible_hflx_lnd', mem%sensible_hflx_lnd,                     &
        & hgrid, surface,                                                                    &
        & t_cf('sensible_hflx_lnd', 'W m-2', 'sensible heat flux at surface (lnd)'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'previous_day_temp_mean', mem%previous_day_temp_mean, &
        & hgrid, surface,                                                          &
        & t_cf('previous_day_temp_mean', '', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_wp )

      CALL jsb_add_var( mem, 'previous_day_temp_min', mem%previous_day_temp_min, &
        & hgrid, surface,                                                        &
        & t_cf('previous_day_temp_min', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=5.0_wp )

      CALL jsb_add_var( mem, 'previous_day_temp_max', mem%previous_day_temp_max, &
        & hgrid, surface,                                                        &
        & t_cf('albedo_nir_previous_day_temp_maxsnow', '', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=15.0_wp )

      CALL jsb_add_var( mem, 'pseudo_soil_temp', mem%pseudo_soil_temp,           &
        & hgrid, surface,                                                        &
        & t_cf('pseudo_soil_temp', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_wp )

      CALL jsb_add_var( mem,'day_temp_sum', mem%day_temp_sum,                    &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_sum', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp )

      CALL jsb_add_var( mem,'day_temp_min', mem%day_temp_min,                    &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_min', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp )

      CALL jsb_add_var( mem,'day_temp_max', mem%day_temp_max,                    &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_max', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._wp)

      CALL jsb_add_var( mem,'N_pseudo_soil_temp', mem%N_pseudo_soil_temp,        &
        & hgrid, surface,                                                        &
        & t_cf('N_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._wp )

      CALL jsb_add_var( mem,'F_pseudo_soil_temp', mem%F_pseudo_soil_temp,        &
        & hgrid, surface,                                                        &
        & t_cf('F_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._wp )

    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      !$acc enter data pcreate(mem%t_lwtr,mem%s_lwtr,mem%qsat_lwtr, &
      !$acc mem%latent_hflx_wtr,mem%sensible_hflx_wtr,mem%fract_lice)

      CALL jsb_add_var( mem, 't_lwtr', mem%t_lwtr,                                     &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 's_lwtr', mem%s_lwtr,                                             &
        & hgrid, surface,                                                                      &
        & t_cf('sfc_dry_static_energy_wtr', 'm2 s-2', 'surface dry static energy (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.FALSE.,                                                                    &
        & initval_r=2.9E5_wp )

      CALL jsb_add_var( mem, 'qsat_lwtr', mem%qsat_lwtr,                                         &
        & hgrid, surface,                                                                        &
        & t_cf('sfc_qsat_lwtr', 'm2 s-2', 'surface specific humidity at saturation (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE.,                                                                      &
        & initval_r=0.0075_wp )

      CALL jsb_add_var( mem, 'latent_hflx_wtr', mem%latent_hflx_wtr,                   &
        & hgrid, surface,                                                              &
        & t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over lake water'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'sensible_hflx_wtr', mem%sensible_hflx_wtr,               &
        & hgrid, surface,                                                              &
        & t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over lake water'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Fraction of ice on lake
      CALL jsb_add_var( mem, 'fract_lice', mem%fract_lice,                             &
        & hgrid, surface,                                                              &
        & t_cf('fract_lice', '', 'fraction of ice on lake'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )
    
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN

        !$acc enter data pcreate(mem%t_lice, & 
        !$acc mem%s_lice,mem%qsat_lice,                                 &
        !$acc mem%latent_hflx_ice,mem%sensible_hflx_ice,mem%depth_lice, &
        !$acc mem%hflx_form_ice,mem%hflx_cond_ice,                      &
        !$acc mem%hflx_melt_lice,mem%hflx_melt_snow)
      
        CALL jsb_add_var( mem, 't_lice', mem%t_lice,                                     &
          & hgrid, surface,                                                              &
          & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=tmelt )
    
        CALL jsb_add_var( mem, 's_lice', mem%s_lice,                                             &
          & hgrid, surface,                                                                      &
          & t_cf('sfc_dry_static_energy_ice', 'm2 s-2', 'surface dry static energy (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
          & prefix, suffix,                                                                      &
          & lrestart=.FALSE.,                                                                    &
          & initval_r=2.9E5_wp )
    
        CALL jsb_add_var( mem, 'qsat_lice', mem%qsat_lice,                                         &
          & hgrid, surface,                                                                        &
          & t_cf('sfc_qsat_lice', 'm2 s-2', 'surface specific humidity at saturation (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix,                                                                        &
          & lrestart=.FALSE.,                                                                      &
          & initval_r=0.0075_wp )
    
        CALL jsb_add_var( mem, 'latent_hflx_ice', mem%latent_hflx_ice,                   &
          & hgrid, surface,                                                              &
          & t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL jsb_add_var( mem, 'sensible_hflx_ice', mem%sensible_hflx_ice,               &
          & hgrid, surface,                                                              &
          & t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over lake ice'),      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL jsb_add_var( mem, 'depth_lice', mem%depth_lice,                             &
          & hgrid, surface,                                                              &
          & t_cf('ice_thickness', 'm', 'thickness of ice on lake'),                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp )
    
        CALL jsb_add_var( mem, 'hflx_form_ice', mem%hflx_form_ice,                       &
          & hgrid, surface,                                                              &
          & t_cf('hflx_form_ice', 'W m-2)', 'heat flux for ice formation on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL jsb_add_var( mem, 'hflx_cond_ice', mem%hflx_cond_ice,                       &
          & hgrid, surface,                                                              &
          & t_cf('hflx_cond_ice', 'W m-2)', 'conductive heat flux on lake ice'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL jsb_add_var( mem, 'hflx_melt_lice', mem%hflx_melt_lice,                     &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_lice', 'W m-2)', 'heat flux from ice melt on lakes'),        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
    
        CALL jsb_add_var( mem, 'hflx_melt_snow', mem%hflx_melt_snow,                     &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_snow', 'W m-2)', 'heat_flux from snow melt on lakes'),       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        END IF

    END IF

  END SUBROUTINE Init_seb_memory

#endif
END MODULE mo_seb_memory_class
