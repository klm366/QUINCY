!> Contains the memory class for the surface energy balance process.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Memory definition of the surface energy balance process
!>
!> The module defines the data structure [[t_seb_memory]] that contains the memory of the surface
!> energy balance (SEB_) process. Each tile of the SEB_ process has its own version of this memory.
!> Which variables are added on a tile depends on its land cover types (land, lake), on
!> `use_tmx` and on the SEB_ configuration ([[mo_seb_config_class:t_seb_config]]), see
!> [[Init_seb_memory]].
!>
!> Sign convention: all energy fluxes are positive downward, i.e. into the surface.
!>
MODULE mo_seb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: lp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d
  USE mo_jsb_physical_constants, ONLY: tmelt
  USE mo_jsb_control,            ONLY: jsbach_runs_standalone

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50   !< Maximum number of 'Add_var' calls in this module

  !> Memory of the SEB_ process.
  !>
  !> _Note for the ford documentation_: The table of components starts with the bookkeeping components
  !> inherited from [[t_jsb_memory]] (ids, owner, links to parent and children, list of registered
  !> variables, ...); the SEB_ variables follow.
  TYPE, EXTENDS(t_jsb_memory) :: t_seb_memory
    !
    ! Common variables
    TYPE(t_jsb_var_real2d) ::   &
      & t,                      & !< Surface temperature; on land the skin temperature if
                                  !< `l_skin_temp`, on lakes the ice-fraction weighted mean of
                                  !< water and ice temperatures                                    [K]
      & t_old,                  & !< Surface temperature at previous timestep                      [K]
      & t_unfilt,               & !< Surface temperature, unfiltered (before Asselin)              [K]
      & t_rad4,                 & !< Surface radiative temperature **4 (TMX only)                  [K4]
      & t_eff4,                 & !< Effective surface temperature **4, linearised around the
                                  !< reference temperature of the energy balance solver;
                                  !< this is the surface temperature used in the atmospheric
                                  !< radiative heating                                             [K4]
      & qsat_star,              & !< Saturation specific humidity at surface (qsat^star, see manual) [kg kg-1]
                                  !< This is the qsat_srf that is used in the vertical diffusion
                                  !< scheme in the atmosphere
      & heat_cap,               & !< Heat capacity of surface layer (land: upper soil layer
                                  !< from SSE_, previous time step; lake: mixed layer / ice)       [J m-2 K-1]
      & latent_hflx,            & !< Latent heat flux                                              [W m-2]
      & sensible_hflx,          & !< Sensible heat flux                                            [W m-2]
      & forc_hflx,              & !< Additional heat flux at surface                               [W m-2]
                                  !< (additional to radiative fluxes and latent and sensible
                                  !<  heat fluxes, i.e. ground heat flux plus remaining latent
                                  !<  energy of phase changes from HYDRO_; zero on lakes)
      & le_phase_change,        & !< Latent energy of phase changes at the land surface within
                                  !< the time step; > 0: freezing (energy released),
                                  !< < 0: melting (energy consumed)                                [J m-2]
      & net_srf_rad_balance,    & !< Residual of all energy fluxes into the surface (net radiation,
                                  !< sensible, latent, forcing and phase change heat flux),
                                  !< i.e. without the change of heat storage                       [W m-2]
      & richardson                !< Richardson number (standalone model)                          []
    !
    ! Additional variables for LAND lct_type
    TYPE(t_jsb_var_real2d) ::   &
      & s_star,                 & !< Skin dry static energy (s^star, see manual)                   [m2 s-2]
                                  !< This is the s that is used in the vertical diffusion
                                  !< scheme in the atmosphere
      & t_unfilt_old,           & !< Skin temperature at previous timestep, unfiltered             [K]
                                  !< (before Asselin)
      & t_filt,                 & !< Skin temperature, filtered (after Asselin)                    [K]
      & t_srf,                  & !< Surface temperature, top soil or snow layer, unfiltered;
                                  !< equal to the skin temperature if not `l_skin_temp`            [K]
      & dt_sk_dt_srf,           & !< Derivative of skin temperature w.r.t. surface temperature     [K K-1]
      & latent_hflx_lnd,        & !< Latent heat flux, land part only (not aggregated with lakes)  [W m-2]
      & sensible_hflx_lnd,      & !< Sensible heat flux, land part only (not aggregated with lakes)[W m-2]
      & previous_day_temp_mean, & !< Mean air temperature of the previous day; updated at the
                                  !< first time step of each new day                               [degC]
      & previous_day_temp_min,  & !< Minimum air temperature of the previous day                   [degC]
      & previous_day_temp_max,  & !< Maximum air temperature of the previous day                   [degC]
      & pseudo_soil_temp,       & !< Pseudo soil temperature: weighted running mean of the air
                                  !< temperature to mimic soil temperature, updated every time
                                  !< step by [[mo_seb_land:calc_pseudo_soil_temp]]                 [degC]
      & day_temp_sum,           & !< Sum of air temperatures since midnight (needed to compute
                                  !< the mean day temperature)                                     [degC]
      & day_temp_min,           & !< Minimum of air temperatures since midnight                    [degC]
      & day_temp_max,           & !< Maximum of air temperatures since midnight                    [degC]
      & N_pseudo_soil_temp,     & !< Normalisation for computing the pseudo soil temperature       []
      & F_pseudo_soil_temp,     & !< Exponential factor for updating the pseudo soil temperature   []
      & skin_conductivity         !< Skin conductivity (bc file, if `l_skin_temp`)                 [W m-2 K-1]
    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & t_lwtr,               & !< Lake surface temperature (water)                                [K]
      & t_lice,               & !< Lake surface temperature (ice)                                  [K]
!!$      & t_lwtr_old,            & !< Lake surface temperature (water) at previous timestep       [K]
!!$      & t_lice_old,            & !< Lake surface temperature (ice) at previous timestep         [K]
      & s_lwtr,               & !< Surface dry static energy over lake water                       [m2 s-2]
      & qsat_lwtr,            & !< Saturation specific humidity at surface over lake water         [kg kg-1]
      & s_lice,               & !< Surface dry static energy over lake ice                         [m2 s-2]
      & qsat_lice,            & !< Saturation specific humidity at surface over lake ice           [kg kg-1]
      & latent_hflx_wtr,      & !< Latent heat flux over lake water                                [W m-2]
      & sensible_hflx_wtr,    & !< Sensible heat flux over lake water                              [W m-2]
      & latent_hflx_ice,      & !< Latent heat flux over lake ice                                  [W m-2]
      & sensible_hflx_ice,    & !< Sensible heat flux over lake ice                                [W m-2]
      & depth_lice,           & !< Thickness of ice on lake                                        [m]
      & fract_lice,           & !< Fraction of lake ice                                            []
      & hflx_form_ice,        & !< Heat flux for ice formation (< 0); kept for the next time step
                                !< if not sufficient to form a new slab of ice                     [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux through ice on lakes                       [W m-2]
      & hflx_melt_lice,       & !< Heat flux from melting ice on lakes                             [W m-2]
      & hflx_melt_snow          !< Heat flux from melting snow on lake ice                         [W m-2]

    ! Diagnostic global land mean variables e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) ::        &
      t_gmean                   !< Global land mean surface temperature (box tile)                 [K]

  CONTAINS
    PROCEDURE :: Init => Init_seb_memory
  END TYPE t_seb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_memory_class'

CONTAINS

  !> Add the SEB_ variables to the memory of a tile.
  !>
  !> Variables common to all tiles are always added; additional variables are added for tiles
  !> containing the land (or, with `use_tmx`, all tiles) and lake land cover types, for the
  !> standalone model (`richardson`), and for the box tile (`t_gmean`).
  !>
  SUBROUTINE Init_seb_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_seb_memory), INTENT(inout), TARGET :: mem !< SEB_ memory of the tile
    CHARACTER(len=*),    INTENT(in)    :: prefix      !< Prefix of variable names (process name)
    CHARACTER(len=*),    INTENT(in)    :: suffix      !< Suffix of variable names (tile name)
    INTEGER,             INTENT(in)    :: lct_ids(:)  !< Ids of the land cover types on the tile
    INTEGER,             INTENT(in)    :: lib_id      !< id of primary lct in lctlib
    INTEGER,             INTENT(in)    :: model_id    !< Model id

    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model   !< This instance of ICON-Land
    TYPE(t_jsb_grid),  POINTER :: hgrid   !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface !< Vertical grid (surface, one layer)
    INTEGER :: table                      !< GRIB table number
    LOGICAL :: lrestart_local             !< add var to restart file, model specific (quincy/jsbach)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_memory'
    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    dsl4jsb_Get_config(SEB_)

    ! set local lrestart switch to TRUE for the MODEL_QUINCY
    lrestart_local = .FALSE.
    IF (model%config%model_scheme == MODEL_QUINCY) THEN
      lrestart_local = .TRUE.
    END IF

    CALL mem%Add_var( 't', mem%t,                                                    &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=280.0_lp )

    IF (model%config%use_tmx) THEN
      CALL mem%Add_var( 't_rad**4', mem%t_rad4,                                      &
        & hgrid, surface,                                                            &
        & t_cf('sfc_temp_rad**4', 'K', 'surface radiative temperature **4'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=280.0_lp**4 )
    END IF

    CALL mem%Add_var( 't_old', mem%t_old,                                            &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart = lrestart_local,                                                   &
      & initval_r=280.0_lp )

    CALL mem%Add_var( 't_unfilt', mem%t_unfilt,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_unfilt', 'K', 'unfiltered surface temperature'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_lp )

    CALL mem%Add_var( 'le_phase_change', mem%le_phase_change,                        &
      & hgrid, surface,                                                              &
      & t_cf('le_phase_change', 'J m-2',                                             &
      &      'Latent energy of phase changes at land surface'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE., output_level=MEDIUM,                                       &
      & initval_r=0.0_lp )

    CALL mem%Add_var( 'heat_cap', mem%heat_cap,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap', 'J m-2 K-1', 'surface layer heat capacity'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=MEDIUM,                                                         &
      & initval_r=0.0_lp )

    CALL mem%Add_var( 'latent_hflx', mem%latent_hflx,                                &
      & hgrid, surface,                                                              &
      & t_cf('latent_hflx', 'W m-2', 'latent heat flux density at surface'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'sensible_hflx', mem%sensible_hflx,                            &
      & hgrid, surface,                                                              &
      & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'forc_hflx', mem%forc_hflx,                                    &
      & hgrid, surface,                                                              &
      & t_cf('forc_hflx', 'W m-2', 'additional heat flux at surface'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'net_srf_rad_balance', mem%net_srf_rad_balance,                &
      & hgrid, surface,                                                              &
      & t_cf('net_srf_rad_balance', 'W m-2',                                         &
      &      'Net surface radiation balance'),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_lp)

    IF (jsbach_runs_standalone()) THEN
      CALL mem%Add_var( 'richardson', mem%richardson,                                &
        & hgrid, surface,                                                            &
        & t_cf('richardson', '-', 'Richardson number'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.TRUE.,                                                           &
        & initval_r=0.0_lp )
    END IF

    IF (model%config%use_tmx .OR. One_of(LAND_TYPE, lct_ids(:)) > 0) THEN
      CALL mem%Add_var( 't_eff**4', mem%t_eff4,                                              &
        & hgrid, surface,                                                                    &
        & t_cf('surface_temperature_eff**4', 'K', 'effective surface temperature**4'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_lp**4 )
    END IF

    ! Additional variables for LAND lct
    IF ( model%config%use_tmx .OR. One_of(LAND_TYPE, lct_ids(:)) > 0) THEN

      CALL mem%Add_var( 's_star', mem%s_star,                                                &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & output_level=FULL,                                                                 &
        & initval_r=2.9E5_lp)

      CALL mem%Add_var( 'qsat_star', mem%qsat_star,                                          &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.TRUE.,                                                                   &
        & output_level=FULL,                                                                 &
        & initval_r=0.0075_lp)

      CALL mem%Add_var( 't_unfilt_old', mem%t_unfilt_old,                                    &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_unfilt_old', 'K', 'unf. surface temperature at previous timestep'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_lp )

      CALL mem%Add_var( 't_filt', mem%t_filt,                                                &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_filt', 'K', 'filtered surface temperature'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=280.0_lp )

      CALL mem%Add_var( 't_srf', mem%t_srf,                                                  &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp', 'K', 'soil or snow surface temperature'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=(    (.NOT. model%config%use_tmx .AND. .NOT. jsbach_runs_standalone())    &
        &           .OR. dsl4jsb_Config(SEB_)%l_skin_temp),                                  &
        & output_level = BASIC,                                                              &
        & initval_r=280.0_lp )

      CALL mem%Add_var( 'dt_sk_dt_srf', mem%dt_sk_dt_srf,                                    &
        & hgrid, surface,                                                                    &
        & t_cf('dt_sk_dt_srf', 'K', 'derivative of skin w.r.t. surface temperature'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=1._lp )

      CALL mem%Add_var( 'latent_hflx_lnd', mem%latent_hflx_lnd,                              &
        & hgrid, surface,                                                                    &
        & t_cf('latent_hflx_lnd', 'W m-2', 'latent heat flux density at surface (lnd)'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'sensible_hflx_lnd', mem%sensible_hflx_lnd,                          &
        & hgrid, surface,                                                                    &
        & t_cf('sensible_hflx_lnd', 'W m-2', 'sensible heat flux at surface (lnd)'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.FALSE.,                                                                  &
        & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'previous_day_temp_mean', mem%previous_day_temp_mean,      &
        & hgrid, surface,                                                          &
        & t_cf('previous_day_temp_mean', '', ''),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_lp )

      CALL mem%Add_var( 'previous_day_temp_min', mem%previous_day_temp_min,      &
        & hgrid, surface,                                                        &
        & t_cf('previous_day_temp_min', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=5.0_lp )

      CALL mem%Add_var( 'previous_day_temp_max', mem%previous_day_temp_max,      &
        & hgrid, surface,                                                        &
        & t_cf('albedo_nir_previous_day_temp_maxsnow', '', ''),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=15.0_lp )

      CALL mem%Add_var( 'pseudo_soil_temp', mem%pseudo_soil_temp,                &
        & hgrid, surface,                                                        &
        & t_cf('pseudo_soil_temp', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10.0_lp )

      CALL mem%Add_var('day_temp_sum', mem%day_temp_sum,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_sum', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._lp )

      CALL mem%Add_var('day_temp_min', mem%day_temp_min,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_min', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._lp )

      CALL mem%Add_var('day_temp_max', mem%day_temp_max,                         &
        & hgrid, surface,                                                        &
        & t_cf('day_temp_max', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0._lp)

      CALL mem%Add_var('N_pseudo_soil_temp', mem%N_pseudo_soil_temp,             &
        & hgrid, surface,                                                        &
        & t_cf('N_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._lp )

      CALL mem%Add_var('F_pseudo_soil_temp', mem%F_pseudo_soil_temp,             &
        & hgrid, surface,                                                        &
        & t_cf('F_pseudo_soil_temp', '', ''),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=10._lp )

      IF (dsl4jsb_Config(SEB_)%l_skin_temp) THEN
        CALL mem%Add_var( 'skin_conductivity', mem%skin_conductivity,              &
          & hgrid, surface,                                                        &
          & t_cf('skin_conductivity', 'W m-2 K-1', 'skin conductivity'),           &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
          & prefix, suffix,                                                        &
          & lrestart=.FALSE., isteptype=TSTEP_CONSTANT,                            &
          & initval_r=0._lp )
      END IF

    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 't_lwtr', mem%t_lwtr,                                          &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_lp )

      CALL mem%Add_var( 's_lwtr', mem%s_lwtr,                                                  &
        & hgrid, surface,                                                                      &
        & t_cf('sfc_dry_static_energy_wtr', 'm2 s-2', 'surface dry static energy (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.FALSE.,                                                                    &
        & initval_r=2.9E5_lp )

      CALL mem%Add_var( 'qsat_lwtr', mem%qsat_lwtr,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('sfc_qsat_lwtr', 'm2 s-2', 'surface specific humidity at saturation (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE.,                                                                      &
        & initval_r=0.0075_lp )

      CALL mem%Add_var( 'latent_hflx_wtr', mem%latent_hflx_wtr,                        &
        & hgrid, surface,                                                              &
        & t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over lake water'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'sensible_hflx_wtr', mem%sensible_hflx_wtr,                    &
        & hgrid, surface,                                                              &
        & t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over lake water'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

      ! Fraction of ice on lake
      CALL mem%Add_var( 'fract_lice', mem%fract_lice,                                  &
        & hgrid, surface,                                                              &
        & t_cf('fract_lice', '', 'fraction of ice on lake'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_lp )

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN

        CALL mem%Add_var( 't_lice', mem%t_lice,                                          &
          & hgrid, surface,                                                              &
          & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=REAL(tmelt,KIND=lp) )

        CALL mem%Add_var( 's_lice', mem%s_lice,                                                  &
          & hgrid, surface,                                                                      &
          & t_cf('sfc_dry_static_energy_ice', 'm2 s-2', 'surface dry static energy (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
          & prefix, suffix,                                                                      &
          & lrestart=.FALSE.,                                                                    &
          & initval_r=2.9E5_lp )

        CALL mem%Add_var( 'qsat_lice', mem%qsat_lice,                                              &
          & hgrid, surface,                                                                        &
          & t_cf('sfc_qsat_lice', 'm2 s-2', 'surface specific humidity at saturation (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix,                                                                        &
          & lrestart=.FALSE.,                                                                      &
          & initval_r=0.0075_lp )

        CALL mem%Add_var( 'latent_hflx_ice', mem%latent_hflx_ice,                        &
          & hgrid, surface,                                                              &
          & t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'sensible_hflx_ice', mem%sensible_hflx_ice,                    &
          & hgrid, surface,                                                              &
          & t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over lake ice'),      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'depth_lice', mem%depth_lice,                                  &
          & hgrid, surface,                                                              &
          & t_cf('ice_thickness', 'm', 'thickness of ice on lake'),                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_lp )

        CALL mem%Add_var( 'hflx_form_ice', mem%hflx_form_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_form_ice', 'W m-2)', 'heat flux for ice formation on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_cond_ice', mem%hflx_cond_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_cond_ice', 'W m-2)', 'conductive heat flux on lake ice'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_melt_lice', mem%hflx_melt_lice,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_lice', 'W m-2)', 'heat flux from ice melt on lakes'),        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_melt_snow', mem%hflx_melt_snow,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_snow', 'W m-2)', 'heat_flux from snow melt on lakes'),       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_lp, l_aggregate_all=.TRUE. )

        END IF

    END IF

    ! Diagnostic global carbon sums for experiment monitoring
    !
    IF (TRIM(suffix) == 'box') THEN

      CALL mem%Add_var( 't_gmean', mem%t_gmean,                                 &
        & hgrid, surface,                                                       &
        & t_cf('t_gmean', 'K', 'Global land mean surface temperature'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                       &
        & output_groups=['mon'], lrestart=.FALSE., initval_r=280.0_lp )

    END IF

  END SUBROUTINE Init_seb_memory

#endif
END MODULE mo_seb_memory_class
