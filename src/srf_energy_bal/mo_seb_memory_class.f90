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
MODULE mo_seb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_seb_memory
    !
    ! Common variables
    TYPE(t_jsb_var_real2d) :: &
      & t,                    & !< Surface temperature                                             [K]
      & t_old,                & !< Surface temperature at previous timestep                        [K]
      & heat_cap,             & !< Heat capacity of surface layer                                  [J m-2 K-1]
      & heat_cap_old,         & !< Heat capacity of surface layer at previous timestep             [J m-2 K-1]
      & latent_hflx,          & !< Latent heat flux                                                [W m-2]
      & sensible_hflx,        & !< Sensible heat flux                                              [W m-2]
      & forc_hflx               !< Additional heat flux at surface                                 [W m-2]
                                !! (additional to radiative fluxes and latent and sensible
                                !!  heat fluxes, e.g. ground heat flux)
    !
    ! Additional variables for LAND lct_type
    TYPE(t_jsb_var_real2d) :: &
      & s_star,               & !< Surface dry static energy (s^star, see manual)                  [m2 s-2]
                                !! This is the s that is used in the vertical diffusion
                                !! scheme in the atmosphere
      & qsat_star,            & !< Saturation specific humidity at surface (qsat^star, see manual) [kg kg-1]
                                !! This is the qsat_srf that is used in the vertical diffusion
                                !! scheme in the atmosphere
      & t_unfilt,             & !< Surface temperature, unfiltered (before Asselin)                [K]
      & t_unfilt_old,         & !< Surface temperature at previous timestep, unfiltered            [K]
                                !! (before Asselin)
      & t_eff,                & !< Effective surface temperature (see manual)                      [K]
                                !! This is the surface temperature used in the atmospheric
                                !! radiative heating
      & fact_q_air,           & !< Factor for lowest atm. level in formulation of surface fluxes
                                !! of specific humidity and dry static energy to account for
                                !! evapotranspiration                                              []
      & fact_qsat_srf,        & !< Similar factor for surface values of specific humidity
                                !! and dry static energy
      & fact_qsat_trans_srf     !< Similar factor for transpiration
    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & t_wtr,                & !< Lake surface temperature (water)                                [K]
      & t_ice,                & !< Lake surface temperature (ice)                                  [K]
!!$      & t_wtr_old,            & !< Lake surface temperature (water) at previous timestep           [K]
!!$      & t_ice_old,            & !< Lake surface temperature (ice) at previous timestep             [K]
      & depth_lice,           & !< Thickness of ice on lake                                        [m]
      & fract_lice,           & !< Fraction of lake ice                                            []
      & hflx_form_ice,        & !< Heat flux for ice formation                                     [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux for ice on lakes
      & hflx_melt_ice,        & !< Heat flux from melting ice on lakes
      & hflx_melt_snow          !< Heat flux from melting snow on lake ice

  CONTAINS
    PROCEDURE :: Init => Init_seb_memory
  END TYPE t_seb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_memory_class'

CONTAINS

  SUBROUTINE Init_seb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_seb_memory), INTENT(inout) :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_seb_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    IF (.NOT. ASSOCIATED(mem%t%ptr)) THEN
      ! Common variables
      CALL jsb_add_var( mem, 't', mem%t,                                               &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_old', mem%t_old,                                       &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 'heat_cap', mem%heat_cap,                                 &
        & hgrid, surface,                                                              &
        & t_cf('sfc_heat_cap', 'J m-2 K-1)', 'surface layer heat capacity'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'heat_cap_old', mem%heat_cap_old,                         &
        & hgrid, surface,                                                              &
        & t_cf('sfc_heat_cap_old', 'J m-2 K-1)', 'srf layer heat cap. prev timestep'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'latent_hflx', mem%latent_hflx,                           &
        & hgrid, surface,                                                              &
        & t_cf('latent_hflx', 'W m-2', 'latent heat flux density at surface'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sensible_hflx', mem%sensible_hflx,                       &
        & hgrid, surface,                                                              &
        & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 'forc_hflx', mem%forc_hflx,                               &
        & hgrid, surface,                                                              &
        & t_cf('forc_hflx', 'W m-2)', 'additional heat flux at surface'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )
    END IF

    ! Additional variables for LAND lct
    IF (One_of(LAND_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%s_star%ptr)) THEN

      CALL jsb_add_var( mem, 's_sfc', mem%s_star,                                            &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=2.9E5_wp )

      CALL jsb_add_var( mem, 'qsat_sfc', mem%qsat_star,                                      &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=0.0075_wp )

      CALL jsb_add_var( mem, 't_unfilt', mem%t_unfilt,                                       &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_unfilt', 'K', 'unfiltered surface temperature'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_unfilt_old', mem%t_unfilt_old,                               &
        & hgrid, surface,                                                                    &
        & t_cf('sfc_temp_unfilt_old', 'K', 'unf. surface temperature at previous timestep'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_sfc_eff', mem%t_eff,                                         &
        & hgrid, surface,                                                                    &
        & t_cf('surface_temperature_eff', 'K', 'effective surface temperature'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 'fact_q_air', mem%fact_q_air,                                   &
        & hgrid, surface,                                                                    &
        & t_cf('fact_q_air', '', ''),                                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'fact_qsat_srf', mem%fact_qsat_srf,                             &
        & hgrid, surface,                                                                    &
        & t_cf('fact_qsat_srf', '', ''),                                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'fact_qsat_trans_srf', mem%fact_qsat_trans_srf,                 &
        & hgrid, surface,                                                                    &
        & t_cf('fact_qsat_trans_srf', '', ''),                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & initval_r=0.0_wp )

    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%t_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 't_wtr', mem%t_wtr,                                       &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL jsb_add_var( mem, 't_ice', mem%t_ice,                                       &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=tmelt )

      CALL jsb_add_var( mem, 'depth_lice', mem%depth_lice,                             &
        & hgrid, surface,                                                              &
        & t_cf('ice_thickness', 'm', 'thickness of ice on lake'),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      ! Fraction of ice on lake, rel. to lake fraction
      CALL jsb_add_var( mem, 'fract_lice', mem%fract_lice,                             &
        & hgrid, surface,                                                              &
        & t_cf('fract_lice', '', 'fraction of ice on lake'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'hflx_form_ice', mem%hflx_form_ice,                       &
        & hgrid, surface,                                                              &
        & t_cf('hflx_form_ice', 'W m-2)', 'heat flux for ice formation on lake ice'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'hflx_cond_ice', mem%hflx_cond_ice,                       &
        & hgrid, surface,                                                              &
        & t_cf('hflx_cond_ice', 'W m-2)', 'conductive heat flux on lake ice'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'hflx_melt_ice', mem%hflx_melt_ice,                       &
        & hgrid, surface,                                                              &
        & t_cf('hflx_melt_ice', 'W m-2)', 'heat flux from ice melt on lakes'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'hflx_melt_snow', mem%hflx_melt_snow,                     &
        & hgrid, surface,                                                              &
        & t_cf('hflx_melt_snow', 'W m-2)', 'heat_flux from snow melt on lakes'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

    END IF

  END SUBROUTINE Init_seb_memory

#endif
END MODULE mo_seb_memory_class
