!>
!! @brief Contains the memory class for the surface.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_srf_memory, max_no_of_vars
!!$            copy_srf_mem, add_srf_mem

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_srf_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      land         ,    & !< Land fraction: 1: land/lake/glacier, 0: sea
      coast_land   ,    & !< mask for coastal cells of the land
      coast_ocean  ,    & !< mask for coastal cells of the ocean
      ocean        ,    & !< mask for ocean cells, 1: ocean
      glacier      ,    & !< Glacier fraction
      lake         ,    & !< Lake fraction
      oro_stddev   ,    & !< Standard deviation of orography
      rough_m_srf  ,    & !< Surface roughness length for momentum
      rough_h_srf  ,    & !< Surface roughness length for heat
      albedo_srf   ,    & !< (Background) albedo of tile
      alb_vis_srf  ,    & !< Albedo in the visible range
      alb_nir_srf         !< Albedo in the NIR range
    !
    ! Fractions
    TYPE(t_jsb_var_real2d) :: &
      fract,    &  !< Cover fraction of tile (relative to parent tile)
      fract_max    !< Maximum potential cover fraction for tile (rel. to parent)
    !
    ! Energy state
    TYPE(t_jsb_var_real2d) :: &
      t_srf            , &
      t_srf_old        , &
      t_srf_unfilt     , &
      t_srf_unfilt_old , &
      t_eff_srf        , &
      s_srf            , &  !< Dry static energy at surface (new timestep)
!!$      s_srf_old        , &  !< Dry static energy at surface (old timestep)
      qsat_srf         , &
!!$      qsat_srf_old     , &
      tte_corr         , &  !< Correction term of tte for snow melt
      fact_q_air       , &  !< cair in echam6
      fact_qsat_srf    , &  !< csat in echam6
      fact_qsat_trans_srf   !< csat_transpiration in echam6/jsbach
    !
    ! Water State
    TYPE(t_jsb_var_real2d) :: &
      wfract_srf    ,   & !< Wet (skin reservoir) fraction of tile []
      wsr_srf       ,   & !< Water content in skin reservoir [m]
      sfract_srf    ,   & !< Snow fraction of tile []
      wsn_srf       ,   & !< Water content of snow reservoir [m water equivalent]
      wsn_glacier         !< Water content of snow on glacier [m water equivalent]
    !
    ! Energy fluxes
    TYPE(t_jsb_var_real2d) :: &
      latent_hflx      ,   & !< latent heat flux
      sensible_hflx    ,   & !< sensible heat flux
      sw_srf_net       ,   & !< net shortwave radiation at surface
      swvis_srf_net    ,   & !< net shortwave radiation at surface in visible range
      swnir_srf_net    ,   & !< net shortwave radiation at surface in NIR range
      lw_srf_net       ,   & !< net longwave radiation at surface
      rad_srf_net            !< net radiation at surface
    !
    ! Water fluxes
    TYPE(t_jsb_var_real2d) :: &
      evapopot         ,   & !< Potential evaporation
      evapotrans       ,   & !< Evapotranspiration
      trans                  !< Transpiration
  CONTAINS
    PROCEDURE :: Init => Init_srf_memory
  END TYPE t_srf_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
!!$  INTERFACE copy_srf_mem
!!$    MODULE PROCEDURE copy_t_srf_energy_state
!!$  END INTERFACE copy_srf_mem
!!$
!!$  INTERFACE add_srf_mem
!!$    MODULE PROCEDURE add_t_srf_energy_state
!!$  END INTERFACE add_srf_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_memory_class'

CONTAINS

  SUBROUTINE Init_srf_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, ZAXIS_SURFACE, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_srf_memory), INTENT(inout) :: mem
    CHARACTER(len=*),   INTENT(in)     :: prefix
    CHARACTER(len=*),   INTENT(in)     :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_srf_memory'

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

      CALL jsb_add_var( mem, prefix//'land_frac'//suffix, mem%land,                     &
        & hgrid, surface,                                                               &
        & t_cf('land_mask', '', 'land fraction', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'ocean_frac'//suffix, mem%ocean,                   &
        & hgrid, surface,                                                               &
        & t_cf('ocean_mask', '', 'ocean fraction', nc_bits),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'lake_frac'//suffix, mem%lake,                     &
        & hgrid, surface,                                                               &
        & t_cf('lake_mask', '', 'lake fraction', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'coast_land'//suffix, mem%coast_land,     &
        & hgrid, surface,                                                               &
        & t_cf('coast_land', '', 'coastal grid cells of the land', nc_bits),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'coast_ocean'//suffix, mem%coast_ocean,   &
        & hgrid, surface,                                                               &
        & t_cf('coast_ocean', '', 'coastal grid cells of the ocean ', nc_bits),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'glacier'//suffix, mem%glacier,           &
        & hgrid, surface,                                                               &
        & t_cf('glacier_mask', '', 'glacier fraction', nc_bits),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'oro_stddev'//suffix, mem%oro_stddev,                         &
        & hgrid, surface,                                                           &
        & t_cf('orography_standard_deviation', '', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'rough_m_srf'//suffix, mem%rough_m_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('roughness_length_momentum', 'm', '', nc_bits),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'rough_h_srf'//suffix, mem%rough_h_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('roughness_length_heat', 'm', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'albedo_srf'//suffix, mem%albedo_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_background', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'alb_vis_srf'//suffix, mem%alb_vis_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_visible', '', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'alb_nir_srf'//suffix, mem%alb_nir_srf,                       &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_nir', '', '', nc_bits),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_srf'//suffix, mem%t_srf,                                   &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature', 'K', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_srf_old'//suffix, mem%t_srf_old,                           &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_old', 'K', '', nc_bits),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_srf_unfilt'//suffix, mem%t_srf_unfilt,                     &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_unfiltered', 'K', '', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_srf_unfilt_old'//suffix, mem%t_srf_unfilt_old,             &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_unfiltered_old', 'K', '', nc_bits),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_eff_srf'//suffix, mem%t_eff_srf,                           &
        & hgrid, surface,                                                           &
        & t_cf('surface_temperature_eff', 'K', '', nc_bits),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( mem, prefix//'t_eff_srf'//suffix, mem%t_eff_srf,                           &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('surface_temperature_eff', 'K', '', nc_bits),                                 &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'s_srf'//suffix, mem%s_srf,                                   &
        & hgrid, surface,                                                           &
        & t_cf('dry_static_energy', '', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( mem, prefix//'s_srf_old'//suffix, mem%s_srf_old,                           &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('dry_static_energy_old', '', '', nc_bits),                                    &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'qsat_srf'//suffix, mem%qsat_srf,                             &
        & hgrid, surface,                                                           &
        & t_cf('sat_specific_humidity', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

!!$      CALL jsb_add_var( mem, prefix//'qsat_srf_old'//suffix, mem%qsat_srf_old,                     &
!!$        & hgrid, surface,                                                           &
!!$        & t_cf('sat_specific_humidity_old', '', '', nc_bits),                                &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!!$        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'tte_corr'//suffix, mem%tte_corr,                             &
        & hgrid, surface,                                                           &
        & t_cf('tte_correction_snowmelt', '', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'fact_q_air'//suffix, mem%fact_q_air,                         &
        & hgrid, surface,                                                           &
        & t_cf('factor_air_moisture', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'fact_qsat_srf'//suffix, mem%fact_qsat_srf,                   &
        & hgrid, surface,                                                           &
        & t_cf('factor_qsat_srf', '', '', nc_bits),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'fact_qsat_trans_srf'//suffix, mem%fact_qsat_trans_srf,       &
        & hgrid, surface,                                                           &
        & t_cf('factor_qsat_trans_srf', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'sw_srf_net'//suffix, mem%sw_srf_net,                        &
        & hgrid, surface,                                                                                 &
        & t_cf('sw_srf_net', '', 'Net surface shortwave radiation', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'swvis_srf_net'//suffix, mem%swvis_srf_net,                  &
        & hgrid, surface,                                                                                 &
        & t_cf('swvis_srf_net', '', 'Net surface shortwave radiation in VIS', nc_bits),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'swnir_srf_net'//suffix, mem%swnir_srf_net,                  &
        & hgrid, surface,                                                                                 &
        & t_cf('swnir_srf_net', '', 'Net surface shortwave radiation in NIR', nc_bits),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'lw_srf_net'//suffix, mem%lw_srf_net,                        &
        & hgrid, surface,                                                                                 &
        & t_cf('lw_srf_net', '', 'Net surface longwave radiation', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'rad_srf_net'//suffix, mem%rad_srf_net,                      &
        & hgrid, surface,                                                                                 &
        & t_cf('rad_srf_net', '', 'Net surface radiation', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'latent_hflx'//suffix, mem%latent_hflx,                      &
        & hgrid, surface,                                                                                 &
        & t_cf('latent_hflx', '', 'Latent heat flux', nc_bits),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'sensible_hflx'//suffix, mem%sensible_hflx,                  &
        & hgrid, surface,                                                                                 &
        & t_cf('sensible_hflx', '', 'Sensible heat flux', nc_bits),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wfract_srf'//suffix, mem%wfract_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_wet_fraction', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsr_srf'//suffix, mem%wsr_srf,                               &
        & hgrid, surface,                                                           &
        & t_cf('surface_skin_reservoir', 'm', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'sfract_srf'//suffix, mem%sfract_srf,                         &
        & hgrid, surface,                                                           &
        & t_cf('surface_snow_fraction', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsn_srf'//suffix, mem%wsn_srf,                               &
        & hgrid, surface,                                                           &
        & t_cf('surface_snow_water_equivalent', 'm', '', nc_bits),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsn_glacier'//suffix, mem%wsn_glacier,                       &
        & hgrid, surface,                                                           &
        & t_cf('glacier_snow_water_equivalent', 'm', '', nc_bits),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'evapopot'//suffix, mem%evapopot,                                     &
        & hgrid, surface,                                                                  &
        & t_cf('surface_evaporation_potential', '', 'Potential evaporation from surface', nc_bits), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'evapotrans'//suffix, mem%evapotrans,                                 &
        & hgrid, surface,                                                                  &
        & t_cf('surface_evapotranspiration', '', 'Evapotranspiration from surface', nc_bits),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'trans'//suffix, mem%trans,                                           &
        & hgrid, surface,                                                                  &
        & t_cf('surface_transpiration', '', 'Transpiration from surface', nc_bits),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & initval_r=0.0_wp )

  END SUBROUTINE Init_srf_memory

!!$  SUBROUTINE copy_t_srf_energy_state(in, out, iblk, ics, ice)
!!$
!!$    TYPE(t_srf_energy_state), INTENT(in)    :: in
!!$    TYPE(t_srf_energy_state), INTENT(inout) :: out
!!$    INTEGER,                  INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_srf         (ics:ice,    iblk) = in%t_srf         (ics:ice,    iblk)
!!$    out%t_srf_old     (ics:ice,    iblk) = in%t_srf_old     (ics:ice,    iblk)
!!$    out%t_srf_unfilt  (ics:ice,    iblk) = in%t_srf_unfilt  (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in%t_eff_srf     (ics:ice,    iblk)
!!$    out%s_srf         (ics:ice,    iblk) = in%s_srf         (ics:ice,    iblk)
!!$    !out%s_srf_old     (ics:ice,    iblk) = in%s_srf_old     (ics:ice,    iblk)
!!$    out%qsat_srf      (ics:ice,    iblk) = in%qsat_srf      (ics:ice,    iblk)
!!$    !out%qsat_srf_old  (ics:ice,    iblk) = in%qsat_srf_old  (ics:ice,    iblk)
!!$    out%tte_corr      (ics:ice,    iblk) = in%tte_corr      (ics:ice,    iblk)
!!$    out%fact_q_air    (ics:ice,    iblk) = in%fact_q_air    (ics:ice,    iblk)
!!$    out%fact_qsat_srf (ics:ice,    iblk) = in%fact_qsat_srf (ics:ice,    iblk)
!!$    out%fact_qsat_trans_srf (ics:ice,    iblk) = in%fact_qsat_trans_srf (ics:ice,    iblk)
!!$
!!$  END SUBROUTINE copy_t_srf_energy_state
!!$
!!$  SUBROUTINE add_t_srf_energy_state(in1, in2, out, iblk, ics, ice)
!!$
!!$    TYPE(t_srf_energy_state), INTENT(in)    :: in1, in2
!!$    TYPE(t_srf_energy_state), INTENT(inout) :: out
!!$    INTEGER,                  INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_srf         (ics:ice,    iblk) = in1%t_srf         (ics:ice,    iblk) + in2%t_srf         (ics:ice,    iblk)
!!$    out%t_srf_old     (ics:ice,    iblk) = in1%t_srf_old     (ics:ice,    iblk) + in2%t_srf_old     (ics:ice,    iblk)
!!$    out%t_srf_unfilt  (ics:ice,    iblk) = in1%t_srf_unfilt  (ics:ice,    iblk) + in2%t_srf_unfilt  (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in1%t_eff_srf     (ics:ice,    iblk) + in2%t_eff_srf     (ics:ice,    iblk)
!!$    out%s_srf         (ics:ice,    iblk) = in1%s_srf         (ics:ice,    iblk) + in2%s_srf         (ics:ice,    iblk)
!!$    !out%s_srf_old     (ics:ice,    iblk) = in1%s_srf_old     (ics:ice,    iblk) + in2%s_srf_old     (ics:ice,    iblk)
!!$    out%qsat_srf      (ics:ice,    iblk) = in1%qsat_srf      (ics:ice,    iblk) + in2%qsat_srf      (ics:ice,    iblk)
!!$    !out%qsat_srf_old  (ics:ice,    iblk) = in1%qsat_srf_old  (ics:ice,    iblk) + in2%qsat_srf_old  (ics:ice,    iblk)
!!$    out%tte_corr      (ics:ice,    iblk) = in1%tte_corr      (ics:ice,    iblk) + in2%tte_corr      (ics:ice,    iblk)
!!$    out%fact_q_air    (ics:ice,    iblk) = in1%fact_q_air    (ics:ice,    iblk) + in2%fact_q_air    (ics:ice,    iblk)
!!$    out%fact_qsat_srf (ics:ice,    iblk) = in1%fact_qsat_srf (ics:ice,    iblk) + in2%fact_qsat_srf (ics:ice,    iblk)
!!$    out%fact_qsat_trans_srf (ics:ice,  iblk) = in1%fact_qsat_trans_srf (ics:ice,  iblk) + in2%fact_qsat_trans_srf (ics:ice,  iblk)
!!$
!!$  END SUBROUTINE add_t_srf_energy_state

#endif
END MODULE mo_srf_memory_class
