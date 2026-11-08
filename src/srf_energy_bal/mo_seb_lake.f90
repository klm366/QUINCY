!> Contains the routines for the surface energy balance on LAKE lct_type.
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
!>#### Surface energy balance on lakes
!>
!> Lakes are represented by a mixed layer of water of constant depth (`lake_mixed_layer_depth`)
!> and, if `l_ice_on_lakes`, by a layer of ice with snow on top. A lake tile is either ice-free or
!> completely covered by ice (`fract_lice` is 0 or 1).
!>
!> - Open water: the water temperature `t_lwtr` changes with the sum of net radiation, sensible and
!>   latent heat fluxes. It cannot drop below the melting point; the excess cooling is accumulated
!>   in `hflx_form_ice` until it suffices to form a slab of ice of thickness `lake_min_ice_depth`
!>   ([[calc_lake_wtr_temp]]).
!> - Ice: the ice thickness changes with the conductive heat flux through the ice, melting at the
!>   surface and sublimation/deposition ([[calc_lice_thickness]]). The temperature of a thin
!>   surface layer of the ice (with the snow on top) is computed from the surface fluxes and the
!>   conductive heat flux to the ice bottom, which is at the melting point; surface heat beyond the
!>   melting point melts the snow first and then the ice ([[calc_lake_ice_temp]]).
!>
!> In contrast to the land surface, the lake energy balance is solved explicitly, with the sensible
!> and latent heat fluxes of the previous time step. The new turbulent fluxes are computed in
!> [[update_surface_fluxes_lake]].
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"

MODULE mo_seb_lake
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: lp
  USE mo_exception, ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: debug_on, jsbach_runs_standalone, acc_stream

  dsl4jsb_Use_processes SEB_, TURB_, RAD_, HYDRO_, A2L_
  dsl4jsb_Use_config(SEB_)

  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_lake, update_surface_fluxes_lake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  !> Compute the new surface temperature on lakes
  !>
  !> The routine computes the new water temperature, the ice thickness and fraction, and the ice
  !> temperature with the heat fluxes from melting snow and ice. From these, it derives the
  !> surface state of the lake tile: the ice-fraction weighted temperature `t`, saturation
  !> specific humidity, dry static energy and heat capacity, and the snow on the lake.
  !> In the standalone model, the surface drag and Richtmyer-Morton coefficients are computed from
  !> the atmospheric forcing ([[mo_phy_schemes:update_drag]]), with the same values for water and ice.
  !>
  SUBROUTINE update_surface_energy_lake(tile, options)

    USE mo_phy_schemes,            ONLY: qsat_water, qsat_ice, qsat_mixed, surface_dry_static_energy, update_drag
    USe mo_jsb_physical_constants, ONLY: cpd, cvd
    USE mo_jsb4_forcing,           ONLY: forcing_options

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t,                    & !< Surface temperature: ice-fraction weighted mean of water and ice [K]
      & s_star,               & !< Surface dry static energy [m2 s-2]
      & t_rad4,               & !< Surface radiative temperature **4 [K4]
      & t_eff4,               & !< Effective surface temperature **4 [K4]
      & t_unfilt,             & !< Surface temperature, unfiltered (before Asselin) [K]
      & t_lwtr,               & !< Lake surface temperature (water) [K]
      & t_lice,               & !< Lake surface temperature (ice) [K]
      & qsat_lwtr,            & !< Saturation specific humidity over lake water [kg kg-1]
      & qsat_lice,            & !< Saturation specific humidity over lake ice [kg kg-1]
      & qsat_star,            & !< Saturation specific humidity at the surface [kg kg-1]
      & s_lwtr,               & !< Surface dry static energy over lake water [m2 s-2]
      & s_lice,               & !< Surface dry static energy over lake ice [m2 s-2]
      & heat_cap,             & !< Heat capacity of the lake surface layer (water or ice) [J m-2 K-1]
      & forc_hflx,            & !< Heat flux into the surface layer from below, zero on lakes [W m-2]
      & press_srf,            & !< Surface pressure (A2L_) [Pa]
      & snow,                 & !< Snow fall rate (A2L_) [kg m-2 s-1]
      & latent_hflx_wtr,      & !< Latent heat flux over lake water [W m-2]
      & sensible_hflx_wtr,    & !< Sensible heat flux over lake water [W m-2]
      & latent_hflx_ice,      & !< Latent heat flux over lake ice [W m-2]
      & sensible_hflx_ice,    & !< Sensible heat flux over lake ice [W m-2]
      & hflx_form_ice,        & !< Heat flux for ice formation on lakes [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux through lake ice [W m-2]
      & hflx_melt_lice,       & !< Heat flux from melting lake ice [W m-2]
      & hflx_melt_snow,       & !< Heat flux from melting snow on lake ice [W m-2]
      & fract_lice,           & !< Fraction of lake ice []
      & depth_lice,           & !< Thickness of ice on lake [m]
      & fract_snow_lice,      & !< Snow fraction on lake ice (HYDRO_) []
      & weq_snow_lice,        & !< Water content of snow on lake ice (HYDRO_) [m water equivalent]
      & weq_snow,             & !< Snow amount on the lake (HYDRO_, set here) [m water equivalent]
      & rad_net_lwtr,         & !< Net radiation over lake water (RAD_) [W m-2]
      & rad_net_lice,         & !< Net radiation over lake ice (RAD_) [W m-2]
      & t_air,                & !< Air temperature at lowest atmospheric level (A2L_) [K]
      & q_air,                & !< Specific humidity, lowest atm. level (A2L_) [kg kg-1]
      & wind_air,             & !< Wind speed at lowest atmospheric level (A2L_) [m s-1]
      & fact_q_air,           & !< Humidity scaling factor for q_air (TURB_) []
      & fact_qsat_srf,        & !< Humidity scaling factor for qsat at surface (TURB_) []
      & rough_h,              & !< Surface roughness length for heat (TURB_) [m]
      & rough_m,              & !< Surface roughness length for momentum (TURB_) [m]
      & drag_wtr,             & !< Surface drag over lake water (A2L_, standalone)
      & t_acoef_wtr,          & !< Richtmyer-Morton coef. A, dry static energy, water []
      & t_bcoef_wtr,          & !< Richtmyer-Morton coef. B, dry static energy, water [m2 s-2]
      & q_acoef_wtr,          & !< Richtmyer-Morton coef. A, specific humidity, water []
      & q_bcoef_wtr,          & !< Richtmyer-Morton coef. B, specific humidity, water [kg kg-1]
      & drag_ice,             & !< Surface drag over lake ice (A2L_, standalone)
      & t_acoef_ice,          & !< Richtmyer-Morton coef. A, dry static energy, ice []
      & t_bcoef_ice,          & !< Richtmyer-Morton coef. B, dry static energy, ice [m2 s-2]
      & q_acoef_ice,          & !< Richtmyer-Morton coef. A, specific humidity, ice []
      & q_bcoef_ice             !< Richtmyer-Morton coef. B, specific humidity, ice [kg kg-1]

    ! Local variables
    REAL(lp), DIMENSION(options%nc) :: &
      & hcap_wtr,           & !< Heat capacity of the lake mixed layer [J m-2 K-1]
      & zero,               & !< Zero heat flux [W m-2]
      & hflx_tmp,           & !< Sum of net radiation and turbulent heat fluxes over water [W m-2]
      & pch_wtr,            & !< Heat transfer coefficient C_h over water (standalone) []
      & ddrag_wtr             !< Derivative of drag over water w.r.t. surface temperature (standalone)

    REAL(lp), ALLOCATABLE :: &
      hcap_ice(:)             !< Heat capacity of the ice surface layer [J m-2 K-1]

    LOGICAL ::                       &
      & l_ice(options%nc),           & !< Indicator for ice cover
      & l_ice_formation(options%nc), & !< Indicator for formation of a new slab of ice
      & l_ice_melting(options%nc)      !< Indicator for complete melting of ice
    LOGICAL :: &
      & config_l_ice_on_lakes, &       !< Ice on lakes (l_ice_on_lakes)
      & use_tmx,               &       !< Coupling to the atmosphere with TMX
      & jsb_standalone                 !< ICON-Land runs standalone

    INTEGER :: &
      & iblk, &              !< Current block index
      & ics,  &              !< Index of first cell in chunk
      & ice,  &              !< Index of last cell in chunk
      & nc,   &              !< Number of cells in chunk
      & ic                   !< Cell index
    INTEGER  :: model_scheme !< Model scheme (MODEL_JSBACH or MODEL_QUINCY)
    REAL(lp) :: dtime        !< Time step length [s]
    REAL(lp) ::         &
      & cpd_or_cvd,     &    !< Heat capacity of dry air (cpd; with TMX: cvd) [J kg-1 K-1]
      & lhflx_tmp,      &    !< Sum of net radiation and turbulent heat fluxes over ice [W m-2]
      & ml_depth,       &    !< Depth of lake mixed layer [m]
      & min_ice_depth,  &    !< Minimum ice thickness to start ice formation [m]
      & min_fract_lice, &    !< Currently not used
      & min_hcap_wtr         !< Currently not used

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_lake'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(SEB_) .OR. .NOT. tile%contains_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx
    model_scheme = model%config%model_scheme

    IF (use_tmx) THEN
      cpd_or_cvd = cvd
    ELSE
      cpd_or_cvd = cpd
    END IF

    dsl4jsb_Get_config(SEB_)

    ml_depth      = dsl4jsb_Config(SEB_)%lake_mixed_layer_depth
    min_ice_depth = dsl4jsb_Config(SEB_)%lake_min_ice_depth
    config_l_ice_on_lakes = dsl4jsb_Config(SEB_)%l_ice_on_lakes

    jsb_standalone = jsbach_runs_standalone()

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(SEB_, t)              ! out
    IF (use_tmx) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_, s_star)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, t_rad4)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, t_eff4)       ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, qsat_star)    ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(SEB_, t_unfilt)       ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr)   ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)         ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lwtr)      ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, s_lwtr)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, heat_cap)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, forc_hflx)      ! out

    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)        ! out
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_wtr)   ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_wtr) ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice)        ! OUT
    IF (jsb_standalone) THEN
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_q_air)       ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_srf)    ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, rough_h)          ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, rough_m)          ! in
      dsl4jsb_Get_var2D_onChunk(A2L_, drag_wtr)          ! OUT if standalone
      dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef_wtr)       ! OUT if standalone
      dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef_wtr)       ! OUT if standalone
      dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef_wtr)       ! OUT if standalone
      dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef_wtr)       ! OUT if standalone
    END IF

    IF (config_l_ice_on_lakes) THEN
      dsl4jsb_Get_var2D_onChunk(A2L_, snow)              ! in
      dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)            ! inout
      dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lice)         ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, s_lice)            ! out
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_form_ice)     ! inout
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_cond_ice)
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_lice)
      dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_snow)

      dsl4jsb_Get_var2D_onChunk(SEB_, depth_lice)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_lice) ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_lice)   ! inout

      dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_ice)   ! in
      dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_ice) ! in

      dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice)      ! in

      IF (jsb_standalone) THEN
        dsl4jsb_Get_var2D_onChunk(A2L_, drag_ice)          ! OUT if standalone
        dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef_ice)       ! OUT if standalone
        dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef_ice)       ! OUT if standalone
        dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef_ice)       ! OUT if standalone
        dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef_ice)       ! OUT if standalone
      END IF

      ALLOCATE(hcap_ice(nc))
      !$ACC ENTER DATA CREATE(hcap_ice) ASYNC(acc_stream)
    ELSE
      ALLOCATE(hflx_form_ice(nc))
      !$ACC ENTER DATA CREATE(hflx_form_ice) ASYNC(acc_stream)
    END IF

    !$ACC DATA CREATE(hcap_wtr, zero, l_ice, l_ice_formation, l_ice_melting, hflx_tmp) ASYNC(acc_stream)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      zero            (ic) = 0._lp
      hcap_wtr        (ic) = 0._lp
      l_ice_formation (ic) = .FALSE.
      l_ice_melting   (ic) = .FALSE.
      hflx_tmp(ic) = rad_net_lwtr(ic) + latent_hflx_wtr(ic) + sensible_hflx_wtr(ic)
      IF (config_l_ice_on_lakes) THEN
        hcap_ice      (ic) = 0._lp
        l_ice         (ic) = fract_lice(ic) > 0.5_lp
      ELSE
        fract_lice    (ic) = 0._lp
        l_ice         (ic) = .FALSE.
        hflx_form_ice (ic) = 0._lp
      END IF
      forc_hflx       (ic) = 0._lp ! Always zero for lake
    END DO
    !$ACC END PARALLEL LOOP

    ! Compute surface drag and exchange coefficients for the offline model
    IF (jsb_standalone) THEN
      !$ACC WAIT
      !$ACC DATA &
      !$ACC   CREATE(ddrag_wtr, pch_wtr)

      ! Update drag and exchange coefficients for lake surface based on external forcing data
      CALL update_drag( &
      ! INTENT in
      & nc, dtime, &
      & model_scheme, &
      & t_air(:), press_srf(:), q_air(:), wind_air(:), &
      & t(:), fact_q_air(:), fact_qsat_srf(:), rough_h(:), rough_m(:), &
      & forcing_options(tile%owner_model_id)%heightWind, forcing_options(tile%owner_model_id)%heightHumidity, &
      & dsl4jsb_Config(SEB_)%coef_ril_tm1, dsl4jsb_Config(SEB_)%coef_ril_t, dsl4jsb_Config(SEB_)%coef_ril_tp1, &
      ! INTENT out
      & drag_wtr(:), ddrag_wtr(:), t_acoef_wtr(:), t_bcoef_wtr(:), q_acoef_wtr(:), q_bcoef_wtr(:), pch_wtr(:))

      IF (config_l_ice_on_lakes) THEN
        ! As input variables for the lake tile refer to the frozen and unfrozen fractions, we cannot
        ! provide separate exchange coefficients for water and ice. Both are set to identical values.
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          drag_ice(ic)    = drag_wtr(ic)
          t_acoef_ice(ic) = t_acoef_wtr(ic)
          t_bcoef_ice(ic) = t_bcoef_wtr(ic)
          q_acoef_ice(ic) = q_acoef_wtr(ic)
          q_bcoef_ice(ic) = q_bcoef_wtr(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      !$ACC END DATA
    END IF

    ! Calculate new water temperature and residual heat flux for ice formation over open water!
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      CALL calc_lake_wtr_temp( &
        & t_lwtr(ic),                                     & !  inout; Surface temperature over open water
        & hflx_tmp(ic),                                   & !  in;    Sum of all radiative and turbulent surface heat fluxes
        & hflx_form_ice(ic),                              & !  inout; Residual heat flux for ice formation (realized or unrealized)
        & hcap_wtr(ic),                                   & !  out;   Heat capacity of lake over open water
        & l_ice_formation(ic),                            & !  out;   Indicator for ice formation
        & l_ice_melting(ic),                              & !  inout; Indicator for ice melting = .FALSE.
        & l_ice(ic),                                      & !  in;    Indicator for ice cover
        & dtime,                                          & !  in;    Time step
        & ml_depth,                                       & !  in;    Mixed-layer depth
        & min_ice_depth)                                    !  in;    Minimum ice thickness
    END DO
    !$ACC END PARALLEL LOOP

    ! TODO How to put this on GPU? The following block is wrong!

      ! !$ACC KERNELS
    ! min_fract_lice = MINVAL(fract_lice)
    ! !$ACC END KERNELS

    ! !$ACC KERNELS
    ! min_hcap_wtr   = MINVAL(hcap_wtr)
    ! !$ACC END KERNELS

    ! IF (min_fract_lice <= 0.5_lp .AND. min_hcap_wtr < EPSILON(1._lp)) THEN
    !   CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    ! END IF

#ifndef _OPENACC
    IF (ANY(fract_lice(:) <= 0.5_lp .AND. hcap_wtr(:) < EPSILON(1._lp))) THEN
      CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    END IF
#endif

    IF (config_l_ice_on_lakes) THEN

      ! Calculate new ice thickness
      ! Tendency = freezing + conductive heat flux - ice melting + sublimation/deposition
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        CALL calc_lice_thickness(                           &
          & depth_lice(ic),                                 & !  inout; Ice thickness [m]
          & fract_lice(ic),                                 & !  inout; Indicator for ice cover
          & hflx_form_ice(ic),                              & !  inout; Heat flux for (complete) ice formation (l_ice_formation = T)
          & latent_hflx_ice(ic),                            & !  in;    Latent heat flux for ice
          & hflx_cond_ice(ic),                              & !  in;    Conductive heat flux for ice
          & hflx_melt_lice(ic),                             & !  inout; Heat flux from unrealized or complete ice melting
          & l_ice_formation(ic),                            & !  in;    Indicator for ice formation
          & l_ice_melting(ic),                              & !  inout; Indicator for ice melting
          & fract_snow_lice(ic),                            & !  in;    Fraction of snow on ice
          & min_ice_depth,                                  & !  in;    Mininum ice thickness
          & dtime                                           & !  in;    Time step
          & )
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate new water temperature for the case of complete ice melting
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
        DO ic=1,nc
          CALL calc_lake_wtr_temp( &
            & t_lwtr(ic),                                     & !  inout; Surface temperature over open water
            & zero(ic),                                       & !  in;    Sum of radiative and turbulent surface heat fluxes (not used)
            & hflx_melt_lice(ic),                             & !  inout; Residual heat flux from unrealized ice melting
            & hcap_wtr(ic),                                   & !  inout; Heat capacity of lake over open water (not changed)
            & l_ice_formation(ic),                            & !  inout; Indicator for ice formation (not changed)
            & l_ice_melting(ic),                              & !  in;    Indicator for ice melting = .TRUE.
            & l_ice(ic),                                      & !  in;    Indicator for ice cover
            & dtime,                                          & !  in;    Time step
            & ml_depth,                                       & !  in;    Mixed-layer depth
            & min_ice_depth)                                    !  in;    Minimum ice thickness
        END DO
      !$ACC END PARALLEL LOOP

      ! Calculate new ice temperature and heat fluxes from melting of snow and ice

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
      !$ACC   PRIVATE(lhflx_tmp)
      DO ic=1,nc
        lhflx_tmp = rad_net_lice(ic) + latent_hflx_ice(ic) + sensible_hflx_ice(ic)

        CALL calc_lake_ice_temp(                            &
          & t_lice(ic),                                     & ! inout
          & hcap_ice(ic),                                   & ! out
          & hflx_cond_ice(ic),                              & ! out
          & hflx_melt_snow(ic),                             & ! out
          & hflx_melt_lice(ic),                             & ! inout
          & lhflx_tmp,                                      & !  in;    Sum of all radiative and turbulent surface heat fluxes
          & latent_hflx_ice(ic),                            & !  in;
          & snow(ic),                                       & !  in; Snow fall
          & depth_lice(ic),                                 &
          & weq_snow_lice(ic),                              &
          & fract_snow_lice(ic),                            &
          & min_ice_depth,                                  &
          & dtime)
      END DO
      !$ACC END PARALLEL LOOP

    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic=1,nc
      IF (use_tmx) THEN
        qsat_lwtr(ic)  = qsat_water(t_lwtr(ic), press_srf(ic))
      ELSE
        qsat_lwtr(ic)  = qsat_mixed(t_lwtr(ic), press_srf(ic))
      END IF
      s_lwtr(ic)     = surface_dry_static_energy(t_lwtr(ic), qsat_lwtr(ic), cpd_or_cvd, jsb_standalone)
    END DO
    !$ACC END PARALLEL LOOP

    IF (config_l_ice_on_lakes) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        IF (use_tmx) THEN
          qsat_lice(ic) = qsat_ice(t_lice(ic), press_srf(ic))
        ELSE
          qsat_lice(ic) = qsat_mixed(t_lice(ic), press_srf(ic))
        END IF
        t(ic)         = (1._lp - fract_lice(ic)) * t_lwtr(ic) + fract_lice(ic) * t_lice(ic)
        t_unfilt(ic)  = t(ic)
        s_lice(ic)    = surface_dry_static_energy(t_lice(ic), qsat_lice(ic), cpd_or_cvd, jsb_standalone)
        IF (use_tmx) THEN
          s_star(ic) = (1._lp - fract_lice(ic)) * s_lwtr(ic)    + fract_lice(ic) * s_lice(ic)
          t_rad4(ic) = (1._lp - fract_lice(ic)) * t_lwtr(ic)**4 + fract_lice(ic) * t_lice(ic)**4
          t_eff4(ic) = t_rad4(ic)
          qsat_star(ic) = (1._lp - fract_lice(ic)) * qsat_lwtr(ic) + fract_lice(ic) * qsat_lice(ic)
        END IF
        heat_cap(ic)  = (1._lp - fract_lice(ic)) * hcap_wtr(ic) + fract_lice(ic) * hcap_ice(ic)
        weq_snow(ic)  = fract_lice(ic) * weq_snow_lice(ic)
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC EXIT DATA DELETE(hcap_ice) ASYNC(acc_stream)
      DEALLOCATE(hcap_ice)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        t(ic) = t_lwtr(ic)
        t_unfilt(ic) = t(ic)
        IF (use_tmx) THEN
          s_star(ic) = s_lwtr(ic)
          t_rad4(ic) = t_lwtr(ic)**4
          t_eff4(ic) = t_rad4(ic)
          qsat_star(ic) = qsat_lwtr(ic)
        END IF
        heat_cap(ic) = hcap_wtr(ic)
        weq_snow(ic) = 0._lp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC EXIT DATA DELETE(hflx_form_ice) ASYNC(acc_stream)
      DEALLOCATE(hflx_form_ice)
    END IF

    !$ACC END DATA

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

  !> Calculate new surface water temperature
  !>
  !> If the lake is ice-free, the water temperature of the mixed layer changes with the surface heat
  !> flux and the residual heat flux for ice formation from the previous time step. If the water
  !> would cool below the melting point, it is kept at the melting point and the excess cooling is
  !> returned in `hflx_resid`; a new slab of ice is formed (`l_ice_formation`) once it suffices to
  !> freeze a layer of thickness `min_ice_depth`. After complete melting of the ice
  !> (`l_ice_melting`), the residual heat from melting heats the water.
  !>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lake_wtr_temp(t_lwtr, hflx, hflx_resid, hcap, l_ice_formation, l_ice_melting, &
    &                                          l_ice, dtime, ml_depth, min_ice_depth)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, clw, rhoi, alf

    REAL(lp), INTENT(inout) :: t_lwtr          !< Surface temperature over open water [K]
    REAL(lp), INTENT(in)    :: hflx            !< Sum of all radiative and turbulent surface heat fluxes
                                               !< [W m-2]
    REAL(lp), INTENT(inout) :: hflx_resid      !< Heat flux for ice formation or melting [W m-2]
                                               !< l_ice_melting = .TRUE.  : heat flux from complete ice melting (input)
                                               !<                           zero on output
                                               !<               = .FALSE. : heat flux from unrealized ice formation from
                                               !<                           previous time step (input),
                                               !<                           heat flux from realized or unrealized ice formation
                                               !<                           in this time step (output)
    REAL(lp), INTENT(inout) :: hcap            !< Heat capacity of lake over open water [J m-2 K-1]
    LOGICAL,  INTENT(inout) :: l_ice_formation !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting   !< Indicator for (complete) ice melting
    LOGICAL,  INTENT(in)    :: l_ice           !< Indicator for ice cover

    REAL(lp), INTENT(in)    :: dtime            !< Time step [s]
    REAL(lp), INTENT(in)    :: ml_depth         !< Mixed-layer depth [m]
    REAL(lp), INTENT(in)    :: min_ice_depth    !< Mininum ice thickness [m]

    REAL(lp) :: &
      & t_star, &   !< Preliminary water temperature [K]
      & t_freeze    !< Water temperature below which a new slab of ice is formed [K]

    ! Lake tile is either not covered by ice or completely covered.
    IF (l_ice) RETURN

    IF (l_ice_melting) THEN
      t_lwtr = t_lwtr + dtime * hflx_resid / hcap
      hflx_resid = 0._lp
      l_ice_melting = .FALSE.
    ELSE
      ! Heat capacity of lake over open water
      ! Note: This is currently constant but in the future, different mixed-layer depths might be used
      hcap = rhoh2o * clw * ml_depth

      ! Water temperature below which a new slab of ice is formed
      t_freeze = tmelt - min_ice_depth * rhoi * alf / hcap

      ! Preliminary water temperature
      t_star = t_lwtr + dtime * (hflx + hflx_resid) / hcap

      IF (t_star >= tmelt) THEN              ! No ice formation
        t_lwtr = t_star
        hflx_resid = 0._lp
        l_ice_formation = .FALSE.
      ELSE                                   ! Realized or unrealized ice formation
        t_lwtr = tmelt
        hflx_resid = (t_star - tmelt) * hcap / dtime      ! < 0
        l_ice_formation = t_star <= t_freeze
        ! l_ice_formation = .TRUE.  : hflx_resid will be used to form a new slab of ice
        !                 = .FALSE. : cooling is not large enough to form a new slab of ice and
        !                             hflx_resid will be used as forcing flux in the next time
        !                             step to compute water temperature t_lwtr
      END IF

    END IF

  END SUBROUTINE calc_lake_wtr_temp

  !> Calculate new ice temperature, snow depth and heat fluxes
  !>
  !> The snow depth on the ice changes with snowfall and sublimation/deposition. The temperature of
  !> a thin surface layer (ice of thickness `hi_min` and the snow on top) is computed implicitly
  !> from the surface heat flux and the conductive heat flux through the ice and snow to the ice
  !> bottom at the melting point, using an effective thickness hi_eff = hi + (ki/ks)(rhoh2o/rhos) hs.
  !> Heat that would warm the surface above the melting point melts the snow first and then the
  !> ice. Finally, the conductive heat flux used for the change of ice thickness in the next time
  !> step is diagnosed. For ice thinner than `hi_min`, the ice temperature is set to the melting
  !> point and the snow is removed.
  !>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lake_ice_temp(t_lice, hcap_ice, hflx_cond_ice, hflx_melt_snow, hflx_melt_lice, &
    &                                          hflx, lhflx, snowfall, hi, hs, fs, hi_min, dtime)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, rhoi, rhoi, rhos, als, alf, ki, ks, ci, cs !, clw

    REAL(lp), INTENT(inout) :: t_lice         !< Lake ice surface temperature [K]
    REAL(lp), INTENT(out)   :: hcap_ice       !< Heat capacity of thin ice surface layer [J m-2 K-1]
    REAL(lp), INTENT(out)   :: hflx_cond_ice  !< Conductive heat flux through ice and snow [W m-2]
    REAL(lp), INTENT(out)   :: hflx_melt_snow !< Heat flux used for melting snow [W m-2]
    REAL(lp), INTENT(inout) :: hflx_melt_lice !< Heat flux from melting of ice [W m-2]
                                              !< (in:  residual heat flux from unrealized ice melting,
                                              !<  out: heat flux from melting of ice)
    REAL(lp), INTENT(in)    :: hflx           !< Sum of all surface heat fluxes (net radiation + sensible
                                              !< + latent) [W m-2]
    REAL(lp), INTENT(in)    :: lhflx          !< Latent surface heat flux for ice [W m-2]
    REAL(lp), INTENT(in)    :: snowfall       !< Snow fall                [kg m-2 s-1]
    REAL(lp), INTENT(in)    :: hi             !< Ice thickness [m]
    REAL(lp), INTENT(inout) :: hs             !< Snow depth [m water equivalent]
    REAL(lp), INTENT(in)    :: fs             !< Fraction of snow on ice []
    REAL(lp), INTENT(in)    :: hi_min         !< Minimum ice thickness, also thickness of the ice
                                              !< surface layer [m]
    REAL(lp), INTENT(in)    :: dtime          !< Time step [s]

    REAL(lp) :: hi_eff       !< Effective ice thickness including snow cover [m]
    REAL(lp) :: hflx_melt    !< Excess heat available for melting snow and ice [W m-2]
    REAL(lp) :: delta_hs     !< Change of snow depth due to melting [m water equivalent]
    REAL(lp) :: delta_t_lice !< Change of ice temperature due to snow melt [K]

    hcap_ice = rhoi * ci * hi_min + rhoh2o * cs * hs

    hi_eff = hi + (ki / ks) * (rhoh2o / rhos) * hs

    IF (hi < hi_min) THEN
      t_lice = tmelt
      hflx_cond_ice = 0._lp
      hflx_melt_snow = 0._lp
      hflx_melt_lice = 0._lp
      hs = 0._lp
      RETURN
    END IF

    ! Preliminary snow depth
    hs = hs + (dtime / rhoh2o) * ( snowfall           & ! Snowfall
      &                          + fs * lhflx / als   & ! Sublimation (< 0) resp. deposition (> 0)
      &                          )
    hs = MAX(hs, 0._lp)

    ! Preliminary ice temperature
    t_lice = ( hflx + hflx_melt_lice + hcap_ice * t_lice / dtime + ki * tmelt / hi_eff ) / ( hcap_ice / dtime + ki / hi_eff)

    hflx_melt_snow = 0._lp
    hflx_melt_lice  = 0._lp

    IF (t_lice > tmelt .AND. hs > 0._lp) THEN        ! Melting for snow first
      ! Excess heat which is available for melting of snow (first) and ice (as soon as snow is  melted away completely)
      hflx_melt = ( hcap_ice / dtime + ki /hi_eff ) * (t_lice - tmelt)
      ! Change in snow depth due to melting
      delta_hs = MIN(hs, dtime * hflx_melt / alf / rhoh2o)             ! positive
      ! Change in ice temperature resulting from snow melt
      hflx_melt_snow = delta_hs * rhoh2o *alf / dtime                  ! Heat that was used for melting the snow
      delta_t_lice = hflx_melt_snow / ( hcap_ice/dtime + ki/hi_eff)     ! positive
      t_lice = t_lice - delta_t_lice
      ! Update snow depth, eff. ice thickness and heat capacity of thin ice layer (including snow)
      hs = hs - delta_hs
      hi_eff = hi + (ki / ks) * (rhoh2o / rhos) * hs
      hcap_ice = rhoi * ci * hi_min + rhoh2o * cs * hs
    END IF
    IF (t_lice > tmelt) THEN                         ! Melting of ice after complete melting of snow
      hflx_melt_lice = ( hcap_ice/dtime + ki/hi_eff) * (t_lice - tmelt)
      t_lice = tmelt
    END IF

    ! Diagnose conductive heat flux for changing ice thickness in next time step
    hflx_cond_ice = - ki * (t_lice - tmelt) / hi_eff

  END SUBROUTINE calc_lake_ice_temp

  !> Calculation of ice depth and fraction
  !>
  !> A new slab of ice is formed from the accumulated heat flux for ice formation. Existing ice
  !> grows or shrinks with the conductive heat flux, melting at the surface and, on the snow-free
  !> part, sublimation/deposition. If the ice becomes thinner than `hi_min`, it is either kept at
  !> `hi_min` and the missing melt energy is passed to the next time step (unrealized melting) or,
  !> if it melts completely, removed and the residual heat is passed to the water.
  !>
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_lice_thickness( &
    & hi,                                        &
    & fi,                                        &
    & hflx_form_ice,                             &
    & lhflx,                                     &
    & hflx_cond_ice,                             &
    & hflx_melt_lice,                            &
    & l_ice_formation,                           &
    & l_ice_melting,                             &
    & fs,                                        &
    & hi_min,                                    &
    & dtime                                      &
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: rhoi, alf, als

    REAL(lp), INTENT(inout) :: hi              !< Ice thickness [m]
    REAL(lp), INTENT(inout) :: fi              !< Ice fraction  []
    REAL(lp), INTENT(inout) :: hflx_form_ice   !< Residual heat flux from ice formation (complete or
                                               !< unrealized) [W m-2]
    REAL(lp), INTENT(in)    :: lhflx           !< Latent heat flux of ice  [W m-2 = J m-2 s-1]
    REAL(lp), INTENT(in)    :: hflx_cond_ice   !< Conductive heat flux of ice [W m-2 = J m-2 s-1]
    REAL(lp), INTENT(inout) :: hflx_melt_lice  !< Heat flux from melting of ice [W m-2 = J m-2 s-1]
                                               !< in:  flux from unrealized ice melting from previous time step
                                               !< out: flux from unrealized or complete ice melting during
                                               !< this time step
    LOGICAL,  INTENT(in)    :: l_ice_formation !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting   !< Indicator for complete ice melting
    REAL(lp), INTENT(in)    :: fs              !< Fraction of snow on ice []

    REAL(lp), INTENT(in)    :: hi_min          !< Minimum ice thickness [m]
    REAL(lp), INTENT(in)    :: dtime           !< Time step [s]

    IF (l_ice_formation) THEN                              ! Formation of a new slab of ice

      hi = - hflx_form_ice * dtime / (rhoi * alf)
      fi = 1._lp
      hflx_form_ice = 0._lp

    ELSE IF (hi >= hi_min) THEN                            ! Ice already present

      hi = hi + (dtime / rhoi) * (   (hflx_cond_ice - hflx_melt_lice) / alf    &
        &                          + (1._lp - fs) * lhflx / als               &
        &                        )

      IF (hi >= hi_min) THEN
        hflx_melt_lice = 0._lp
        fi = 1._lp
        l_ice_melting = .FALSE.
      ELSE IF (hi <= 0._lp) THEN                              ! Complete melting
        hflx_melt_lice = - rhoi * alf * hi / dtime             ! >= 0
        hi = 0._lp
        fi = 0._lp
        l_ice_melting = .TRUE.
      ELSE                                                    ! Unrealized melting
        hflx_melt_lice = rhoi * alf * (hi_min - hi) / dtime    ! > 0
        hi = hi_min
        fi = 1._lp
        l_ice_melting = .FALSE.
      END IF

    END IF

  END SUBROUTINE calc_lice_thickness

  !> Compute the sensible and latent heat fluxes on lakes
  !>
  !> The fluxes are computed separately for water and ice from the respective surface dry static
  !> energy and evaporation (HYDRO_), and then weighted with the ice fraction. Without TMX, the dry
  !> static energy of the lowest atmospheric level is obtained by back-substitution with the
  !> Richtmyer-Morton coefficients for water and ice; with TMX, its old value `t_bcoef` and the
  !> exchange coefficient `ch` of TURB_ are used for both water and ice.
  !>
  SUBROUTINE update_surface_fluxes_lake(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef ! q_effective, surface_dry_static_energy
    USE mo_jsb_physical_constants, ONLY: alv, als

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options !< Additional run-time parameters

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: sensible_hflx     !< Sensible heat flux of the lake (water and ice) [W m-2]
    dsl4jsb_Real2D_onChunk :: latent_hflx       !< Latent heat flux of the lake (water and ice) [W m-2]
    dsl4jsb_Real2D_onChunk :: s_lwtr            !< Surface dry static energy over lake water [m2 s-2]
    dsl4jsb_Real2D_onChunk :: t_acoef           !< Richtmyer-Morton coef. A, dry static energy (TMX) []
    dsl4jsb_Real2D_onChunk :: t_bcoef           !< Richtmyer-Morton coef. B, dry static energy (TMX) [m2 s-2]
    dsl4jsb_Real2D_onChunk :: t_acoef_wtr       !< Richtmyer-Morton coef. A, dry static energy, water []
    dsl4jsb_Real2D_onChunk :: t_bcoef_wtr       !< Richtmyer-Morton coef. B, dry static energy, water [m2 s-2]
    dsl4jsb_Real2D_onChunk :: drag_wtr          !< Surface drag over lake water (A2L_, without TMX)
    dsl4jsb_Real2D_onChunk :: latent_hflx_wtr   !< Latent heat flux over lake water [W m-2]
    dsl4jsb_Real2D_onChunk :: sensible_hflx_wtr !< Sensible heat flux over lake water [W m-2]
    dsl4jsb_Real2D_onChunk :: evapo_wtr         !< Evaporation from lake water (HYDRO_) [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: s_lice            !< Surface dry static energy over lake ice [m2 s-2]
    dsl4jsb_Real2D_onChunk :: t_acoef_ice       !< Richtmyer-Morton coef. A, dry static energy, ice []
    dsl4jsb_Real2D_onChunk :: t_bcoef_ice       !< Richtmyer-Morton coef. B, dry static energy, ice [m2 s-2]
    dsl4jsb_Real2D_onChunk :: drag_ice          !< Surface drag over lake ice (A2L_, without TMX)
    dsl4jsb_Real2D_onChunk :: latent_hflx_ice   !< Latent heat flux over lake ice [W m-2]
    dsl4jsb_Real2D_onChunk :: sensible_hflx_ice !< Sensible heat flux over lake ice [W m-2]
    dsl4jsb_Real2D_onChunk :: evapo_ice         !< Evaporation from lake ice (HYDRO_) [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: fract_lice        !< Fraction of lake ice []
    dsl4jsb_Real2D_onChunk :: ch                !< Transfer coefficient for heat (TURB_, TMX) [kg m-2 s-1]

    REAL(lp) ::                  &
      & s_air,     &  !< Dry static energy at lowest atmospheric level [m2 s-2]
      & heat_tcoef    !< Heat transfer coefficient (rho*C_h*|v|) [kg m-2 s-1]

    INTEGER :: &
      & iblk, &         !< Current block index
      & ics,  &         !< Index of first cell in chunk
      & ice,  &         !< Index of last cell in chunk
      & nc,   &         !< Number of cells in chunk
      & ic              !< Cell index
    REAL(lp) :: &
      & dtime, &        !< Time step length [s]
      & alpha           !< Implicitness factor of the atmosphere-surface coupling []
    LOGICAL  :: use_tmx !< Coupling to the atmosphere with TMX

    TYPE(t_jsb_model), POINTER :: model !< This instance of ICON-Land

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes_lake'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    alpha   = options%alpha

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model => Get_model(tile%owner_model_id)
    use_tmx = model%config%use_tmx

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    IF (use_tmx) THEN
      dsl4jsb_Get_memory(TURB_)
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)       ! IN
      dsl4jsb_Get_var2D_onChunk(TURB_, ch)               ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef_wtr)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef_wtr)       ! IN
      dsl4jsb_Get_var2D_onChunk(A2L_,   drag_wtr)        ! in
    END IF

    dsl4jsb_Get_var2D_onChunk(SEB_,   s_lwtr)            ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx)     ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_wtr) ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_wtr)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_wtr)         ! in

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      IF (.NOT. use_tmx) THEN
        dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef_ice)       ! IN
        dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef_ice)       ! IN
        dsl4jsb_Get_var2D_onChunk(A2L_,   drag_ice)      ! in
      END IF
      dsl4jsb_Get_var2D_onChunk(SEB_,   s_lice)            ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_ice) ! out
      dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_ice)   ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)         ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)        ! in
    END IF

    ! ================================================================================================================================
    ! Surface fluxes for lake water

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(s_air, heat_tcoef)
    DO ic=1,nc

      IF (use_tmx) THEN
        s_air = t_bcoef(ic)  ! Old dry static energy at lowest atmospheric level
        heat_tcoef = ch(ic)  ! Transfer coefficient; TODO: distinguish between wtr and ice?
        sensible_hflx_wtr(ic) = heat_tcoef * (s_air - s_lwtr(ic))      ! Sensible energy flux
      ELSE
        s_air = t_acoef_wtr(ic) * s_lwtr(ic) + t_bcoef_wtr(ic)    ! New dry static energy at lowest atmospheric level by backsubstitution
        heat_tcoef = heat_transfer_coef(drag_wtr(ic), dtime, alpha)  ! Transfer coefficient
        sensible_hflx_wtr(ic) = heat_tcoef * (s_air - s_lwtr(ic))      ! Sensible heat flux
      END IF
      latent_hflx_wtr(ic) = alv * evapo_wtr(ic) ! Latent heat flux

    END DO
    !$ACC END PARALLEL LOOP

    ! ================================================================================================================================
    ! Surface fluxes for lake ice

    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(acc_stream) GANG VECTOR PRIVATE(s_air, heat_tcoef)
      DO ic=1,nc

        IF (use_tmx) THEN
          s_air = t_bcoef(ic) ! Old dry static energy at lowest atmospheric level
          heat_tcoef = ch(ic) ! Transfer coefficient; TODO: distinguish between wtr and ice?
          sensible_hflx_ice(ic) = heat_tcoef * (s_air - s_lice(ic))      ! Sensible energy flux
        ELSE
          s_air = t_acoef_ice(ic) * s_lice(ic) + t_bcoef_ice(ic)    ! New dry static energy at lowest atmospheric level by backsubstitution
          heat_tcoef = heat_transfer_coef(drag_ice(ic), dtime, alpha)  ! Transfer coefficient
          sensible_hflx_ice(ic) = heat_tcoef * (s_air - s_lice(ic))      ! Sensible heat flux
        END IF
        latent_hflx_ice(ic) = als * evapo_ice(ic) ! Latent heat flux

        ! ================================================================================================================================
        ! Surface fluxes for lake (wtr and ice aggregated)
        sensible_hflx(ic) = (1._lp - fract_lice(ic)) * sensible_hflx_wtr(ic) + fract_lice(ic) * sensible_hflx_ice(ic)
        latent_hflx  (ic) = (1._lp - fract_lice(ic)) * latent_hflx_wtr  (ic) + fract_lice(ic) * latent_hflx_ice  (ic)

      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
      DO ic=1,nc
        sensible_hflx(ic) = sensible_hflx_wtr(ic)
        latent_hflx  (ic) = latent_hflx_wtr  (ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_lake

#endif
END MODULE mo_seb_lake
