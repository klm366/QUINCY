!>
!! @brief Contains the routines for the soil and snow energy processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_process
#ifndef __NO_JSBACH__

  USE mo_jsb_physical_constants, ONLY: rhoh2o
  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: relative_humidity_soil, init_soil_temperature, &
    & calc_vol_heat_capacity, calc_thermal_conductivity, calc_soil_temperature, calc_snow_temperature, &
    & calc_soiltemp_old, Get_liquid_max

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_process'

CONTAINS

  !!
  ! @brief      Calculates the volumetric heat capacity of soil (de Vries, 1963)
  !
  ! @param [in]  vol_heat_cap_soil_dry_mineral  Heat capacity of dry mineral soil solid
  ! @param [in]  porosity                       Volumetric soil porosity
  ! @param [in]  fract_organic                  Volume fraction (rel. to 1-porosity) of organic soil
  ! @param [in]  fract_water                    Volume fraction of liquid water
  ! @param [in]  fract_ice                      Volume fraction of frozen water
  ! @param [in]  layer_depth                    (dynamic) Layer depth above bedrock
  ! @param [in]  dz                             (fixed) soil layer depth
  !
  ! @return     The volumetric heat capacity of soil.
  !!
  ELEMENTAL PURE FUNCTION calc_vol_heat_capacity( &
    & vol_heat_cap_soil_dry_mineral, &
    & porosity,              &
    & fract_organic,         &
    & fract_water,           &
    & fract_ice,             &
    & layer_depth,           &
    & dz                     &
    & ) RESULT(vol_heat_cap)

    USE mo_sse_constants, ONLY: vol_hcap_water, vol_hcap_ice, vol_hcap_org, vol_hcap_air, vol_hcap_bedrock

    REAL(wp), INTENT(in) :: porosity, vol_heat_cap_soil_dry_mineral, fract_organic, fract_water, fract_ice, layer_depth, dz
    REAL(wp)             :: vol_heat_cap

    IF (layer_depth > 0._wp) THEN  ! above bedrock
      vol_heat_cap = (1._wp - porosity)                                            & ! vol. fraction of solid (mineral + organic)
        &              * ( (1._wp - fract_organic) * vol_heat_cap_soil_dry_mineral & ! heat cap of mineral soil
        &                  +        fract_organic  * vol_hcap_org                  & ! heat cap of organic soil
        &                ) &
        &            + fract_water * vol_hcap_water                                & ! heat cap of liquid water
        &            + fract_ice   * vol_hcap_ice                                  & ! heat cap of frozen water
        ! @todo : porosity - fract_water - fract_ice should not become lower than zero.
        &            + MAX((porosity - fract_water - fract_ice),0._wp) * vol_hcap_air  ! heat cap of remaining free pores
      vol_heat_cap = ( vol_heat_cap * layer_depth + vol_hcap_bedrock * (dz - layer_depth) ) / dz ! soil and bedrock mixture
    ELSE
      vol_heat_cap = vol_hcap_bedrock
    END IF

  END FUNCTION calc_vol_heat_capacity

  !!
  ! @brief      Calculates the soil thermal conductivity after Johansen 1977.
  !
  ! @param      heat_cond_soil_mineral  The heat condition soil mineral
  ! @param      water                   The water
  ! @param      ice                     The ice
  ! @param      porosity                Volumentric porosity of soil, including soil mineral and organic part
  ! @param      vol_field_cap           The volume field capability
  ! @param      fract_organic           The fract organic
  ! @param      layer_depth             The layer depth
  ! @param      dz                      { parameter_description }
  ! @param      heat_cond               The heat condition
  !
  ! @return     The thermal conductivity.
  !!
  ELEMENTAL PURE FUNCTION calc_thermal_conductivity( &
    & heat_cond_soil_mineral,   &
    & water, ice,               &
    & porosity, vol_field_cap,  &
    & fract_organic,            &
    & layer_depth,              &
    & dz                        &
    & ) RESULT(heat_cond)

    USE mo_jsb_physical_constants, ONLY: dens_soil
    USE mo_sse_constants,          ONLY: hcond_org, hcond_dry_org, hcond_water, hcond_ice, hcond_bedrock
    USE mo_hydro_constants,        ONLY: vol_porosity_org

    REAL(wp), INTENT(in) :: heat_cond_soil_mineral, water, ice, porosity, vol_field_cap, fract_organic, layer_depth, dz
    REAL(wp)             :: heat_cond

    REAL(wp) :: &
      & dens_bulk,      & !< Bulk density of dry soil (including pores)
      & porosity_mineral, & !< Volumentric porosity of mineral soil part
      & hcond_dry_min,  & !< Dry thermal conductivity soil (including mineral soil and pores)
      & hcond_dry,      & !< Dry thermal conductivity of soil (including mineral and organic soil and pores)
      & vpor_sat,       & !< Effective porosity for saturated soil computation (pores remaining if soil is saturated)
      & saturation,     & !< Degree of saturation
      & kersten,        & !< Kersten number
      & hcond_soil,     & !< Thermal conductivity of soil solids (mineral and organic)
      & hcond_sat,      & !< Thermal conductivity for saturated soil
      & fract_water,    & !< Effective fractional water volume for saturated soil computation
      & fract_ice         !< Effective fractional ice volume for saturated soil computation

    IF (layer_depth > 0._wp) THEN  ! above bedrock
      ! Thermal conductivity of dry soil
      porosity_mineral = (porosity - fract_organic * vol_porosity_org) / (1._wp - fract_organic)
      dens_bulk = dens_soil * (1._wp - porosity_mineral)
      hcond_dry_min = (0.135_wp * dens_bulk + 64.7_wp) / (dens_soil - 0.947_wp * dens_bulk)
      hcond_dry = (1._wp - fract_organic) * hcond_dry_min + fract_organic * hcond_dry_org

      IF (water + ice <= EPSILON(1._wp)) THEN
        heat_cond = hcond_dry
        RETURN
      END IF

      ! Thermal conductivity of soil solids
      hcond_soil = (1._wp - fract_organic) * heat_cond_soil_mineral + fract_organic * hcond_org

      ! Thermal conductivity of saturated soil
      saturation = (water + ice) / layer_depth / porosity  ! fractional volume of water/ice rel. to pore volume

      ! ! --- Alternate formulation ---
      ! ! @todo: in jsbach3, frozen soil is checked by t_soil <= tmelt, however, the computations in calc_soil_temperature
      ! ! are such that some ice can remain in a layer even if t_soil > tmelt.
      ! IF (ice > 0._wp) THEN  ! Frozen case
      !   kersten = saturation
      !   fract_water = water / layer_depth
      !   fract_ice   = porosity - fract_water
      !   hcond_sat =   hcond_soil      ** (1._wp-porosity) &
      !     &         * hcond_water     ** fract_water      &
      !     &         * hcond_ice       ** fract_ice
      ! ELSE                   ! Unfrozen case
      !   ! ! This is the case for fine unfrozen soil
      !   ! IF (saturation > 0.1_wp) THEN
      !   !   kersten = LOG10(saturation) + 1._wp
      !   ! ELSE
      !   !   kersten = 0._wp
      !   ! END IF
      !   ! This is the case for coarse unfrozen soil
      !   IF (saturation > 0.05_wp) THEN
      !     kersten = 0.7_wp * LOG10(saturation) + 1._wp
      !   ELSE
      !     kersten = 0._wp
      !   END IF
      !   fract_water = porosity
      !   hcond_sat =   hcond_soil      ** (1._wp-porosity) &
      !     &         * hcond_water     ** porosity
      ! END IF

      ! @todo: in jsbach3, frozen soil is checked by t_soil <= tmelt, however, the computations in calc_soil_temperature
      ! are such that some ice can remain in a layer even if t_soil > tmelt.
      IF (ice > 0._wp) THEN  ! Frozen case
        kersten = saturation
      ELSE                   ! Unfrozen case
        ! This is the case for fine unfrozen soil
        IF (saturation > 0.1_wp) THEN
          kersten = LOG10(saturation) + 1._wp
        ELSE
          kersten = 0._wp
        END IF
        ! ! This is the case for coarse unfrozen soil
        ! IF (saturation > 0.05_wp) THEN
        !   kersten = 0.7_wp * LOG10(saturation) + 1._wp
        ! ELSE
        !   kersten = 0._wp
        ! END IF
      END IF
      ! @todo: vol_field_cap should probably be replaced by MIN(vol_field_cap, liquid_max/layer_depth) if l_supercool=T
      vpor_sat = porosity - vol_field_cap  ! Empty pores remaining if soil is saturated
      fract_water = water / (1._wp - vpor_sat) / layer_depth ! Water fraction rel. to completely saturated soil
      fract_ice   =   ice / (1._wp - vpor_sat) / layer_depth ! Ice fraction rel. to completely saturated soil
      hcond_sat =   hcond_soil      ** (1._wp-fract_water-fract_ice) &
        &         * hcond_water     ** fract_water                   &
        &         * hcond_ice       ** fract_ice

      heat_cond = (hcond_sat - hcond_dry) * kersten + hcond_dry

      ! Finally, consider soil and bedrock mixture
      heat_cond = ( heat_cond * layer_depth + hcond_bedrock * (dz - layer_depth) ) / dz

      heat_cond = MAX(heat_cond, hcond_dry)

    ELSE
      heat_cond = hcond_bedrock
    END IF

  END FUNCTION calc_thermal_conductivity

  REAL(wp) ELEMENTAL FUNCTION relative_humidity_soil(ws, fc)

    USE mo_math_constants, ONLY: pi

    REAL(wp), INTENT(in) :: &
      ws, & !< moisture content [m], e.g. of uppermost soil layer
      fc    !< field capacity (maximum moisture content) [m], e.g. of uppermost soil layer

    REAL(wp) :: moisture

    moisture = MIN(ws, fc)
    IF (moisture > 0._wp) THEN
      relative_humidity_soil = 0.5_wp * (1._wp - COS(moisture * pi / fc)) ! relative_humidity increases with higher moisture and lower fc
                                                                     ! Note, fc has to be  >= than moisture!
    ELSE
      relative_humidity_soil = 0._wp
    END IF

  END FUNCTION relative_humidity_soil

  !>
  !
  ! Initialize soil and surface temperature
  !
  ! Method:
  !
  ! Starting from the tslclim field temperatures are set in relation
  ! to the depth of the soil layer and position of the initial
  ! day in the annual cycle.
  !
  ! tsl is at 0.07 m
  ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m (s. soiltemp)
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 2000, original source
  ! L. Dumenil, MPI, June 1989, original source of *initemp*
  !             which is now part of this routine
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! E. Roeckner, MPI, Sep 2002, initialization of mixed layer ocean
  ! T. Raddatz, MPI, Mai 2005, adaption to JSBACH
  !
  SUBROUTINE init_soil_temperature(cmid, tclim, tsoil, tsurf)

    USE mo_jsb_time,       ONLY: t_datetime, get_time_start, get_year, &
                                 get_year_length, get_month_length, get_day_length, get_year_day
    USE mo_math_constants, ONLY: pi
    USE mo_exception,      ONLY: finish

    REAL(wp), INTENT(in)  :: cmid(:)     !< Depths of soil layer mid points
    REAL(wp), INTENT(in)  :: tclim(:,:)  !< Climatological monthly surface temperature
    REAL(wp), INTENT(out) :: tsoil(:,:)  !< Soil temperatures
    REAL(wp), INTENT(out) :: tsurf(:)    !< Surface temperature

    INTEGER               :: nsoil, ndim
    REAL(wp)              :: zkap, zsqrt
    REAL(wp), ALLOCATABLE :: zrange(:), zdmax(:)
    INTEGER,  ALLOCATABLE :: jmax(:)
    TYPE(t_datetime)      :: inidate
    INTEGER               :: iniyear
    INTEGER               :: im, i, jl
    REAL(wp)              :: year_len, day_in_year
    REAL(wp)              :: month_len(12), month_mid(12)

    REAL(wp), PARAMETER :: cmid_offset = 0.07_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_temperature'

    ndim  = SIZE(tsoil,1)
    nsoil = SIZE(tsoil,2)

    ALLOCATE(zrange(ndim), jmax(ndim), zdmax(ndim))

    zkap = 7.5e-7_wp

    inidate = get_time_start()
    iniyear = get_year(inidate)

    year_len    = REAL(get_year_length(iniyear), wp)
    day_in_year = INT(get_year_day(inidate))

    DO im=1,12
      month_len(im) = REAL(get_month_length(iniyear, im), wp)
      month_mid(im) = month_len(im) * 0.5_wp
    END DO

    zsqrt = SQRT(zkap * year_len * get_day_length() / pi)

    ! Calendar month of maximum surface temperature
    jmax(:) = MAXLOC(tclim(:,:), DIM=2)

    ! Difference between months of maximum and minimum surface temperature
    zrange(:) = MAXVAL(tclim(:,:), DIM=2) - MINVAL(tclim(:,:), DIM=2)

    ! Calendar day of maximum surface temperature
    zdmax(:) = 0._wp
    DO jl=1,ndim
      DO im=1,jmax(jl)
        zdmax(jl) = zdmax(jl) + month_len(im)
      END DO
      zdmax(jl) = zdmax(jl) - month_mid(jmax(jl))
    END DO

    ! Initialize soil temperatures
    DO i=1,nsoil
      tsoil(:,i) = SUM(tclim(:,:), DIM=2) / 12._wp &
                 & + 0.5_wp * zrange(:) * EXP(-(cmid(i)-cmid_offset) / zsqrt) &
                 &   * COS(2._wp * pi * (day_in_year - zdmax(:)) / year_len - (cmid(i)-cmid_offset) / zsqrt)
    END DO

    ! Set surface temperature to uppermost soil temperature
    tsurf(:) = tsoil(:,1)

    DEALLOCATE(zrange, jmax, zdmax)

  END SUBROUTINE init_soil_temperature

  !
  ! @brief      Calculates the soil temperature.
  !
  ! @param      nidx          The nidx
  ! @param      nsoil         The nsoil
  ! @param      dz            { parameter_description }
  ! @param      mids          The cmid
  ! @param      delta_time    The delta time
  ! @param      lstart        The lstart
  ! @param      is_glacier
  ! @param [in] l_freeze
  ! @param [in] l_supercool
  ! @param      t_srf_unfilt  The t srf unfilt
  ! @param      vol_heat_cap  The volume heat capability
  ! @param      heat_cond     The heat condition
  ! @param [in]      ws_max
  ! @param [in]      matrix_pot
  ! @param [in]      bclapp
  ! @param [in,out]  t_soil        The t soil
  ! @param [in,out]  t_soil_acoef  The t soil acoef
  ! @param [in,out]  t_soil_bcoef  The t soil bcoef
  ! @param [in,out]  ws
  ! @param [in,out]  ice
  ! @param [out]     hcap_grnd     The hcap grnd
  ! @param [out]     grnd_hflx     The grnd hflx
  ! @param [out]     frozen
  ! @param [out]     melt
  ! @param [out]     thaw_depth
  !
  ! @return     { description_of_the_return_value } !
  !
  SUBROUTINE calc_soil_temperature(       &
    & nidx, nsoil, dz, mids,                         &
    & delta_time, lstart, &
    & l_freeze, l_supercool, &
    & t_srf_unfilt,              &
    & vol_heat_cap, heat_cond,          &
    & t_soil, t_soil_acoef, t_soil_bcoef,   &
    & hcap_grnd, grnd_hflx, &
    & thaw_depth, &
    ! optional
    & ws_max, matrix_pot, bclapp, ws, ice, frozen, melt &
    & )
    !
    !
    !   AUTHOR:  FREDERIC HOURDIN     30/01/92
    !
    !            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
    !            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
    !
    !            J.-P. SCHULZ   MPI - OCTOBER 1997 :
    !               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
    !               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
    !               ECHAM4 GCM.
    !            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
    !            U.Schlese DKRZ - February 2000  new soil temperatures
    !            L Kornblueh, MPI, January 2003, removed MERGE
    !
    !            Adapted to JSBACH by Thomas Raddatz, Mai 2004
    !
    !   OBJECTIVE:  COMPUTATION OF:
    !               THE GROUND TEMPERATURE EVOLUTION
    !               THE GROUND SPECIFIC HEAT "CAPCAL"
    !               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
    !
    !
    !   METHOD:  IMPLICIT TIME INTEGRATION
    !
    !   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
    !           T(K+1) = C(K) + D(K)*T(K)  (1)
    !   THE COEFFICIENTS C (=t_soil_acoef) AND D (=t_soil_bcoef) ARE COMPUTED AT THE
    !   T-DT TIME-STEP.
    !   ROUTINE STRUCTURE:
    !   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
    !   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
    !     PROFILE FOR THE T+DT TIME-STEP
    !   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
    !     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
    !            FDIFF = A + B TS(T+DT)
    !     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
    !            WITH F0 = A + B (TS(T))
    !                 CAPCAL = B*DT
    !
    !
    !     ------------------------------------------------------------------
    !
    !     MODIFIED BY ALTUG EKICI
    !
    !     INCLUDED:
    !     -PHASE CHANGE
    !     -HEAT TRANSFER PARAMETERS AS A FUNCTION OF SOIL WATER AND ICE
    !     -SUPERCOOLED WATER EQUATION
    !     -THAW DEPTH CALCULATION
    !     -MULTI LAYERED SNOW SCHEME and KEEPING 5 SOIL LAYERS AT ALL TIMES
    !     -ORGANIC LAYER ABOVE THE SOIL
    !
    !----------------------------------------------------------------------------------------------

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, rhoh2o, rhoi

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    ::  nidx                     !! length of the vector
    INTEGER,  INTENT(in)    ::  nsoil                    !! number of soil layers
    REAL(wp), INTENT(in)    ::  delta_time               !! Time step
    LOGICAL,  INTENT(in)    ::  lstart                   !! if T, start of experiment
    LOGICAL,  INTENT(in)    ::  l_freeze, l_supercool
    REAL(wp), INTENT(in)    ::  dz(nsoil)                !! soil layer thickness [m]
    REAL(wp), INTENT(in)    ::  mids(nsoil)              !! depth of mids of soil layers [m]
    REAL(wp), INTENT(in)    ::  t_srf_unfilt(nidx)       !! surface temperature at top of soil [K]
    REAL(wp), INTENT(in)    ::  vol_heat_cap(nidx,nsoil) !! soil heat capacity [J/m^3K]
    REAL(wp), INTENT(in)    ::  heat_cond(nidx,nsoil)    !! soil thermal conductivity [J/m/s/K]

    REAL(wp), INTENT(inout) ::  t_soil(nidx,nsoil)       !! soil temperature [K]
    REAL(wp), INTENT(inout) ::  t_soil_acoef(nidx,nsoil) !! soil A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) ::  t_soil_bcoef(nidx,nsoil) !! soil B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   ::  hcap_grnd(nidx)          !! heat capacity of the ground [J m-2 K-1]
    REAL(wp), INTENT(out)   ::  grnd_hflx(nidx)          !! ground heat flux [J m-2 s-1]
    REAL(wp), INTENT(out)   ::  thaw_depth(nidx)         !! Thawing depth [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  ws_max(nidx,nsoil)       !! Maximum water storage [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  matrix_pot(nidx,nsoil)   !! Matrix potential
    REAL(wp), INTENT(in), OPTIONAL    ::  bclapp(nidx,nsoil)       !! Clapp and Hornberger parameter
    REAL(wp), INTENT(inout), OPTIONAL ::  ws(nidx,nsoil)           !! Water storage [m]
    REAL(wp), INTENT(inout), OPTIONAL ::  ice(nidx,nsoil)          !! Ice storage [m]
    REAL(wp), INTENT(out), OPTIONAL   ::  frozen(nidx,nsoil)       !! Water flux from freezing soil water [kg m-2 s-1]
    REAL(wp), INTENT(out), OPTIONAL  ::  melt(nidx,nsoil)         !! Water flux from melting soil ice    [kg m-2 s-1]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: jl, jk
    REAL(wp) :: z1(nidx)
    REAL(wp) :: zd1(nsoil)
    REAL(wp) :: zdz1(nidx,nsoil), zdz2(nidx,nsoil)
    REAL(wp) :: heat_cap(nidx,nsoil), liquid_max(nidx,nsoil), freezemelt(nidx,nsoil)

    !-----------------------------------------------------------------------------------------------
    !  Computation of useful constants

    DO jk = 1,nsoil-1
       zd1(jk) = 1._wp / (mids(jk+1) - mids(jk))
    END DO

    heat_cap(:,:) = SPREAD(dz(:), DIM=1, ncopies=nidx) * vol_heat_cap(:,:)

    !-----------------------------------------------------------------------------------------------
    !  Computation of soil temperaturs
    !    from surface temperatur and A and B coefficients computed at the previous time step
    !-----------------------------------------------------------------------------------------------

    IF (PRESENT(ws)) THEN
      frozen(:,:)   = 0._wp
      melt  (:,:)   = 0._wp
    END IF
    thaw_depth(:) = 0._wp

    IF (.NOT. lstart) THEN

      ! uppermost layer
      t_soil(:,1) = t_srf_unfilt(:)
      ! deeper layers
      DO jk = 1, nsoil-1
        t_soil(:,jk+1) = t_soil_acoef(:,jk) + t_soil_bcoef(:,jk) * t_soil(:,jk)
      END DO

      ! Phase change calculations
      IF (l_freeze .AND. PRESENT(ws)) THEN

        liquid_max(:,:) = Get_liquid_max( &
          & l_supercool    , &
          & t_soil(:,:)    , &
          & ws_max(:,:)    , &
          & matrix_pot(:,:), &
          & bclapp(:,:)      &
          & )

        WHERE (t_soil(:,:) < tmelt .AND. ws(:,:) > liquid_max(:,:))
          frozen(:,:) = MIN(ws(:,:) - liquid_max(:,:), &
            &               heat_cap(:,:) * (tmelt - t_soil(:,:)) / (alf * rhoh2o) &
            &              )
          ws (:,:) = ws (:,:) - frozen(:,:)
          ice(:,:) = ice(:,:) + frozen(:,:)
          t_soil(:,:) = t_soil(:,:) + frozen(:,:) * alf * rhoh2o / heat_cap(:,:)
        ELSE WHERE (t_soil(:,:) > tmelt .AND. ice(:,:) > 0._wp)
          melt(:,:) = MIN(ice(:,:), heat_cap(:,:) * (t_soil(:,:) - tmelt) / (alf * rhoi))
          ws (:,:) = ws (:,:) + melt(:,:)
          ice(:,:) = ice(:,:) - melt(:,:)
          t_soil(:,:) = t_soil(:,:) - melt(:,:) * alf * rhoi / heat_cap(:,:)
        END WHERE

        frozen(:,:) = frozen(:,:) * rhoh2o / delta_time
        melt  (:,:) = melt  (:,:) * rhoi   / delta_time

      END IF

      ! Thawing depth
      thaw_depth(:) = mids(nsoil)   ! Initialize with maximum thawing depth
      WHERE (t_soil(:,1) <= tmelt)  ! Top soil layer is frozen
        thaw_depth(:) = 0._wp
      END WHERE
      DO jl=1,nidx
        DO jk=2,nsoil
          IF  (       thaw_depth(jl)  == mids(nsoil) & ! still the inital value
            &   .AND. t_soil(jl,jk-1) >  tmelt &
            &   .AND. t_soil(jl,jk  ) <= tmelt &
            & ) THEN
            thaw_depth(jl) = mids(jk-1) + (mids(jk)-mids(jk-1)) * (t_soil(jl,jk-1) - tmelt) / (t_soil(jl,jk-1) - t_soil(jl,jk))
          END IF
        END DO
      END DO

    END IF


    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------

    DO jk = 1, nsoil
      zdz2(:,jk) = heat_cap(:,jk) / delta_time
    END DO

    DO jk = 1, nsoil-1
      zdz1(:,jk) = zd1(jk) * heat_cond(:,jk)
    END DO

    ! lowest soil layer
    z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
    t_soil_acoef(:,nsoil-1) = zdz2(:,nsoil) * t_soil(:,nsoil) / z1(:)
    t_soil_bcoef(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)

    ! soil layers above
    DO jk = nsoil-1, 2, -1
      z1(:) = 1._wp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._wp - t_soil_bcoef(:,jk)))
      t_soil_acoef(:,jk-1) = (t_soil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * t_soil_acoef(:,jk)) * z1(:)
      t_soil_bcoef(:,jk-1) = zdz1(:,jk-1) * z1(:)
    END DO

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------

    grnd_hflx(:) = zdz1(:,1) * (t_soil_acoef(:,1) + (t_soil_bcoef(:,1) - 1._wp) * t_soil(:,1))
    hcap_grnd(:) = (zdz2(:,1) * delta_time + delta_time * (1._wp - t_soil_bcoef(:,1)) * zdz1(:,1))

  END SUBROUTINE calc_soil_temperature

  SUBROUTINE calc_snow_temperature(       &
    & nidx, nsnow,                          &
    & dz_soil, &
    & delta_time, lstart,                   &
    & dz, itop_old, itop, &
    & t_srf_unfilt,              &
    & vol_heat_cap, heat_cond,          &
    & vol_heat_cap_soil, heat_cond_soil, &
    & t_soil, t_soil_acoef, t_soil_bcoef, &
    & t_snow, t_snow_acoef, t_snow_bcoef,   &
    & t_soil_top, hcap_grnd, grnd_hflx)
    !
    !
    !   AUTHOR:  FREDERIC HOURDIN     30/01/92
    !
    !            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
    !            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
    !
    !            J.-P. SCHULZ   MPI - OCTOBER 1997 :
    !               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
    !               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
    !               ECHAM4 GCM.
    !            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
    !            U.Schlese DKRZ - February 2000  new soil temperatures
    !            L Kornblueh, MPI, January 2003, removed MERGE
    !
    !            Adapted to JSBACH by Thomas Raddatz, Mai 2004
    !
    !   OBJECTIVE:  COMPUTATION OF:
    !               THE GROUND TEMPERATURE EVOLUTION
    !               THE GROUND SPECIFIC HEAT "CAPCAL"
    !               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
    !
    !
    !   METHOD:  IMPLICIT TIME INTEGRATION
    !
    !   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
    !           T(K+1) = C(K) + D(K)*T(K)  (1)
    !   THE COEFFICIENTS C (=t_soil_acoef) AND D (=t_soil_bcoef) ARE COMPUTED AT THE
    !   T-DT TIME-STEP.
    !   ROUTINE STRUCTURE:
    !   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
    !   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
    !     PROFILE FOR THE T+DT TIME-STEP
    !   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
    !     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
    !            FDIFF = A + B TS(T+DT)
    !     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
    !            WITH F0 = A + B (TS(T))
    !                 CAPCAL = B*DT
    !
    !
    !     ------------------------------------------------------------------
    !
    !     MODIFIED BY ALTUG EKICI
    !
    !     INCLUDED:
    !     -PHASE CHANGE
    !     -HEAT TRANSFER PARAMETERS AS A FUNCTION OF SOIL WATER AND ICE
    !     -SUPERCOOLED WATER EQUATION
    !     -THAW DEPTH CALCULATION
    !     -MULTI LAYERED SNOW SCHEME and KEEPING 5 SOIL LAYERS AT ALL TIMES
    !     -ORGANIC LAYER ABOVE THE SOIL
    !
    !----------------------------------------------------------------------------------------------

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    :: nidx                       !! length of the vector
    INTEGER,  INTENT(in)    :: nsnow                      !! number of (fixed) snow layers
    REAL(wp), INTENT(in)    :: delta_time                 !! Time step
    LOGICAL,  INTENT(in)    :: lstart                     !! if T, start of experiment
    REAL(wp), INTENT(in)    :: dz(nidx,nsnow)             !! (dynamic) snow layer thickness [m]
    REAL(wp), INTENT(in)    :: dz_soil(2)                 !! Thickness of two uppermost soil layers [m]
    INTEGER,  INTENT(in)    :: itop(nidx), itop_old(nidx) !! Top snow layer for previous and new time step
    REAL(wp), INTENT(in)    :: t_srf_unfilt(nidx)         !! surface temperature at top of snow [K]
    REAL(wp), INTENT(in)    :: vol_heat_cap(nidx)         !! snow volumetric heat capacity [J/m^3K]
    REAL(wp), INTENT(in)    :: heat_cond(nidx)            !! snow thermal conductivity [J/m/s/K]
    REAL(wp), INTENT(in)    :: vol_heat_cap_soil(nidx)    !! top soil layer volumetric heat capacity
    REAL(wp), INTENT(in)    :: heat_cond_soil(nidx)       !! top soil layer heat conductivity
    REAL(wp), INTENT(in)    :: t_soil(nidx)               !! top soil layer temperature
    REAL(wp), INTENT(in)    :: t_soil_acoef(nidx)         !! top soil layer A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_soil_bcoef(nidx)         !! top soil layer B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(inout) :: t_snow(nidx,nsnow)         !! snow temperature [K]
    REAL(wp), INTENT(inout) :: t_snow_acoef(nidx,nsnow)   !! snow A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) :: t_snow_bcoef(nidx,nsnow)   !! snow B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   :: t_soil_top(nidx)           !! new temperature of top soil layer
    REAL(wp), INTENT(out)   :: hcap_grnd(nidx)            !! heat capacity of the top snow layer [J m-2 K-1]
    REAL(wp), INTENT(out)   :: grnd_hflx(nidx)            !! ground heat flux [J m-2 s-1]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: jl, jk, itop_new !, itop_prev
    REAL(wp) :: zmid(nidx,nsnow+2)
    REAL(wp) :: z1(nidx)
    REAL(wp) :: zd1(nidx,nsnow+1)
    REAL(wp) :: zdz1(nidx,nsnow+1), zdz2(nidx,nsnow+1)

    grnd_hflx(:) = 0._wp
    hcap_grnd(:) = 0._wp
    t_soil_top(:) = t_srf_unfilt(:)

    !
    !-----------------------------------------------------------------------------------------------
    !  Computation of useful constants
    !
    ! dz and zmid are zero for unused snow layers
    !
    ! mid of snow layers
    zmid(:,1) = dz(:,1) / 2._wp
    DO jk=2,nsnow
      zmid(:,jk) = zmid(:,jk-1) + 0.5_wp * (dz(:,jk-1) + dz(:,jk))
      WHERE (zmid(:,jk) > zmid(:,jk-1))
        zd1(:,jk-1) = 1._wp / (zmid(:,jk) - zmid(:,jk-1))
      END WHERE
    END DO
    ! mid of two uppermost soil layers
    zmid(:,nsnow+1) = zmid(:,nsnow) + 0.5_wp * (dz(:,nsnow) + dz_soil(1))
    zmid(:,nsnow+2) = zmid(:,nsnow+1) + 0.5_wp * (dz_soil(1) + dz_soil(2))
    zd1(:,nsnow) = 1._wp / (zmid(:,nsnow+1) - zmid(:,nsnow))
    zd1(:,nsnow+1) = 1._wp / (zmid(:,nsnow+2) - zmid(:,nsnow+1))

    ! Reset old values accounting for changes in number of snow layers
    !
    DO jl=1,nidx
      ! itop_prev = itop_old(jl)
      itop_new = itop(jl)
      ! t_snow_acoef(jl,1:itop_new-1) = 0._wp
      ! t_snow_bcoef(jl,1:itop_new-1) = 0._wp
      t_snow      (jl,1:itop_new-1) = 0._wp
      ! DO jk=itop_new,nsnow                                     ! At least one snow layer present
      !   IF (itop_new < itop_prev .AND. jk < itop_prev) THEN    ! New snow layer
      !     IF (itop_prev == nsnow+1) THEN                       ! No snow layers present in previous time step
      !       t_snow_acoef(jl,jk) = t_soil_acoef(jl)             ! so we use coefficients and old temperature from top soil layer
      !       t_snow_bcoef(jl,jk) = t_soil_bcoef(jl)
      !       t_snow      (jl,jk) = t_soil      (jl)
      !     ELSE                                                 ! Use values from previous top snow layer
      !       t_snow_acoef(jl,jk) = t_snow_acoef(jl,itop_prev)
      !       t_snow_bcoef(jl,jk) = t_snow_bcoef(jl,itop_prev)
      !       t_snow      (jl,jk) = t_snow      (jl,itop_prev)
      !     END IF
      !   ! ELSE                                                 ! Snow layer was present before
      !   !   ! Nothing to do
      !   END IF
      ! END DO
    END DO
    !
    !-----------------------------------------------------------------------------------------------
    !  Computation of new snow temperatures
    !-----------------------------------------------------------------------------------------------
    !
    IF (.NOT. lstart) THEN
      t_snow(:,:) = 0._wp
      DO jk=1,nsnow
        WHERE (itop(:) <= jk .AND. jk <= itop_old(:))           ! zero+ new snow layers: layers from new top layer to top prev layer
          t_snow(:,jk) = t_srf_unfilt(:)
        ELSE WHERE (itop(:) > itop_old(:) .AND. jk == itop(:))  ! snow present but fewer snow layers : new top layer
          t_snow(:,jk) = t_srf_unfilt(:)
        END WHERE
      END DO
      DO jk=2,nsnow
        WHERE (jk > MAX(itop(:), itop_old(:)))                   ! Snow layers below top new and old layer
          t_snow(:,jk) = t_snow_acoef(:,jk-1) + t_snow_bcoef(:,jk-1) * t_snow(:,jk-1)
        END WHERE
      END DO
      ! New temperature for top soil layer (t_srf_unfilt where there was or is no snow)
      WHERE (MAX(itop(:),itop_old(:)) <= nsnow)
        t_soil_top(:) = t_snow_acoef(:,nsnow) + t_snow_bcoef(:,nsnow) * t_snow(:,nsnow)
      ELSE WHERE
        t_soil_top(:) = t_srf_unfilt(:)
      END WHERE
    END IF

    !
    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------

    DO jl=1,nidx
      itop_new = itop(jl)
      t_snow_acoef(jl,1:itop_new-1) = 0._wp
      t_snow_bcoef(jl,1:itop_new-1) = 0._wp
    END DO

    DO jk = 1,nsnow
      WHERE (jk >= itop(:))
        zdz2(:,jk) = vol_heat_cap(:) * dz(:,jk) / delta_time
        zdz1(:,jk) = zd1(:,jk) * heat_cond(:)
      END WHERE
    END DO
    zdz2(:,nsnow+1) = vol_heat_cap_soil(:) * dz_soil(1) / delta_time
    zdz1(:,nsnow+1) = zd1(:,nsnow+1) * heat_cond_soil(:)

    ! lowest snow layer
    WHERE (itop(:) <= nsnow)  ! Snow present
      z1(:) = 1._wp / (zdz2(:,nsnow+1) + zdz1(:,nsnow) + zdz1(:,nsnow+1) * (1._wp - t_soil_bcoef(:)))
      t_snow_acoef(:,nsnow) = (t_soil(:) * zdz2(:,nsnow+1) + zdz1(:,nsnow+1) * t_soil_acoef(:)) * z1(:)
      t_snow_bcoef(:,nsnow) = zdz1(:,nsnow) * z1(:)
    ELSE WHERE
      t_snow_acoef(:,nsnow) = 0._wp
      t_snow_bcoef(:,nsnow) = 0._wp
    END WHERE
    ! snow layers above up to top snow layer
    DO jk = nsnow, 2, -1
      WHERE (jk > itop(:))
        z1(:) = 1._wp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._wp - t_snow_bcoef(:,jk)))
        t_snow_acoef(:,jk-1) = (t_snow(:,jk) * zdz2(:,jk) + zdz1(:,jk) * t_snow_acoef(:,jk)) * z1(:)
        t_snow_bcoef(:,jk-1) = zdz1(:,jk-1) * z1(:)
      ! ELSE WHERE
      !   t_snow_acoef(:,nsnow) = 0._wp
      !   t_snow_bcoef(:,nsnow) = 0._wp
      END WHERE
    END DO

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !
    DO jk=1,nsnow
      WHERE (itop(:) == jk)  ! Snow present
        grnd_hflx(:) = zdz1(:,jk) * (t_snow_acoef(:,jk) + (t_snow_bcoef(:,jk) - 1._wp) * t_snow(:,jk))
        hcap_grnd(:) = (zdz2(:,jk) * delta_time + delta_time * (1._wp - t_snow_bcoef(:,jk)) * zdz1(:,jk))
      END WHERE
    END DO

  END SUBROUTINE calc_snow_temperature

  SUBROUTINE calc_soiltemp_old(nidx,nsoil,cdel,cmid            &
                           , delta_time                        &
                           , lstart                            &
                           , pts                               &
                           , psn                               &
                           , psocond, prgcgn                   &
                           , pgrndc, pgrndd                    &
                           , ptsoil                            &
                           , pgrndcapc, pgrndhflx              &
                           , ldglac                            &
                           )
  !
  !   AUTHOR:  FREDERIC HOURDIN     30/01/92
  !
  !            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
  !            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
  !
  !            J.-P. SCHULZ   MPI - OCTOBER 1997 :
  !               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
  !               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
  !               ECHAM4 GCM.
  !            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
  !            U.Schlese DKRZ - February 2000  new soil temperatures
  !            L Kornblueh, MPI, January 2003, removed MERGE
  !
  !            Adapted to JSBACH by Thomas Raddatz, Mai 2004
  !            Adapted to ICON  by Thomas Raddatz, Sep 2011
  !   OBJECTIVE:  COMPUTATION OF:
  !               THE GROUND TEMPERATURE EVOLUTION
  !               THE GROUND SPECIFIC HEAT "CAPCAL"
  !               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
  !
  !
  !   METHOD:  IMPLICIT TIME INTEGRATION
  !
  !   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
  !           T(K+1) = C(K) + D(K)*T(K)  (1)
  !   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
  !   T-DT TIME-STEP.
  !   ROUTINE STRUCTURE:
  !   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
  !   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
  !     PROFILE FOR THE T+DT TIME-STEP
  !   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
  !     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
  !            FDIFF = A + B TS(T+DT)
  !     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
  !            WITH F0 = A + B (TS(T))
  !                 CAPCAL = B*DT
  !
  !     ------------------------------------------------------------------
  !
  !   DECLARATIONS:
  !
  !!$ TR USE mo_jsbach_constants   , ONLY: RhoH2O
  !!$ TR USE mo_time_control       , ONLY: lstart
  !
  !-----------------------------------------------------------------------
  !  ARGUMENTS
  !
    INTEGER, Intent(in)  ::  nidx                 !! length of the vector
    INTEGER, Intent(in)  ::  nsoil                !! number of soil layers (fixed to 5)
    REAL(wp), Intent(in) ::  cdel(nsoil)
    REAL(wp), Intent(in) ::  cmid(nsoil)
    REAL(wp), Intent(in)     ::  delta_time       !! time step
    LOGICAL, INTENT(in)  :: lstart
    REAL(wp), Intent(in)     ::  pts(nidx)            !! surface temperature at top of soil [K]
    REAL(wp), Intent(in)     ::  psn(nidx)            !! equivalent snow depth [m water]
    REAL(wp), Intent(in)     ::  psocond(nidx)        !! soil heat conductivity
    REAL(wp), Intent(in)     ::  prgcgn(nidx)         !! soil heat capacity [J/m^3K]
    REAL(wp), Intent(inout)  ::  pgrndc(nidx,nsoil)   !!
    REAL(wp), Intent(inout)  ::  pgrndd(nidx,nsoil)   !!
    REAL(wp), Intent(inout)  ::  ptsoil(nidx,nsoil)   !! soil temperature [K]
    REAL(wp), Intent(out)    ::  pgrndcapc(nidx)      !!
    REAL(wp), Intent(out)    ::  pgrndhflx(nidx)      !! ground heat flux
    LOGICAL, Intent(in)  ::  ldglac(nidx)         !! glacier mask
  !
  !     ------------------------------------------------------------------
  !
  !  local Variables
  !
    INTEGER :: jk
    REAL(wp) :: zso_cond(nidx), zso_capa(nidx)
    REAL(wp) :: z1(nidx)
    REAL(wp) :: zd1(nsoil)
    REAL(wp) :: zdz1(nidx,nsoil),   zdz2(nidx,nsoil)
    REAL(wp) :: zkappa(nidx,nsoil), zcapa(nidx,nsoil)
    REAL(wp) :: zsnow_h(nidx), zx1(nidx), zx2(nidx)
    REAL(wp) :: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa
  !
  !     ------------------------------------------------------------------
  !
  !*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
  !
  !*    1.1 SOME CONSTANTS.
  !
    zrici = 2.09e+06_wp                                 !! volumetric heat capacity of ice [j/m**3/k]
    zdifiz = 12.e-07_wp                                 !! temperature diffusivity of ice  [m**2/s]
    zsn_cond = 0.31_wp                                  !! snow thermal conductivity [j/s/m/k]
    zsn_dens = 330.0_wp                                 !! snow density              [kg/m**3]
    zsn_capa = 634500.0_wp                              !! snow  heat capacity   [j/m**3/k]
  !
  !*    1.2 COMPUTING SOME USEFUL CONSTANTS.
  !
    DO jk = 1,nsoil-1
       zd1(jk) = 1._wp / (cmid(jk+1) - cmid(jk))
    END DO
  !
  !*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
  !*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
  !
    WHERE (ldglac(:))
      zso_capa(:) = zrici
      zso_cond(:) = zso_capa(:) * zdifiz
    ELSEWHERE
      zso_capa(:) = prgcgn(:)
      zso_cond(:) = psocond(:)
    END WHERE
  !
  !*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
  !
    DO jk = 1,nsoil
        zkappa(:,jk) = zso_cond(:)
        zcapa(:,jk)  = zso_capa(:)
    END DO

  !
  !   --------------------------------------------------------------
  !   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
  !   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
  !   --------------------------------------------------------------
  !
  !   Upper layer
  !
    IF (.NOT. lstart) THEN

      ptsoil(:,1) = pts(:)
  !
  !   Deeper layers
  !
      DO jk = 1,nsoil-1
        ptsoil(:,jk+1) = pgrndc(:,jk) + pgrndd(:,jk) * ptsoil(:,jk)
      END DO

    END IF
  !
  !   ---------------------------------------------------------------
  !   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
  !   ---------------------------------------------------------------
  !
    zsnow_h(:) = psn(:) * RhoH2O / zsn_dens
  !
  !*       Special treatment for first layer
  !
    WHERE ( zsnow_h(:) > cmid(2) )
      zcapa(:,1) = zsn_capa
      zkappa(:,1) = zsn_cond
    ELSEWHERE( zsnow_h(:) > 0.0_wp .AND. zsnow_h(:) <= cmid(2) )
      zx1 = zsnow_h(:) / cmid(2)
      zx2 = ( cmid(2) - zsnow_h(:)) / cmid(2)
      zcapa(:,1) = zx1 * zsn_capa + zx2 * zso_capa(:)
      zkappa(:,1) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
    ELSEWHERE
      zcapa(:,1) = zso_capa(:)
      zkappa(:,1) = zso_cond(:)
    ENDWHERE
  !
    DO jk = 2, nsoil - 2
       WHERE ( zsnow_h(:) > cmid(jk+1) )
          zcapa(:,jk) = zsn_capa
          zkappa(:,jk) = zsn_cond
       ELSEWHERE ( zsnow_h(:) > cmid(jk) .AND. zsnow_h(:) <= cmid(jk+1) )
          zx1 = (zsnow_h(:) - cmid(jk)) * zd1(jk)
          zx2 = ( cmid(jk+1) - zsnow_h(:)) * zd1(jk)
          zcapa(:,jk) = zx1 * zsn_capa + zx2 * zso_capa(:)
          zkappa(:,jk) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
       ELSEWHERE
          zcapa(:,jk) = zso_capa(:)
          zkappa(:,jk) = zso_cond(:)
       ENDWHERE
    END DO
  !
    DO jk=1,nsoil
      zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
    END DO
  !
    DO jk=1,nsoil-1
      zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
    END DO
  !
    z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
    pgrndc(:,nsoil-1) = zdz2(:,nsoil) * ptsoil(:,nsoil) / z1(:)
    pgrndd(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)
  !
    DO jk=nsoil-1,2,-1
      z1(:) = 1._wp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._wp - pgrndd(:,jk)))
      pgrndc(:,jk-1) = (ptsoil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * pgrndc(:,jk)) * z1(:)
      pgrndd(:,jk-1) = zdz1(:,jk-1) * z1(:)
    END DO
  !
  !   ---------------------------------------------------------
  !   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
  !   CALORIFIC CAPACITY OF THE GROUND:
  !   ---------------------------------------------------------
  !
    pgrndhflx(:) = zdz1(:,1) * (pgrndc(:,1) + (pgrndd(:,1) - 1._wp) * ptsoil(:,1))
    pgrndcapc(:) = (zdz2(:,1) * delta_time + delta_time * (1._wp - pgrndd(:,1)) * zdz1(:,1))

  END SUBROUTINE calc_soiltemp_old

  PURE ELEMENTAL FUNCTION Get_liquid_max( &
    & l_supercool, &
    & temp, w_max, matrix_pot, bclapp &
    & ) RESULT(liquid_max)

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, grav

    LOGICAL,  INTENT(in) :: l_supercool
    REAL(wp), INTENT(in) :: &
      & temp, w_max, matrix_pot, bclapp
    REAL(wp) :: liquid_max

    liquid_max = 0._wp

    IF (l_supercool) THEN
      IF (temp < tmelt) THEN
        ! supercooled water equation (Niu&Yang,2006)
        liquid_max =                                                              &
          & w_max                                                     & ! Maximum water storage
          & * ( alf * (tmelt - temp) / MAX(1.e-10_wp, grav * temp * matrix_pot) ) &
          &   ** ( -1._wp / MAX(1._wp, bclapp))
      END IF
    END IF

  END FUNCTION Get_liquid_max

#endif
END MODULE mo_sse_process
