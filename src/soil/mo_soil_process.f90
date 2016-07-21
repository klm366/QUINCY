!>
!! @brief Contains the routines for the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_process
#ifndef __NO_JSBACH__

  USE mo_jsb_physical_constants, ONLY: rhoh2o
  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Get_snow_fract_noforest, relative_humidity, update_surf_down_soil, update_surf_down_glacier, &
    & update_soiltemp, update_surf_up, init_soil_temperature

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_process'

CONTAINS

  !< Roesch et al. 2002, Climate Dynamics
  !
  REAL(wp) ELEMENTAL FUNCTION Get_snow_fract_noforest(snow, orodev)

    USE mo_soil_constants, ONLY: wsn2frac_const, wsn2frac_eps, wsn2frac_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   &
      & orodev

    Get_snow_fract_noforest = wsn2frac_const * TANH(snow * 100._wp) &
                               & * SQRT(snow * 1000._wp / (snow * 1000._wp + wsn2frac_eps + wsn2frac_sigfac * orodev))

  END FUNCTION Get_snow_fract_noforest

  REAL(wp) ELEMENTAL FUNCTION relative_humidity(ws, fc)

    USE mo_math_constants, ONLY: pi

    REAL(wp), INTENT(in) :: &
      ws, & !< moisture content
      fc    !< field capacity (maximum moisture content)

    REAL(wp) :: moisture

    moisture = MIN(ws, fc)
    IF (moisture > 0._wp) THEN
      relative_humidity = 0.5_wp * (1._wp - COS(moisture * pi / fc))
    ELSE
      relative_humidity = 0._wp
    END IF

  END FUNCTION relative_humidity

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
  ! TBD: remove fixed soil depths
  !
  SUBROUTINE init_soil_temperature(tclim, tsoil, tsurf)

    USE mo_jsb_time,       ONLY: t_datetime, get_time_start, get_year, &
                                 get_year_length, get_month_length, get_day_length, get_year_day
    USE mo_math_constants, ONLY: pi
    USE mo_exception,      ONLY: finish

    REAL(wp), INTENT(in)  :: tclim(:,:)
    REAL(wp), INTENT(out) :: tsoil(:,:), tsurf(:)

    INTEGER :: nsoil, ndim
    REAL(wp) :: zkap, zsqrt
    REAL(wp), POINTER :: zrange(:), zdmax(:), zd(:)
    INTEGER, POINTER :: jmax(:)
    TYPE(t_datetime) :: inidate
    INTEGER :: iniyear
    INTEGER :: im, i, jl
    REAL(wp) :: year_len, day_in_year
    REAL(wp) :: month_len(12), month_mid(12)

    ndim  = SIZE(tsoil,1)
    nsoil = SIZE(tsoil,2)
    IF (nsoil /= 5) CALL finish('init_soil_temperature', 'Currently only works with 5-layer soil')

    ALLOCATE(zrange(ndim), jmax(ndim), zdmax(ndim))
    ALLOCATE(zd(nsoil))

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

    ! layer depths
    ! TBD: should not be hardcoded here!
    zd(1) =            0.5_wp * 0.065_wp - 0.07_wp
    zd(2) = 0.065_wp + 0.5_wp * 0.254_wp - 0.07_wp
    zd(3) = 0.065_wp + 0.254_wp + 0.5_wp * 0.913_wp - 0.07_wp
    zd(4) = 0.065_wp + 0.254_wp + 0.913_wp + 0.5_wp * 2.902_wp - 0.07_wp
    zd(5) = 0.065_wp + 0.254_wp + 0.913_wp + 2.902_wp + 0.5_wp * 5.7_wp - 0.07_wp

    ! Initialize soil temperatures
    DO i=1,nsoil
      tsoil(:,i) = SUM(tclim(:,:), DIM=2) / 12._wp &
                 & + 0.5_wp * zrange(:) * EXP(-zd(i) / zsqrt) &
                 &   * COS(2._wp * pi * (day_in_year - zdmax(:)) / year_len - zd(i) / zsqrt)
    END DO

    ! Set surface temperature to uppermost soil temperature
    tsurf(:) = tsoil(:,1)

    DEALLOCATE(zrange, jmax, zdmax, zd)

  END SUBROUTINE init_soil_temperature

  PURE ELEMENTAL SUBROUTINE update_surf_down_soil ( delta_time, lstart  &
     , pqm1, ptsl, pspeed10, ptm1  &
     , pwl, pcvw, pwlmx &
     , psn, pcvs, psnc             &
     , pgrndcapc &
!!$     , paphm1, paphms1, ptte             &     ! ATTENTION!
     , pevapl,     pevapot, pevwsd                    &
     , prain, psnow                            &
!!$     , psnacl, psnmel&
     , pmlres &
!++mgs : add instantaneous snow melt for drydep calculation
     , psmelt, tte_corr   &
!--mgs

     )

    USE mo_jsb_physical_constants, ONLY: alf, cpd, cpvd1, tmelt, g=>grav

    !
    !     ------------------------------------------------------------------
    ! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
    ! pspeed10      Wind speed at 10m height [m/s] ... from 'vdiff'
    ! ptm1          Air temp [K] at timestep t-dt (filtered)
    ! pwl           Water content [m] in skin reservoir
    !               (vegetation and bare soil)
    ! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
    ! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
    !               of vegetation index and leaf area index)
    ! psn           Snow depth [m water equivalent] at the ground
    ! pcvs          Fractional snow cover (function of psn in *physc*)
    ! psnc          Snow depth [m water equivalent] at the canopy
    ! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
    ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
    ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
    ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
    ! prain         Totalrainfall [kg/m**2/s]
    ! psnow         Total snowfall [kg/m**2/s]
    ! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
    ! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
    !
    ! The following local variables represent the respective fluxes
    ! integrated over one timestep (delta_time) and divided by the density of
    ! water (rhoh2o). Units are m water equivalent.
    !
    ! zraind        Total rain
    ! zsnowd        Total snowfall
    ! zevttd        Total evaporation
    ! zevsnd        Sublimation from snow
    ! zevwld        Evaporation from the skin reservoir
    ! zsnmel        Snow/ice melt at land and glacier points
    ! zsncmelt      Snow melt in the canopy
    ! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
    ! pmlres        Residual melt water available for infiltration into the
    !               non-frozen soil after being intercepted by the
    !               skin reservoir
    ! tte_corr      Correction term for temperature tendencies accounting for snow melt
    ! psmelt        instantaneous snow melt rate (for dry deposition)
    !
    !       Rest folgt spaeter ....
    !
    !     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
    !              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
    !              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
    !              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
    !              DURING TIMESTEP DELTA_TIME.
    !
    !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
    !     MODIFIED BY
    !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
    !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
    !     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
    !                                 COUPLING BETWEEN LAND SURFACE
    !                                 AND ATMOSPHERE.
    !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
    !     MODIFIED BY M. ESCH        MPI - APR  1999
    !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
    !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
    !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir
    !                                                for snow changed
    !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
    !
    !     MODIFICATION
    !
    !     PURPOSE
    !
    !     INTERFACE.
    !
    !          *SURF* IS CALLED FROM *PHYSC*.
    !
    !     METHOD.
    !
    !     EXTERNALS.
    !
    !          NONE.
    !
    !     REFERENCE.
    !
    !          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
    !     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
    !
  !!$ TR  USE mo_param_switches ! only lsurf is used
  !!$ TR  USE mo_physc2
  !!$ TR  USE mo_constants
  !!$ TR  USE mo_vegetation ! only cvinter is used
  !!$ TR  USE mo_rad

    REAL(wp), INTENT(in) :: delta_time
    LOGICAL,  INTENT(in) :: lstart
    REAL(wp), INTENT(in) :: pqm1
    REAL(wp), INTENT(in) :: ptsl
    REAL(wp), INTENT(in) :: pspeed10
    REAL(wp), INTENT(in) :: ptm1
    REAL(wp), INTENT(in) :: pcvw
    REAL(wp), INTENT(in) :: pwlmx
    REAL(wp), INTENT(in) :: pcvs
    REAL(wp), INTENT(in) :: pgrndcapc
    REAL(wp), INTENT(in) :: pevapl
    REAL(wp), INTENT(in) :: pevapot
    REAL(wp), INTENT(in) :: prain
    REAL(wp), INTENT(in) :: psnow

    REAL(wp), INTENT(inout) :: pwl
    REAL(wp), INTENT(inout) :: psn
    REAL(wp), INTENT(inout) :: psnc
  !!$  REAL(wp), INTENT(inout) :: psnacl
  !!$  REAL(wp), INTENT(inout) :: psnmel
  !!$  REAL(wp), INTENT(inout) :: progl

    REAL(wp), INTENT(out)   :: pevwsd
    REAL(wp), INTENT(out)   :: tte_corr
    REAL(wp), INTENT(out)   :: pmlres
    REAL(wp), INTENT(out)   :: psmelt
    !
    !  local arrays
    !
    INTEGER :: jl
    REAL(wp) ::                                                             &
         zraind ,       zsnowd,         zevttd        &
         , zevsnd,        zevwld        &
         , zrogl         &
         , zsnmel,        zsn                                 &
         ,        zsncmelt
  !!$       , zdp,
    REAL(wp) :: zlfdcp
    !
    !  local scalars
    !
    REAL(wp) ::                                                             &
         zsmelt, zsnmlt,      &
         zmprcp,    &
         zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
    REAL(wp) :: zdtime, zrcp, zsncfac
    !
    !      Parameters
    !

    REAL(wp), PARAMETER :: cvinter = 0.25_wp  !!$ TR for JSBACH testing
    REAL(wp), PARAMETER :: cwlmax = 2.E-4_wp  !!$ TR for JSBACH testing

    zdtime = delta_time
    zc2=1.87E5_wp
    zc3=1.56E5_wp
    zsncfac=rhoh2o*g/zdtime

    zrogl    = 0._wp
    pmlres   = 0._wp
    zsnmel   = 0._wp
    zsncmelt = 0._wp

    ! Convert water fluxes to [m water equivalent * timestep]
    zraind=prain             *zdtime/rhoh2o                   !rain     [kg m-2 s-1] -> [m]
    zsnowd=psnow             *zdtime/rhoh2o                   !snowfall [kg m-2 s-1] -> [m]
    zevttd=pevapl                        *zdtime/rhoh2o    !evapotransp. [kg m-2 s-1] -> [m]
    zevsnd=pcvs*pevapot              *zdtime/rhoh2o    !pot. sublimation from snow [m]
    zevwld=(1._wp-pcvs)*pcvw*pevapot*zdtime/rhoh2o !pot. evaporation from skinreservoir [m]
    pevwsd=zevttd-zevsnd                               ! [m]

    zrcp=1._wp/(cpd*(1._wp+cpvd1*MAX(0.0_wp,pqm1)))
    zlfdcp=alf*zrcp

    !
    !     ------------------------------------------------------------------
    !
    !*    2.     Budgets of snow (canopy and ground) and glaciers
    !
    !*    2.1    Snow changes in the canopy (interception of snowfall,
    !            sublimation, melting, unloading due to wind)
    !
    zsn=zsnowd+zevsnd                       ! [m]
    zsncmax=MAX(0.0_wp,pwlmx-cwlmax)
    zmprcp=MIN(zsnowd*cvinter,zsncmax-psnc)
    zsncp=psnc+zmprcp
    zsnowd=zsnowd-zmprcp
    psnc=MIN(MAX(0._wp,zsncp+zevsnd),zsncmax)
    zevsnd=zevsnd-(psnc-zsncp)
    zexpt=MAX(0._wp,ptm1+3._wp-tmelt)*zdtime/zc2
    zexpw=pspeed10*zdtime/zc3
    zsncmelt=psnc*(1._wp-EXP(-zexpt))
    psnc=psnc-zsncmelt
    zsncwind=psnc*(1._wp-EXP(-zexpw))
    psnc=psnc-zsncwind
    zsnowd=zsnowd+zsncwind
    ! TBD - Question: Why is tte_corr only done for snowmelt in canopy?
    tte_corr=zsncmelt*zsncfac*zlfdcp ! unit: K kg m-1 s-3 = K Pa s-1 (will be devided by pressure difference in atmosphere)
    !
    !   pwl=pwl+zsncmelt see section 2.5
    !
    !
    !*    2.2    Snowfall and sublimation on land (excluding glaciers)
    !
    psn=psn+zsnowd+zevsnd
    IF (psn.LT.0._wp) THEN
      pevwsd=pevwsd+psn
      psn=0._wp
    END IF
    !
    !*    2.4    Snow and glacier melt
    !
    IF (.NOT. lstart) THEN
      IF (ptsl.GT.tmelt) THEN               !vg  temperature > 0 degC
        IF (psn.GT.0._wp) THEN                         !vg  snow cover at non glacier land:
          zsmelt=pgrndcapc*(ptsl-tmelt)/(alf*rhoh2o)      !vg  potential snow and ice melt [m]
          zsnmel=MIN(psn,zsmelt)                          !vg  actual snow and ice melt [m]
          psn=MAX(psn-zsnmel,0._wp)                   !vg  re-calculation of snow depth
        END IF
      END IF
    END IF
    !
    !*    2.5    Snow budget and meltwater (glacier-free land only)
    !
    pwl=pwl+zsncmelt                  ! Add melt water on canopy to skin reservoir
    zsnmlt=zsnmel+MAX(0._wp,pwl-pwlmx)   ! Excess water on canopy drips to ground and
    pwl=MIN(pwlmx,pwl)                ! skin reservoir on canopy is set to maximum
    pmlres=zsnmlt
    zsnmel=zsnmel+zsncmelt
    zsn=zsn-zsnmel
    !
    ! Accumulate water fluxes
    ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
    !       accumulation - multiply by zdtime
!!$            IF (.NOT. lpglac) psnacl = psnacl  +zsn * rhoh2o
!!$            psnmel = psnmel  +zsnmel * rhoh2o
!!$            IF (lpglac) progl  = progl   +zrogl  *rhoh2o

    psmelt = zsnmel * rhoh2o / zdtime

  END SUBROUTINE update_surf_down_soil

  PURE ELEMENTAL SUBROUTINE update_surf_down_glacier ( delta_time, lstart,  &
    & ptsl,            & !< in
    & pcvs,            & !< in 
    & pgrndcapc,       & !< in
    & pevapl,          & !< in
    & pevapot,         & !< in
    & prain, psnow,    & !< in
    & pgld,            & !< inout
    & pevwsd,          & !< out

!!$     , paphm1, paphms1, ptte             &     ! ATTENTION!
    & palac,           & !< out
    & psmelt,          & !< out
!!$     , psnacl, psnmel, progl &
    & tte_corr &
    & )

    USE mo_jsb_physical_constants, ONLY: alf, cpd, cpvd1, tmelt, g=>grav

  !
  !     ------------------------------------------------------------------
  ! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
  ! pcvs          Fractional snow cover (function of psn in *physc*)
  ! pgld          Glacier depth (including snow) [m water equivalent]
  ! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
  ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
  ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
  ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  ! prain         Totalrainfall [kg/m**2/s]
  ! psnow         Total snowfall [kg/m**2/s]
  ! palac         Precipitation minus sublimation at glacier points [kg/m**2/s]
  ! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
  ! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
  ! progl         Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
  !
  ! The following local variables represent the respective fluxes
  ! integrated over one timestep (delta_time) and divided by the density of
  ! water (rhoh2o). Units are m water equivalent.
  !
  ! zraind        Total rain
  ! zsnowd        Total snowfall
  ! zevttd        Total evaporation
  ! zevsnd        Sublimation from snow
  ! zrogl         Runoff at glacier points (rain and melt, but no calving)
  ! zsnmel        Snow/ice melt at land and glacier points
  ! tte_corr      Correction term for temperature tendencies accounting for snow melt
  ! psmelt        instantaneous snow melt rate (for dry deposition)
  !
  !       Rest folgt spaeter ....
  !
  !     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
  !              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
  !              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
  !              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
  !              DURING TIMESTEP DELTA_TIME.
  !
  !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
  !     MODIFIED BY
  !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
  !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
  !     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
  !                                 COUPLING BETWEEN LAND SURFACE
  !                                 AND ATMOSPHERE.
  !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
  !     MODIFIED BY M. ESCH        MPI - APR  1999
  !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
  !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
  !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir
  !                                                for snow changed
  !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
  !
  !     MODIFICATION
  !
  !     PURPOSE
  !
  !     INTERFACE.
  !
  !          *SURF* IS CALLED FROM *PHYSC*.
  !
  !     METHOD.
  !
  !     EXTERNALS.
  !
  !          NONE.
  !
  !     REFERENCE.
  !
  !          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
  !     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
  !
!!$ TR  USE mo_param_switches ! only lsurf is used
!!$ TR  USE mo_physc2
!!$ TR  USE mo_constants
!!$ TR  USE mo_vegetation ! only cvinter is used
!!$ TR  USE mo_rad

    REAL(wp), INTENT(in) :: delta_time
    LOGICAL,  INTENT(in) :: lstart
    REAL(wp), INTENT(in) :: &
      & ptsl, &
      & pcvs, &
      & pgrndcapc, &
      & pevapl, &
      & pevapot, &
      & prain, &
      & psnow

    REAL(wp), INTENT(inout) :: &
      & pgld
!!$  REAL(wp), INTENT(inout) :: psnacl
!!$  REAL(wp), INTENT(inout) :: psnmel
!!$  REAL(wp), INTENT(inout) :: progl

    REAL(wp), INTENT(out)   :: &
      & pevwsd, &
      & palac, &
      & tte_corr, &
      & psmelt
    !
    !  local arrays
    !
    REAL(wp) ::                                                             &
      & zraind, &
      & zsnowd, &
      & zevttd, &
      & zevsnd, &
      & zrogl,  &
      & zsnmel, &
!!$       , zdp,
      & zsmelt, &
      & zsnmlt
    REAL(wp) :: zdtime
    !
    !      Parameters
    !
    REAL(wp), PARAMETER :: cvinter = 0.25_wp  !!$ TR for JSBACH testing
    REAL(wp), PARAMETER :: cwlmax = 2.E-4_wp  !!$ TR for JSBACH testing

    zdtime = delta_time
    !
    !     ------------------------------------------------------------------
    !
    !*    1.     Convert water fluxes to [m water equivalent * timestep]
    !
    zrogl=0._wp
    zsnmel=0._wp
    zraind=prain          * zdtime/rhoh2o    !rain     [kg m-2 s-1] -> [m]
    zsnowd=psnow          * zdtime/rhoh2o    !snowfall [kg m-2 s-1] -> [m]
    zevttd=pevapl         * zdtime/rhoh2o    !evapotransp. [kg m-2 s-1] -> [m]
    zevsnd=pcvs * pevapot * zdtime/rhoh2o    !pot. sublimation from snow [m]
    !
    !     ------------------------------------------------------------------
    pevwsd = zevttd - zevsnd                               ! [m]
    !
    tte_corr = 0._wp      ! no snow melt correction for glaciers
    !
    !*    2.3    Snowfall and sublimation on glaciers and diagnostics
    !
    pgld = pgld + zsnowd + zevsnd      !vg  glacier depth [m] = glacier depth + snowfall + sublimation of snow
    palac = zraind + zsnowd + zevttd   !vg  P-E on glaciers [m] = rainfall + snowfall + evaporation (incl. sublim.) 
    zrogl = zraind                         !vg  glacier runoff = rainfall [m]
    !
    !*    2.4    Snow and glacier melt
    !
    IF (.NOT. lstart .AND. ptsl > tmelt) THEN
      zsnmel = pgrndcapc * (ptsl - tmelt) / (alf * rhoh2o)  ! snow and ice melt [m]
      pgld   = pgld - zsnmel                                ! adapt glacier depth to snow and ice melt
      zrogl  = zrogl + zsnmel                               ! add melt water to glacier runoff
    END IF
    !
    ! Accumulate water fluxes
    ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
    !       accumulation - multiply by zdtime
!!$            psnmel(jl) = psnmel(jl)  +zsnmel(jl) * rhoh2o
!!$            IF (lpglac(jl)) progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o
    palac = palac * rhoh2o/zdtime  ! [m] -> [kg m-2 s-1]

    psmelt = zsnmel * rhoh2o / zdtime

  END SUBROUTINE update_surf_down_glacier

  SUBROUTINE update_soiltemp(nidx,nsoil                        &
                           , delta_time                        &
                           , lstart                            &
                           , pts                               &
                           , psn                               &
                           , psodif, prgcgn                    &
                           , pgrndc, pgrndd                    &
                           , ptsoil                            &
                           , pgrndcapc, pgrndhflx              &
                           , ldglac                     &
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
    REAL(wp), Intent(in)     ::  delta_time       !! time step
    LOGICAL, INTENT(in)  :: lstart
    REAL(wp), Intent(in)     ::  pts(nidx)            !! surface temperature at top of soil [K]
    REAL(wp), Intent(in)     ::  psn(nidx)            !! equivalent snow depth [m water]
    REAL(wp), Intent(in)     ::  psodif(nidx)         !! soil temperature diffusivity [m^2/s]
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
    REAL(wp) :: cdel(nsoil),cmid(nsoil)
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
    cdel = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)         !! thicknesses of soil layers [m]
    cmid = (/0.0325_wp,0.192_wp,0.7755_wp,2.683_wp,6.984_wp/)       !! depth of mids of soil layers [m]
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
      zso_cond(:) = zso_capa(:) * psodif(:)
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

  END SUBROUTINE update_soiltemp

  SUBROUTINE update_surf_up ( kbdim, delta_time, nsoil &
       , pwl, pcvw, pwlmx &
       , pws, pwsmx &
       , ptsoil, pcvs          &
       , porostd &
       , pevapl,     pevapot, pevwsd                     &
       , prain                           &
       , prunoff, pdrain &
  !vg      , pros_hd,    pdrain_hd                   &
       , pmlres &
       , zevwld &
  !     Soil moisture layers, evap. fluxes over land
       , zdel, wsi, zwsfc, zwpwp           &
       , dzr, dzs, fksat , fmpot , bclapp, vpor     &
       , etrans, redevap                    &
  !!$     , vfc, vpwp &
       , spsi                           &
       , w_soil_overflow                       &
       , jllog                                     &
       , l_tpe_open                      &
       , l_tpe_closed                    &
       )
  !
  !     ------------------------------------------------------------------
  ! pwl           Water content [m] in skin reservoir
  !               (vegetation and bare soil)
  ! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
  !               of vegetation index and leaf area index)
  ! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
  ! pws           Soil water content [m]
  ! pwsmx         Water holding capacity [m] of the soil
  ! ptsoil        Temperature [K] in upper soil layer
  ! pcvs          Fractional snow cover (function of psn in *physc*)
  ! porostd       Subgrid standard diviation [m] used in runoff scheme
  ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
  ! etrans        Transpiration [kg/m**2/s]
  ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
  ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  !               BUT in this routine it is [m], such as zevwld.
  ! prain         Totalrainfall [kg/m**2/s]
  ! prunoff       Runoff [kg/m**2/s] at non-glacier points (not including drainage)
  ! pdrain        Drainage at non-glacier points [kg/m**2/s]
  !vg ! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
  !vg ! pdrain_hd     Drainage for HD-Model [m]
  ! 
  !
  ! The following local variables represent the respective fluxes
  ! integrated over one timestep (delta_time) and divided by the density of
  ! water (rhoh2o). Units are m water equivalent.
  !
  ! zraind        Total rain [m]
  ! zevwld        Evaporation from the skin reservoir [m]
  ! zros          Runoff (not including drainage) at non-glacier points [m]
  ! zdrain        Drainage at non-glacier points [m]
  ! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
  ! pmlres        Residual melt water available for infiltration into the
  !               non-frozen soil after being intercepted by the
  !               skin reservoir
  !
  !  ******** 5 layer scheme
  !     ***  kdeep = Number of soil layers = 5
  !     ***   zdel = Thickness of the 5 soil layers [m]
  !     ***  dzrsi = Rooted depth per soil layer (until rooting depth DZR) [m]
  !     ***   dzsi = Depth that can be saturated with water per soil layer
  !     ***          (until bedrock DZS) [m]
  !     ***    dzr = Rooting depth in [m]
  !     ***    dzs = Soil depth derived from textures in [m]
  !     ***  fksat = Saturated hydraulic conductivity: unit [m/s]
  !     ***  fmpot = Matrix potential [m]
  !     *** bclapp = Exponent B in Clapp and Hornberger
  !     ***   vpor = Volumetric soil porosity [m/m]   following Cosby et al.
  !     ***    vfc = Volumetric soil field capacity [m/m]
  !     ***   vpwp = Volumetric soil wilting point [m/m]
  !     ***   spsi = Soil Pore size distribution index
  !     ***  zwsat = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
  !     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
  !     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS dzs, WO DER BEDROCK BEGINNT.
  !     ***  zwsfc = Field capacity of soil layer I (analog zu ZWSAT) [m]
  !     ***  zwpwp = Soil wilting point of soil layer I (analog to ZWSAT) [m]
  !     ** aebsoil = Bare soil evaporation in [m], accumulated over TWODT/2
  !     ** aetrans = Transpiration in [m] accumulated over one time step
  !     ***          (calculated from etrans which is in [kg/m**2/s])
  !     *** redevap  = Diagnostic desired reduction of evapotranspiration
  !     ***            due to limited storage in soil layers [kg/m**2/s]
  !     *** wsi(I) = Soil moisture content in I. soil layer [m]
  !     ***          Is updated with the fluxes calculated in SOILHYD
  !     *** wsim1m(i) = Soil Moisture of previous time step in i. Soil layer [m]
  !     *** ZDRAIN = Drainage: Volumen in time step: unit [m]
  !
  !     *** jllog          relative Index of grid box for log output between 1 and kbdim
  !
  !       The rest is to be added later on ....
  !
  !     *SURF* - Updates land values of temperature, moisture and snow.
  !              Calculate fluxes of total rain, total snow and evapo-
  !              ration from the three reservoirs (SN, WS, WL)
  !              convert fluxes (kg/m**2*s) to changes of water levels (m)
  !              during timestep delta_time.
  !
  !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
  !     MODIFIED BY
  !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
  !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
  !                 J.-P. SCHULZ   MPI - 1997 : Implementation of implicit coupling between
  !                                             land surface and atmosphere.
  !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
  !     MODIFIED BY M. ESCH        MPI - APR  1999
  !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
  !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
  !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir
  !                                                for snow changed
  !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
  !     MODIFIED BY ??? R. SCHNUR?        JSBACH Version
  !
  !     MODIFIED BY S. HAGEMANN    MPI - February  2009
  !     *** Implementation of 5 layer soil hydrology scheme
  !
  !     MODIFICATION
  !
  !     PURPOSE
  !
  !     INTERFACE.
  !
  !          *SURF* is called from *PHYSC*.
  !
  !     METHOD.
  !
  !     EXTERNALS.
  !
  !          NONE.
  !
  !     REFERENCE.
  !
  !          See soil processes' part of the model's documentation for
  !     details about the mathematics of this routine.
  !

    USE mo_jsb_physical_constants, ONLY: tmelt, g=>grav

    INTEGER,  INTENT(in)    :: kbdim
    REAL(wp), INTENT(in)    :: delta_time
    INTEGER,  INTENT(in)    :: nsoil ! also used as switch for 5 layer calculation
    REAL(wp), INTENT(in)    :: pcvw(kbdim)
    REAL(wp), INTENT(in)    :: pwlmx(kbdim)
    REAL(wp), INTENT(in)    :: pwsmx(kbdim)
    REAL(wp), INTENT(in)    :: ptsoil(kbdim)
    REAL(wp), INTENT(in)    :: pcvs(kbdim)
    REAL(wp), INTENT(in)    :: porostd(kbdim)
    REAL(wp), INTENT(in)    :: pevapl(kbdim)      !! not used
    REAL(wp), INTENT(in)    :: pevapot(kbdim)
    REAL(wp), INTENT(in)    :: prain(kbdim)
    REAL(wp), INTENT(in)    :: pmlres(kbdim)

    REAL(wp), INTENT(inout) :: pwl(kbdim)
    REAL(wp), INTENT(inout) :: pws(kbdim)
    REAL(wp), INTENT(inout) :: pevwsd(kbdim)
    REAL(wp), INTENT(inout) :: prunoff(kbdim)
    REAL(wp), INTENT(inout) :: pdrain(kbdim)

  !vg   REAL(wp), INTENT(inout) :: pros_hd(kbdim)
  !vg   REAL(wp), INTENT(inout) :: pdrain_hd(kbdim)
    REAL(wp), INTENT(out)   :: zevwld(kbdim) ! From optional part
  !
  ! Soil moisture layers, evap. fluxes over land.
    REAL(wp), INTENT(in)   , optional :: zdel(nsoil)
    REAL(wp), INTENT(inout), optional :: wsi(kbdim,nsoil)
    REAL(wp), INTENT(in), optional :: zwsfc(kbdim,nsoil)
    REAL(wp), INTENT(in), optional :: zwpwp(kbdim,nsoil)
    REAL(wp), INTENT(in)   , optional :: dzr(kbdim)
    REAL(wp), INTENT(in)   , optional :: dzs(kbdim)
    REAL(wp), INTENT(in)   , optional :: fksat(kbdim)
    REAL(wp), INTENT(in)   , optional :: fmpot(kbdim)
    REAL(wp), INTENT(in)   , optional :: bclapp(kbdim)
    REAL(wp), INTENT(in)   , optional :: vpor(kbdim)
    REAL(wp), INTENT(in)   , optional :: etrans(kbdim)
    REAL(wp), INTENT(out)  , optional :: redevap(kbdim) ! reduction of evapotranspiration from SOILCHANGE
  !!$  REAL(wp), INTENT(in)   , optional :: vfc(kbdim)
  !!$  REAL(wp), INTENT(in)   , optional :: vpwp(kbdim)
    REAL(wp), INTENT(in)   , optional :: spsi(kbdim)
    REAL(wp), INTENT(inout), optional :: w_soil_overflow(kbdim) ! excess of water in the terraplanet setup (l_tpe_closed = true)

    INTEGER,  INTENT(inout), optional :: jllog      ! gridbox number for logoutput in SOILHYD
    LOGICAL,  INTENT(in),    OPTIONAL :: l_tpe_open   ! Switch for terraplanet setup, open water cycle
    LOGICAL,  INTENT(in),    OPTIONAL :: l_tpe_closed ! Switch for terraplanet setup, closed water cycle
  !
  !  local variables
  !
    INTEGER :: jl,jk
    REAL(wp) ::                                                             &
         zraind (kbdim)      &
       , zros(kbdim),          zdrain(kbdim)

    REAL(wp) ::                                                             &
         zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt,              &
         zmprcp, zwlp, zwptr, zwdtr, zwslim, zconw2, zconw3, zconw4,    &
         zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zprfl,           &
         zlysic, zwsup, zsncp, zexpt, zexpw, zsncmax, zsncwind
    REAL(wp) :: zdtime, zrcp, zsncfac
    ! variables for terraplanet
    REAL(wp) :: remoist
    LOGICAL  :: ltpe_open, ltpe_closed

    ! Variables necessary for evaporation fluxes and 5 soil layers
    REAL(wp), PARAMETER :: cvinter = 0.25  !!$ TR for JSBACH testing
    REAL(wp), PARAMETER :: zeps  = 1.e-10_wp
    REAL(wp), PARAMETER :: zfak1 = 1._wp           ! =1., Dummy for SOILCHANGE (not used)
    INTEGER,  PARAMETER :: kb1   = 1               ! = 1, Digit for 1st field value for transfer to subroutine
    INTEGER,  PARAMETER :: ilog  = 0               ! Switch for logoutput in SOILHYD
  !!!!  PARAMETER (jllog=918, ilog=3)        ! No logoutput for ilog=0

    REAL(wp) ::                                                         &
         dzrsi(kbdim, nsoil), dzsi(kbdim, nsoil),                       &
         zwsat(kbdim,nsoil),                        &
  !!$       zwpwp(kbdim,nsoil),                                            &
         wsim1m(kbdim, nsoil), aebsoil(kbdim), aetrans(kbdim)

    REAL(wp) :: zinfil(kbdim)    ! Infiltration for transfer to SOILCHANGE - was previously a local scalar
    INTEGER  :: isch             ! Switch to control SOILCHANGE

  !
  !      Parameters
  !
    zdtime  = delta_time
    zorvari = 100._wp
    ! TBD:
!#ifndef STANDALONE
!    zorvars = 1000._dp*64._dp/ngl
!#else
!    zorvars = 1000._dp*64._dp/nlat
!    IF (nlat == 2) zorvars = 1000._dp*64._dp/96  !vg: for site level jsbach: T63 value
!#endif
    zorvars = 1000._wp*64._wp/180._wp
    zdrmin  = 0.001_wp/(3600._wp*1000._wp)
    zdrmax  = 0.1_wp/(3600._wp*1000._wp)
    zdrexp  = 1.5_wp

    IF (PRESENT(l_tpe_open)) THEN
      ltpe_open = l_tpe_open
    ELSE
      ltpe_open = .false.
    ENDIF
    IF (PRESENT(l_tpe_closed)) THEN
      ltpe_closed = l_tpe_closed
    ELSE
      ltpe_closed = .false.
    ENDIF

  !     ------------------------------------------------------------------
  !
  !*    1.     Convert water fluxes to [m water equivalent * timestep]

    zros   = 0._wp
    zdrain = 0._wp
    zraind = prain * zdtime / rhoh2o   !vg: prain [kg/m2/s] -> zraind [m]
    zevwld = 0._wp
    DO jl=1,kbdim
  !!$      zros(jl) = 0._wp
  !!$      zdrain(jl) = 0._wp
  !!$      zraind(jl) = prain(jl) * zdtime/rhoh2o
        !vg: potential evaporation from skin reservoir
        !vg:  = (1-snow_fract_noforest) * skinres_fract * evap_pot  [kg/m2/s] -> [m]
        zevwld(jl) = (1._wp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
    ENDDO

  ! Initialization of working arrays for 5 soil layer scheme
    IF (nsoil == 5) THEN
      CALL SoilDef(kb1, kbdim, kbdim, nsoil,               &
                   dzr, dzs, vpor &
  !!$                 , vfc, vpwp &
                   , zdel, &
                   dzrsi, dzsi, zwsat &
  !!$                , zwsfc, zwpwp &
                   )

  !   *** Unit conversion for Transpiration [kg/m**2/s] --> [m]
      aetrans(1:kbdim) = etrans(1:kbdim)  * zdtime / rhoh2o
      redevap(1:kbdim) = 0._wp

    ENDIF

  !     ------------------------------------------------------------------
  !*  4.     Water budget
  !
  !*  4.1    Skin reservoir (vegetation and bare soil)
  !
      DO jl=1,kbdim
  !
  !*      4.1.1  Interception of rain
  !
          zmprcp     = MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))  !vg  amount of water going into interception reservoir
          zwlp       = pwl(jl)    + zmprcp                        !vg  update interception reservoir
          zraind(jl) = zraind(jl) - zmprcp                        !vg  remaining water reaching the soil 
  !
  !*      4.1.2  Evaporation or dew collection
  !
          pwl(jl)    = MIN(MAX(0._wp,zwlp+zevwld(jl)),pwlmx(jl))  !vg  update interception reservoir with respect to pot. evaporation
          pevwsd(jl) = pevwsd(jl)-(pwl(jl)-zwlp)                  !vg  reduce amount of evaporation accordingly -> evapor. from soil
  !HAG
  !       Get sure that ESKIN comprises only amounts that really change the skin reservoir
          zevwld(jl) = pwl(jl)-zwlp ! This is supposedly the correct way     !vg  calculate actual evaporation from pot. evaporation
  !        if (nsoil==5) zevwld(jl) = pwl(jl)-zwlp ! Use this for compatibility with older jsbach revisions
      ENDDO

  !   Calculation of Bare Soil Evaporation
  !   Unit of pevwsd is [m]
      IF (nsoil == 5) aebsoil(1:kbdim) = pevwsd(1:kbdim) - aetrans(1:kbdim)  !vg  evap_baresoil = evapotrans_soil - transpiration
  !
  !*  4.2    Soil reservoir
  !
      DO jl=1,kbdim

        zinfil(jl) = 0._wp
          zwslim = 0.05_wp*pwsmx(jl)

          zroeff = MAX(0._wp, porostd(jl)-zorvari)   &
                       /(porostd(jl)+zorvars)
          zbws   = MAX(MIN(zroeff,0.5_wp),0.01_wp)
          zb1    = 1._wp+zbws
          zbm    = 1._wp/zb1
          zconw1 = pwsmx(jl)*zb1
          zlyeps = 0._wp
          zvol   = 0._wp

  !*    4.2.1  Surface runoff, infiltration and evaporation from soil

          ! Positive Enoskin_nosnow is added to throughfall if the old bucket is used

  !         ******* Positive EBSoil and Etrans are added to throughfall separately
  !         ******* Negative fluxes change soil moisture in SOILCHANGE.
            zprfl = pmlres(jl)+zraind(jl)                           !vg  water reaching the soil = meltwater + rain reaching the soil
            IF (aebsoil(jl) >= 0.0_wp) zprfl = zprfl + aebsoil(jl)  !vg  add bare soil evaporation
            IF (aetrans(jl) >= 0.0_wp) zprfl = zprfl + aetrans(jl)  !vg  add transpiration
          
          IF (ltpe_closed) THEN   !vg  only used with terra planet: runoff swiched off
            zros(jl)   = 0._wp
            zinfil(jl) = zprfl
          ELSE
            IF (ptsoil(jl) < tmelt) THEN
              zros(jl) = zprfl            !vg  no water infiltrated if uppermost soil layer is frozen -> runoff
            ELSE
              IF (zprfl > 0._wp .AND. pws(jl) > zwslim) THEN   !vg  water reaches soil and soil water content above 5%-limit
                IF (pws(jl) > pwsmx(jl)) THEN                  !vg  soil water content above water holding capacity
                  zlyeps = pws(jl)-pwsmx(jl)                   !vg     excess water
                ELSE
                  zlyeps = 0._wp                               !vg     no excess water
                END IF
                zlysic     = MIN(1._wp,(pws(jl)-zlyeps)/pwsmx(jl))        !vg  relative soil moisture (MIN() needed??)
                zvol       = (1._wp-zlysic)**zbm-zprfl/zconw1             !vg  ??
                zros(jl)   = zprfl-(pwsmx(jl)-pws(jl))                    !vg  runoff = remaining water, when pwsmx reached
                IF (zvol > 0._wp) zros(jl) = zros(jl)+pwsmx(jl)*zvol**zb1 !vg  ??
                zros(jl)   = MAX(zros(jl),0._wp)                          !vg  no negative runoff (code would be clearer, if excess 
                                                                          !vg       water was considered two lines above)
                zinfil(jl) = zprfl-zros(jl)                               !vg  water not going into runoff is infiltrated
              ELSE
                zros(jl)   = 0._wp   !vg  no runoff as no water reaches the soil or soil is too dry (below 5%-limit)
                zinfil(jl) = zprfl
              END IF
            END IF
          END IF
  !
  !*    4.2.2  Drainage and total runoff
  !
      ENDDO
  !
        !  Change of soil moisture in 5 layers
        !  change belongs to one time step = zdtime

        ! Soil moisture of previous time step: wsim1m
        wsim1m(:,:) = wsi(:,:)

        ! Zero initialization for drainage field before SOILCHANGE
        zdrain(1:kbdim) = 0._wp

        ! Soilchange expects aetrans and aebsoil in [m]
        isch = 1     ! Changes by evaporation fluxes
                     ! Soilchange only uses negative evaporation fluxes as positive fluxes,
                     ! e.g. by dew formation, are part of the throughfall zprfl.
        CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
              isch, zdtime, zfak1,                    &
              wsi , wsim1m,                           &
              dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
              zinfil, zdrain, aetrans, aebsoil,       &
              redevap)
        ! Unit of zinfil is [m] but soilchange expects m/s
        ! For numerical reasons, half of the water is infiltrated
        ! before SOILHYD, half if it after
        zinfil(1:kbdim) = zinfil(1:kbdim) / zdtime / 2._wp

        isch = 2     ! Changes by infiltration

        CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
              isch, zdtime, zfak1,                    &
              wsi , wsim1m,                           &
              dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
              zinfil, zdrain, aetrans, aebsoil,       &
              redevap, w_soil_overflow, ltpe_closed)

        ! Drainage calculation [m] with 5 layer soil scheme
        CALL SoilHyd(kb1, kbdim, kbdim, nsoil,      &
             zdtime, ilog, jllog,                   &
             wsi, dzsi, zwsat, zwsfc,               &
             fksat, vpor, bclapp, fmpot,            &
             zdrain, spsi, zwpwp, zdel,             &
             w_soil_overflow, ltpe_closed)

        ! Second half of the water is infiltrated after SOILHYD
        isch = 2     ! Changes by infiltration
        CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
              isch, zdtime, zfak1,                    &
              wsi , wsim1m,                           &
              dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
              zinfil, zdrain, aetrans, aebsoil,       &
              redevap, w_soil_overflow, ltpe_closed)

        ! Does overflow occur due to infiltration in the layers? --> added to drainage
        DO jk=1, nsoil
          DO jl=kb1, kbdim
            IF (dzsi(jl,1) > 0) THEN                                     !vg  there is soil (above bedrock)
              IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN                        !vg  soil moisture above field capacity
                IF (ltpe_closed) THEN
                  ! modification for TPE
                  ! put drainage in w_soil_overflow pool
                  w_soil_overflow(jl) = w_soil_overflow(jl) + wsi(jl,jk) - zwsfc(jl,jk)
                ELSE
                  zdrain(jl) = zdrain(jl) + wsi(jl,jk) - zwsfc(jl,jk)    !vg  excess moisture goes to drainage
                ENDIF
                wsi(jl,jk) = zwsfc(jl,jk)                                !vg  soil moisture reduced accordingly
              ELSE ! wsi <= zwsfc
                IF (ltpe_closed) THEN
                  ! modification for TPE
                  ! refill soil moisture from w_soil_overflow pool
                  remoist = MIN(w_soil_overflow(jl), zwsfc(jl,jk) - wsi(jl,jk))
                  wsi(jl,jk) = wsi(jl,jk) + remoist
                  w_soil_overflow(jl) = w_soil_overflow(jl) - remoist
                ENDIF
              ENDIF
              IF (ltpe_open) THEN
                ! modification for TPE
                ! keep the two lowest soil layers always very wet 
                IF (jk == nsoil - 1 .OR. jk == nsoil) THEN
                  wsi(jl,jk) = MAX(wsi(jl,jk),zwsfc(jl,jk) * 0.9_wp)
                END IF
              END IF
            ENDIF
          ENDDO
        ENDDO

        ! WS must be re-derived from wsi after SOILHYD
        pws(1:kbdim) = 0.0_wp
        DO jk=1, nsoil
          DO jl=kb1, kbdim
              IF (dzrsi(jl,jk) >= zdel(jk)) THEN
                pws(jl) = pws(jl) + wsi(jl,jk)
              ELSE IF (dzrsi(jl,jk) > 0._wp) THEN
                pws(jl) = pws(jl) + wsi(jl,jk) * dzrsi(jl,jk)/dzsi(jl,jk)
              ENDIF
          ENDDO
        ENDDO

      DO jl=1,kbdim
          !  [m] -> [kg/m2/s]
          prunoff(jl) = zros(jl)   *rhoh2o/zdtime
          pdrain(jl)  = zdrain(jl) *rhoh2o/zdtime
      END DO

  END SUBROUTINE update_surf_up

  SUBROUTINE soilhyd(kidia, kfdia, klon, kdeep, zdiagt, ilog, jllog,  &
      wsi, dzsi, zwsat, zwsfc, fksat, vpor, bclapp, fmpot,  &
      zdrain, spsi, zwpwp, zdel, w_soil_overflow, ltpe_closed )
  !     **************************************************************************

  !     ******* SOIL HYDROLOGY - Routine to be build into REMO
  !     *******   to calculate soil hydrology for 5 soil layers
  !     ******* SOIL HYDROLOGY = percolation and vertical diffusion of water
  !     ******* Type of diffusiont:  Richtmeyr/Morton

  !     *** PROGRAMMING AND DEVELOPMENT: STEFAN HAGEMANN
  !     ****
  !     *** VS. 1.0 - Test version for soil processes without run-off and others
  !              percolation dependent of IPERC type
  !              REMO calculations in [m], but output is in [mm]/output interval.
  !     ***
  !     *** VS. 1.1 - Version for REMO
  !     ***           Units must be made consistent.

  !     *** VS. 1.2 - Limitation to IPERC=2, IDIFF=2
  !     ***           --> alternatives removed

  !     *** VS. 1.3 - Last loop splitted and loop order changed to get longest loop inside

  !     *** VS. 1.4 - September 2007 - HAG
  !     ***           Loop from Numerical Recipes changed to get longest loop inside.
  !     ***               Therefore were ZTRI(KDEEP) and ZDC changed into fields
  !     ***               ZTRI(KLON,KDEEP) and ZDC(KLON,KDEEP).
  !     ***           Inconsistency for percolation overflow of layers eliminated
  !     ***               now overflow exists if WSI > ZWSFC (instead of ZWSAT as used before)
  !     ***           Test for 1E-10. was changed due to the calculation accuracy of the NEC
  !     ***               at WSMX-1.E-10 < ZEPS (1E-15)
  !     ***
  !     ***
  !     *** Fluxes are always calculated to be valid two time steps ahead.
  !     *** Therefore the fluxes must be divided by 2 for accumulation.
  !     *** The calculations themselves uses volumes, not fluxes.
  !     *** The flux is calculated so that later the unit of the flux is the accumulated volume
  !     *** per time step.
  !     ***
  !     *** VS. 1.5 - December 2007 - HAG
  !     ***           Use of ZWSFC instead of ZWSAT for percolation

  !     ***           Juli 2008 - S. HAGEMANN (R. PODZUN)
  !     ***           Removal of all occurances of ZVFC in soilhyd.f and surf.f
  !     ***           UNSAFE! Use may cause problems. on nec no obvious errors were caused.

  !     *** VS. 1.6 - March 2009 - HAG
  !     ***           Introduction of maximum percolation rate and solving the percolation
  !     ***           equation with midpoint method instead of euler method.

  !     *** VS. 1.7 - March-May 2009 - HAG
  !     ***           Implementation of van genuchetn equation for percolation.
  !     ***           preliminary using iperc to allow switching to previuos
  !     ***               Clapp & Hornberger method.
  !     ***           Allow ECHAM drainage for lowest layer above bedrock (IDRECH).
  !     ***           Allow reduction of hydraulic conductivity with depth for percolation.

  !     *** VS. 1.8 - October 2009 - HAG
  !     ***           Percolation overflow is not allowed anymore from layers below the surface
  !     ***               layer. This has led to a too strong drying
  !     ***               of soil layers in wet areas, especially from the second lowest layer
  !     ***               above the bedrock. Now it assumed that water piles upwards from
  !     ***               SATURATEDsoil5layer.f90 layers, very similar to a rising water table.
  !     ***               Whether this may actually be used to describe a water table, has to be investigated.
  !     ***               But this can be also done offline.

  !     *** VS. 1.9 - June 2010 - HAG
  !     ***            Initialization of zdrain with zero must happen in calling routine before
  !     ***                soilchange is called.
  !
  !     *** vs. 1.10 - Jay 2011 - HAG
  !     ***            Code converted using TO_F90 by Alan Miller
  !     ***            plus REAL --> REAL(wp) and fine tuning afterwards
  !
  !     SOILHYD
  !     ***  FKSAT = Saturated hydraulic conductivity: unit [m/s] -> Must be checked in ECMGBDT.F
  !     ***  FMPOT = Matrix potential [m]
  !     *** BCLAPP = Exponent B in Clapp and Hornberger
  !     ***   VPOR = Volumetric soil porosity [m/m]
  !     ***   SPSI = Soil pore size distribution index used in Van Genuchten

  !     *** ZDIAGT = Time step: unit [s]   ##### org. was [day]

  !     *** ZDRAIN = Drainage: Volumen during time step: unit [m]
  !     ***  ZWSAT = Maximum storage in soil layer I (porosity) [m]
  !     ***          Important for the take up of the infiltration in the soil layers
  !     ***          Defined down to the depth of the soil - that is where the bedrock begins.
  !     ***  ZWSFC = Field capacity for soil layer I (analog to ZWSAT) [m]
  !     ***  ZWPWP = Texture depending wilting point of the soil layer I (analog to ZWSAT) [m]
  !     ***
  !     ***  ZPERC = Downward flux from soil layer I: First in m/day, then in [m].
  !     ***          First because of gravitational drainage.
  !     ***  ZDIFF = Diffusion of water between soil layer i and layer i+1 [m^2/day]
  !     ***   ZMVG = Exponent M in Van Genuchten scheme = SPSI / (SPSI+1)
  !     ***   ZDEL = Thickness of the 5 soil layers [m]
  !     *** ZKREDU = Reduction factor of hydraulic conductivity with depth for percolation
  !     ***  IREDU = Reduction of conductivity with depth yes/no (1/0)
  !     ***
  !     ***  NSOIL - Number of soil layers
  !     *** WSI(I) = Soil moisture content in I. soil layer [m]
  !     ***          Updated through the fluxes calculated in SOILHYD
  !     ***  ILOG =  1 --> Log output to StdOut with keword SOILFLUX for grid cell JLLOG
  !     ***                5 * DIFF. FLUXES, 5 PERC. FLUXES, ZDRAIN
  !     ***          2 --> Log output to StdOut with keyword SOILHYD1 for grid cell JLLOG

  !     ***  IPERC = Type of percolation
  !     ***          1  ECHAM4 drainage (removed, not available)
  !     ***          2  Clapp & Hornberger
  !     ***          3  Van Genuchten methode 2, after these
  !     ***
  !     ***  IDRECH = ECHAM drainage
  !     ***          0   No, drainage is percolation from lowest layer above bedrock.
  !     ***          1   Yes, for lowest layer above bedrock (replaces percolation)
  !     ***          2   For all layers above bedrock (in addition to percolation, except
  !     ***              for lowest layer above bedrock where it replaces percolation.

  INTEGER,  INTENT(IN)                     :: kidia
  INTEGER,  INTENT(IN)                     :: kfdia
  INTEGER,  INTENT(IN)                     :: klon
  INTEGER,  INTENT(IN)                     :: kdeep
  REAL(wp), INTENT(IN)                     :: zdiagt
  INTEGER,  INTENT(IN)                     :: ilog
  INTEGER,  INTENT(INOUT)                  :: jllog
  REAL(wp), INTENT(INOUT)                  :: wsi(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: dzsi(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: zwsat(klon,kdeep)
  REAL(wp), INTENT(IN)                     :: zwsfc(klon,kdeep)
  REAL(wp), INTENT(IN)                     :: fksat(klon)
  REAL(wp), INTENT(IN)                     :: vpor(klon)
  REAL(wp), INTENT(IN)                     :: bclapp(klon)
  REAL(wp), INTENT(IN)                     :: fmpot(klon)
  REAL(wp), INTENT(INOUT)                  :: zdrain(klon)
  REAL(wp), INTENT(IN)                     :: spsi(klon)
  REAL(wp), INTENT(IN)                     :: zwpwp(klon,kdeep)
  REAL(wp), INTENT(IN)                     :: zdel(kdeep)
  REAL(wp), INTENT(INOUT), OPTIONAL        :: w_soil_overflow(klon)
  LOGICAL,  INTENT(IN),    OPTIONAL        :: ltpe_closed

  REAL(wp), PARAMETER :: zeps   = 1.e-15_wp
  INTEGER,  PARAMETER :: iperc  = 3
  INTEGER,  PARAMETER :: idrech = 2
  INTEGER,  PARAMETER :: iredu  = 0
  REAL(wp), PARAMETER :: credu  = 1.0_wp

  REAL(wp) :: zdum, zdum2, zdum3
  REAL(wp) :: zperc(klon,kdeep), zdum4, zctsat
  REAL(wp) :: zmvg(klon), zwrel
  REAL(wp) :: zdiflog(10), wslog(kdeep)
  !     *** ECHAM-Drainage-Parameter
  REAL(wp) :: zdrmin, zdrmax, zdrexp, zwdmin, zwdcri

  !     *** Diffusion
  REAL(wp) :: zda(klon,kdeep), zdb(klon,kdeep), zdiff(klon,kdeep-1)
  REAL(wp) :: ztri(klon,kdeep), zdc(klon,kdeep)

  REAL(wp) :: zpmax, zp2         ! Maximum percolation rate
  REAL(wp) :: zkredu(kdeep)      ! KSAT reduction factor with depth

  !     *** Dummy variables, running indices, etc.
  INTEGER  :: i, jl

  LOGICAL  :: l_tpe_closed

    IF (PRESENT(ltpe_closed)) THEN
      l_tpe_closed = ltpe_closed
    ELSE
      l_tpe_closed = .false.
    END IF

    zwdmin = 0.05_wp
    zwdcri = 0.9_wp
    zdrmin = 0.001_wp / (3600._wp*1000._wp)
    zdrmax = 0.1_wp / (3600._wp*1000._wp)
    zdrexp = 1.5_wp
    IF (iredu == 1) THEN
      zkredu(1) = (1._wp + EXP(-credu*zdel(1)) )/2._wp
      DO i=2, kdeep
        zkredu(i) = (EXP(-credu*SUM(zdel(1:i-1))) +  &
            EXP(-credu*SUM(zdel(1:i))) )/2._wp
      END DO
    ELSE
      zkredu(1:kdeep) = 1._wp
    END IF

  !     ***** Initialization of local fields
    zperc(kidia:kfdia,1:kdeep) = 0._wp
    zmvg(kidia:kfdia) = spsi(kidia:kfdia) / (spsi(kidia:kfdia)+1._wp)

  !     *** Percolation = Gravitational drainage - calculated from bottom upward.
  !     *** Maximum percoplation due to numeric instability at large time steps
  !     *** and small spatial distances.
    IF (iperc == 2) THEN          ! Clapp & Hornberger-formulation
      DO i=1, kdeep
        DO jl=kidia, kfdia
  !cc        IF (DZSI(JL,I).GT.0) THEN
          IF (zwsfc(jl,i)-1.e-10_wp > zeps) THEN

  !         *** MAXIMUM PERCOLATION RATE ZPMAX
            zdum = ( zwsfc(jl,i) -  &
                fksat(jl)*zkredu(i) * zdiagt - zwpwp(jl,i) ) /  &
                ( zwsfc(jl,i) - zwpwp(jl,i))
            IF (zdum > zeps) THEN
              zp2 = fksat(jl)*zkredu(i) * zdum ** (2._wp*bclapp(jl)+3._wp)
            ELSE
              zp2=0._wp
            END IF
            zpmax  = MAX( (fksat(jl)*zkredu(i) + zp2)/2._wp , 0._wp)
            zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._wp)  &
                / (fksat(jl)*zkredu(i)) /zdiagt , 1._wp)
            zpmax  = zctsat * fksat(jl)*zkredu(i) + (1._wp-zctsat) * zpmax

            zwrel  = MAX(wsi(jl,i)-zwpwp(jl,i), 0._wp) / (zwsfc(jl,i)-zwpwp(jl,i))
            zwrel  = MIN(zwrel, 1._wp)

            IF (zwrel-1.e-10_wp > zeps) THEN
              zp2   = fksat(jl)*zkredu(i) * zdiagt / 2._wp * zwrel ** (2._wp*bclapp(jl)+3._wp)
              zdum2 = MAX( MAX(wsi(jl,i)-zp2,0._wp)  -zwpwp(jl,i), 0._wp)/  &
                  (zwsfc(jl,i)-zwpwp(jl,i))
              zdum2 = MIN(zdum2, 1._wp)
              zperc(jl,i) = fksat(jl)*zkredu(i) * zdum2 ** (2._wp*bclapp(jl)+3._wp)
              zperc(jl,i) = MIN(zperc(jl,i) , zpmax )
            ELSE
              zperc(jl,i) = 0._wp
            END IF
          END IF
        END DO
      END DO
    ELSE IF (iperc == 3) THEN          ! Van Genuchten-formulation
      DO i = 1, kdeep
        DO jl=kidia, kfdia
          IF (zwsfc(jl,i)-1.e-10_wp > zeps) THEN

  !         *** MAXIMUM PERCOLATION RATE ZPMAX
            zdum = MAX(zwsfc(jl,i)-zwpwp(jl,i) -  &
                fksat(jl)*zkredu(i)*zdiagt, 0._wp) / (zwsfc(jl,i)-zwpwp(jl,i))
            IF (zdum > zeps) THEN
              zp2 = fksat(jl)*zkredu(i) * SQRT(zdum) * ( 1._wp - ( 1._wp -  &
                  zdum**(1._wp/zmvg(jl)) )**zmvg(jl) )**2._wp
            ELSE
              zp2=0._wp
            END IF
            zpmax = MAX( (fksat(jl)*zkredu(i) + zp2)/2._wp , 0._wp)
            zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._wp)  &
                / (fksat(jl)*zkredu(i)) /zdiagt , 1._wp)
            zpmax = zctsat * fksat(jl)*zkredu(i) + (1._wp-zctsat) * zpmax

  !         *** Midpoint-Method

            zwrel = MAX(wsi(jl,i)-zwpwp(jl,i), 0._wp) / (zwsfc(jl,i)-zwpwp(jl,i))
            zwrel = MIN(zwrel, 1._wp)

            zdum2 = fksat(jl)*zkredu(i) * SQRT(zwrel) * ( 1._wp - ( 1._wp -  &
                zwrel**(1._wp/zmvg(jl)) )**zmvg(jl) )**2._wp
            zdum3 = wsi(jl,i) - MIN(zdum2,zpmax) * zdiagt / 2._wp

            zdum4 = MAX(zdum3-zwpwp(jl,i), 0._wp) / (zwsfc(jl,i)-zwpwp(jl,i))
            zdum4 = MIN(zdum4, 1._wp)

            zperc(jl,i) = fksat(jl)*zkredu(i) * SQRT(zdum4) * ( 1._wp -  &
                ( 1._wp - zdum4**(1._wp/zmvg(jl)) )**zmvg(jl) )**2._wp

            zperc(jl,i) = MIN(zperc(jl,i) , zpmax)
          ELSE
            zperc(jl,i) = 0._wp
          END IF
        END DO
      END DO
    END IF

  ! *** Using ECHAM4 drainage calculation? for the lowest layer above bedrock
  ! *** this replaces the percolation value.

    IF (idrech == 1 .OR. idrech == 2) THEN
      DO i = 1, kdeep-1
        DO jl=kidia, kfdia
          IF (dzsi(jl,1) > 0) THEN
            IF (dzsi(jl,i) > 0) THEN
              IF (dzsi(jl,i+1) <= 0) THEN           ! *** Lowest layer above bedrock

  !cc              ZDUM = ZWSFC(JL,I) * ZWDMIN
                zdum = zwpwp(jl,i)
                zdum2 = zwsfc(jl,i) * zwdcri
                IF (wsi(jl,i) <= zdum) THEN
                  zperc(jl,i) = 0._wp
                ELSE
                  zdum3 = zdrmin * (wsi(jl,i)-zdum) / (zwsfc(jl,i)-zdum)
                  IF (wsi(jl,i) > zdum2) THEN
                    zdum4 = (zdrmax-zdrmin)*( (wsi(jl,i)-zdum2)  &
                        / (zwsat(jl,i) - zdum2) )**zdrexp
                    zdum3 = zdum3 + MIN(zdum4, (wsi(jl,i)-zdum2)/zdiagt)
                  END IF
                  zperc(jl,i) = MIN(zdum3, (wsi(jl,i)-zdum)/zdiagt)
                END IF
              ELSE IF (idrech == 2) THEN            ! *** Other layers for IDRECH=2

  !             *** Drainage is calculated with percolation removed from WSI
                zp2 = wsi(jl,i) - zperc(jl,i) * zdiagt
                zdum = zwpwp(jl,i)
                zdum2 = zwsfc(jl,i) * zwdcri
                IF (zp2 > zdum) THEN         ! Drainage occurs only above wilting point
                  zdum3 = zdrmin * (zp2-zdum) / (zwsfc(jl,i)-zdum)
                  IF (zp2 > zdum2) THEN
                    zdum4 = (zdrmax-zdrmin)*( (zp2-zdum2)  &
                        / (zwsat(jl,i) - zdum2) )**zdrexp
                    zdum3 = zdum3 + MIN(zdum4, (zp2-zdum2)/zdiagt)
                  END IF

  !               *** Except for the lowest layer above bedrock, drainage is removed from
  !               *** soil moisture layer wsi before diffusion occurs
  !               *** here, drainage is used as accumulated value over ZDIAGT
                  zdum4 = MIN(zdum3 * zdiagt, (zp2-zdum))
                  IF (.NOT. l_tpe_closed) THEN ! no drainage for TPE case
                    wsi(jl,i) = wsi(jl,i) - zdum4
                    zdrain(jl) = zdrain(jl) + zdum4
                  END IF
                END IF
              END IF                                 ! Lowest layer / IDRECH=2
            END IF
          END IF
        END DO
      END DO

  !   *** Lowest layer (=5)
      DO jl=kidia, kfdia
        IF (dzsi(jl,kdeep) > 0) THEN
  !cc          ZDUM = ZWSFC(JL,KDEEP) * ZWDMIN
          zdum = zwpwp(jl,kdeep)
          zdum2 = zwsfc(jl,kdeep) * zwdcri
          IF (wsi(jl,kdeep) <= zdum) THEN
            zperc(jl,kdeep) = 0._wp
          ELSE
            zdum3 = zdrmin * (wsi(jl,kdeep)-zdum) / (zwsfc(jl,kdeep)-zdum)
            IF (wsi(jl,kdeep) > zdum2) THEN
              zdum4 = (zdrmax-zdrmin)*( (wsi(jl,kdeep)-zdum2)  &
                  / (zwsat(jl,kdeep) - zdum2) )**zdrexp
              zdum3= zdum3 + MIN(zdum4, (wsi(jl,kdeep)-zdum2)/zdiagt)
            END IF
            zperc(jl,kdeep) = MIN(zdum3, (wsi(jl,kdeep)-zdum)/zdiagt)
          END IF
        END IF
      END DO
    END IF

  ! *** Original unit of percolation fluxes (=unit of FKSAT) was [m/day]
  ! *** Now it's accumulatied [m] over ZDIAGE (unit of ZDIAGT=[s]
  ! *** --> Take time step into account : * ZDIAGT
    zperc(kidia:kfdia,1:kdeep) = zperc(kidia:kfdia,1:kdeep) * zdiagt

  ! *** Diffusion of water *******************************************

  ! *** Numerical recipes/Richtmyer & -Morton diffusion  *************

    DO i = 1, kdeep-1
      DO jl=kidia, kfdia
        IF (dzsi(jl,i) > 0) THEN

  !       ***  Soil moisture diffusivity [m^2/day]
  !       ***  method c: diffusivity of both levels is averaged
  !       ***            to calculate diffusivity between layers
          IF (wsi(jl,i+1) > 0) THEN
            zdum = wsi(jl,i+1)/dzsi(jl,i+1)

  !         *** Calculating the diffusivity of layer i+1
            zdum4 = zdum / vpor(jl)
            IF (zdum4-1.e-10_wp > zeps) THEN
              zdum2 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                  zdum * zdum4 ** (bclapp(jl)+3.)
            ELSE
              zdum2 = 0._wp
            END IF
          ELSE
            zdum2=0._wp
          END IF
          IF (wsi(jl,i) > 0._wp) THEN
            zdum = wsi(jl,i)/dzsi(jl,i)

  !         *** Calculating the diffusivity of layer i
            zdum4 = zdum / vpor(jl)
            IF (zdum4-1.e-10_wp > zeps) THEN
              zdum3 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                  zdum *  zdum4 ** (bclapp(jl)+3._wp)
            ELSE
              zdum3 = 0._wp
            END IF
          ELSE
            zdum3=0._wp
          END IF

  !       *** Calculating the diffusivity at the bottom of layer i
          zdiff(jl,i) = (zdum2+zdum3)*0.5_wp
        END IF
      END DO
    END DO

  ! *** Calculation of diffusion coefficients
  ! *** Deepest layer KDEEP = 5
    DO jl=kidia, kfdia
      zda(jl, kdeep) = 0._wp
      IF (dzsi(jl,kdeep) > 0._wp) THEN

  !     *** [ZDIAGT] = s
        zdb(jl, kdeep) = zdiff(jl,kdeep-1) *zdiagt /dzsi(jl,kdeep)  &
            / (dzsi(jl,kdeep)+dzsi(jl,kdeep-1)) * 2._wp
      ELSE
        zdb(jl, kdeep) = 0._wp
      END IF
    END DO

  ! *** Layers i = 4,3,2
    DO i=kdeep-1, 2,-1
      DO jl=kidia, kfdia
        IF (dzsi(jl,i) > 0) THEN
          IF (dzsi(jl,i+1) > 0._wp) THEN
            zda(jl,i) = zdiff(jl,i) * zdiagt / dzsi(jl,i) /  &
                (dzsi(jl,i)+dzsi(jl,i+1)) * 2._wp
          ELSE
            zda(jl,i) = 0._wp
          END IF
          zdb(jl,i) = zdiff(jl,i-1) * zdiagt / dzsi(jl,i) /  &
              (dzsi(jl,i)+dzsi(jl,i-1)) * 2._wp
        ELSE
          zda(jl, i) = 0._wp
          zdb(jl, i) = 0._wp
        END IF
      END DO
    END DO

  ! *** Layer 1
    DO jl=kidia, kfdia
      zdb(jl, 1) = 0._wp
      IF (dzsi(jl,1) > 0) THEN
        IF (dzsi(jl,2) > 0) THEN
          zda(jl,1) = zdiff(jl,1) * zdiagt / dzsi(jl,1) /  &
              (dzsi(jl,1)+dzsi(jl,2)) * 2._wp
        ELSE
          zda(jl,1) = 0._wp
        END IF
      ELSE
        zda(jl, 1) = 0._wp
      END IF
    END DO

  ! *** Temporary storage for log output
    IF (ilog == 1 .OR. ilog == 3) wslog(1:kdeep) = wsi(jllog, 1:kdeep)

  ! *** Routine TRIDIAG from Numerical Recipes, p. 43
  ! *** -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WSI(T)=RI, WSI(T+1)=UI
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        zdc(jl,1) = zda(jl,1)+zdb(jl,1)+1._wp
        wsi(jl,1) = wsi(jl,1) / zdc(jl,1)
      END IF
    END DO
  ! *** Decomposition and forward substitution
    DO i=2, kdeep
      DO jl=kidia, kfdia
        IF (dzsi(jl,1) > 0) THEN

          IF (dzsi(jl,i) > 0) THEN
            ztri(jl,i) = -zda(jl,i-1) / zdc(jl,i-1)
            zdc(jl,i) = zda(jl,i)+zdb(jl,i)+1._wp + zdb(jl,i)*ztri(jl,i)
            wsi(jl,i) = (wsi(jl,i)/dzsi(jl,i) + zdb(jl,i)  &
                * wsi(jl,i-1)/dzsi(jl,i-1) ) / zdc(jl,i) * dzsi(jl,i)
          END IF
        END IF
      END DO
    END DO

    DO i=kdeep-1,1,-1       ! Backsubstitution
      DO jl=kidia, kfdia
        IF (dzsi(jl,1) > 0) THEN

          IF (dzsi(jl,i+1) > 0) THEN
            wsi(jl,i) = wsi(jl,i) - ztri(jl,i+1)*wsi(jl,i+1)*  &
                (dzsi(jl,i)/dzsi(jl,i+1))
          END IF
        END IF
      END DO
    END DO

  ! *** Log output of fluxes ?
    IF (ilog == 1 .OR. ilog == 3) THEN
      DO i=1, kdeep
        zdiflog(i) = wslog(i) - wsi(jllog, i)
      END DO
    END IF

  ! *** If layer overflow occurs after diffusion, it is added to drainage.
    DO i = 1, kdeep
      DO jl=kidia, kfdia
        IF (dzsi(jl,i) > 0) THEN
          IF (wsi(jl,i) > zwsfc(jl,i)) THEN
            IF (l_tpe_closed) THEN
              ! modification for TPE
              ! redirect overflow from drainage to w_soil_overflow pool
              w_soil_overflow(jl) = w_soil_overflow(jl) + wsi(jl,i) - zwsfc(jl,i)
            ELSE
              zdrain(jl) = zdrain(jl) + wsi(jl,i) - zwsfc(jl,i)
            END IF
            wsi(jl,i) = zwsfc(jl,i)
          END IF
        END IF
      END DO
    END DO

    IF (ilog == 2 .OR. ilog == 3) THEN
      WRITE(*,'(A10,1X, 11(E16.9E2, 1X) )')  &
          "SOILHYD1: ", fksat(jllog)*1000._wp, wsi(jllog,1:kdeep)*1000._wp,  &
          zwsfc(jllog,1:kdeep)*1000._wp
    END IF

  ! ******* NEW SOIL LAYER SOIL MOISTURE CONTENTS *********************

  ! *** Percolation per time step must be smaller than the usable soil
  ! *** layer contents.
    DO i = 1, kdeep
      DO jl=kidia, kfdia
        IF (dzsi(jl,i) > 0 .AND. zperc(jl,i) > 0) THEN
          zdum = zwdmin*zwsfc(jl,i)
          IF (wsi(jl,i)-zperc(jl,i) < zdum) THEN
            IF (wsi(jl,i) > zdum) THEN
              zperc(jl,i) = wsi(jl,i) - zdum
              wsi(jl,i) = zdum
            ELSE
              zperc(jl,i) = 0._wp
            END IF
          ELSE
            wsi(jl,i) = wsi(jl,i) - zperc(jl,i)
          END IF
        END IF
      END DO
    END DO

  ! *** Percolation added to downward layer or drainage
  ! *** if lowest layer or bedrock is reached.
    DO i = 1, kdeep-1
      DO jl=kidia, kfdia
        IF (dzsi(jl,1) > 0) THEN
          IF (dzsi(jl,i+1) > 0) THEN
            wsi(jl,i+1) = wsi(jl,i+1) + zperc(jl,i)
          ELSE
            IF (l_tpe_closed) THEN
              ! modification for TPE
              ! redirect percolation from drainage to overflow pool
              w_soil_overflow(jl) = w_soil_overflow(jl) + zperc(jl,i)
            ELSE
              zdrain(jl) = zdrain(jl) + zperc(jl,i)
            END IF
            zperc(jl,i) = 0._wp
          END IF
        END IF
      END DO
    END DO

    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        IF (l_tpe_closed) THEN
          ! modification for TPE
          ! redirect percolation from drainage to overflow pool
          w_soil_overflow(jl) = w_soil_overflow(jl) + zperc(jl,kdeep)
        ELSE
          zdrain(jl) = zdrain(jl) + zperc(jl,kdeep)
        END IF
      END IF
    END DO

  ! *** Percolation overflow is not allowed from layers below the surface layer.
  ! *** Instead it is assumed that water piles upwards from saturated layers.
    DO i = kdeep,2, -1
      DO jl=kidia, kfdia
        IF (dzsi(jl,i) > 0) THEN
          IF (wsi(jl,i) > zwsfc(jl,i)) THEN
            zperc(jl,i-1) = MAX( zperc(jl,i-1) - wsi(jl,i) + zwsfc(jl,i), 0._wp)
            wsi(jl,i-1) = wsi(jl,i-1) + wsi(jl,i) - zwsfc(jl,i)
            wsi(jl,i) = zwsfc(jl,i)
          END IF
        END IF
      END DO
    END DO

  ! *** Overflow from the surface layer enters drainage. this should not happen,
  ! *** so this is just for security.
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        IF (wsi(jl,1) > zwsfc(jl,1)) THEN
          IF (l_tpe_closed) THEN
            ! modification for TPE
            ! redirect w_soil_overflow from drainage to overflow pool
            w_soil_overflow(jl) = w_soil_overflow(jl) + wsi(jl,1) - zwsfc(jl,1)
          ELSE
            zdrain(jl) = zdrain(jl) + wsi(jl,1) - zwsfc(jl,1)
          END IF
          wsi(jl,1) = zwsfc(jl,1)
        END IF
      END IF
    END DO

  ! *** Log output of the fluxes ?
    IF (ilog == 1 .OR. ilog == 3) THEN
      DO i=1, kdeep
        zdiflog(kdeep+i) = zperc(jllog,i)
      END DO
      WRITE(*, '(A10,1X,11(E15.6E3,1X) )') &
         "SOILFLUX",zdiflog(1:10)*1000._wp, zdrain(jllog)*1000._wp
    END IF

  ! *** Log output of WSI!
    IF (ilog == 2 .OR. ilog == 3) THEN
      WRITE(*,'(A10,1X,10(E16.9E2,1X) )' ) &
         "SOILHYD2: ", wsi(jllog,1:kdeep)*1000._wp, dzsi(jllog,1:kdeep)*1000._wp
    END IF

  END SUBROUTINE soilhyd

  SUBROUTINE soildef (kidia, kfdia, klon, kdeep,  &
      dzr, dzs, vpor &
  !!$    , vfc, vpwp &
      , zdel, dzrsi, dzsi, zwsat &
  !!$   , zwsfc, zwpwp &
      )
  !     ***********************************************************************

  !     *** Routine defining various fields for the 5-layer soil scheme.

  !     ******** Version 1.0 - Juli 2006
  !              Programming and development: Stefan Hagemann
  !     ******** Version 1.1 - Juni 2007
  !              Optimization : Ralf Podzun
  !     ******** Version 1.2 - September 2007
  !              A loop optimized : Stefan Hagemann
  !     ******** Version 1.3 - Oktober 2007
  !              Minor technical correction without effect: Stefan Hagemann
  !     ******** Version 1.4 - January 2008
  !              Minor technical correction without effect: Stefan Hagemann
  !     ******** Version 1.5 - February 2008 - Stefan Hagemann
  !              Minor technical correction that may effect only 1-2 gridboxes
  !              where DZR=0. (near lakes)
  !     ******** Version 1.6 - March 2009 - Stefan Hagemann
  !              Volumetric soil field capacity vfc and wilting point vpwp are given
  !              via input list and do not need to be computed in soildef.
  !              new output: ZWPWP, WSMX eliminated.

  !     ***   ZDEL = Thickness of the 5 soil layers [m]
  !     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR) [m]
  !     ***   DZSI = Depth that can be saturated with water per soil layer
  !     ***          (until bedrock DZS) [m]
  !     ***    DZR = Rooting depth in [m]
  !     ***    DZS = Soil depth derived from textures in [m]
  !     ***   VPOR = Volumetric porosity [m water / m soil]
  !     ***    VFC = Volumetric field capacity [m water / m soil]
  !     ***   VPWP = Volumetric permanent wilting point [m water / m soil]
  !     ***  KDEEP = Number of soil layers
  !     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
  !     ***          Important for the uptake of infiltration in the soil layers.
  !     ***          Defined to the depth of the soil, that is, where the bedrock begins.
  !     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]
  !     ***  ZWPWP = Texture depending wilting point of soil layer I (analog to ZWSAT) [m]
  !
  !
  INTEGER,  INTENT(IN)                     :: kidia
  INTEGER,  INTENT(IN)                     :: kfdia
  INTEGER,  INTENT(IN)                     :: klon
  INTEGER,  INTENT(IN)                     :: kdeep

  REAL(wp), INTENT(IN)                     :: dzr(klon)
  REAL(wp), INTENT(IN)                     :: dzs(klon)
  REAL(wp), INTENT(IN)                     :: vpor(klon)
  !!$REAL(wp), INTENT(IN)                     :: vfc(klon)
  !!$REAL(wp), INTENT(IN)                     :: vpwp(klon)
  REAL(wp), INTENT(IN)                     :: zdel(kdeep)

  REAL(wp), INTENT(OUT)                    :: dzrsi(klon, kdeep)
  REAL(wp), INTENT(OUT)                    :: dzsi(klon, kdeep)
  REAL(wp), INTENT(OUT)                    :: zwsat(klon,kdeep)
  !!$REAL(wp), INTENT(OUT)                    :: zwsfc(klon,kdeep)
  !!$REAL(wp), INTENT(OUT)                    :: zwpwp(klon,kdeep)


  REAL(wp), PARAMETER :: zeps = 1.e-15_wp

  LOGICAL  :: loflag1(klon), loflag2(klon)

  INTEGER  :: i, jl
  REAL(wp) :: zdels(kdeep), zdelsm1(kdeep)

  ! ************ Field initialization -
    DO  jl=1,klon
      loflag1(jl) = .true.
      loflag2(jl) = .true.
    END DO

    dzsi (:,:) = 0._wp
    dzrsi(:,:) = 0._wp
  !  CALL setra(dzsi,klon*kdeep,0._wp)
  !  CALL setra(dzrsi,klon*kdeep,0._wp)

  !CDIR NOVECTOR
    DO  i=2,kdeep
      zdels(i)   = sum(zdel(1:i))
      zdelsm1(i) = sum(zdel(1:i-1))
    END DO

    zdels(1)   = zdel(1)
    zdelsm1(1) = 0._wp

    DO  i=1,kdeep
      DO  jl=kidia, kfdia
        IF (loflag1(jl)) THEN
          IF (dzr(jl) > zdels(i)) THEN
            dzrsi(jl,i) = zdel(i)
          ELSE
            dzrsi(jl,i) = dzr(jl) - zdelsm1(i)
            loflag1(jl) = .false.
          END IF
        END IF
      END DO
    END DO

  ! *** Saturatable depth
    DO  i=1,kdeep
      DO  jl=kidia, kfdia
        IF (loflag2(jl)) THEN
          IF (dzs(jl) > zdels(i)) THEN
            dzsi(jl,i) = zdel(i)
          ELSE
            dzsi(jl,i)  = dzs(jl) - zdelsm1(i)
            loflag2(jl) = .false.
          END IF
        END IF
      END DO
    END DO

    DO  i=1,kdeep
      DO  jl=1,klon
        zwsat(jl,i) = vpor(jl) * dzsi(jl,i)
  !!$      zwsfc(jl,i) = vfc(jl)  * dzsi(jl,i)
  !!$      zwpwp(jl,i) = vpwp(jl) * dzsi(jl,i)
      END DO
    END DO

  END SUBROUTINE soildef

  SUBROUTINE soilchange (kidia, kfdia, klon, kdeep,  isch ,  &
      dt, zdiags, &
  !         *** - Soil moisture layers,
      wsi   , wsim1m, dzr, dzrsi  , dzsi, zwsat, zwsfc, &
  !         ***   Vertical water fluxes  -- evap. fluxes over land
      finfil, drain, aetrans, aebsoil,                  &
  !         ***   Reduction of evaporation if necessary
      redevap, w_soil_overflow, ltpe_closed)

  !     Intermediate routine for changing 5 soil moisture layers due to
  !     evaporation fluxes, infiltration and drainage that are calculated
  !     with the old bucket scheme. The soil moisture changes are calculated
  !     for one time step, not two as previously done for ws using the
  !     3 time step scheme

  !     *** VS. 1.0 -- Juli 2006 -- Stefan Hagemann

  !     *** VS 1.1 August 2006 -- Stefan Hagemann
  !     *** Security check moved to ISCH =1 loop

  !     *** VS 1.2 -- January 2007 -- Stefan Hagemann
  !     *** IWORK renamed to ISCH .

  !     *** VS 1.3 -- September 2007 -- Stefan Hagemann
  !     *** Small error in transpiration loop eliminated that occured if transpiration
  !     ***   was taken from a layer where WSI < wilting point. The error lead to the
  !     ***   the result that the WSI of the corresponding layer was set to the wilting point.
  !     ***   The error only causes problems in dry region gridboxes with small WSMX that are
  !     ***   close to the bare soil evaporation water depth of 10 cm, and probably
  !     ***   only if small time steps are used, e.g. at 0.25 degree.
  !     ***   resolution.
  !     *** Introduction of zeps as small increment to be used instead of zero in some
  !     ***   comparison-if statements. small flux: 0.01 mm/day = 1.16e-10 m/s
  !     ***   ZEPS = 1e-15 << against a small flux.
  !     *** If WSI < WWILT for one layer, transpiration amount will be reduced
  !     ***   instead of re-distribution to other layers.
  !     ***   Thus, AETRANS is changed in SOILCHANGE, but also EVAP and EVAPL
  !     ***   must be changed.
  !     *** Bare soil evaporation will only be taken from the most upper layer.
  !     ***   overshooting amounts will be disregarded and AEBSOIL will be reduced.
  !     ***   thus, AEBSOIL is changed in SOILCHANGE, but also EVAP and EVAPL
  !     ***   must be changed.
  !     ***  --> New array submitted to routine: REDEVAP
  !     ***
  !     *** VS 1.4 -- March 2008 -- Stefan Hagemann
  !     *** Reduction of evaporation is taken out again as there was no feedback
  !     *** (technical) to the atmosphere. An implementation of this feedback would
  !     *** require more complex changes to VDIFF and comprehensive testing. In
  !     *** my opinion a reduction of transpiration in VDIFF would be desirable.
  !     *** Field REDEVAP is kept for diagnostic purposes.
  !     ***
  !     *** VS 1.5 - October 2009 Stefan Hagemann
  !     ***   Small error in soil moisture change for REDEVAP not zero was corrected.
  !     ***   This has led to errors in dry areas.
  !     ***   Infiltration fills the layers from above until field capacity. Initially
  !     ***   it was filled up until saturation but this caused (likely) unreal(wp)istic
  !     ***   high layer overflow values after diffusion(and thus drainage) in dry areas
  !     ***   as congo.
  !     ***
  !     *** VS 1.6 - February 2010 Stefan Hagemann
  !     *** Adding infiltration security check if soil for some reason
  !     *** cannot take it all

  !     ******** Input fields
  !     ***AEBSOIL = Bare soil evaporation in [m], accumulated over TWODT/2
  !     ***AETRANS = Transpiration in [m], accumulated over TWODT/2
  !     ***  DRAIN = Drainage in [m], accumulated over TWODT
  !     ***          --> Conversion to TWODT/2 needed
  !     *** FINFIL = Infiltration in [m/s] --> multiplication with DT required.

  !     ***    DZR = Rooting depth in [m]
  !     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR)
  !     ***   DZSI = Depth that can be saturated with water per soil layer
  !     ***          (until bedrock DZS)
  !     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
  !     ***          Important for the uptake of infiltration in the soil layers.
  !     ***          Defined to the depth of the soil, that is, where the bedrock begins.
  !     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]

  !     ***    DT = time step, usually [s]
  !     *** ISCH  = Switch which values will change the soil moisture
  !     ***         0 = No change: set WSI = WSIM1M
  !     ***         1 = change by evaporation fluxes
  !     ***         2 = change by infiltration
  !     ***         3 = change by infiltration and drainage
  !     ---------------------

  INTEGER,  INTENT(IN)                     :: kidia
  INTEGER,  INTENT(IN)                     :: kfdia
  INTEGER,  INTENT(IN)                     :: klon
  INTEGER,  INTENT(IN)                     :: kdeep
  INTEGER,  INTENT(IN)                     :: isch
  REAL(wp), INTENT(IN)                     :: dt
  REAL(wp), INTENT(IN)                     :: zdiags
  REAL(wp), INTENT(INOUT)                  :: wsi(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: wsim1m(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: dzr(klon)
  REAL(wp), INTENT(IN)                     :: dzrsi(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: dzsi(klon, kdeep)
  REAL(wp), INTENT(IN)                     :: zwsat(klon,kdeep)
  REAL(wp), INTENT(IN)                     :: zwsfc(klon,kdeep)
  REAL(wp), INTENT(IN)                     :: finfil(klon)
  REAL(wp), INTENT(INOUT)                  :: drain(klon)
  REAL(wp), INTENT(IN)                     :: aetrans(klon)
  REAL(wp), INTENT(IN)                     :: aebsoil(klon)
  REAL(wp), INTENT(INOUT)                  :: redevap(klon)
  REAL(wp), INTENT(INOUT), OPTIONAL        :: w_soil_overflow(klon)
  LOGICAL,  INTENT(IN),    OPTIONAL        :: ltpe_closed

  REAL(wp), PARAMETER :: zeps = 1.e-15_wp

  INTEGER :: jl, jk

  REAL(wp) :: zwdmin, zwilt
  REAL(wp) :: zrebsoil(klon), zretrans(klon)

  REAL(wp) :: zdum, zdum2, zrel
  REAL(wp) :: zover(klon)

  LOGICAL  :: l_tpe_closed

    IF (PRESENT(ltpe_closed)) THEN
      l_tpe_closed = ltpe_closed
    ELSE
      l_tpe_closed = .false.
    END IF

  ! *** Constants
    zwdmin = 0.05_wp
    zwilt  = 0.35_wp ! Read WAVA later on

  ! ******* No change = set WSI to WSIM1M
    IF (isch == 0) THEN
      wsi(kidia:kfdia,1:kdeep) = wsim1m(kidia:kfdia,1:kdeep)
      RETURN
    END IF

  ! ******* Change by evaporation fluxes
    IF (isch == 1) THEN
      zrebsoil(kidia:kfdia) = 0._wp
      zretrans(kidia:kfdia) = 0._wp
      zover(kidia:kfdia) = aebsoil(kidia:kfdia)   !vg bare soil evaporation [m]

  !   *** Bare soil evaporation loop
      DO  jl=kidia, kfdia
  !        *** TAU is packed in SURF in precipitation     !vg  dew added to rain in update_surf_up
  !        ***     --> Is contained in the infiltration!
          IF (zover(jl) < -zeps) THEN             !vg  considerable evaporation (why not < 0. ?)
            wsi(jl,1) = wsi(jl,1) + zover(jl)     !vg  take from soil water of uppermost layer 
            IF (wsi(jl,1) < 0._wp) THEN           !vg  if water of uppermost layer gets negative
              zover(jl) = wsi(jl,1)
              wsi(jl,1) = 0._wp
            ELSE
              zover(jl)=0._wp
            END IF
          END IF
      END DO

  !       *** Reduction of AEBSOIL necessary? --> diagnostic
      DO jl=kidia, kfdia
          IF (zover(jl) < -zeps) THEN             !vg  considerable amount of water could not be taken from first level
            zrebsoil(jl) = zover(jl)              !vg  reduce evaporation accordingly
  !!!            AEBSOIL(JL) = MIN(AEBSOIL(JL) - ZOVER(JL), 0.)
          END IF
      END DO

  !   *** transpiration loop
      zover(kidia:kfdia) = 0._wp
      DO  jk=1, kdeep
        DO  jl=kidia, kfdia
            IF (aetrans(jl) < -zeps) THEN         !vg  considerable transpiration (why not < 0. ?)

  !           *** Distribute transpiration over rooting depth
  !           *** for this, calculate rooted water depth
              IF (dzr(jl) <= zeps) THEN                                  !vg  rooting depth almost zero

  !             *** Special case: sealed grid cell
                IF (jk == 1) zretrans(jl) = zretrans(jl) + aetrans(jl)   !vg    no transpiration possible, add to reduction term

              ELSE IF (dzrsi(jl,jk) > 0) THEN                            !vg  roots reach into the soil layer
  !             *** Note that transpiration should only access water until the
  !             *** rooting depth --> pay regard where the relative term of the
  !             *** rooted depth within the layer has to be applied
                zrel = MIN(dzrsi(jl,jk)/dzsi(jl,jk), 1._wp)                     !vg  rooting depth relative to soil layer depth    
                                                                                !vg     (why should that become greater than 1 ?)
                zdum2 = zwsfc(jl,jk)*zwilt * zrel                               !vg  amount of water in root zone, not available for 
                                                                                !vg     transpiration
                IF (wsi(jl,jk)*zrel >= zdum2) THEN                              !vg  there is water in root zone available
                  zdum = wsi(jl,jk) * zrel + aetrans(jl) * dzrsi(jl,jk)/dzr(jl) !vg    amount of soil water after pot. transpiration
                                                                                !vg     (root depth in layer rel. to total root depth
                                                                                !vg      for portion of transpiration from the layer)
                  IF (zdum < zdum2) THEN                                        !vg    more water transpired than available
                    zover(jl) = zdum - zdum2                                    !vg      water deficit
                    wsi(jl,jk) = zdum2 + wsi(jl,jk) * MAX(1._wp-zrel, 0._wp)    !vg      soil water = amount in root zone not avail.
                                                                                !vg        for transpiration + water below root zone
                  ELSE                                                          !vg    enough water available
                    wsi(jl,jk) = zdum + wsi(jl,jk) * MAX(1._wp-zrel, 0._wp)     !vg      soil water = amount in root zone after tran-
                                                                                !vg        spiration + water below root zone
                    zover(jl)=0._wp                                             !vg      no deficit
                  END IF
                ELSE                                                            !vg  no water in root zone available for transpiration
                  zover(jl) = aetrans(jl) * dzrsi(jl,jk)/dzr(jl)                !vg    add pot. transpiration of layer to deficit 
                END IF
                zretrans(jl) = zretrans(jl) + zover(jl)                         !vg  add deficit to reduction term
              END IF
            END IF
        END DO
      END DO

  !   *** Adding ZRETRANS and ZREBSOIL to REDEVAP
      DO jl=kidia, kfdia
  !CC          IF (zretrans(JL).LT.-zeps) THEN
  !CC            aetrans(JL) = MIN(aetrans(JL) - redevap(JL), 0.)
  !CC          ENDIF
          redevap(jl) = redevap(jl) + zrebsoil(jl) + zretrans(jl)
      END DO

  !   *** If redevap not zero, soil moisture must be changed            !vg  why again???
  !   *** to close the land surface water balance.
  !   *** Reduction will be relatively distributed to the wetness
  !   *** of each layer.
      DO  jl=kidia, kfdia
          IF (redevap(jl) < -zeps) THEN
            zdum2 = SUM(wsi(jl,:))
            IF (zdum2 > 0._wp) THEN
              DO jk=1, kdeep
                zdum = wsi(jl,jk) + redevap(jl) * wsi(jl,jk)/zdum2
                wsi(jl,jk) = MAX(zdum, 0._wp)
              END DO
            END IF
          END IF
      END DO
    END IF

  ! ******* Change by infiltration
    IF (isch == 2 .OR. isch == 3) THEN
      zover(kidia:kfdia) = finfil(kidia:kfdia)*dt

      DO  jk = 1, kdeep
        DO  jl=kidia, kfdia

  !         *** infiltration-loop for 5 layer routine
            wsi(jl,jk) = wsi(jl,jk) + zover(jl)
            IF (l_tpe_closed) THEN
              IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
                w_soil_overflow(jl) = w_soil_overflow(jl) + wsi(jl,jk) - zwsfc(jl,jk)
                wsi(jl,jk) = zwsfc(jl,jk)
              END IF
              zover(jl) = 0._wp
            ELSE
              IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
                zover(jl)  = wsi(jl,jk) - zwsfc(jl,jk)
                wsi(jl,jk) = zwsfc(jl,jk)
              ELSE
                zover(jl) = 0._wp
              END IF
            END IF
        END DO
      END DO

  !   *** If infiltration cannot fully enter the soil
      DO jl=kidia, kfdia
          IF (zover(jl) > 0) drain(jl) = drain(jl)+zover(jl)
      END DO
    END IF

  ! ******* Change by drainage (unit m in SURF), used only for testing
    IF (isch == 3) THEN
      DO  jl=kidia, kfdia
        IF (drain(jl) > 0._wp) THEN
          zover(jl) = -drain(jl) * zdiags
          DO jk=kdeep,1,-1
  !         *** Start saturated(?) depth
            zdum = zwsfc(jl,jk)*zwdmin
            IF (zwsat(jl,jk) > 0 .AND. wsi(jl,jk)-zdum > 0) THEN
              wsi(jl,jk) =  wsi(jl,jk) + zover(jl)
              IF (wsi(jl,jk)-zdum < 0) THEN
                zover(jl) = wsi(jl,jk)-zdum
                wsi(jl,jk) = zdum
              ELSE
                zover(jl) =0._wp
                EXIT      ! Exit soil layer loop
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF

  ! ****** Avoid negative values - shouldn't happen
    wsi(kidia:kfdia,1:kdeep) = MAX( wsi(kidia:kfdia,1:kdeep), 0._wp)

  END SUBROUTINE soilchange

#endif
END MODULE mo_soil_process
