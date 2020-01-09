!
! This file was automatically generated from JSBACH source file:
!    src/srf_energy_bal/mo_seb_lake.f90
! with md5 hash: 9f138d4dcf8dd4cb0b15b39fa5419855
!
!>
!! @brief Contains the routines for the surface energy balance on LAKE lct_type.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-22)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_seb_lake
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_process_class, ONLY: seb_, rad_, hydro_, a2l_
  USE mo_seb_config_class, ONLY: t_seb_config

  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_rad_memory_class, ONLY: t_rad_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_lake, update_surface_fluxes_lake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  SUBROUTINE update_surface_energy_lake(tile, options)

    USE mo_phy_schemes,            ONLY: qsat, surface_dry_static_energy

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & t,                    &
      & t_unfilt,             &
      & t_lwtr,               &
      & t_lice,               &
      & qsat_lwtr,            &
      & qsat_lice,            &
      & s_lwtr,               &
      & s_lice,               &
      & heat_cap,             &
      & press_srf,            &
      & snow,                 &
      & latent_hflx_wtr,      &
      & sensible_hflx_wtr,    &
      & latent_hflx_ice,      &
      & sensible_hflx_ice,    &
      & hflx_form_ice,        &
      & hflx_cond_ice,        &
      & hflx_melt_lice,       &
      & hflx_melt_snow,       &
      & fract_lice,           &
      & depth_lice,           &
      & fract_snow_lice,      &
      & w_snow_lice,          &
      & w_snow,               &
      & rad_net_lwtr,         &
      & rad_net_lice

    ! Local variables
    REAL(wp) ::                   &
      & hcap_wtr(options%nc),       &
      & zero(options%nc),           &
      & lhflx_tmp(options%nc),      &
      & hflx_tmp(options%nc)
    REAL(wp), ALLOCATABLE :: &
      hcap_ice(:)

    LOGICAL ::                       &
      & l_ice(options%nc),           & ! Indicator for ice cover
      & l_ice_formation(options%nc), &
      & l_ice_melting(options%nc)
    LOGICAL :: config_l_ice_on_lakes

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime
    REAL(wp) ::   &
      & ml_depth, &      !< Depth of lake mixed layer
      & min_ice_depth    !< Minimum ice thickness to start ice formation

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_lake'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(SEB_) .OR. .NOT. tile%contains_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT

    ml_depth      = seb__conf%lake_mixed_layer_depth
    min_ice_depth = seb__conf%lake_min_ice_depth
    config_l_ice_on_lakes = seb__conf%l_ice_on_lakes

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(rad_)%p); TYPE IS (t_rad_memory); rad__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT

    t => seb__mem%t%ptr(ics:ice,iblk)              ! out
    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)       ! out
    rad_net_lwtr => rad__mem%rad_net_lwtr%ptr(ics:ice,iblk)   ! in
    t_lwtr => seb__mem%t_lwtr%ptr(ics:ice,iblk)         ! inout
    qsat_lwtr => seb__mem%qsat_lwtr%ptr(ics:ice,iblk)      ! out
    s_lwtr => seb__mem%s_lwtr%ptr(ics:ice,iblk)         ! out
    heat_cap => seb__mem%heat_cap%ptr(ics:ice,iblk)       ! out

    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)         ! in
    w_snow => hydro__mem%w_snow%ptr(ics:ice,iblk)          ! out
    latent_hflx_wtr => seb__mem%latent_hflx_wtr%ptr(ics:ice,iblk)   ! in
    sensible_hflx_wtr => seb__mem%sensible_hflx_wtr%ptr(ics:ice,iblk) ! in
    fract_lice => seb__mem%fract_lice%ptr(ics:ice,iblk)        ! OUT

    IF (config_l_ice_on_lakes) THEN
      snow => a2l__mem%snow%ptr(ics:ice,iblk)              ! in
      t_lice => seb__mem%t_lice%ptr(ics:ice,iblk)            ! inout
      qsat_lice => seb__mem%qsat_lice%ptr(ics:ice,iblk)         ! out
      s_lice => seb__mem%s_lice%ptr(ics:ice,iblk)            ! out
      hflx_form_ice => seb__mem%hflx_form_ice%ptr(ics:ice,iblk)     ! inout
      hflx_cond_ice => seb__mem%hflx_cond_ice%ptr(ics:ice,iblk)
      hflx_melt_lice => seb__mem%hflx_melt_lice%ptr(ics:ice,iblk)
      hflx_melt_snow => seb__mem%hflx_melt_snow%ptr(ics:ice,iblk)

      depth_lice => seb__mem%depth_lice%ptr(ics:ice,iblk)        ! OUT
      fract_snow_lice => hydro__mem%fract_snow_lice%ptr(ics:ice,iblk) ! in
      w_snow_lice => hydro__mem%w_snow_lice%ptr(ics:ice,iblk)     ! inout

      latent_hflx_ice => seb__mem%latent_hflx_ice%ptr(ics:ice,iblk)   ! in
      sensible_hflx_ice => seb__mem%sensible_hflx_ice%ptr(ics:ice,iblk) ! in

      rad_net_lice => rad__mem%rad_net_lice%ptr(ics:ice,iblk)      ! in

      ALLOCATE(hcap_ice(nc))
      !$claw acc enter data create(hcap_ice)
    ELSE
      ALLOCATE(hflx_form_ice(nc))
      !$claw acc enter data create(hflx_form_ice)
    END IF

    !$claw acc data create(hcap_wtr,zero,l_ice, & 
    !$claw acc l_ice_formation,l_ice_melting,lhflx_tmp,hflx_tmp)

    !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_0
    zero            (:) = 0._wp
    hcap_wtr        (:) = 0._wp
    l_ice_formation (:) = .FALSE.
    l_ice_melting   (:) = .FALSE.
    hflx_tmp(:) = rad_net_lwtr(:) + latent_hflx_wtr(:) + sensible_hflx_wtr(:) 
    !$claw end expand
    IF (config_l_ice_on_lakes) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_1
      hcap_ice        (:) = 0._wp
      l_ice           (:) = fract_lice(:) > 0.5_wp
      !$claw end expand
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_2
      fract_lice      (:) = 0._wp
      l_ice           (:) = .FALSE.
      hflx_form_ice   (:) = 0._wp
      !$claw end expand
    END IF

    ! Calculate new water temperature and residual heat flux for ice formation over open water!
    !$claw sca forward update savepoint=ser_jsb_seb_lake_calc_lake_wtr_temp_1
    CALL calc_lake_wtr_temp( &
      & t_lwtr(:),                                      & !< inout; Surface temperature over open water
      & hflx_tmp(:),                                    & !< in;    Sum of all radiative and turbulent surface heat fluxes
      & hflx_form_ice(:),                               & !< inout; Residual heat flux for ice formation (realized or unrealized)
      & hcap_wtr(:),                                    & !< out;   Heat capacity of lake over open water
      & l_ice_formation(:),                             & !< out;   Indicator for ice formation
      & l_ice_melting(:),                               & !< inout; Indicator for ice melting = .FALSE.
      & l_ice(:),                                       & !< in;    Indicator for ice cover
      & dtime,                                          & !< in;    Time step
      & ml_depth,                                       & !< in;    Mixed-layer depth
      & min_ice_depth)                                    !< in;    Minimum ice thickness

    !$claw acc update host(fract_lice(:),hcap_wtr(:))
    IF (ANY(fract_lice(:) <= 0.5_wp .AND. hcap_wtr(:) < EPSILON(1._wp))) THEN
      CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    END IF

    IF (config_l_ice_on_lakes) THEN

      ! Calculate new ice thickness
      ! Tendency = freezing + conductive heat flux - ice melting + sublimation/deposition
      !$claw sca forward update savepoint=ser_jsb_seb_lake_calc_lice_thickness_1
      CALL calc_lice_thickness(                           &
        & depth_lice(:),                                  & !< inout; Ice thickness [m]
        & fract_lice(:),                                  & !< inout; Indicator for ice cover
        & hflx_form_ice(:),                               & !< inout; Heat flux for (complete) ice formation (l_ice_formation = T)
        & latent_hflx_ice(:),                             & !< in;    Latent heat flux for ice
        & hflx_cond_ice(:),                               & !< in;    Conductive heat flux for ice
        & hflx_melt_lice(:),                              & !< inout; Heat flux from unrealized or complete ice melting
        & l_ice_formation(:),                             & !< in;    Indicator for ice formation
        & l_ice_melting(:),                               & !< inout; Indicator for ice melting
        & fract_snow_lice(:),                             & !< in;    Fraction of snow on ice
        & min_ice_depth,                                  & !< in;    Mininum ice thickness
        & dtime                                           & !< in;    Time step
        & )
    
      ! Calculate new water temperature for the case of complete ice melting
      !$claw sca forward update savepoint=ser_jsb_seb_lake_calc_lake_wtr_temp_2
      CALL calc_lake_wtr_temp( &
        & t_lwtr(:),                                      & !< inout; Surface temperature over open water
        & zero(:),                                        & !< in;    Sum of radiative and turbulent surface heat fluxes (not used)
        & hflx_melt_lice(:),                              & !< inout; Residual heat flux from unrealized ice melting
        & hcap_wtr(:),                                    & !< inout; Heat capacity of lake over open water (not changed)
        & l_ice_formation(:),                             & !< inout; Indicator for ice formation (not changed)
        & l_ice_melting(:),                               & !< in;    Indicator for ice melting = .TRUE.
        & l_ice(:),                                       & !< in;    Indicator for ice cover
        & dtime,                                          & !< in;    Time step
        & ml_depth,                                       & !< in;    Mixed-layer depth
        & min_ice_depth)                                    !< in;    Minimum ice thickness
    
      ! Calculate new ice temperature and heat fluxes from melting of snow and ice

      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_3
      lhflx_tmp(:) = rad_net_lice(:) + latent_hflx_ice(:) + sensible_hflx_ice(:)

      !$claw sca forward update savepoint=ser_jsb_seb_lake_calc_lake_ice_temp_1
      CALL calc_lake_ice_temp(                            &
        & t_lice(:),                                      & ! inout
        & hcap_ice(:),                                    & ! out
        & hflx_cond_ice(:),                               & ! out
        & hflx_melt_snow(:),                              & ! out
        & hflx_melt_lice(:),                              & ! inout
        & lhflx_tmp(:),                                   & !< in;    Sum of all radiative and turbulent surface heat fluxes
        & latent_hflx_ice(:),                             & !< in;
        & snow(:),                                        & !< in; Snow fall
        & depth_lice(:),                                  &
        & w_snow_lice(:),                                 &
        & fract_snow_lice(:),                             &
        & min_ice_depth,                                  &
        & dtime)

    END IF

    !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_4
    qsat_lwtr(:)  = qsat(t_lwtr(:), press_srf(:))
    s_lwtr(:)     = surface_dry_static_energy(t_lwtr(:), qsat_lwtr(:))
    !$claw end expand

    IF (config_l_ice_on_lakes) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_5
      qsat_lice(:) = qsat(t_lice(:), press_srf(:))
      s_lice(:)    = surface_dry_static_energy(t_lice(:), qsat_lice(:))
      t(:) = (1._wp - fract_lice(:)) * t_lwtr   (:) + fract_lice(:) * t_lice   (:)
      heat_cap(:) = (1._wp - fract_lice(:)) * hcap_wtr(:) + fract_lice(:) * hcap_ice(:)
      w_snow(:) = fract_lice(:) * w_snow_lice(:)
      t_unfilt(:) = t(:)
      !$claw end expand

      !$claw acc exit data delete(hcap_ice)
      DEALLOCATE(hcap_ice)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_energy_lake_expand_6
      t(:) = t_lwtr(:)
      heat_cap(:) = hcap_wtr(:)
      w_snow(:) = 0._wp
      t_unfilt(:) = t(:)
      !$claw end expand

      !$claw acc exit data delete(hflx_form_ice)
      DEALLOCATE(hflx_form_ice)
    END IF

    !$claw acc end data

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

  ! Calculate new surface water temperature
  ELEMENTAL PURE SUBROUTINE calc_lake_wtr_temp(t_lwtr, hflx, hflx_resid, hcap, l_ice_formation, l_ice_melting, &
    &                                          l_ice, dtime, ml_depth, min_ice_depth)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, clw, rhoi, alf

    !$claw model-data
    REAL(wp), INTENT(inout) :: t_lwtr            !< Surface temperature over open water
    REAL(wp), INTENT(in)    :: hflx             !< Sum of all radiative and turbulent surface heat fluxes
    REAL(wp), INTENT(inout) :: hflx_resid       !< Heat flux for ice formation or melting
                                                !! l_ice_melting = .TRUE.  : heat flux from complete ice melting (input)
                                                !!                           zero on output
                                                !!               = .FALSE. : heat flux from unrealized ice formation from
                                                !!                           previous time step (input),
                                                !!                           heat flux from realized or unrealized ice formation
                                                !!                           in this time step (output)
    REAL(wp), INTENT(inout) :: hcap             !< Heat capacity of lake over open water
    LOGICAL,  INTENT(inout) :: l_ice_formation  !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting    !< Indicator for (complete) ice melting
    LOGICAL,  INTENT(in)    :: l_ice            !< Indicator for ice cover
    !$claw end model-data

    REAL(wp), INTENT(in)    :: dtime            !< Time step
    REAL(wp), INTENT(in)    :: ml_depth         !< Mixed-layer depth
    REAL(wp), INTENT(in)    :: min_ice_depth    !< Mininum ice thickness

    REAL(wp) :: &
      & t_star, &
      & t_freeze

    !$claw sca

    ! Lake tile is either not covered by ice or completely covered.
    IF (l_ice) RETURN

    IF (l_ice_melting) THEN
      t_lwtr = t_lwtr + dtime * hflx_resid / hcap
      hflx_resid = 0._wp
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
        hflx_resid = 0._wp
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

  ! Calculate new ice temperature, snow depth and heat fluxes
  ELEMENTAL PURE SUBROUTINE calc_lake_ice_temp(t_lice, hcap_ice, hflx_cond_ice, hflx_melt_snow, hflx_melt_lice, &
    &                                          hflx, lhflx, snowfall, hi, hs, fs, hi_min, dtime)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, rhoi, rhoi, rhos, als, alf, ki, ks, ci, cs !, clw

    !$claw model-data
    REAL(wp), INTENT(inout) :: t_lice
    REAL(wp), INTENT(out)   :: hcap_ice       !< Heat capacity of thin ice surface layer
    REAL(wp), INTENT(out)   :: hflx_cond_ice
    REAL(wp), INTENT(out)   :: hflx_melt_snow
    REAL(wp), INTENT(inout) :: hflx_melt_lice  !< in:  Residual heat flux from unrealized ice melting
                                              !! out: Heat flux from melting of ice
    REAL(wp), INTENT(in)    :: hflx           !< Sum of all surface heat fluxes (net radiation + sensible + latent)
    REAL(wp), INTENT(in)    :: lhflx          !< Latent surface heat flux for ice
    REAL(wp), INTENT(in)    :: snowfall       !< Snow fall                [kg m-2 s-1]
    REAL(wp), INTENT(in)    :: hi             !< Ice thickness [m]
    REAL(wp), INTENT(inout) :: hs             !< Snow depth [m water equivalent]
    REAL(wp), INTENT(in)    :: fs             !< Fraction of snow on ice
    !$claw end model-data
    REAL(wp), INTENT(in)    :: hi_min
    REAL(wp), INTENT(in)    :: dtime

    REAL(wp) :: hi_eff     !< Effective ice thickness including snow cover
    REAL(wp) :: hflx_melt  !< Excess heat available for melting snow and ice
    REAL(wp) :: delta_hs
    REAL(wp) :: delta_t_lice

    !$claw sca 

    hcap_ice = rhoi * ci * hi_min + rhoh2o * cs * hs

    hi_eff = hi + (ki / ks) * (rhoh2o / rhos) * hs

    IF (hi < hi_min) THEN
      t_lice = tmelt
      hflx_cond_ice = 0._wp
      hflx_melt_snow = 0._wp
      hflx_melt_lice = 0._wp
      hs = 0._wp
      RETURN
    END IF

    ! Preliminary snow depth
    hs = hs + (dtime / rhoh2o) * ( snowfall           & ! Snowfall
      &                          + fs * lhflx / als   & ! Sublimation (< 0) resp. deposition (> 0)
      &                          )
    hs = MAX(hs, 0._wp)

    ! Preliminary ice temperature
    t_lice = ( hflx + hflx_melt_lice + hcap_ice * t_lice / dtime + ki * tmelt / hi_eff ) / ( hcap_ice / dtime + ki / hi_eff)

    hflx_melt_snow = 0._wp
    hflx_melt_lice  = 0._wp

    IF (t_lice > tmelt .AND. hs > 0._wp) THEN        ! Melting for snow first
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

  ! Calculation of ice depth and fraction
  PURE ELEMENTAL SUBROUTINE calc_lice_thickness( &
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

    USE mo_jsb_physical_constants, ONLY: rhoi, alf, als

    !$claw model-data
    REAL(wp), INTENT(inout) :: hi              !< Ice thickness [m]
    REAL(wp), INTENT(inout) :: fi              !< Ice fraction  []
    REAL(wp), INTENT(inout) :: hflx_form_ice   !< Residual heat flux from ice formation (complete or unrealized) [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(in)    :: lhflx           !< Latent heat flux of ice  [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(in)    :: hflx_cond_ice   !< Conductive heat flux of ice [W m-2 = J m-2 s-1]
    REAL(wp), INTENT(inout) :: hflx_melt_lice   !< Heat flux from melting of ice [W m-2 = J m-2 s-1]
                                               !! in:  flux from unrealized ice melting from previous time step
                                               !! out: flux from unrealized or complete ice melting during this time step
    LOGICAL,  INTENT(in)    :: l_ice_formation !< Indicator for ice formation
    LOGICAL,  INTENT(inout) :: l_ice_melting   !< Indicator for complete ice melting
    REAL(wp), INTENT(in)    :: fs              !< Fraction of snow on ice
    !$claw end model-data

    REAL(wp), INTENT(in)    :: hi_min          !< Minimum ice thickness [m]
    REAL(wp), INTENT(in)    :: dtime

    !$claw sca

    IF (l_ice_formation) THEN                              ! Formation of a new slab of ice

      hi = - hflx_form_ice * dtime / (rhoi * alf)
      fi = 1._wp
      hflx_form_ice = 0._wp

    ELSE IF (hi >= hi_min) THEN                            ! Ice already present

      hi = hi + (dtime / rhoi) * (   (hflx_cond_ice - hflx_melt_lice) / alf    &
        &                          + (1._wp - fs) * lhflx / als               &
        &                        )

      IF (hi >= hi_min) THEN
        hflx_melt_lice = 0._wp
        fi = 1._wp
        l_ice_melting = .FALSE.
      ELSE IF (hi <= 0._wp) THEN                              ! Complete melting
        hflx_melt_lice = - rhoi * alf * hi / dtime             ! >= 0
        hi = 0._wp
        fi = 0._wp
        l_ice_melting = .TRUE.
      ELSE                                                    ! Unrealized melting
        hflx_melt_lice = rhoi * alf * (hi_min - hi) / dtime    ! > 0
        hi = hi_min
        fi = 1._wp
        l_ice_melting = .FALSE.
      END IF

    END IF

  END SUBROUTINE calc_lice_thickness

  SUBROUTINE update_surface_fluxes_lake(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef ! q_effective, surface_dry_static_energy
    USE mo_jsb_physical_constants, ONLY: alv, als

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: sensible_hflx
    REAL(wp), POINTER, DIMENSION(:) :: latent_hflx
    REAL(wp), POINTER, DIMENSION(:) :: s_lwtr
    REAL(wp), POINTER, DIMENSION(:) :: t_acoef_wtr
    REAL(wp), POINTER, DIMENSION(:) :: t_bcoef_wtr
    REAL(wp), POINTER, DIMENSION(:) :: drag_wtr
    REAL(wp), POINTER, DIMENSION(:) :: latent_hflx_wtr
    REAL(wp), POINTER, DIMENSION(:) :: sensible_hflx_wtr
    REAL(wp), POINTER, DIMENSION(:) :: evapo_wtr
    REAL(wp), POINTER, DIMENSION(:) :: s_lice
    REAL(wp), POINTER, DIMENSION(:) :: t_acoef_ice
    REAL(wp), POINTER, DIMENSION(:) :: t_bcoef_ice
    REAL(wp), POINTER, DIMENSION(:) :: drag_ice
    REAL(wp), POINTER, DIMENSION(:) :: latent_hflx_ice
    REAL(wp), POINTER, DIMENSION(:) :: sensible_hflx_ice
    REAL(wp), POINTER, DIMENSION(:) :: evapo_ice
    REAL(wp), POINTER, DIMENSION(:) :: fract_lice

    ! Locally allocated vectors
    !
    ! REAL(wp) ::                  &
    !   & s_air      (options%nc), &  !< Dry static energy at lowest atmospheric level
    !   & heat_tcoef (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)
    REAL(wp), ALLOCATABLE ::                  &
      & s_air      (:), &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef (:)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes_lake'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    t_acoef_wtr => a2l__mem%t_acoef_wtr%ptr(ics:ice,iblk)       ! IN
    t_bcoef_wtr => a2l__mem%t_bcoef_wtr%ptr(ics:ice,iblk)       ! IN
    drag_wtr => a2l__mem%drag_wtr%ptr(ics:ice,iblk)          ! IN
    drag_ice => a2l__mem%drag_ice%ptr(ics:ice,iblk)          ! IN

    s_lwtr => seb__mem%s_lwtr%ptr(ics:ice,iblk)            ! in
    sensible_hflx => seb__mem%sensible_hflx%ptr(ics:ice,iblk)     ! out
    sensible_hflx_wtr => seb__mem%sensible_hflx_wtr%ptr(ics:ice,iblk) ! out
    latent_hflx => seb__mem%latent_hflx%ptr(ics:ice,iblk)       ! out
    latent_hflx_wtr => seb__mem%latent_hflx_wtr%ptr(ics:ice,iblk)   ! out
    evapo_wtr => hydro__mem%evapo_wtr%ptr(ics:ice,iblk)         ! in

    IF (seb__conf%l_ice_on_lakes) THEN
      t_acoef_ice => a2l__mem%t_acoef_ice%ptr(ics:ice,iblk)       ! IN
      t_bcoef_ice => a2l__mem%t_bcoef_ice%ptr(ics:ice,iblk)       ! IN
      s_lice => seb__mem%s_lice%ptr(ics:ice,iblk)            ! in
      sensible_hflx_ice => seb__mem%sensible_hflx_ice%ptr(ics:ice,iblk) ! out
      latent_hflx_ice => seb__mem%latent_hflx_ice%ptr(ics:ice,iblk)   ! out
      evapo_ice => hydro__mem%evapo_ice%ptr(ics:ice,iblk)         ! in
      fract_lice => seb__mem%fract_lice%ptr(ics:ice,iblk)        ! in
    END IF
    ! ---------------------------
    ! Go

    ! TODO which one to use? The first one gave an error
    !$cla$w acc data create(s_air, heat_tcoef)
    ALLOCATE(s_air(nc), heat_tcoef(nc))
    !$claw acc enter data create(s_air,heat_tcoef)

    ! ================================================================================================================================
    ! Surface fluxes for lake water

    ! Compute new dry static energy at lowest atmospheric level by back-substitution
    !
    !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_fluxes_lake_expand_0
    s_air(:) = t_acoef_wtr(:) * s_lwtr  (:) + t_bcoef_wtr(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_wtr(:), steplen)

    ! Compute sensible heat flux
    !
    sensible_hflx_wtr(:) = heat_tcoef(:) * (s_air(:) - s_lwtr(:))

    ! Compute latent heat flux
    !
    latent_hflx_wtr(:) = alv * evapo_wtr(:)
    !$claw end expand 

    ! ================================================================================================================================
    ! Surface fluxes for lake ice

    IF (seb__conf%l_ice_on_lakes) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_fluxes_lake_expand_1
      !
      ! Compute new dry static energy at lowest atmospheric level by back-substitution
      !
      s_air(:) = t_acoef_ice(:) * s_lice  (:) + t_bcoef_ice(:)

      ! Transfer coefficient
      heat_tcoef(:) = heat_transfer_coef(drag_ice(:), steplen)

      ! Compute sensible heat flux
      !
      sensible_hflx_ice(:) = heat_tcoef(:) * (s_air(:) - s_lice(:))

      ! Compute latent heat flux
      !
      latent_hflx_ice(:) = als * evapo_ice(:)

      ! ================================================================================================================================
      ! Surface fluxes for lake (wtr and ice aggregated)

      sensible_hflx(:) = (1._wp - fract_lice(:)) * sensible_hflx_wtr(:) + fract_lice(:) * sensible_hflx_ice(:)
      latent_hflx  (:) = (1._wp - fract_lice(:)) * latent_hflx_wtr  (:) + fract_lice(:) * latent_hflx_ice  (:)
      !$claw end expand 
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_lake_update_surface_fluxes_lake_expand_2
      sensible_hflx(:) = sensible_hflx_wtr(:)
      latent_hflx  (:) = latent_hflx_wtr  (:)
      !$claw end expand 
    END IF

    !$cla$w acc end data
    !$claw acc exit data delete(s_air,heat_tcoef)
    DEALLOCATE(s_air, heat_tcoef)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_lake

#endif
END MODULE mo_seb_lake
