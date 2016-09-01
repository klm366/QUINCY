!>
!! @brief Contains the routines for the surface energy balance on LAKE lct_type.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-22)
!!
MODULE mo_seb_lake
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_control,        ONLY: debug_on

  dsl4jsb_Use_processes SEB_, RAD_, HYDRO_, ATM2LAND_
  dsl4jsb_Use_config(SEB_)

  dsl4jsb_Use_memory(ATM2LAND_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_lake, update_surface_fluxes_lake

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_lake'

CONTAINS

  SUBROUTINE update_surface_energy_lake(tile, options)

    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_lct_class,          ONLY: LAKE_TYPE, Contains_lct
    USE mo_jsb_physical_constants, ONLY: tmelt

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(ATM2LAND_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t, &
      & t_lwtr, &
      & t_lice, &
      & qsat_lwtr, &
      & qsat_lice, &
      & s_wtr, &
      & s_ice, &
      & heat_cap, &
      & press_srf, &
      & snow, &
      & latent_hflx_wtr, &
      & sensible_hflx_wtr, &
      & latent_hflx_ice, &
      & sensible_hflx_ice, &
      & hflx_form_ice, &
      & hflx_cond_ice, &
      & hflx_melt_lice, &
      & hflx_melt_snow, &
      & rad_srf_net, &
      & fract_lice, &
      & depth_lice, &
      & fract_snow_lice, &
      & w_snow_lice, &
      & w_snow, &
      & rad_net_lwtr, &
      & rad_net_lice

    ! Local variables
    REAL(wp) :: &
      hcap_wtr(options%nc), &
      hcap_ice(options%nc), &
      t_freeze(options%nc), &
      zero(options%nc)
    LOGICAL :: &
      & l_ice(options%nc),           & ! Indicator for ice cover
      & l_ice_formation(options%nc), &
      & l_ice_melting(options%nc)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime
    REAL(wp) :: &
      & ml_depth, &      !< Depth of lake mixed layer
      & min_ice_depth    !< Minimum ice thickness to start ice formation
    LOGICAL  :: lstart

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_lake'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_active(SEB_) .OR. .NOT. Contains_lct(tile%lcts, LAKE_TYPE)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)

    ml_depth      = dsl4jsb_Config(SEB_)%lake_mixed_layer_depth
    min_ice_depth = dsl4jsb_Config(SEB_)%lake_min_ice_depth

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(SEB_, t)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lwtr)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lice)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, s_wtr)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, s_ice)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, heat_cap)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, hflx_form_ice) ! inout
    dsl4jsb_Get_var2D_onChunk(SEB_, hflx_cond_ice)
    dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_lice)
    dsl4jsb_Get_var2D_onChunk(SEB_, hflx_melt_snow)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, press_srf)     ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, snow)          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_wtr)    ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_wtr)  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_ice)    ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_ice)  ! in

    dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice) ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_, depth_lice) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_lice)    ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow_lice) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow)      ! out

    dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr)         ! in
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice)         ! in

    zero            (:) = 0._wp
    hcap_wtr        (:) = 0._wp
    hcap_ice        (:) = 0._wp
    l_ice_formation (:) = .FALSE.
    l_ice_melting   (:) = .FALSE.
    l_ice           (:) = fract_lice(:) > 0.5_wp

    ! Calculate new water temperature and residual heat flux for ice formation over open water!
    CALL calc_lake_wtr_temp( &
      & t_lwtr(:),                                      & !< inout; Surface temperature over open water
      & rad_net_lwtr(:) +                               &
      &   latent_hflx_wtr(:)+sensible_hflx_wtr(:),      & !< in;    Sum of all radiative and turbulent surface heat fluxes
      & hflx_form_ice(:),                               & !< inout; Residual heat flux for ice formation (realized or unrealized)
      & hcap_wtr(:),                                    & !< inout; Heat capacity of lake over open water
      & l_ice_formation(:),                             & !< out;   Indicator for ice formation
      & l_ice_melting(:),                               & !< inout; Indicator for ice melting = .FALSE.
      & l_ice(:),                                       & !< in;    Indicator for ice cover
      & dtime,                                          & !< in;    Time step
      & ml_depth,                                       & !< in;    Mixed-layer depth
      & min_ice_depth)                                    !< in;    Minimum ice thickness

    IF (ANY(fract_lice(:) <= 0.5_wp .AND. hcap_wtr(:) < EPSILON(1._wp))) THEN
      CALL finish(TRIM(routine),'hcap_wtr can"t be zero')
    END IF

    ! Calculate new ice thickness
    ! Tendency = freezing + conductive heat flux - ice melting + sublimation/deposition
    CALL calc_lice_thickness(                           &
      & depth_lice(:),                                  & !< inout; Ice thickness [m]
      & fract_lice(:),                                  & !< inout; Indicator for ice cover
      & hflx_form_ice(:),                               & !< inout; Heat flux for (complete) ice formation (l_ice_formation = T)
      & latent_hflx_ice(:),                             & !< in;    Latent heat flux for ice
      & hflx_cond_ice(:),                               & !< in;    Conductive heat flux for ice
      & hflx_melt_lice(:),                               & !< inout; Heat flux from unrealized or complete ice melting
      & l_ice_formation(:),                             & !< in;    Indicator for ice formation
      & l_ice_melting(:),                               & !< inout; Indicator for ice melting
      & fract_snow_lice(:),                             & !< in;    Fraction of snow on ice
      & min_ice_depth,                                  & !< in;    Mininum ice thickness
      & dtime                                           & !< in;    Time step
      & )

    ! Calculate new water temperature for the case of complete ice melting
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
    CALL calc_lake_ice_temp( &
      & t_lice(:), &
      & hcap_ice, &
      & hflx_cond_ice, &
      & hflx_melt_snow, &
      & hflx_melt_lice, &
      & rad_net_lice(:) +                               &
      &   latent_hflx_ice(:)+sensible_hflx_ice(:),      & !< in;    Sum of all radiative and turbulent surface heat fluxes
      & latent_hflx_ice(:),                             & !< in;
      & snow(:),                                        & !< in; Snow fall
      & depth_lice, &
      & w_snow_lice, &
      & fract_snow_lice, &
      & min_ice_depth, &
      & dtime)

    qsat_lwtr(:)  = qsat(t_lwtr(:), press_srf(:))
    qsat_lice(:) = qsat(t_lice(:), press_srf(:))
    s_wtr(:)     = surface_dry_static_energy(t_lwtr(:), qsat_lwtr(:))
    s_ice(:)     = surface_dry_static_energy(t_lice(:), qsat_lice(:))

    t(:) = (1._wp - fract_lice(:)) * t_lwtr   (:) + fract_lice(:) * t_lice   (:)
    heat_cap(:) = (1._wp - fract_lice(:)) * hcap_wtr(:) + fract_lice(:) * hcap_ice(:)
    w_snow(:) = fract_lice(:) * w_snow_lice(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_lake

  ! Calculate new surface water temperature
  ELEMENTAL PURE SUBROUTINE calc_lake_wtr_temp(t_lwtr, hflx, hflx_resid, hcap, l_ice_formation, l_ice_melting, &
    &                                          l_ice, dtime, ml_depth, min_ice_depth)

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, clw, rhoi, alf

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
    REAL(wp), INTENT(in)    :: dtime            !< Time step
    REAL(wp), INTENT(in)    :: ml_depth         !< Mixed-layer depth
    REAL(wp), INTENT(in)    :: min_ice_depth    !< Mininum ice thickness

    REAL(wp) :: &
      & t_star, &
      & t_freeze

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

    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, rhoi, clw, rhoi, rhos, als, alf, ki, ks, ci, cs
    
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
    REAL(wp), INTENT(in)    :: hi_min
    REAL(wp), INTENT(in)    :: dtime

    REAL(wp) :: hi_eff     !< Effective ice thickness including snow cover
    REAL(wp) :: hflx_melt  !< Excess heat available for melting snow and ice
    REAL(wp) :: delta_hs
    REAL(wp) :: delta_t_lice

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
    REAL(wp), INTENT(in)    :: hi_min          !< Minimum ice thickness [m]
    REAL(wp), INTENT(in)    :: dtime

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

    USE mo_phy_schemes,            ONLY: q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: alv, als

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(ATM2LAND_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & evapopot,       &
      & evapotranspiration
    dsl4jsb_Real2D_onChunk :: qsat_lwtr
    dsl4jsb_Real2D_onChunk :: s_wtr
    dsl4jsb_Real2D_onChunk :: t_acoef_wtr
    dsl4jsb_Real2D_onChunk :: t_bcoef_wtr
    dsl4jsb_Real2D_onChunk :: q_acoef_wtr
    dsl4jsb_Real2D_onChunk :: q_bcoef_wtr
    dsl4jsb_Real2D_onChunk :: drag_wtr
    dsl4jsb_Real2D_onChunk :: latent_hflx_wtr
    dsl4jsb_Real2D_onChunk :: sensible_hflx_wtr
    dsl4jsb_Real2D_onChunk :: evapo_wtr
    dsl4jsb_Real2D_onChunk :: qsat_lice
    dsl4jsb_Real2D_onChunk :: s_ice
    dsl4jsb_Real2D_onChunk :: t_acoef_ice
    dsl4jsb_Real2D_onChunk :: t_bcoef_ice
    dsl4jsb_Real2D_onChunk :: q_acoef_ice
    dsl4jsb_Real2D_onChunk :: q_bcoef_ice
    dsl4jsb_Real2D_onChunk :: drag_ice
    dsl4jsb_Real2D_onChunk :: latent_hflx_ice
    dsl4jsb_Real2D_onChunk :: sensible_hflx_ice
    dsl4jsb_Real2D_onChunk :: evapo_ice
    dsl4jsb_Real2D_onChunk :: fract_lice

    ! Locally allocated vectors
    !
    REAL(wp) ::                &
      & s_air          (options%nc),   &  !< Dry static energy at lowest atmospheric level
      & q_air          (options%nc),   &  !< Humidity at lowest atmospheric level
      heat_tcoef       (options%nc)       !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_fluxes_lake'

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
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(ATM2LAND_)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_acoef_wtr)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_bcoef_wtr)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef_wtr)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef_wtr)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_acoef_ice)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_bcoef_ice)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef_ice)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef_ice)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_wtr)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_ice)          ! IN

    dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_lwtr)         ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      s_wtr)             ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_lice)         ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      s_ice)             ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      fract_lice)        ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx_wtr) ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx_wtr)   ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx_ice) ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx_ice)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_wtr)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_ice)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)          ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotranspiration) ! out
    
    ! ---------------------------
    ! Go
    
    ! ================================================================================================================================
    ! Surface fluxes for lake water

    ! Compute new dry static energy and moisture at lowest atmospheric level by back-substitution
    !
    s_air(:) = t_acoef_wtr(:) * s_wtr   (:) + t_bcoef_wtr(:)
    q_air(:) = q_acoef_wtr(:) * qsat_lwtr(:) + q_bcoef_wtr(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_wtr(:), steplen)

    ! Compute moisture fluxes
    !
    evapo_wtr (:) = heat_tcoef(:) * (q_air(:) - qsat_lwtr(:))       ! Potential evaporation

    ! Compute sensible heat flux
    !
    sensible_hflx_wtr(:) = heat_tcoef(:) * (s_air(:) - s_wtr(:))

    ! Compute latent heat flux
    !
    latent_hflx_wtr(:) = alv * evapo_wtr(:)

    ! ================================================================================================================================
    ! Surface fluxes for lake ice

    ! Compute new dry static energy and moisture at lowest atmospheric level by back-substitution
    !
    s_air(:) = t_acoef_ice(:) * s_ice   (:) + t_bcoef_ice(:)
    q_air(:) = q_acoef_ice(:) * qsat_lice(:) + q_bcoef_ice(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_ice(:), steplen)

    ! Compute moisture fluxes
    !
    evapo_ice (:) = heat_tcoef(:) * (q_air(:) - qsat_lice(:))       ! Potential evaporation

    ! Compute sensible heat flux
    !
    sensible_hflx_ice(:) = heat_tcoef(:) * (s_air(:) - s_ice(:))

    ! Compute latent heat flux
    !
    latent_hflx_ice(:) = als * evapo_ice(:)

    evapopot          (:) = (1._wp - fract_lice(:)) * evapo_wtr   (:) + fract_lice(:) * evapo_ice   (:)
    evapotranspiration(:) = evapopot(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_lake

#endif
END MODULE mo_seb_lake
