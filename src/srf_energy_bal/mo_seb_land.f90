!>
!! @brief Contains the routines for the surface energy balance on LAND lct_type.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version in mo_srf_interface.f90 by Reiner Schnur (2013-04-04)
!! Separated out the surface energy balance for LAND part by Reiner Schnur (2016-03-22)
!!
MODULE mo_seb_land
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, Contains_lct
  USE mo_jsb_control,        ONLY: debug_on

  dsl4jsb_Use_processes SEB_, RAD_, TURB_, ATM2LAND_, LAND2ATM_, SOIL_, HYDRO_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(HYDRO_)

  dsl4jsb_Use_memory(ATM2LAND_)
  ! dsl4jsb_Use_memory(LAND2ATM_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(SOIL_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_land, update_asselin_land, update_surface_fluxes_land

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_land'

CONTAINS
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_surface_energy_land(tile, options)

    !USE mo_seb_process,            ONLY: calc_previous_day_variables, calc_pseudo_soil_temp
    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: tmelt, tpfac2, tpfac3
    USE mo_jsb_time,               ONLY: get_asselin_coef, is_newday, timeStep_in_days, time_steps_per_day

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(ATM2LAND_)
    ! dsl4jsb_Def_memory(LAND2ATM_)

    LOGICAL                :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: day_temp_min
    dsl4jsb_Real2D_onChunk :: day_temp_max
    dsl4jsb_Real2D_onChunk :: day_temp_sum
    dsl4jsb_Real2D_onChunk :: previous_day_temp_mean    
    dsl4jsb_Real2D_onChunk :: previous_day_NPP_rate
    dsl4jsb_Real2D_onChunk :: previous_day_temp_min
    dsl4jsb_Real2D_onChunk :: previous_day_temp_max
    dsl4jsb_Real2D_onChunk :: F_pseudo_soil_temp
    dsl4jsb_Real2D_onChunk :: N_pseudo_soil_temp
    dsl4jsb_Real2D_onChunk :: pseudo_soil_temp
    
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: t_old
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: t_unfilt_old
    dsl4jsb_Real2D_onChunk :: t_eff
    dsl4jsb_Real2D_onChunk :: qsat_star
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: heat_cap
    dsl4jsb_Real2D_onChunk :: forc_hflx
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: w_snow
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc)  :: &
      fract,            &
      qsat_srf_old,     &
      s_old,            &
      dQdT,             & !< Sensitivity of saturated surface specific humidity to temperature
      t2s_conv,         & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
      heat_tcoef,       & !< Heat transfer coefficient (rho*C_h*|v|)
      t_air_in_Celcius  
      
    LOGICAL :: &
      is_glacier(options%nc)

    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: steplen !, dtime
    REAL(wp) :: eps
    INTEGER  :: ilct

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen
    !dtime   = options%dtime
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(ATM2LAND_)
    ! dsl4jsb_Get_memory(LAND2ATM_)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, press_srf)      ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_srf)       ! IN
    
    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_mean) ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,      day_temp_sum)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      day_temp_min)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      day_temp_max)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_min)  ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      previous_day_temp_max)  ! OUT
    
    dsl4jsb_Get_var2D_onChunk(SEB_,      N_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,      F_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,      pseudo_soil_temp)       ! INOUT

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt_old)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_eff)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      qsat_star)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      s_star)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      heat_cap)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,      forc_hflx)      ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,     fact_q_air)     ! in
    dsl4jsb_Get_var2D_onChunk(TURB_,     fact_qsat_srf)  ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      rad_srf_net)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)     ! in

    dsl4jsb_Get_var2D_onChunk(SOIL_,     grnd_hflx)      ! IN
    dsl4jsb_Get_var2D_onChunk(SOIL_,     hcap_grnd)      ! IN


    ! In the moment the calls of calc_previous_day_variables and calc_pseudo_soil_temp are only
    ! necessary for the phenology process
    t_air_in_Celcius = t_air - tmelt  ! convert Kelvin in Celcius
    ! Sum up of the previous day temperatures
    CALL calc_previous_day_variables(is_newday(),               & ! Input
                                     t_air_in_Celcius(:),       & ! Input
                                     !dtime,                     & ! Input
                                     time_steps_per_day(),      & ! Input
                                     previous_day_temp_mean(:), & ! Output (for summer- and evergreen)
                                     day_temp_sum(:),           & ! InOut
                                     previous_day_temp_min(:),  & ! Output (for crop)
                                     day_temp_min(:),           & ! InOut
                                     previous_day_temp_max(:),  & ! Output (for crop)
                                     day_temp_max(:) )            ! InOut

    ! Update of pseudo-soil temperature for each time step
    CALL calc_pseudo_soil_temp(t_air_in_Celcius(:)  , & ! Input
                               N_pseudo_soil_temp(:), & ! Input
                               F_pseudo_soil_temp(:), & ! Input
                               pseudo_soil_temp(:)  )   ! InOut (for summer-, evergreen and crop)



    
    ! Save old surface temperature and saturation spec. humidity
!!$    t_srf_unfilt_old(:) = t_srf_unfilt(:)
!!$    t_srf_old   (:) = t_srf(:)
    ! Temporary fix for temperatures over ocean getting too warm or too cold, leading to crash in convect_tables
    t_unfilt_old(:) = MAX(MIN(t_unfilt(:), 380._wp), 100._wp)
    t_old(:)        = MAX(MIN(t(:),        380._wp), 100._wp)
    qsat_srf_old(:) = qsat(t_old(:), press_srf(:))

    ! The heat capacity of the surface layer is currently taken as the heat capacity of the upper soil layer
    heat_cap(:) = hcap_grnd(:)

    ! The last term of the rhs of the energy balance equation is the conductive heat flux from the ground below
    forc_hflx(:) = grnd_hflx(:)

    IF (ALL(heat_cap(:) < EPSILON(1._wp))) THEN            ! Start of experiment
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    ! Update surface temperature
    !
    ! 1. Compute sensitivity of surface saturation specific humidity to temperature
    dQdT(:) = ( qsat(t(:) + 0.001_wp, press_srf(:)) - qsat_srf_old(:) ) * 1000._wp

    ! 2. Old dry static energy
    s_old(:) = surface_dry_static_energy( &
      & t_old(:), &
      & q_effective(qsat_srf_old(:), q_air(:), fact_qsat_srf(:), fact_q_air(:)))

    t2s_conv(:) = s_old(:) / t_old(:)  !

    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute new surface temperature and saturation humidity
    IF (lstart) THEN
      s_star(:) = s_old(:)
    ELSE
      CALL surface_temp_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), forc_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        fract_snow(:), heat_cap(:), &
        & s_star(:)) ! < the only output
    END IF

    ! TBD - hack for oceans boxes
    CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))
    WHERE (fract(:) < EPSILON(1._wp))
      s_star(:) = s_old(:)
    END WHERE

    ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    t_unfilt(:) = ( tpfac2 * s_star(:) + tpfac3 * s_old(:) ) / t2s_conv(:)

    ! TBD: currently, fraction of glacier has to be either zero or one!
    ! This is because seb does not run on glacier tile but on land tile ...
    IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN ! R: IF-Abfrage nicht unbedingt notwendig...
      DO ilct=1,SIZE(tile%lcts)
        IF (tile%lcts(ilct)%id == GLACIER_TYPE) THEN
          is_glacier(:) = tile%lcts(ilct)%fract(ics:ice,iblk) > 0._wp
          EXIT
        END IF
      END DO
    ELSE
      is_glacier(:) = .FALSE.
    END IF

    ! Correct new dry static energy for snowmelt
    WHERE (w_snow(:) > dsl4jsb_Config(HYDRO_)%wsn_critical .OR. is_glacier(:))
      s_star(:) = MIN(s_star(:) / t2s_conv(:), tmelt) * t2s_conv(:)
    END WHERE

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    IF (lstart) THEN
      qsat_star(:) = qsat_srf_old(:)
    ELSE
      qsat_star(:) = qsat_srf_old(:) + dQdT(:) * (s_star(:) - s_old(:)) / t2s_conv(:)
    END IF

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    t_eff(:) = (t_old(:)**3 * (4._wp * s_star(:)/t2s_conv(:) - 3._wp * t_old(:))) ** 0.25_wp

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_land

  SUBROUTINE surface_temp_implicit(klon,      & ! nc
                                   steplen,   & ! steplen
                                   pcp,       & ! t2s_conv
                                   pescoe,    & ! t_acoef
                                   pfscoe,    & ! t_bcoef
                                   peqcoe,    & ! q_acoef
                                   pfqcoe,    & ! q_bcoef
                                   psold,     & ! s_old
                                   pqsold,    & ! qsat_srf_old
                                   pdqsold,   & ! dQdT
                                   pnetrad,   & ! rad_srf_net
                                   pgrdfl,    & ! forc_hflx
                                   pcfh,      & ! heat_tcoef
                                   pcair,     & ! fact_q_air
                                   pcsat,     & ! fact_qsat_srf
                                   pfracsu,   & ! fract_snow
                                   pgrdcap,   & ! heat_cap
                                   psnew)       ! s_star. The only OUT.

    USE mo_jsb_physical_constants, ONLY: &
      & cpd, &
      & cpvd1, &
      & zemiss_def, &
      & stbo, &
      & alv,  &
      & als,  &
      & cvdifts


    INTEGER,  INTENT(in)    :: klon
    REAL(wp), INTENT(in)    :: steplen
    REAL(wp), INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL(wp), INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
    REAL(wp), INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
    REAL(wp), INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL(wp), INTENT(in)    :: pgrdcap(klon)
    REAL(wp), INTENT(out)   :: psnew(klon)
    !
    REAL(wp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
    REAL(wp) :: pdt
    REAL(wp) :: pc16

    pdt     = cvdifts * steplen
    pc16    = cpd * cpvd1

    !
    zicp(:) = 1._wp / pcp(:)
    !
    zca(:)    = als * pfracsu(:) +  alv * (pcair(:) - pfracsu(:))
    zcs(:)    = als * pfracsu(:) +  alv * (pcsat(:) - pfracsu(:))
    !
    zcolin(:) = pgrdcap(:)*zicp(:) +                                                              &
              &      pdt * (zicp(:) * 4._wp * zemiss_def * stbo * ((zicp(:) * psold(:))**3)  -    &
              &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:)) *                                    &
              &      zicp(:) * pdqsold(:))
    !
    zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
    !
    zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:) *                              &
              &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:)) + pgrdfl(:))
    !
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))

  END SUBROUTINE surface_temp_implicit
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_asselin_land(tile, options)

    USE mo_jsb_time,               ONLY: get_asselin_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t,           &
      & t_filt,      &
      & t_old,       &
      & t_unfilt,    &
      & t_unfilt_old

    ! Locally allocated variables
    !
    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: eps

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_asselin_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_filt)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)       ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt_old)   ! in

    ! Asselin time filter, if applicable
    eps = get_asselin_coef()
    IF (eps > 0._wp) THEN
      t_filt(:) = t_unfilt_old(:) + eps * (t_old(:) - 2._wp * t_unfilt_old(:) + t_unfilt(:))
    ELSE
      t_filt(:) = t_unfilt(:)
    END IF
    t(:) = t_filt(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_asselin_land

  SUBROUTINE update_surface_fluxes_land(tile, options)

    USE mo_phy_schemes,            ONLY: q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: alv, als

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(ATM2LAND_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: latent_hflx
    dsl4jsb_Real2D_onChunk :: sensible_hflx
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapotranspiration
    dsl4jsb_Real2D_onChunk :: fract_snow

    ! Locally allocated vectors
    !
    REAL(wp) ::                &
      & s_air          (options%nc),   &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef     (options%nc)       !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_fluxes_land'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(ATM2LAND_)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_acoef)             ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_bcoef)             ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_srf)            ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      s_star)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotranspiration)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)          ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx)         ! out

    ! Compute new dry static energy at lowest atmospheric level by back-substitution
    !
    s_air(:) = t_acoef(:) * s_star   (:) + t_bcoef(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute sensible heat flux
    !
    sensible_hflx(:) = heat_tcoef(:) * (s_air(:) - s_star(:))

    ! Compute latent heat flux
    !
    latent_hflx(:) = alv * evapotranspiration(:) + (als - alv) * fract_snow(:) * evapopot(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_land


    ! Updates the mean day temperature ("previous_day_temp_mean(:)") from air temperature and the mean day NPP-rate
  ! (field previous_day_NPP(:,:)
  ELEMENTAL SUBROUTINE calc_previous_day_variables(is_newday,              & ! Input
                                                   t_air_in_Celcius,       & ! Input
                                                   !dtime,                  & ! Input
                                                   time_steps_per_day,     & ! Input
                                                   previous_day_temp_mean, & ! Output
                                                   day_temp_sum,           & ! InOut
                                                   previous_day_temp_min,  & ! Output
                                                   day_temp_min,           & ! InOut
                                                   previous_day_temp_max,  & ! Output
                                                   day_temp_max )            ! InOut

    !-----------------------------------------------------------------------
    !  DECLARATIONS

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    LOGICAL,   intent(in)    :: is_newday

    REAL(wp),  intent(in)    :: t_air_in_Celcius!,  & ! air temperature at current time step in lowest atmospheric layer in Celcius
                                !dtime

    INTEGER,  intent(in)     :: time_steps_per_day

    REAL(wp),  intent(out)   :: previous_day_temp_mean, &
                                previous_day_temp_min,  &
                                previous_day_temp_max

    REAL(wp),  intent(inout) :: day_temp_sum,           &
                                day_temp_min,           &
                                day_temp_max


    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES


    !-----------------------------------------------------------------------
    ! CONTENT

    ! --- update mean day values

    ! Note that updating is done globally at the same time step, i.e. for different longitudes the updating happens at different
    ! local times.

    IF(is_newday) THEN  ! day has changed --> recompute mean, min, max day temperature of previous day
                        ! and reinitialize day_temp_sum(), day_temp_min(), day_temp_max()
       previous_day_temp_mean = day_temp_sum/time_steps_per_day
       day_temp_sum = t_air_in_Celcius

       previous_day_temp_min = day_temp_min
       day_temp_min = t_air_in_Celcius

       previous_day_temp_max = day_temp_max
       day_temp_max = t_air_in_Celcius

    ELSE  ! day has not changed
       day_temp_sum = day_temp_sum  + t_air_in_Celcius

       IF (day_temp_min > t_air_in_Celcius)   day_temp_min = t_air_in_Celcius
       IF (day_temp_max < t_air_in_Celcius)   day_temp_max = t_air_in_Celcius
    END IF

  END SUBROUTINE calc_previous_day_variables


  ! --- update_pseudo_soil_temp() --------------------------------------------------------------------------------------------------
  !
  ! This routine computes a weighted running mean of the air temperature, which is interpreted here as a pseudo soil temperature:
  ! (1)   T_ps(t) = N^(-1) * SUM(t'=-infty,t) T(t')*exp(-(t-t')*delta/tau_soil),
  ! where "T(t)" is the air temperature at time step with number "t", "delta" the length of the time step (in days) and
  ! "tau_soil" is the characteristic time for loosing the memory of temperature in the soil (also in days; this is a tuning
  ! parameter! The normalization "N" is
  ! (2)   N = SUM(t'=-infty,t) exp(-(t-t')*delta/tau_soil) = 1/(1 - exp(-delta/tau_soil)).
  ! This normalization constant (called "N_pseudo_soil_temp") is computed during initialization of this phenology module.
  ! Computation of T_ps(t) is performed iteratively: it follows from (1)
  !                    T(t+1)          delta
  ! (3)   T_ps(t+1) =  ------ + exp(- --------) * T_ps(t).
  !                      N            tau_soil
  ! The exponential factor (called F_pseudo_soil_temp) is computed during initialization of this phenology module.
  !
  ! Technically the only effect of this routine is an update of the field "pseudo_soil_temp(:)". The routine has to be called
  ! once every time step for every grid point.
  !
  ! Remark: Instead of air-temperature one could try to use the bottom temperature to compute the pseudo_soil_temp.
  !
  ELEMENTAL SUBROUTINE calc_pseudo_soil_temp(t_air_in_Celcius,   & ! Input
                                             N_pseudo_soil_temp, & ! Input
                                             F_pseudo_soil_temp, & ! Input
                                             pseudo_soil_temp    & ! InOut
                                             )

    REAL(wp),  intent(in)    ::                     &
                                t_air_in_Celcius,   &
                                N_pseudo_soil_temp, &
                                F_pseudo_soil_temp

    REAL(wp),  intent(inout) :: pseudo_soil_temp


    pseudo_soil_temp= t_air_in_Celcius / N_pseudo_soil_temp  +  F_pseudo_soil_temp * pseudo_soil_temp

  END SUBROUTINE calc_pseudo_soil_temp

#endif
END MODULE mo_seb_land
