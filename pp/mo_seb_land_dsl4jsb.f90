!
! This file was automatically generated from:
!    /users/rschnur/icon-les-dev-staging.pgi_gpu/externals/jsbach/src/srf_energy_bal/mo_seb_land.f90
!    with md5 hash: 27dba4c2180bcae046e27dbe5002be8f
!
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
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_seb_land
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, Contains_lct
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start

  USE mo_jsb_process_class, ONLY: seb_, rad_, turb_, a2l_, sse_, hydro_
  !USE mo_seb_config_class, ONLY: t_seb_config
  USE mo_sse_config_class, ONLY: t_sse_config

  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_seb_memory_class, ONLY: t_seb_memory
  USE mo_rad_memory_class, ONLY: t_rad_memory
  USE mo_turb_memory_class, ONLY: t_turb_memory
  USE mo_sse_memory_class, ONLY: t_sse_memory
  USE mo_hydro_memory_class, ONLY: t_hydro_memory

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
    USE mo_sse_process,            ONLY: Get_liquid_max
    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, tpfac2, tpfac3
    USE mo_sse_constants,          ONLY: snow_depth_min
    USE mo_jsb_time,               ONLY: is_newday, time_steps_per_day ! get_asselin_coef, timeStep_in_days

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    !TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    TYPE(t_turb_memory), POINTER :: turb__mem
    TYPE(t_sse_memory), POINTER :: sse__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    LOGICAL                :: lstart

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: t_air
    REAL(wp), POINTER, DIMENSION(:) :: day_temp_min
    REAL(wp), POINTER, DIMENSION(:) :: day_temp_max
    REAL(wp), POINTER, DIMENSION(:) :: day_temp_sum
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_mean
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_min
    REAL(wp), POINTER, DIMENSION(:) :: previous_day_temp_max
    REAL(wp), POINTER, DIMENSION(:) :: f_pseudo_soil_temp
    REAL(wp), POINTER, DIMENSION(:) :: n_pseudo_soil_temp
    REAL(wp), POINTER, DIMENSION(:) :: pseudo_soil_temp

    REAL(wp), POINTER, DIMENSION(:) :: t
    REAL(wp), POINTER, DIMENSION(:) :: t_old
    REAL(wp), POINTER, DIMENSION(:) :: t_unfilt
    REAL(wp), POINTER, DIMENSION(:) :: t_unfilt_old
    REAL(wp), POINTER, DIMENSION(:) :: t_eff
    REAL(wp), POINTER, DIMENSION(:) :: qsat_star
    REAL(wp), POINTER, DIMENSION(:) :: s_star
    REAL(wp), POINTER, DIMENSION(:) :: heat_cap
    REAL(wp), POINTER, DIMENSION(:) :: forc_hflx
    REAL(wp), POINTER, DIMENSION(:) :: press_srf
    REAL(wp), POINTER, DIMENSION(:) :: q_air
    REAL(wp), POINTER, DIMENSION(:) :: fact_qsat_srf
    REAL(wp), POINTER, DIMENSION(:) :: fact_q_air
    REAL(wp), POINTER, DIMENSION(:) :: drag_srf
    REAL(wp), POINTER, DIMENSION(:) :: rad_srf_net
    REAL(wp), POINTER, DIMENSION(:) :: grnd_hflx
    REAL(wp), POINTER, DIMENSION(:) :: hcap_grnd
    REAL(wp), POINTER, DIMENSION(:) :: w_snow
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow
    REAL(wp), POINTER, DIMENSION(:) :: snow_soil_dens
    REAL(wp), POINTER, DIMENSION(:,:) :: soil_depth_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_soil_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: w_ice_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: vol_porosity_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: matrix_pot_sl
    REAL(wp), POINTER, DIMENSION(:,:) :: bclapp_sl
    REAL(wp), POINTER, DIMENSION(:) :: t_acoef
    REAL(wp), POINTER, DIMENSION(:) :: t_bcoef
    REAL(wp), POINTER, DIMENSION(:) :: q_acoef
    REAL(wp), POINTER, DIMENSION(:) :: q_bcoef

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc)  :: &
     & fract,                           &
     & qsat_srf_old,                    &
     & s_old,                           &
     & dQdT,                            & !< Sensitivity of saturated surface specific humidity to temperature
     & t2s_conv,                        & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
     & heat_tcoef,                      & !< Heat transfer coefficient (rho*C_h*|v|)
     & t_star,                          &
     & liquid_max,                      &
     & t_air_in_Celcius

    LOGICAL, DIMENSION(options%nc) :: &
     & is_glacier,                    &
     & has_snow_layer

    REAL(wp) ::    & ! Avoid compute in param
     & soil_depth_sl_vol_porosity_sl_tmp(size(liquid_max,dim=1))
    REAL(wp) :: w_soil_critical_tmp

    INTEGER :: iblk, ics, ice, nc, i
    REAL(wp) :: steplen
    INTEGER  :: ilct

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(rad_)%p); TYPE IS (t_rad_memory); rad__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(turb_)%p); TYPE IS (t_turb_memory); turb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(sse_)%p); TYPE IS (t_sse_memory); sse__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    t_air => a2l__mem%t_air%ptr(ics:ice,iblk)          ! IN
    t_acoef => a2l__mem%t_acoef%ptr(ics:ice,iblk)        ! IN
    t_bcoef => a2l__mem%t_bcoef%ptr(ics:ice,iblk)        ! IN
    q_acoef => a2l__mem%q_acoef%ptr(ics:ice,iblk)        ! IN
    q_bcoef => a2l__mem%q_bcoef%ptr(ics:ice,iblk)        ! IN
    press_srf => a2l__mem%press_srf%ptr(ics:ice,iblk)      ! IN
    q_air => a2l__mem%q_air%ptr(ics:ice,iblk)          ! IN
    drag_srf => a2l__mem%drag_srf%ptr(ics:ice,iblk)       ! IN

    previous_day_temp_mean => seb__mem%previous_day_temp_mean%ptr(ics:ice,iblk) ! IN
    day_temp_sum => seb__mem%day_temp_sum%ptr(ics:ice,iblk)           ! INOUT
    day_temp_min => seb__mem%day_temp_min%ptr(ics:ice,iblk)           ! INOUT
    day_temp_max => seb__mem%day_temp_max%ptr(ics:ice,iblk)           ! INOUT
    previous_day_temp_min => seb__mem%previous_day_temp_min%ptr(ics:ice,iblk)  ! OUT
    previous_day_temp_max => seb__mem%previous_day_temp_max%ptr(ics:ice,iblk)  ! OUT

    n_pseudo_soil_temp => seb__mem%n_pseudo_soil_temp%ptr(ics:ice,iblk)     ! IN
    f_pseudo_soil_temp => seb__mem%f_pseudo_soil_temp%ptr(ics:ice,iblk)     ! IN
    pseudo_soil_temp => seb__mem%pseudo_soil_temp%ptr(ics:ice,iblk)       ! INOUT

    t => seb__mem%t%ptr(ics:ice,iblk)              ! in
    t_old => seb__mem%t_old%ptr(ics:ice,iblk)          ! OUT
    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)       ! OUT
    t_unfilt_old => seb__mem%t_unfilt_old%ptr(ics:ice,iblk)   ! OUT
    t_eff => seb__mem%t_eff%ptr(ics:ice,iblk)          ! OUT
    qsat_star => seb__mem%qsat_star%ptr(ics:ice,iblk)      ! OUT
    s_star => seb__mem%s_star%ptr(ics:ice,iblk)         ! OUT
    heat_cap => seb__mem%heat_cap%ptr(ics:ice,iblk)       ! OUT
    forc_hflx => seb__mem%forc_hflx%ptr(ics:ice,iblk)      ! OUT
    fact_q_air => turb__mem%fact_q_air%ptr(ics:ice,iblk)     ! in
    fact_qsat_srf => turb__mem%fact_qsat_srf%ptr(ics:ice,iblk)  ! in
    rad_srf_net => rad__mem%rad_srf_net%ptr(ics:ice,iblk)    ! in
    w_snow => hydro__mem%w_snow%ptr(ics:ice,iblk)         ! in
    fract_snow => hydro__mem%fract_snow%ptr(ics:ice,iblk)     ! in
    snow_soil_dens => hydro__mem%snow_soil_dens%ptr(ics:ice,iblk) ! in
    soil_depth_sl => hydro__mem%soil_depth_sl%ptr(ics:ice,:,iblk)  ! in
    w_soil_sl => hydro__mem%w_soil_sl%ptr(ics:ice,:,iblk)      ! in
    w_ice_sl => hydro__mem%w_ice_sl%ptr(ics:ice,:,iblk)       ! in
    vol_porosity_sl => hydro__mem%vol_porosity_sl%ptr(ics:ice,:,iblk) ! in
    matrix_pot_sl => hydro__mem%matrix_pot_sl%ptr(ics:ice,:,iblk)  ! in
    bclapp_sl => hydro__mem%bclapp_sl%ptr(ics:ice,:,iblk)      ! in

    grnd_hflx => sse__mem%grnd_hflx%ptr(ics:ice,iblk)      ! IN
    hcap_grnd => sse__mem%hcap_grnd%ptr(ics:ice,iblk)      ! IN

    !$claw acc data create(fract,qsat_srf_old,s_old,dQdT,t2s_conv,heat_tcoef, &
    !$claw acc t_star,liquid_max,t_air_in_Celcius,is_glacier,has_snow_layer,  & 
    !$claw acc soil_depth_sl_vol_porosity_sl_tmp)

    ! @todo: currently, fraction of glacier has to be either zero or one!
    !$claw acc data present(is_glacier)
    IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN
      DO ilct=1,SIZE(tile%lcts)
        IF (tile%lcts(ilct)%id == GLACIER_TYPE) THEN
          is_glacier(:) = tile%lcts(ilct)%fract(ics:ice,iblk) > 0._wp
          EXIT
        END IF
      END DO
    ELSE
      is_glacier(:) = .FALSE.
    END IF
    !$claw acc update device(is_glacier) ! 
    !$claw acc end data

    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_0_0 tile_name=tile%name
    soil_depth_sl_vol_porosity_sl_tmp(:) = soil_depth_sl(:,1) * vol_porosity_sl(:,1)

    ! In the moment, the calls of calc_previous_day_variables and calc_pseudo_soil_temp are only
    ! necessary for the phenology process
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_0 tile_name=tile%name
    t_air_in_Celcius(:) = t_air(:) - tmelt  ! convert Kelvin in Celcius
    
    ! Sum up of the previous day temperatures
    !$claw sca forward update savepoint=ser_jsb_seb_land_calc_previous_day_variables_1 tile_name=tile%name
    CALL calc_previous_day_variables(is_newday(),               & ! Input
                                     t_air_in_Celcius(:),       & ! Input
                                     !dtime,                     & ! Input
                                     time_steps_per_day(),      & ! Input
                                     previous_day_temp_mean(:), & ! InOut (for summer- and evergreen)
                                     day_temp_sum(:),           & ! InOut
                                     previous_day_temp_min(:),  & ! InOut (for crop)
                                     day_temp_min(:),           & ! InOut
                                     previous_day_temp_max(:),  & ! InOut (for crop)
                                     day_temp_max(:) )            ! InOut

    ! Update of pseudo-soil temperature for each time step
    !$claw sca forward update savepoint=ser_jsb_seb_land_calc_pseudo_soil_temp_1 tile_name=tile%name
    CALL calc_pseudo_soil_temp(t_air_in_Celcius(:)  , & ! Input
                               N_pseudo_soil_temp(:), & ! Input
                               F_pseudo_soil_temp(:), & ! Input
                               pseudo_soil_temp(:)  )   ! InOut (for summer-, evergreen and crop)

    ! Save old surface temperature and saturation spec. humidity
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_1 tile_name=tile%name
    t_unfilt_old(:) = t_unfilt(:)
    t_old   (:) = t(:)
    !$claw end expand

    ! Temporary fix for temperatures over ocean getting too warm or too cold, leading to crash in convect_tables
    ! t_unfilt_old(:) = MAX(MIN(t_unfilt(:), 380._wp), 100._wp)
    ! t_old(:)        = MAX(MIN(t(:),        380._wp), 100._wp)
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_1_1 tile_name=tile%name
    qsat_srf_old(:) = qsat(t_old(:), press_srf(:))

    ! The heat capacity of the surface layer is currently taken as the heat capacity of the upper soil layer
    heat_cap(:) = hcap_grnd(:)

    ! The last term of the rhs of the energy balance equation is the conductive heat flux from the ground below
    forc_hflx(:) = grnd_hflx(:)
    !$claw end expand

    IF (is_time_experiment_start()) THEN            ! Start of experiment
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    ! Update surface temperature
    !
    ! 1. Compute sensitivity of surface saturation specific humidity to temperature
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_2 tile_name=tile%name
    dQdT(:) = ( qsat(t(:) + 0.001_wp, press_srf(:)) - qsat_srf_old(:) ) * 1000._wp

    ! 2. Old dry static energy
    s_old(:) = surface_dry_static_energy( &
      & t_old(:), &
      & q_effective(qsat_srf_old(:), q_air(:), fact_qsat_srf(:), fact_q_air(:)))

    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)
    !$claw end expand

    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_2_1 tile_name=tile%name
    t2s_conv(:) = s_old(:) / t_old(:)  !

    ! Compute new surface temperature and saturation humidity
    IF (lstart) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_3 tile_name=tile%name
      s_star(:) = s_old(:)
    ELSE
      !$claw sca forward update savepoint=ser_jsb_seb_land_surface_temp_implicit_1 tile_name=tile%name
      CALL surface_temp_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), forc_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        fract_snow(:), heat_cap(:), &
        & s_star(:))
    END IF

    ! @todo - hack for oceans boxes (land tile is direct child of root tile)
    !$claw acc data present(fract)
    CALL tile%parent_tile%Get_fraction(ics, ice, iblk, fract=fract(:))
    !$claw acc update device(fract)
    !$claw acc end data

    !$claw acc data present(fract,s_star,t2s_conv)
    !$claw acc parallel loop gang vector
    DO i=1, size(s_star,dim=1)
      IF (fract(i) < EPSILON(1._wp)) THEN
        ! s_star(:) = s_old(:)
        s_star(i) = 280._wp * t2s_conv(i)
      END IF
    END DO
    !$claw acc end data

    ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_4 tile_name=tile%name
    t_unfilt(:) = ( tpfac2 * s_star(:) + tpfac3 * s_old(:) ) / t2s_conv(:)

    ! Correct new dry static energy for melting of snow on the surface and (if l_freeze=.TRUE.) melting/freezing of
    ! ice/water in the top soil layer that will happen during this time step
    ! @todo: shouldn't this also consider wether there is actually enough energy to melt water/ice or freeze water?
    t_star(:) = s_star(:) / t2s_conv(:)
    !$claw end expand

    w_soil_critical_tmp = sse__conf%w_soil_critical
    IF (sse__conf%l_snow) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_5 tile_name=tile%name
      has_snow_layer(:) = w_snow(:) * rhoh2o / snow_soil_dens(:) > snow_depth_min
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_6 tile_name=tile%name
      has_snow_layer(:) = w_snow(:) > w_soil_critical_tmp
    END IF

    ! Actual snowmelt will be computed later in the hydrology process (update_surface_hydrology)
    !$claw acc data present(has_snow_layer,is_glacier,t_star)
    !$claw acc parallel loop gang vector
    DO i=1,size(t_star,dim=1)
      IF (has_snow_layer(i) .OR. is_glacier(i)) THEN
        ! Snow melt where glacier or snow layer
        t_star(i) = MIN(t_star(i), tmelt)
      END IF
    END DO
    !$claw acc end data

    IF (sse__conf%l_freeze) THEN

      IF (sse__conf%l_supercool) THEN
        !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_6_1 tile_name=tile%name
        liquid_max(:) = Get_liquid_max( &
          & t_star(:)                                   , &
          & soil_depth_sl_vol_porosity_sl_tmp(:)        , & ! Maximum water storage
          & matrix_pot_sl  (:,1)                        , &
          & bclapp_sl      (:,1)                          &
          & )
        !$claw end expand
      ELSE
        liquid_max(:) = 0._wp
      END IF

      !$claw acc data present(has_snow_layer,is_glacier,w_soil_sl,liquid_max,t_star,w_ice_sl)
      !$claw acc parallel loop gang vector
      DO i=1,size(t_star, dim=1)  
        IF (.NOT. has_snow_layer(i) .AND. .NOT. is_glacier(i)) THEN
          IF (w_soil_sl(i,1) > MAX(liquid_max(i), w_soil_critical_tmp)) THEN
            ! Freezing of water in uppermost soil layer where no snow layers above surface and considerable amount of moisture
            t_star(i) = MAX(t_star(i), tmelt)
          END IF
          IF (w_ice_sl(i,1) > w_soil_critical_tmp ) THEN
            ! Ice melt in uppermost soil layer where no snow layers above surface and considerable amount of ice
            t_star(i) = MIN(t_star(i), tmelt)
          END IF
        END IF
      END DO
      !$claw acc end data
    END IF

    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_7 tile_name=tile%name
    s_star(:) = t_star(:) * t2s_conv(:)

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    IF (lstart) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_8 tile_name=tile%name
      qsat_star(:) = qsat_srf_old(:)
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_9 tile_name=tile%name
      qsat_star(:) = qsat_srf_old(:) + dQdT(:) * (s_star(:) - s_old(:)) / t2s_conv(:)
    END IF

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_10 tile_name=tile%name
    t_eff(:) = (t_old(:)**3 * (4._wp * t_star(:) - 3._wp * t_old(:))) ** 0.25_wp

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

    !$claw acc end data

  END SUBROUTINE update_surface_energy_land 

  SUBROUTINE surface_temp_implicit(klon, steplen,             &
    &            pcp,                                         &
    &            pescoe, pfscoe, peqcoe, pfqcoe,              &
    &            psold, pqsold, pdqsold,                      &
    &            pnetrad, pgrdfl,                             &
    &            pcfh, pcair, pcsat, pfracsu, pgrdcap,        &
    &            psnew)

    USE mo_jsb_physical_constants, ONLY: &
      !& cpd, &
      !& cpvd1, &
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
    INTEGER :: i
    !REAL(wp) :: pc16

    !$claw acc data & 
    !$claw acc create(zcolin,zcohfl,zcoind,zicp,zca,zcs)
    pdt     = cvdifts * steplen
    !pc16    = cpd * cpvd1

    !
    !$claw expand parallel update savepoint=ser_jsb_seb_land_surface_temp_implicit_1
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
    !$claw end expand

    !$claw acc end data
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
    !TYPE(t_seb_config), POINTER :: seb__conf
    TYPE(t_seb_memory), POINTER :: seb__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: &
      & t,           &
      & t_filt,      &
      & t_old,       &
      & t_unfilt,    &
      & t_unfilt_old

    ! Locally allocated variables
    !
    INTEGER :: iblk, ics, ice
    REAL(wp) :: eps

    !TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    !SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT

    t => seb__mem%t%ptr(ics:ice,iblk)              ! out
    t_filt => seb__mem%t_filt%ptr(ics:ice,iblk)         ! out
    t_old => seb__mem%t_old%ptr(ics:ice,iblk)          ! in
    t_unfilt => seb__mem%t_unfilt%ptr(ics:ice,iblk)       ! in
    t_unfilt_old => seb__mem%t_unfilt_old%ptr(ics:ice,iblk)   ! in

    ! Asselin time filter, if applicable
    eps = get_asselin_coef()
    IF (eps > 0._wp) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_asselin_land_expand_0 tile_name=tile%name
      t_filt(:) = t_unfilt_old(:) + eps * (t_old(:) - 2._wp * t_unfilt_old(:) + t_unfilt(:))
    ELSE
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_asselin_land_expand_1 tile_name=tile%name
      t_filt(:) = t_unfilt(:)
    END IF
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_asselin_land_expand_2 tile_name=tile%name
    t(:) = t_filt(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_asselin_land

  SUBROUTINE update_surface_fluxes_land(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: alv, als

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_seb_memory), POINTER :: seb__mem
    TYPE(t_hydro_memory), POINTER :: hydro__mem
    TYPE(t_a2l_memory), POINTER :: a2l__mem

    ! Pointers to variables in memory
    REAL(wp), POINTER, DIMENSION(:) :: s_star
    REAL(wp), POINTER, DIMENSION(:) :: t_acoef
    REAL(wp), POINTER, DIMENSION(:) :: t_bcoef
    REAL(wp), POINTER, DIMENSION(:) :: drag_srf
    REAL(wp), POINTER, DIMENSION(:) :: latent_hflx
    REAL(wp), POINTER, DIMENSION(:) :: sensible_hflx
    REAL(wp), POINTER, DIMENSION(:) :: latent_hflx_lnd
    REAL(wp), POINTER, DIMENSION(:) :: sensible_hflx_lnd
    REAL(wp), POINTER, DIMENSION(:) :: evapopot
    REAL(wp), POINTER, DIMENSION(:) :: evapotrans
    REAL(wp), POINTER, DIMENSION(:) :: fract_snow

    ! Locally allocated vectors
    !
    ! REAL(wp) ::                      &
    !   & s_air          (options%nc), &  !< Dry static energy at lowest atmospheric level
    !   & heat_tcoef     (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)
    REAL(wp), ALLOCATABLE ::                      &
      & s_air          (:), &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef     (:)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    !TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes_land'

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SEB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(hydro_)%p); TYPE IS (t_hydro_memory); hydro__mem => mm ; END SELECT
    SELECT TYPE (mm => tile%mem(a2l_)%p); TYPE IS (t_a2l_memory); a2l__mem => mm ; END SELECT

    t_acoef => a2l__mem%t_acoef%ptr(ics:ice,iblk)             ! in
    t_bcoef => a2l__mem%t_bcoef%ptr(ics:ice,iblk)             ! in
    drag_srf => a2l__mem%drag_srf%ptr(ics:ice,iblk)            ! in

    s_star => seb__mem%s_star%ptr(ics:ice,iblk)              ! in
    evapotrans => hydro__mem%evapotrans%ptr(ics:ice,iblk)  ! in
    evapopot => hydro__mem%evapopot%ptr(ics:ice,iblk)            ! in
    fract_snow => hydro__mem%fract_snow%ptr(ics:ice,iblk)          ! in

    sensible_hflx => seb__mem%sensible_hflx%ptr(ics:ice,iblk)       ! out
    latent_hflx => seb__mem%latent_hflx%ptr(ics:ice,iblk)         ! out
    sensible_hflx_lnd => seb__mem%sensible_hflx_lnd%ptr(ics:ice,iblk)   ! out
    latent_hflx_lnd => seb__mem%latent_hflx_lnd%ptr(ics:ice,iblk)     ! out

    ! Compute new dry static energy at lowest atmospheric level by back-substitution
    !
    !$cla$w acc data create(s_air, heat_tcoef)
    ALLOCATE(s_air(nc), heat_tcoef(nc))
    !$claw acc enter data create(s_air, heat_tcoef)

    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_fluxes_land_expand_0 tile_name=tile%name
    s_air(:) = t_acoef(:) * s_star   (:) + t_bcoef(:)

    ! Transfer coefficient
    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute sensible heat flux
    !
    sensible_hflx(:) = heat_tcoef(:) * (s_air(:) - s_star(:))

    ! Compute latent heat flux
    !
    latent_hflx(:) = alv * evapotrans(:) + (als - alv) * fract_snow(:) * evapopot(:)

    ! These two variables are not aggregated together with lake and are the fluxes given back to the atmosphere
    sensible_hflx_lnd(:) = sensible_hflx(:)
    latent_hflx_lnd  (:) = latent_hflx  (:)
    !$claw end expand

    !$cla$w acc end data
    !$claw acc exit data delete(s_air, heat_tcoef)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_fluxes_land


    ! Updates the mean day temperature ("previous_day_temp_mean(:)") from air temperature and the mean day NPP-rate
  ! (field previous_day_NPP(:,:)
  ELEMENTAL SUBROUTINE calc_previous_day_variables( &
    & is_newday,                                    & ! Input
    & t_air_in_Celcius,                             & ! Input
    ! & dtime,                                        & ! Input
    & time_steps_per_day,                           & ! Input
    & previous_day_temp_mean,                       & ! InOut
    & day_temp_sum,                                 & ! InOut
    & previous_day_temp_min,                        & ! InOut
    & day_temp_min,                                 & ! InOut
    & previous_day_temp_max,                        & ! InOut
    &  day_temp_max )                                 ! InOut

    !-----------------------------------------------------------------------
    !  DECLARATIONS

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    LOGICAL,   intent(in)    :: is_newday

    !$claw model-data
    REAL(wp),  intent(in)    :: t_air_in_Celcius!,  & ! air temperature at current time step in lowest atmospheric layer in Celcius
                                !dtime
    !$claw end model-data

    INTEGER,  intent(in)     :: time_steps_per_day

    !$claw model-data
    REAL(wp),  intent(inout) :: previous_day_temp_mean, & ! "in" because if it is not calculated new, it should remain
                                previous_day_temp_min,  & ! as before. Without "in" it would be NaN in the output.
                                previous_day_temp_max

    REAL(wp),  intent(inout) :: day_temp_sum,           &
                                day_temp_min,           &
                                day_temp_max
    !$claw end model-data

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES


    !-----------------------------------------------------------------------
    ! CONTENT

    ! --- update mean day values

    ! Note that updating is done globally at the same time step, i.e. for different longitudes the updating happens at different
    ! local times.

    !$claw sca                            

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
  ELEMENTAL SUBROUTINE calc_pseudo_soil_temp( &
    & t_air_in_Celcius,                       & ! Input
    & N_pseudo_soil_temp,                     & ! Input
    & F_pseudo_soil_temp,                     & ! Input
    & pseudo_soil_temp                        & ! InOut
    & )

    !$claw model-data
    REAL(wp), intent(in) :: &
     & t_air_in_Celcius,    &
     & N_pseudo_soil_temp,  &
     & F_pseudo_soil_temp

    REAL(wp),  intent(inout) :: pseudo_soil_temp
    !$claw end model-data

    !$claw sca

    pseudo_soil_temp= t_air_in_Celcius / N_pseudo_soil_temp  +  F_pseudo_soil_temp * pseudo_soil_temp

  END SUBROUTINE calc_pseudo_soil_temp

#endif
END MODULE mo_seb_land
