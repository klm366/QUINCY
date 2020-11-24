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
  USE mo_exception, ONLY: message, message_text, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, Contains_lct
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start

  dsl4jsb_Use_processes SEB_, RAD_, TURB_, A2L_, SSE_, HYDRO_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)

  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(SSE_)
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
    USE mo_sse_process,            ONLY: Get_liquid_max
    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: tmelt, rhoh2o, tpfac2, tpfac3, cvdifts
    USE mo_sse_constants,          ONLY: snow_depth_min
    USE mo_jsb_time,               ONLY: is_newday, timesteps_per_day, get_date_components ! get_asselin_coef, timeStep_in_days
    USE mo_jsb_grid,               ONLY: Get_grid
    USE mo_time_control,           ONLY: current_date


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    LOGICAL                :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: day_temp_min
    dsl4jsb_Real2D_onChunk :: day_temp_max
    dsl4jsb_Real2D_onChunk :: day_temp_sum
    dsl4jsb_Real2D_onChunk :: previous_day_temp_mean
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
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_ice_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real3D_onChunk :: t_soil
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef
    dsl4jsb_Real2D_onChunk :: lat

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

    REAL(wp), DIMENSION(options%nc)  :: &
     & soil_stor,                       &
     & t_soil_top

    LOGICAL, DIMENSION(options%nc) :: &
     & is_glacier,                    &
     & has_snow_layer

    REAL(wp) ::    & ! Avoid compute in param
     & soil_depth_sl_vol_porosity_sl_tmp(size(liquid_max,dim=1))
    REAL(wp) :: w_soil_critical_tmp

    INTEGER  :: iblk, ics, ice, nc, i
    REAL(wp) :: steplen
    INTEGER  :: tmp_timesteps_per_day
    LOGICAL  :: l_is_newday
    LOGICAL  :: l_is_clas
    LOGICAL  :: l_is_skin
    INTEGER  :: ilct
    INTEGER  :: yr, mo, dy, hr, mn, se

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    lat   => grid%lat(ics:ice, iblk)

    l_is_newday = is_newday(options%current_datetime, options%dtime)
    tmp_timesteps_per_day = timesteps_per_day(options%dtime)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   press_srf)      ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(A2L_,   drag_srf)       ! IN

    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_mean) ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_sum)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_min)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   day_temp_max)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_min)  ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   previous_day_temp_max)  ! OUT

    dsl4jsb_Get_var2D_onChunk(SEB_,   N_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   F_pseudo_soil_temp)     ! IN
    dsl4jsb_Get_var2D_onChunk(SEB_,   pseudo_soil_temp)       ! INOUT

    dsl4jsb_Get_var2D_onChunk(SEB_,   t)              ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_old)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt_old)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_eff)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   qsat_star)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   s_star)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   heat_cap)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB_,   forc_hflx)      ! OUT
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_q_air)     ! in
    dsl4jsb_Get_var2D_onChunk(TURB_,  fact_qsat_srf)  ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   rad_srf_net)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens) ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, w_soil_sl)      ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, w_ice_sl)       ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl) ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)      ! in

    dsl4jsb_Get_var3D_onChunk(SSE_,   t_soil)      ! in

    dsl4jsb_Get_var2D_onChunk(SSE_,   grnd_hflx)      ! IN
    dsl4jsb_Get_var2D_onChunk(SSE_,   hcap_grnd)      ! IN

    !$claw acc data create(fract,qsat_srf_old,s_old,dQdT,t2s_conv,heat_tcoef, &
    !$claw acc t_star,liquid_max,t_air_in_Celcius,is_glacier,has_snow_layer,  &
    !$claw acc soil_depth_sl_vol_porosity_sl_tmp)


    l_is_clas = dsl4jsb_Config(SEB_)%l_scheme_canopy_clas
    l_is_skin = dsl4jsb_Config(SEB_)%l_scheme_canopy_skin


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

    IF (is_time_experiment_start(options%current_datetime)) THEN            ! Start of experiment
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_0_0 tile_name=tile%name
    soil_depth_sl_vol_porosity_sl_tmp(:) = soil_depth_sl(:,1) * vol_porosity_sl(:,1)

    ! In the moment, the calls of calc_previous_day_variables and calc_pseudo_soil_temp are only
    ! necessary for the phenology process
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_0 tile_name=tile%name
    t_air_in_Celcius(:) = t_air(:) - tmelt  ! convert Kelvin in Celcius

    ! Sum up of the previous day temperatures
    !$claw sca forward update savepoint=ser_jsb_seb_land_calc_previous_day_variables_1 tile_name=tile%name
    CALL calc_previous_day_variables(l_is_newday,               & ! Input
                                     lstart,                    & ! Input
                                     t_air_in_Celcius(:),       & ! Input
                                     !dtime,                     & ! Input
                                     tmp_timesteps_per_day,     & ! Input
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

    open(555, file='output.txt')
    CALL get_date_components (current_date, yr, mo, dy, hr, mn, se)
    999 FORMAT(I0.2, I0.2, I0.2, 1X, I0.2, I0.2, I0.2, 1X, A4, 1X, ES15.7E3)
    997 FORMAT(I0.2, I0.2, I0.2, 1X, I0.2, I0.2, I0.2, 1X, A4, 1X, L1)

    t_soil_top(:) = t_soil(:,1)

    ! Compute new surface temperature and saturation humidity
    IF (lstart) THEN
      !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_3 tile_name=tile%name
      s_star(:) = s_old(:)
    ELSE

      IF (l_is_clas) THEN
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

      IF (l_is_skin) THEN


        CALL surface_temp_implicit_skin( &

          ! Input
            nc               &  ! Number of land grid points
          , steplen          &  ! Length of timestep (2 * delta t)
          , t2s_conv(:)      &  ! Conversion factor for humidity from dry static energy
          , t_acoef(:)       &  ! Linear   Richtmyer Morton coefficient for dry static energy (As)
          , t_bcoef(:)       &  ! Additive Richtmyer Morton coefficient for dry static energy (Bs)
          , q_acoef(:)       &  ! Linear   Richtmyer Morton coefficient for humidity (Aq)
          , q_bcoef(:)       &  ! Additive Richtmyer Morton coefficient for humidity (Bq)
          , s_old(:)         &  ! Dry static energy of the surface (old timestep)
          , qsat_srf_old(:)  &  ! Saturated specific humidity of the surface (old timestep)
          , dQdT(:)          &  ! Sensitivity of surface saturation specific humidity to temperature
          , rad_srf_net(:)   &  ! Net radiation
          , heat_tcoef(:)    &  ! Drag coefficient for heat
          , fact_q_air(:)    &  ! See (*)
          , fact_qsat_srf(:) &  ! Similar to fact_q_air (formerly c_sat)
          , fract_snow(:)    &  ! Snow fraction
          , heat_cap(:)      &  ! Heat capacity of uppermost soil layer (only diagnostic for SkIn+)
          , press_srf(:)     &  ! Surface pressure
          , t_soil_top(:)    &  ! Soil temperature of top/uppermost soil layer


        ! Output
          , forc_hflx(:)     &  ! Ground heat flux
          , s_star(:)        )  ! Dry static energy of the surface (new timestep)

        ! (*) Dimensionless factor for lowest atm. level in formulation of surface fluxes of specific humidity and
        ! dry static energy to account for evapotranspiration including snow, wetskin and vegetation cover fractions,
        ! stomatal and atmospheric resistance, as well as the specific humidity of air and the saturated specific humidity
        ! of the surface (formerly c_air)

      END IF


    END IF

    soil_stor(:)  = heat_cap(:) / t2s_conv(:)*(s_star(:) - s_old(:)) / (steplen*cvdifts)

      IF (yr == 1999  .AND. mo == 10 .AND. (dy == 23 .OR. dy == 24 .OR. dy == 25 .OR. dy == 26)) THEN
        DO i=1,size(rad_srf_net, dim=1)
          IF (lat(i) > 37._wp) THEN
  			    write(555,999) yr,mo,dy, hr,mn,se, 'Rnet', rad_srf_net(i)
            write(555,999) yr,mo,dy, hr,mn,se, 'Gflx', grnd_hflx(i)
            write(555,999) yr,mo,dy, hr,mn,se, 'Ssoi', soil_stor(i)
            write(555,999) yr,mo,dy, hr,mn,se, 'Tsoi', t_soil_top(i)
            write(555,999) yr,mo,dy, hr,mn,se, 'lati', lat(i)
            write(555,997) yr,mo,dy, hr,mn,se, 'Clas', l_is_clas
            write(555,997) yr,mo,dy, hr,mn,se, 'Skin', l_is_skin

          END IF
        END DO
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

    w_soil_critical_tmp = dsl4jsb_Config(SSE_)%w_soil_critical
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
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

    IF (dsl4jsb_Config(SSE_)%l_freeze) THEN

      IF (dsl4jsb_Config(SSE_)%l_supercool) THEN
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
#ifndef _OPENACC
    IF (ANY(4._wp * t_star(:) - 3._wp * t_old(:) <= 0._wp)) THEN
       i = MINLOC(4._wp * t_star(:) - 3._wp * t_old(:), DIM=1)
       WRITE (message_text,*) 'Instability: Extreme temperature difference from one time step to the next at ',     &
            & '(', grid%lon(i,iblk), ';', grid%lat(i,iblk), '): t_star: ', t_star(i), 'K,  t_old: ' , t_old(i), 'K'
       CALL finish(TRIM(routine), message_text)
    END IF
#endif
    !$claw expand parallel update savepoint=ser_jsb_seb_land_update_surface_energy_land_expand_10 tile_name=tile%name
    t_eff(:) = (t_old(:)**3 * (4._wp * t_star(:) - 3._wp * t_old(:))) ** 0.25_wp

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

    !$claw acc end data

  END SUBROUTINE update_surface_energy_land
  !
  ! ================================================================================================================================
  !
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
  SUBROUTINE surface_temp_implicit_skin( &

    ! Input
      nc            &  ! Number of land grid points
    , steplen       &  ! Length of timestep (2 * delta t)
    , t2s_conv      &  ! Conversion factor for humidity from dry static energy
    , t_acoef       &  ! Linear   Richtmyer Morton coefficient for dry static energy (As)
    , t_bcoef       &  ! Additive Richtmyer Morton coefficient for dry static energy (Bs)
    , q_acoef       &  ! Linear   Richtmyer Morton coefficient for humidity (Aq)
    , q_bcoef       &  ! Additive Richtmyer Morton coefficient for humidity (Bq)
    , s_old         &  ! Dry static energy of the surface (old timestep)
    , qsat_srf_old  &  ! Saturated specific humidity of the surface (old timestep)
    , dQdT          &  ! Sensitivity of surface saturation specific humidity to temperature
    , rad_srf_net   &  ! Net radiation
    , heat_tcoef    &  ! Drag coefficient for heat
    , fact_q_air    &  ! See (*)
    , fact_qsat_srf &  ! Similar to fact_q_air (formerly c_sat)
    , fract_snow    &  ! Snow fraction
    , heat_cap      &  ! Heat capacity of uppermost soil layer (only diagnostic for SkIn+)
    , press_srf     &  ! Surface pressure
    , t_soil_top    &  ! Soil temperature of top/uppermost soil layer



  ! Output
    , forc_hflx     &  ! Ground heat flux
    , s_star        )  ! Dry static energy of the surface (new timestep)

  ! (*) Dimensionless factor for lowest atm. level in formulation of surface fluxes of specific humidity and
  ! dry static energy to account for evapotranspiration including snow, wetskin and vegetation cover fractions,
  ! stomatal and atmospheric resistance, as well as the specific humidity of air and the saturated specific humidity
  ! of the surface (formerly c_air)

    USE mo_phy_schemes,            ONLY: qsat
    USE mo_jsb_physical_constants, ONLY: &
      !& cpd, &
      !& cpvd1, &
      & zemiss_def, &
      & stbo, &
      & alv,  &
      & als,  &
      & cvdifts

  ! Input
    INTEGER,  INTENT(in)    :: nc
    REAL(wp), INTENT(in)    :: steplen
    REAL(wp), INTENT(in)    :: t2s_conv(nc)
    REAL(wp), INTENT(in)    :: t_bcoef(nc)
    REAL(wp), INTENT(in)    :: t_acoef(nc)
    REAL(wp), INTENT(in)    :: q_bcoef(nc)
    REAL(wp), INTENT(in)    :: q_acoef(nc)
    REAL(wp), INTENT(in)    :: s_old(nc)
    REAL(wp), INTENT(in)    :: qsat_srf_old(nc)
    REAL(wp), INTENT(in)    :: dQdT(nc)
    REAL(wp), INTENT(in)    :: rad_srf_net(nc)
    REAL(wp), INTENT(in)    :: heat_tcoef(nc)
    REAL(wp), INTENT(in)    :: fact_q_air(nc)
    REAL(wp), INTENT(in)    :: fact_qsat_srf(nc)
    REAL(wp), INTENT(in)    :: fract_snow(nc)
    REAL(wp), INTENT(in)    :: heat_cap(nc)
    REAL(wp), INTENT(in)    :: press_srf(nc)
    REAL(wp), INTENT(in)    :: t_soil_top(nc)
  ! Output
    REAL(wp), INTENT(out)   :: forc_hflx(nc)
    REAL(wp), INTENT(out)   :: s_star(nc)



    REAL(wp) :: zdt ! time step (later needed for CanoStor in SkIn+)
    INTEGER  :: i

    ! Old local variables needed for classic SEB
    REAL(wp), DIMENSION(nc) :: &
        zcolin			&	! factors of linearized terms in energy balance
      , zcohfl			&	! factors from sensible and ground heat flux in energy balance
      , zcoind			&	! factors of independent terms in energy balance
      , zicp			  &	! inverse specific heat capacity of humid air
      , zca				  &	! cair, snow fraction, Lv and Lsub
      , zcs					  ! csat, snow fraction, Lv and Lsub

    ! Iteration variables
    REAL(wp), DIMENSION(nc) :: &
        zsend			  &	! final dry static energy
    	, zTdiffend		&	! final temperature difference between two iteration steps (for potential diag.)
    	, zsdiffend		&	! final dry static energy difference between two iteration steps (for potential diag.)
      , zTdiff			&	! temperature difference between two iteration steps (for potential diag.)
      , zsdiff				! dry static energy difference between two iteration steps

    REAL(wp), DIMENSION(nc) :: &
    	  T_old			  &	! surface temperature (old time step)
      , zTold			  &	! surface temperature (old iteration step)
    	, zTnew			  &	! new surface temperature
      , zsold			  &	! surface dry static energy (old iteration step)
  	  , zqsold			&	! saturated surface specific humidity (old iteration step)
  	  , zdqsold			&	! derivative saturated surface specific humidity (old iteration step)
  	  , znetrad			&	! net radiation (iteration)
      , zradcon			&	! constant part of net radiation (SWnet + LWdown)
      , zhtc				  ! stable/unstable heat transfer coefficient (skin conductivity)


    ! Parameters needed for heat transfer coefficient
    REAL(wp) :: 		&
    	  zslope			&	! slope of hyberbolic tangent
    	, zaxSym				! axial symmetry of hyberbolic tangent

    !
    REAL(wp) :: 		&
        zTeps			  & ! Small temperature difference used for derivative of saturated specific humidity
      , zpressfrac    ! Fraction between surface and lowest atmospheric level pressure (should be read dynamically)

    ! Iteration parameters
    REAL(wp) :: 		&
        zlam			  &	! iteration weighting factor
      , zlammin			&	! minimal iteration weighting factor
      , zlamdec			&	! decrease of iteration weighting factor
      , zsdelta				! stop criterium

    INTEGER,  DIMENSION(nc) :: &
        zfinaloop	    ! maximal number of iterations
    LOGICAL,  DIMENSION(nc) :: &
        literate		  ! grid points which have to be iterated

    INTEGER ::			&
    	  zloop			  &	! iteration loop
    	, zendloop		&	! final iteration loop
    	, zlamfreq			! number of iteration loops per decrease of iteration weighting factor

    ! Input/Output Variables needed for SkIn are set to dummy values in a first approach
    REAL(wp), DIMENSION(nc) :: &
    	  pspress     &
      , pTsoil


    ! Input/Output Variables needed for SkIn are set to dummy values in a first approach
      pspress(:) = 101300._wp
      pTsoil(:)  = s_old(:) * zicp(:)



    ! Heat transfer coefficient (Skin conductivity from ECMWF Part IV
    ! Physical Processes, Chapter 8 Surface Parametrization)
      zhtc(:) = 35._wp 		! [W/m^2/K]
      zslope  = -2._wp	 	! slope of hyperbolic tangent transition (this value can be adjusted)
      zaxSym  =  0.5_wp  	! axial symmetry of hyperbolic tangent (this value has to be 0.5)

    ! Small temperature difference used for derivative of saturated specific humidity
      zTeps   = 0.001_wp

    ! Fraction between surface and lowest atmospheric level pressure (should be read dynamically)
      zpressfrac   = 0.99615_wp

    ! Iteration options
      zloop       = 0
      zendloop    = 50
      zsdelta  	  = 100._wp
      literate(:) = .TRUE.

    ! Weighting factor to ensure convergence of iteration
      zlam       = 2._wp/3._wp
      zlammin    = 0.001_wp
      zlamdec    = 0.5_wp
      zlamfreq   = 10

    ! Starting/Missing values of iteration
      zsend(:)     = 0._wp
      zTdiffend(:) = 9999._wp
      zfinaloop(:) = 9999




      !-----------------------------------------------------------------------------------------
      ! Update SurfaceTemp (no iteration)
      !-----------------------------------------------------------------------------------------

      !$claw acc data &
      !$claw acc create(zcolin,zcohfl,zcoind,zicp,zca,zcs)
      ! Overimplicit time stepping: cvdifts = 3/2, steplen = 2 * delta_time
        zdt     = cvdifts * steplen

      !
      !$claw expand parallel update savepoint=ser_jsb_seb_land_surface_temp_implicit_skin_1

      ! Inverse specific heat capacity 1/cp
        zicp(:)     = 1._wp / t2s_conv(:)

      ! Factors used for latent heat flux (including atmospheric/stomatal resistance in cair/csat
        zca(:)      = als*fract_snow(:) +  alv*( fact_q_air(:) - fract_snow(:) )
        zcs(:)      = als*fract_snow(:) +  alv*( fact_qsat_srf(:) - fract_snow(:) )


      ! Surface temperature (old time step)
        T_old(:)    = s_old(:) * zicp(:)


      ! Net radiation without upward longwave radiation (SWnet + LWdown) so that it is
      ! independent from surface temperature (dry static energy)
        zradcon(:)  = rad_srf_net(:) + zemiss_def * stbo * T_old(:)**4

      ! Old time step values as starting points for iteration
        zTold(:)    = T_old(:)
        zsold(:)    = s_old(:)
        zqsold(:)   = qsat_srf_old(:)
        zdqsold(:)  = dQdT(:)
        znetrad(:)  = rad_srf_net(:)


    !----------------------------------------------------------------------------------------
    ! While loop
    !-----------------------------------------------------------------------------------------

      DO WHILE ( zloop <= zendloop)

    	  zloop = zloop + 1

    	!------------------------------------------------------------------------------------
    	! Update temperature dependent variables
    	!------------------------------------------------------------------------------------

    	  IF (zloop > 1) THEN
    	  ! Update dry static energy, temperature, humidity and net radiation
    	  	zsold(:)   = s_star(:)
    	  	zTold(:)   = zsold(:) * zicp(:)
    	  	zqsold(:)  = qsat(zTold(:), pspress(:))
    	  	zdqsold(:) = ( qsat(zTold(:)+zTeps, pspress(:)) - zqsold(:) ) / zTeps
    	  	znetrad(:) = zradcon(:) - zemiss_def * stbo * zTold(:)**4
    	  END IF

    	! Heat Transfer Coefficient (Hyperbolic tangent ensures smooth transition from stable to unstable)
    	  !zhtc(:) = zaxSym*(phtcu(:) - phtcs(:)) * tanh(-zslope*zril(:)) + zaxSym*(phtcu(:) + phtcs(:))


    	!------------------------------------------------------------------------------------
    	! Update surface temperature (iteration)
    	!------------------------------------------------------------------------------------




    	  zcohfl(:)     = zhtc(:)*zicp(:) - heat_tcoef(:) * (t_acoef(:) - 1._wp)

    	  zcolin(:)     = 4._wp * zicp(:) * zemiss_def * stbo * ( (zTold(:))**3 )  &
    					        - heat_tcoef(:) * ( zca(:)*q_acoef(:) - zcs(:) )                  &
                      * zicp(:) * zdqsold(:)


    	  zcoind(:)     = znetrad(:) + zhtc(:) * pTsoil(:)    &
    	  			        + heat_tcoef(:)*( t_bcoef(:) + (zca(:) * q_acoef(:) - zcs(:)) * zqsold(:) + zca(:) * q_bcoef(:) )

    	  s_star(:)      = (zcolin(:) * zsold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))

    	! Iteration filter
    	  s_star(:)      = zsold(:) + zlam * (s_star(:) - zsold(:))


    	  zTnew(:)      = s_star(:) * zicp(:)

    	  zsdiff(:)     = ABS(s_star(:) - zsold(:))
    	  zTdiff(:)     = ABS(zTnew(:) - zTold(:))


    	!------------------------------------------------------------------------------------
    	! Check iteration criteria
    	!------------------------------------------------------------------------------------

    	! Check which grid points are already converged and save their final values
    	  WHERE (zsdiff(:) < zsdelta .AND. literate(:))
    		literate(:)  = .FALSE.
    		zsend(:)     = s_star(:)
    		zTdiffend(:) = zTdiff(:)
    		zsdiffend(:) = zsdiff(:)
    		zfinaloop(:) = zloop
    	  END WHERE

    	! Decrease iteration weighting factor at certain iteration steps
    	  IF(MOD(zloop,zlamfreq) == 0) THEN
    		zlam = MAX(zlamdec * zlam, zlammin)
    	  END IF

    	! If the maximal number of iteration loops is reached, save the current values
    	  IF(zloop == zendloop) THEN
    		WRITE(*,*) 'Iteration did not converge after ', zloop ,' loops'
    		WHERE( literate(:) )
    		  literate(:)  = .FALSE.
    		  zsend(:)     = s_star(:)
    		  zTdiffend(:) = zTdiff(:)
    		  zsdiffend(:) = zsdiff(:)
    		  zfinaloop(:) = zloop
    		END WHERE
    	  END IF

      	! End iteration prematurely if all points are converged
    	  IF ( .NOT. ANY(literate(:)) ) THEN
    		EXIT
    	  END IF

      END DO
    !----------------------------------------------------------------------------------------
    ! End of while loop
    !-----------------------------------------------------------------------------------------

    ! Rewrite output variable with saved final values
      s_star(:)   = zsend(:)
      zTnew(:)   = s_star(:) * zicp(:)


    ! Ground heat flux
      forc_hflx(:)  = zhtc(:)*(pTsoil(:) - zTnew(:))



    !$claw end expand

    !$claw acc end data

  END SUBROUTINE surface_temp_implicit_skin
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_asselin_land(tile, options)

    USE mo_jsb_time,               ONLY: get_asselin_coef, get_date_components
    USE mo_time_control,           ONLY: current_date
    USE mo_jsb_grid,               ONLY: Get_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    !dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & t,           &
      & t_filt,      &
      & t_old,       &
      & t_unfilt,    &
      & t_unfilt_old, &
      & lat


    ! Locally allocated variables
    !
    INTEGER :: iblk, ics, ice, i
    INTEGER :: yr, mo, dy, hr, mn, se
    REAL(wp) :: eps

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    !dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_filt)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_old)          ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)       ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt_old)   ! in

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    lat   => grid%lat(ics:ice, iblk)

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

    CALL get_date_components (current_date, yr, mo, dy, hr, mn, se)
    999 FORMAT(I0.2, I0.2, I0.2, 1X,  I0.2, I0.2, I0.2, 1X, A4, 1X, ES15.7E3)

      IF (yr == 1999  .AND. mo == 10 .AND. (dy == 23 .OR. dy == 24 .OR. dy == 25 .OR. dy == 26)) THEN
        DO i=1,size(t_filt, dim=1)
          IF (lat(i) > 37._wp) THEN
            write(555,999) yr,mo,dy, hr,mn,se, 'Tsfc', t_filt(i)
          END IF
        END DO
      END IF

    t(:) = t_filt(:)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_asselin_land

  SUBROUTINE update_surface_fluxes_land(tile, options)

    USE mo_phy_schemes,            ONLY: heat_transfer_coef
    USE mo_jsb_physical_constants, ONLY: alv, als
    USE mo_jsb_time,               ONLY: get_date_components
    USE mo_time_control,           ONLY: current_date
    USE mo_jsb_grid,               ONLY: Get_grid

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: latent_hflx
    dsl4jsb_Real2D_onChunk :: sensible_hflx
    dsl4jsb_Real2D_onChunk :: latent_hflx_lnd
    dsl4jsb_Real2D_onChunk :: sensible_hflx_lnd
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: lat

    ! Locally allocated vectors
    !
    ! REAL(wp) ::                      &
    !   & s_air          (options%nc), &  !< Dry static energy at lowest atmospheric level
    !   & heat_tcoef     (options%nc)     !< Heat transfer coefficient (rho*C_h*|v|)
    REAL(wp), ALLOCATABLE ::                      &
      & s_air          (:), &  !< Dry static energy at lowest atmospheric level
      & heat_tcoef     (:)     !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc, i
    INTEGER  :: yr, mo, dy, hr, mn, se
    REAL(wp) :: steplen

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

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
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_acoef)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_bcoef)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   drag_srf)            ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,   s_star)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)          ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx)         ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   sensible_hflx_lnd)   ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   latent_hflx_lnd)     ! out

    ! Compute new dry static energy at lowest atmospheric level by back-substitution
    !
    !$cla$w acc data create(s_air, heat_tcoef)
    ALLOCATE(s_air(nc), heat_tcoef(nc))
    !$claw acc enter data create(s_air, heat_tcoef)

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    lat   => grid%lat(ics:ice, iblk)

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

    CALL get_date_components (current_date, yr, mo, dy, hr, mn, se)
    999 FORMAT(I0.2, I0.2, I0.2, 1X,  I0.2, I0.2, I0.2, 1X, A4, 1X, ES15.7E3)

      IF (yr == 1999  .AND. mo == 10 .AND. (dy == 23 .OR. dy == 24 .OR. dy == 25 .OR. dy == 26)) THEN
        DO i=1,size(latent_hflx, dim=1)
          IF (lat(i) > 37._wp) THEN
            write(555,999) yr,mo,dy, hr,mn,se, 'Lflx', latent_hflx(i)
            write(555,999) yr,mo,dy, hr,mn,se, 'Hflx', sensible_hflx(i)
          END IF
        END DO

      END IF

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
    & l_start,                                      & ! Input
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
    LOGICAL,   intent(in)    :: l_start

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

    IF (.NOT. l_start .AND. is_newday) THEN  ! day has changed --> recompute mean, min, max day temperature of previous day
                                             ! and reinitialize day_temp_sum(), day_temp_min(), day_temp_max()
       previous_day_temp_mean = day_temp_sum/time_steps_per_day
       day_temp_sum = t_air_in_Celcius

       previous_day_temp_min = day_temp_min
       day_temp_min = t_air_in_Celcius

       previous_day_temp_max = day_temp_max
       day_temp_max = t_air_in_Celcius

    ELSE  ! day has not changed or start of experiment (day_temp_sum is initialized with zero!)
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
