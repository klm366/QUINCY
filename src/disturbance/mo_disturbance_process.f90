!>
!! @brief Contains the routines for the disturb processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! 2018-04-25, StW: Changed routine names for consistency with namelists
!!
MODULE mo_disturb_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_exception,        ONLY: message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: burned_fractt_jsb3, broken_woody_fract_jsb3, get_relative_humidity_air
  PUBLIC :: calc_fire_spitfire_fuel, calc_fire_spitfire_fuel_moist, calc_fire_spitfire_fire_danger_index
  PUBLIC :: calc_fire_spitfire_fire_geometry, calc_fire_spitfire_fire_density,calc_fire_spitfire_fuel_consumption
  PUBLIC :: calc_vegetation_height, calc_fire_spitfire_fire_intensity, calc_fire_spitfire_fire_mortality
  PUBLIC :: calc_fire_spitfire_pool_consumption, fire_spitfire_relocate_carbon

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  PURE ELEMENTAL SUBROUTINE burned_fractt_jsb3( &
    & fire_rel_hum_threshold,   & ! in
    & fire_litter_threshold,    & ! in
    & fire_minimum,             & ! in
    & fire_tau,                 & ! in
    & q_rel_air_climbuf,        & ! in
    & fuel,                     & ! in
    & burned_fract              & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in) ::     &
      & fire_rel_hum_threshold, &
      & fire_litter_threshold,  &
      & fire_minimum,           &
      & fire_tau,               &
      & q_rel_air_climbuf,      &
      & fuel

    ! Output Arguments
    REAL(wp), INTENT(out) :: burned_fract

    ! Locals
    REAL(wp) :: delta_time_yr
                !fire_rel_hum_threshold, & ! maximal relative humidity for fire
                !fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                !fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                !fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                !fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                !fire_tau_grass,         & ! return period of fire for grass PFT [year] at 0% relative humidity

    delta_time_yr = 1._wp / 365._wp

    burned_fract = fire_minimum * delta_time_yr

    IF (       q_rel_air_climbuf < fire_rel_hum_threshold   &
         & .AND. fuel > fire_litter_threshold               &
       ) THEN
      burned_fract = burned_fract                                                                        &
        &           + (fire_rel_hum_threshold - q_rel_air_climbuf) / (fire_tau * fire_rel_hum_threshold) &
        &             * delta_time_yr
    END IF

  END SUBROUTINE burned_fractt_jsb3

  ! Calculate windbreak for woody types
  PURE ELEMENTAL SUBROUTINE broken_woody_fract_jsb3( &
    & wind_threshold,                                & ! in
    & wind_damage_scale,                             & ! in
    & cover_fract_pot,                               & ! in
    & prev_day_max_wind_10m,                         & ! in
    & max_wind_10m,                                  & ! in
    & damaged_fract                                  & ! out
    & )

    ! Arguments
    REAL(wp), INTENT(in) ::    &
      & wind_threshold,        &
      & wind_damage_scale,     &
      & cover_fract_pot,       &
      & prev_day_max_wind_10m, &
      & max_wind_10m

    REAL(wp), INTENT(out) ::   &
      & damaged_fract

    ! Locals
    REAL(wp) :: delta_time_yr

    delta_time_yr = 1._wp / 365._wp

    IF (      cover_fract_pot > EPSILON(1._wp)                      &
        .AND. prev_day_max_wind_10m > max_wind_10m * wind_threshold &
       ) THEN
      damaged_fract = MIN( 1._wp - EPSILON(1._wp) / cover_fract_pot,                                          &
        &                 delta_time_yr * wind_damage_scale * prev_day_max_wind_10m ** 3._wp / max_wind_10m   &
        &               )
    ENDIF

  END SUBROUTINE broken_woody_fract_jsb3

  FUNCTION get_relative_humidity_air(nc, air_moisture, air_temperature, air_press) RESULT(relative_humidity_air)

    USE mo_physical_constants_iface, ONLY: rd  ! gas constant for dry air [J/(K*kg)]
    USE mo_physical_constants_iface, ONLY: rv  ! gas constant for water vapor [J/(K*kg)]
    !USE mo_jsb4_forcing,             ONLY: sat_specific_humidity ! Function

    INTEGER, INTENT(in)  :: nc
    REAL(wp), INTENT(in) :: &
      air_moisture(:),      & ! Specific humidity of air [kg kg-1]
      air_temperature(:),   & ! Temperature of air [K]
      air_press(:)            ! Pressure of air [Pa]

    REAL(wp) :: air_qsat(SIZE(air_moisture)), &
                relative_humidity_air(SIZE(air_moisture))

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_relative_humidity_air'

    CALL sat_specific_humidity(nc,               & ! in
                               air_temperature,  & ! in
                               air_press,        & ! in
                               air_qsat          & ! out
                              )

        IF (ANY( air_qsat(:) > 0.99_wp*HUGE(1._wp) )) THEN
           CALL finish(TRIM(routine), 'lookup table overflow.')
        END IF

    relative_humidity_air(:) = 100._wp * air_moisture(:)                                       &
                               * ( (1._wp - rd / rv)  * air_qsat(:)  + (rd / rv) )             &
                               / ( ((1._wp - rd / rv) * air_moisture(:) + (rd / rv)) * air_qsat(:) )

  END FUNCTION get_relative_humidity_air

  ! R: TBD: I copied this subroutine from MO_JSB4_FORCING, because
  !    this mo is not compiled with icon/echam.
  !    Reiner has to decide if we change config or the place of this procedure.
  !    I made some changes because the dimensions (vector lenghts) were unnecessary
  !    cryptic. We could even skip "nc" if we would change the DO loop below...
  SUBROUTINE sat_specific_humidity(nc, temp, pressure, qsat)
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
    USE mo_exception,              ONLY: message_text
    USE mo_echam_convect_tables,   ONLY: tlucua, jptlucu1, jptlucu2
    USE mo_jsb_physical_constants, ONLY: rvd1

    INTEGER,  INTENT(IN) :: nc
    REAL(wp), INTENT(IN) :: temp(nc)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(nc)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(nc)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sat_specific_humidity'

    INTEGER  :: i
    INTEGER  :: it(nc)
    REAL(wp) :: tluc(nc)

    !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
    DO i = 1,nc
      IF ((temp(i) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)) .OR. (temp(i) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp))) THEN
        WRITE (message_text,*) &
          & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
          & '- temperature was: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
    ENDDO

    it = NINT(temp*1000._wp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
      qsat = tluc / (pressure - rvd1*tluc)
    ELSEWHERE
      qsat = HUGE(1._wp)
    END WHERE

  END SUBROUTINE sat_specific_humidity
  ! --------------------------------------------------------------------------------------------------------
  !   StW: Below here comes the process routines needed for the Spitfire fire model (Thonicke et al. 2010)
  ! --------------------------------------------------------------------------------------------------------

  ! Calculates available fuel on tile
  SUBROUTINE calc_fire_spitfire_fuel(nc,woody_pft,C4_pft,fuel_dens_tile,                                &
                                     fract_green_active,SurfArea2Vol,                                   &
                                     veg_fract_correction,c_green,c_reserve,c_wood,                     &
                                     c_acid_ag1,c_water_ag1,c_ethanol_ag1,c_nonsoluble_ag1,             &
                                     c_acid_ag2,c_water_ag2,c_ethanol_ag2,c_nonsoluble_ag2,             &
                                     fract_litter_wood_new,                                             &
                                     fract_1hr_wood,fract_10hr_wood,fract_100hr_wood,fract_1000hr_wood, &
                                     fuel,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,                    &
                                     c_green_g,livegrass,Fuel1To100hr,net_fuel,char_net_fuel,           &
                                     fuel_SurfArea2Vol,ratio_live_fuel,ratio_dead_fuel,                 &
                                     char_dens_fuel_avg)
  USE mo_jsb_physical_constants, ONLY: ratio_C_2_biomass
  USE mo_disturbance_constants,  ONLY: fract_1hr_wood_new, fract_10hr_wood_new, fract_100hr_wood_new, fract_1000hr_wood_new, &
                                       fbd10hr_2_1hr, fbd100hr_2_1hr, fbd_C3_livegrass, fbd_C4_livegrass, MINER_TOT
  USE mo_carbon_constants,       ONLY: fract_green_aboveGround, molarMassC_kg

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fuel'
  INTEGER,  INTENT(in)     :: nc                                                                   ! #pts in vectors
  LOGICAL,  INTENT(in)     :: woody_pft, C4_pft                                                    ! Parameters from lctlib
  REAL(wp), INTENT(in)     :: fuel_dens_tile                                                       ! Parameter from lctlib
  REAL(wp), INTENT(in)     :: fract_green_active, SurfArea2Vol(3)                                   ! Parameters from namelists
  REAL(wp), INTENT(in)     :: veg_fract_correction(:)                                              ! m2 canopy -> m2 veg area cnv
  REAL(wp), INTENT(in)     :: c_green(:), c_reserve(:), c_wood(:)                                  ! Living C-pools
  REAL(wp), INTENT(in)     :: c_acid_ag1(:), c_water_ag1(:), c_ethanol_ag1(:), c_nonsoluble_ag1(:) ! Green litter pools
  REAL(wp), INTENT(in)     :: c_acid_ag2(:), c_water_ag2(:), c_ethanol_ag2(:), c_nonsoluble_ag2(:) ! Woody litter pools
  REAL(wp), INTENT(in)     :: fract_litter_wood_new(:)                                              ! Fract of litter put since last
  REAL(wp), INTENT(inout)  :: fract_1hr_wood(:), fract_10hr_wood(:), fract_100hr_wood(:), fract_1000hr_wood(:) !Fract diff fuel clss
  REAL(wp), INTENT(out)    :: fuel(:), fuel_1hr(:), fuel_10hr(:), fuel_100hr(:), fuel_1000hr(:)    ! Fuel in diff fuel classes
  REAL(wp), INTENT(out)    :: c_green_g(:), livegrass(:), Fuel1To100hr(:), net_fuel(:), char_net_fuel(:), &
                              fuel_SurfArea2Vol(:), ratio_live_fuel(:), ratio_dead_fuel(:),               &
                              char_dens_fuel_avg(:)                                                ! Various fuel properties

  REAL(wp) :: canopy2activeveg(nc)
  REAL(wp) :: c_litter_green_ag_g(nc)
  REAL(wp) :: c_litter_wood_ag_g(nc)
  REAL(wp) :: ratio_C3_livegrass(nc), ratio_C4_livegrass(nc)
  REAL(wp) :: ratio_fbd(nc)
  REAL(wp) :: dens_fuel_avg(nc), dens_livegrass_avg(nc)

    ! Tile fraction of potentially burning area
    ! StW: Doesn't make any sense in "all separate" approach: 1 should be cover_fract(tile)/sum(cover_fract(all_tile_with_fire))
    canopy2activeveg(:) = veg_fract_correction(:) * 1._wp

    ! update fuel fractions
    fract_1hr_wood   (:) = fract_1hr_wood   (:) * (1._wp-fract_litter_wood_new(:)) + fract_1hr_wood_new    * fract_litter_wood_new(:)
    fract_10hr_wood  (:) = fract_10hr_wood  (:) * (1._wp-fract_litter_wood_new(:)) + fract_10hr_wood_new   * fract_litter_wood_new(:)
    fract_100hr_wood (:) = fract_100hr_wood (:) * (1._wp-fract_litter_wood_new(:)) + fract_100hr_wood_new  * fract_litter_wood_new(:)
    fract_1000hr_wood(:) = fract_1000hr_wood(:) * (1._wp-fract_litter_wood_new(:)) + fract_1000hr_wood_new * fract_litter_wood_new(:)

    ! Get total available fuel
    c_litter_green_ag_g(:) = molarMassC_kg * 1000._wp * (c_acid_ag1(:) + c_water_ag1(:) + c_ethanol_ag1(:) + c_nonsoluble_ag1(:))
    c_litter_wood_ag_g (:) = molarMassC_kg * 1000._wp * (c_acid_ag2(:) + c_water_ag2(:) + c_ethanol_ag2(:) + c_nonsoluble_ag2(:))
    c_green_g          (:) = molarMassC_kg * 1000._wp * c_green(:)

    ! Fuel in different fuel classes
    IF (woody_pft) THEN
      fuel_1hr   (:) = ( c_litter_green_ag_g(:) * fract_green_active  + &
                        fract_1hr_wood   (:) * c_litter_wood_ag_g(:)) * canopy2activeveg(:) / ratio_C_2_biomass
      fuel_10hr  (:) =  fract_10hr_wood  (:) * c_litter_wood_ag_g(:)  * canopy2activeveg(:) / ratio_C_2_biomass
      fuel_100hr (:) =  fract_100hr_wood (:) * c_litter_wood_ag_g(:)  * canopy2activeveg(:) / ratio_C_2_biomass
      fuel_1000hr(:) =  fract_1000hr_wood(:) * c_litter_wood_ag_g(:)  * canopy2activeveg(:) / ratio_C_2_biomass
      livegrass  (:) = 0._wp
    ELSE
      fuel_1hr   (:) = c_litter_green_ag_g(:) * canopy2activeveg(:) / ratio_C_2_biomass
      fuel_10hr  (:) = 0._wp
      fuel_100hr (:) = 0._wp
      fuel_1000hr(:) = 0._wp
      livegrass  (:) = fuel_1hr(:) * fract_green_aboveground
    ENDIF

    ! Aggregation of litter fuel
    fuel         (:) = fuel_1hr(:) + fuel_10hr(:) + fuel_100hr(:) + fuel_1000hr(:)
    Fuel1To100hr (:) = fuel_1hr(:) + fuel_10hr(:) + fuel_100hr(:)
    net_fuel     (:) = (1.0_wp - MINER_TOT) * Fuel1To100hr(:) / 1000.0_wp
    char_net_fuel(:) = (1.0_wp - MINER_TOT) * livegrass   (:) / 1000.0_wp

    ! Surface area to volumne ratio of dead fuel
    WHERE (Fuel1To100hr(:) > EPSILON(1._wp))
      fuel_SurfArea2Vol(:) = (  fuel_1hr  (:) * SurfArea2Vol(1) &
                              + fuel_10hr (:) * SurfArea2Vol(2) &
                              + fuel_100hr(:) * SurfArea2Vol(3) &
                             ) / Fuel1To100hr(:)
    ELSEWHERE
      fuel_SurfArea2Vol(:) = 0.00001_wp
    ENDWHERE

    ! Ratios of living carbon for fuel calculations
    ratio_C3_livegrass(:) = 0.0_wp
    ratio_C4_livegrass(:) = 0.0_wp
    ratio_dead_fuel   (:) = 0.0_wp
    ratio_live_fuel   (:) = 0.0_wp
    ! StW: This doesn't make much sense for "all separated" approach
    IF (C4_pft) THEN
      WHERE (livegrass(:) > EPSILON(1._wp))
        ratio_C4_livegrass(:) = (c_green_g(:) / ratio_C_2_biomass &
                                * canopy2activeveg(:) * fract_green_aboveGround) / livegrass(:)
      END WHERE
    ENDIF
    WHERE (livegrass(:) > EPSILON(1._wp))
      ratio_C3_livegrass(:) = 1._wp - ratio_C4_livegrass(:)
    END WHERE
    WHERE (Fuel1To100hr(:) + livegrass(:) > EPSILON(1._wp))
      ratio_dead_fuel(:) = Fuel1To100hr(:) / (Fuel1To100hr(:) + livegrass(:))
      ratio_live_fuel(:) = livegrass   (:) / (Fuel1To100hr(:) + livegrass(:))
    END WHERE

    ! Fuel densities
    !      influence of livegrass on FBD
    !    Kirsten Thonicke in LPJ code: presence of livegrass additionally lowers FBD, when trees mixed
    !    with grasses fbd for C3 and C4 = 4
    dens_livegrass_avg(:) = fbd_C3_livegrass * ratio_C3_livegrass(:) &
                          + fbd_C4_livegrass * ratio_C4_livegrass(:)

    WHERE (Fuel1To100hr(:) > EPSILON(1._wp))
      ratio_fbd(:) = (  fuel_1hr  (:)                  &
                      + fuel_10hr (:) * fbd10hr_2_1hr  &
                      + fuel_100hr(:) * fbd100hr_2_1hr &
                     ) / Fuel1To100hr(:)
      dens_fuel_avg(:) = dens_fuel_avg(:) + ratio_fbd(:) * fuel_dens_tile
    ELSEWHERE
      dens_fuel_avg(:) = 0._wp
    ENDWHERE
    char_dens_fuel_avg(:) = ratio_dead_fuel(:) * dens_fuel_avg(:)      + &
                            ratio_live_fuel(:) * dens_livegrass_avg(:)

  END SUBROUTINE calc_fire_spitfire_fuel

  ! Calculates the humidity of available fuel
  SUBROUTINE calc_fire_spitfire_fuel_moist(nc,only10,moist_extinction,moisture_scaling,SurfArea2Vol,      &
                                           fuel,fuel_1hr,fuel_10hr,fuel_100hr,Fuel1To100hr,               &
                                           ratio_live_fuel,ratio_dead_fuel,SoilMoisture,NI_acc,           &
                                           char_moistfactor,leaf_moisture,char_alpha_fuel,RelFuelMoisture)
  USE mo_disturbance_constants, ONLY: moistfactor_livegrass

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fuel_moist'
  INTEGER,  INTENT(in)  :: nc
  LOGICAL,  INTENT(in)  :: only10
  REAL(wp), INTENT(in)  :: moist_extinction
  REAL(wp), INTENT(in)  :: moisture_scaling
  REAL(wp), INTENT(in)  :: SurfArea2Vol(3)
  REAL(wp), INTENT(in)  :: fuel(:)
  REAL(wp), INTENT(in)  :: fuel_1hr(:), fuel_10hr(:), fuel_100hr(:), Fuel1To100hr(:)
  REAL(wp), INTENT(in)  :: ratio_live_fuel(:)
  REAL(wp), INTENT(in)  :: ratio_dead_fuel(:)
  REAL(wp), INTENT(in)  :: SoilMoisture(:)
  REAL(wp), INTENT(in)  :: NI_acc(:)
  REAL(wp), INTENT(out) :: char_moistfactor(:), leaf_moisture(:)
  REAL(wp), INTENT(out) :: char_alpha_fuel(:)
  REAL(wp), INTENT(out) :: RelFuelMoisture(:)

  REAL(wp) :: alpha_fuel_1hr, alpha_fuel_10hr, alpha_fuel_100hr
  REAL(wp) :: alpha_fuel(nc), alpha_livegrass(nc)
  REAL(wp) :: moistfactor(nc)

    WHERE (fuel(:) > EPSILON(1._wp))
      ! StW: Should be fuel on tile / total fuel - this version does not make much sense
      moistfactor(:) = fuel(:) / fuel(:) * moist_extinction
    ELSEWHERE
      moistfactor(:) = 0._wp
    ENDWHERE

    ! leaf moisture content of grasses (RelFuelMoisture_lg in spitfire) eq. B2
    leaf_moisture(:) = MAX(0._wp,((10._wp/9._wp) * SoilMoisture(:) - (1._wp/9._wp)))

    char_moistfactor(:) = moistfactor(:) * ratio_dead_fuel(:) + moistfactor_livegrass * ratio_live_fuel(:)

    alpha_fuel_1hr   = SurfArea2Vol(1) / moisture_scaling
    alpha_fuel_10hr  = SurfArea2Vol(2) / moisture_scaling
    alpha_fuel_100hr = SurfArea2Vol(3) / moisture_scaling

    alpha_fuel(:) = 0._wp
    IF (only10) THEN
      WHERE (Fuel1To100hr(:) > EPSILON(1._wp))
        alpha_fuel(:) = ( alpha_fuel_1hr  * fuel_1hr (:) + &
                          alpha_fuel_10hr * fuel_10hr(:)   &
                        ) / (fuel_1hr + fuel_10hr(:))
      ENDWHERE
    ELSE
      WHERE (Fuel1To100hr(:) > EPSILON(1._wp))
        alpha_fuel(:) = (  alpha_fuel_1hr   * fuel_1hr  (:) &
                         + alpha_fuel_10hr  * fuel_10hr (:) &
                         + alpha_fuel_100hr * fuel_100hr(:) &
                        ) / Fuel1To100hr(:)
      ENDWHERE
    ENDIF

    alpha_livegrass(:) = 0.0_wp
    WHERE ((NI_acc(:) > EPSILON(1._wp)) .AND. (leaf_moisture(:) > EPSILON(1._wp)))
      alpha_livegrass(:) = -LOG(leaf_moisture(:)) / NI_acc(:)
    ENDWHERE
    char_alpha_fuel(:) = alpha_fuel(:) * ratio_dead_fuel(:) + alpha_livegrass(:) * ratio_live_fuel(:)

    ! weighted average of the relative moisture content eq. 6 in TH2010
    RelFuelMoisture(:) = EXP(-char_alpha_fuel(:) * NI_acc(:))  ! used for ROS calcs. 1000hr therefore ignored for alpha.

  END SUBROUTINE calc_fire_spitfire_fuel_moist

  ! calculate fire danger index d_fdi eq. 8 in TH2010
  PURE ELEMENTAL SUBROUTINE calc_fire_spitfire_fire_danger_index(Fuel1To100hr,char_moistfactor,char_alpha_fuel,NI_acc,FDI)
  REAL(wp), INTENT(in)  :: Fuel1To100hr
  REAL(wp), INTENT(in)  :: char_moistfactor
  REAL(wp), INTENT(in)  :: char_alpha_fuel
  REAL(wp), INTENT(in)  :: NI_acc
  REAL(wp), INTENT(out) :: FDI

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fire_danger_index'

    FDI = 0.0_wp
    IF ((NI_acc > EPSILON(1._wp)) .AND. (char_moistfactor > 0._wp) .AND. (Fuel1To100hr > EPSILON(1._wp))) &
      FDI = MAX(0.0_wp,(1.0_wp-((1.0_wp/char_moistfactor) * (EXP(-char_alpha_fuel*NI_acc)))))

  END SUBROUTINE calc_fire_spitfire_fire_danger_index

  ! Calculates the geometry (in space and time) of a single fire
  SUBROUTINE calc_fire_spitfire_fire_geometry(nc,lduration_popd,           &
      lwind_speed_limit,wind_max,wind_slope,WindSpeed, population_density, &
      fuel_SurfArea2Vol,RelFuelMoisture,                                   &
      net_fuel,char_dens_fuel_avg,char_net_fuel,char_moistfactor,          &
      wood_fract,grass_fract, FDI,                                         &
      gamma,IR,ROS,FireDuration,Area)
  USE mo_math_constants,        ONLY: pi
  USE mo_disturbance_constants, ONLY: ft_per_m, MINER_DAMP, part_dens, HeatContentFuel
  INTEGER , INTENT(in)  :: nc
  LOGICAL,  INTENT(in)  :: lduration_popd
  LOGICAL,  INTENT(in)  :: lwind_speed_limit
  REAL(wp), INTENT(in)  :: wind_max, wind_slope
  REAL(wp), INTENT(in)  :: WindSpeed(:)
  REAL(wp), INTENT(in)  :: population_density(:)
  REAL(wp), INTENT(in)  :: fuel_SurfArea2Vol(:), RelFuelMoisture(:)
  REAL(wp), INTENT(in)  :: net_fuel(:)
  REAL(wp), INTENT(in)  :: char_dens_fuel_avg(:)
  REAL(wp), INTENT(in)  :: char_net_fuel(:)
  REAL(wp), INTENT(in)  :: char_moistfactor(:)
  REAL(wp), INTENT(in)  :: wood_fract(:)
  REAL(wp), INTENT(in)  :: grass_fract(:)
  REAL(wp), INTENT(in)  :: FDI(:)
  REAL(wp), INTENT(out) :: gamma(:)        ! used for FireResidenceTime for mortality
  REAL(wp), INTENT(out) :: IR(:)           ! reaction intensity
  REAL(wp), INTENT(out) :: ROS(:)          ! rate of spread
  REAL(wp), INTENT(out) :: FireDuration(:) ! Duration of a fire [min]
  REAL(wp), INTENT(out) :: Area(:)         ! Area of a single fire [m^2]

  ! locals
  REAL(wp) :: beta(nc), beta_op(nc), bet(nc)
  REAL(wp) :: q_ig(nc)               ! heat of pre-ignition
  REAL(wp) :: eps(nc)                ! effective heating number
  REAL(wp) :: b(nc), c(nc), e(nc), a(nc)
  REAL(wp) :: phi_wind(nc)
  REAL(wp) :: xi(nc)
  REAL(wp) :: dummy(nc), gamma_max(nc),gamma_aptr(nc)
  REAL(wp) :: mw_weight(nc), moist_damp(nc)
  REAL(wp) :: ROS_b(nc), LengthOfMajAxis(nc)
  REAL(wp) :: WindSpeed_m_min(nc), WindSpeed_ft_min(nc)
  REAL(wp) :: length2breadthratio(nc)

    ! fire duration as a function of Fire danger index and optinal population density
    FireDuration(:) = 241._wp / (1.0_wp + (240._wp * EXP(-11.06_wp * FDI(:)))) * 3.0_wp
    IF (lduration_popd) THEN
      FireDuration(:) = FireDuration(:) * (-LOG10(MIN(MAX(population_density(:),0.01_wp),100._wp)) + 4_wp) * 0.5_wp
    ENDIF

    ! moisture dampening coefficient, see tab. A1
    WHERE (char_moistfactor(:) > 0._wp)
      mw_weight(:) = RelFuelMoisture(:) / char_moistfactor(:)
    ELSEWHERE
      mw_weight(:) = 0._wp
    ENDWHERE

    WHERE ((fuel_SurfArea2Vol(:) > 0.0001_wp) .AND. (char_dens_fuel_avg(:) > EPSILON(1._wp)))

      ! influence of wind speed (eq. A5-A9)
      beta   (:) = char_dens_fuel_avg(:) / part_dens
      beta_op(:) = 0.200395_wp * (fuel_SurfArea2Vol(:)**(-0.8189_wp))
      bet    (:) = beta(:) / beta_op(:)

      a(:) = 8.9033_wp  * (fuel_SurfArea2Vol(:)**(-0.7913_wp))
      b(:) = 0.15988_wp * (fuel_SurfArea2Vol(:)**0.54_wp)
      c(:) = 7.47_wp    * (EXP(-0.8711_wp  * (fuel_SurfArea2Vol(:)**0.55_wp)))
      e(:) = 0.715_wp   * (EXP(-0.01094_wp *  fuel_SurfArea2Vol(:)))

      ! heat of pre-ignition eq. A4 in TH2010
      q_ig(:) = 581.0_wp + 2594.0_wp * RelFuelMoisture(:)

      ! effective heating number eq. A3 in  TH2010
      eps(:) = EXP(-4.528_wp/fuel_SurfArea2Vol(:))

      ! propagating flux ratio xi Eq. A2 in TH2010
      xi(:) = EXP((0.792_wp + 3.7597_wp * (fuel_SurfArea2Vol(:)**0.5_wp)) * (beta(:) + 0.1_wp)) / &
                  (192._wp  + 7.9095_wp *  fuel_SurfArea2Vol(:))

      ! optimum reaction velocity table A1 in TH2010
      dummy     (:) = EXP(a(:) * (1._wp - bet(:)))
      gamma_max (:) = 1._wp / (0.0591_wp + 2.926_wp * (fuel_SurfArea2Vol(:)**(-1.5_wp))) ! maximum reaction velocity
      gamma_aptr(:) = gamma_max(:) * (bet(:)**a(:)) * dummy(:)               ! optimum reaction velocity

      moist_damp(:) = MAX(0._wp,(  1._wp - (2.59_wp * mw_weight(:)) &
                                 + (5.11_wp * (mw_weight(:)**2._wp)) - (3.52_wp * (mw_weight(:)**3._wp))))

      ! reaction intensity IR eq. A1 in TH2010
      IR(:) = gamma_aptr(:) * char_net_fuel(:) * HeatContentFuel * moist_damp(:) * MINER_DAMP

      ! for use in postfire mortality (computation of FireResidenceTime)
      gamma(:) = gamma_aptr(:) * moist_damp(:) * MINER_DAMP * char_net_fuel(:)

    ELSEWHERE
      IR   (:) = 0._wp
      gamma(:) = 0._wp
      bet  (:) = 0._wp
      q_ig (:) = 0._wp
      eps  (:) = 0._wp
    ENDWHERE

    ! wind speed in m/min needed instead of m/s and reduction with respect to grass_fpc (0.6) and wood_fpc (0.4) (TH2010 p.1997)
    WindSpeed_m_min(:) = (  (wood_fract (:) * 0.4_wp) &
                          + (grass_fract(:) * 0.6_wp) &
                         ) * WindSpeed(:) * 60._wp
    WindSpeed_ft_min(:) = WindSpeed_m_min(:) * ft_per_m ! Windspeed in ft/min

    ! length to breadth ratio of burned ellipse according to presence of grasses and forests eq. 12 and 13 in TH2010
    ! in paper units of wind in eq 12,13 is m/min in code factor 0.06 included --> km/h, in original reference km/h is correct
    WHERE (WindSpeed_m_min(:) < 16.67_wp)
      Length2BreadthRatio(:) = 1.0_wp
    ELSEWHERE
      Length2BreadthRatio(:) = MAX(MIN(                                                                        &
        wood_fract (:) * (1.0_wp+(8.729_wp*((1.0_wp-(EXP(-0.03_wp*0.06_wp*WindSpeed_m_min(:))))**2.155_wp))) + &
        grass_fract(:) * (1.1_wp+((0.06_wp*WindSpeed_m_min(:))**0.0464_wp)),                                   &
        8._wp),0.001_wp)
    ENDWHERE

    phi_wind(:) = 0._wp
    WHERE (bet(:) > 0._wp)
      phi_wind(:) = c(:) * (WindSpeed_ft_min(:)**b(:)) * (bet(:)**(-e(:)))
    ENDWHERE

    ! wind speed limitation
    IF (lwind_speed_limit) THEN
      WHERE ((IR(:) > 0._wp) .AND. (WindSpeed_ft_min(:) > wind_max))
        phi_wind(:) = c(:)*(MAX(wind_max*(1._wp + wind_slope) - wind_slope*WindSpeed_ft_min(:),0._wp)**b(:))*(bet(:)**(-e(:)))
      ENDWHERE
    ENDIF

    ! rate of spread eq. 9in TH2010
    ROS(:) = 0._wp
    WHERE ((char_dens_fuel_avg(:) > 0._wp) .AND. (eps(:) > 0._wp) .AND. (q_ig(:) > 0._wp))
      ROS(:) = (IR(:) * xi(:) * (1._wp + phi_wind(:))) / (char_dens_fuel_avg(:) * eps(:) * q_ig(:))
    ENDWHERE

    ! backward rate of spread
    ROS_b(:) = ROS(:) * EXP(-0.012_wp * WindSpeed_m_min(:)) ! windspeed in m/min and averaged for grass and trees
    WHERE (ROS_b(:) < 0.05_wp)
      ROS_b(:) = 0.0_wp
    ENDWHERE

    ! HIRSCH (1996)
    WHERE (net_fuel(:) > EPSILON(1._wp))
      LengthOfMajAxis(:) = (ROS_b(:) + ROS(:)) * FireDuration(:) ! in min
    ELSEWHERE
      LengthOfMajAxis(:) = 0._wp
      FireDuration   (:) = 0._wp
      ROS            (:) = 0._wp
    ENDWHERE

    ! Area of a fire
    area(:) = 0._wp
    WHERE ((Length2BreadthRatio(:) > EPSILON(1._wp)) .AND. (ROS(:) > EPSILON(1._wp)))
      area(:) = ((pi / (4._wp * Length2BreadthRatio(:))) * ((LengthOfMajAxis(:))**2._wp)) / 10000._wp ! burned area per ha
    ENDWHERE

  END SUBROUTINE calc_fire_spitfire_fire_geometry

  ! Calculate fire density [#/ha/day]
  SUBROUTINE calc_fire_spitfire_fire_density(ign_para,human_para,const_ignitions,human_ignitions,lightning_c2g, &
                                             has_fire,FDI,fire_dens)
  REAL(wp), INTENT(in)  :: ign_para
  REAL(wp), INTENT(in)  :: human_para
  REAL(wp), INTENT(in)  :: const_ignitions
  REAL(wp), INTENT(in)  :: human_ignitions(:)
  REAL(wp), INTENT(in)  :: lightning_c2g(:)
  LOGICAL,  INTENT(in)  :: has_fire(:)
  REAL(wp), INTENT(in)  :: FDI(:)
  REAL(wp), INTENT(out) :: fire_dens(:)

    ! only a fraction of cloud to ground flashes ignites a fire
    IF (const_ignitions < 0._wp) THEN
      WHERE (has_fire(:))
        ! probability of human ignitions +probability of lightning ignition per ha
        fire_dens(:) = ign_para * 0.01 * FDI(:) * (human_para * human_ignitions(:) +  &
                                                   0.04_wp    * lightning_c2g  (:))
      ELSEWHERE  !not enough fuel
        fire_dens(:) = 0.0_wp
      ENDWHERE
    ELSE
      WHERE (has_fire(:))
        ! constant probability of ignitions per ha
        fire_dens(:) = FDI(:) * const_ignitions * 0.01_wp
      ELSEWHERE  !not enough fuel
        fire_dens(:) = 0.0_wp
      ENDWHERE
    ENDIF

  END SUBROUTINE calc_fire_spitfire_fire_density

  PURE ELEMENTAL SUBROUTINE calc_fire_spitfire_fuel_consumption(NI_acc, burned_fract,       &
                              fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr, fuel_1hr_total, &
                              RelFuelMoisture, leaf_moisture, livegrass, char_moistfactor,  &
                              fuel_consum, fuel_consum1000perAreaBurned,                    &
                              CombustCompleteness_lg,                                       &
                              CombustCompleteness_1hr, CombustCompleteness_10hr,            &
                              CombustCompleteness_100hr, CombustCompleteness_1000hr)
  USE mo_disturbance_constants, ONLY: MINER_TOT, alpha_1hr

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fuel_consumption'
  REAL(wp), INTENT(in)  :: NI_acc
  REAL(wp), INTENT(in)  :: RelFuelMoisture
  REAL(wp), INTENT(in)  :: fuel_1hr, fuel_10hr, fuel_100hr, fuel_1000hr
  REAL(wp), INTENT(in)  :: fuel_1hr_total                ! StW: Should be from from veg tile
  REAL(wp), INTENT(in)  :: burned_fract                  ! fraction burned per pft, grid cell
  REAL(wp), INTENT(in)  :: leaf_moisture
  REAL(wp), INTENT(in)  :: livegrass
  REAL(wp), INTENT(in)  :: char_moistfactor              ! extinction moist, weighted with rel amount of fuel per tile
  REAL(wp), INTENT(out) :: fuel_consum                                           ! fuel consumption per vegetated area
  REAL(wp), INTENT(out) :: fuel_consum1000perAreaBurned                          ! fuel consumption per area burned
  REAL(wp), INTENT(out) :: CombustCompleteness_lg
  REAL(wp), INTENT(out) :: CombustCompleteness_1hr,CombustCompleteness_10hr
  REAL(wp), INTENT(out) :: CombustCompleteness_100hr, CombustCompleteness_1000hr

  ! local variables:
  REAL(wp) :: fc_1hr, fc_10hr
  REAL(wp) :: fc_100hr, fc_1000hr
  REAL(wp) :: moist_lg_d1hr, moist_1hr, moist_10_100hr
  REAL(wp) :: RelFuelMoisture_1hr

    ! initialize variables
    fuel_consum1000perAreaBurned = 0._wp
    fuel_consum                  = 0._wp
    CombustCompleteness_lg       = 0._wp
    CombustCompleteness_1hr      = 0._wp
    CombustCompleteness_10hr     = 0._wp
    CombustCompleteness_100hr    = 0._wp
    CombustCompleteness_1000hr   = 0._wp
    fc_1hr                       = 0._wp
    fc_10hr                      = 0._wp
    fc_100hr                     = 0._wp
    fc_1000hr                    = 0._wp

    RelFuelMoisture_1hr = EXP(-alpha_1hr * NI_acc)
    moist_lg_d1hr = 0._wp
    IF (fuel_1hr_total > 0._wp) THEN
      ! GL 3/2015 modification of moisture averaging, now different from SPITFIRE documentation but more reasonable
      moist_lg_d1hr = (RelFuelMoisture_1hr * fuel_1hr_total + leaf_moisture * livegrass) &
                       / (fuel_1hr_total + livegrass)
    ENDIF
    moist_1hr      = 1._wp
    moist_10_100hr = 1._wp
    IF (char_moistfactor > 0._wp) THEN
      moist_1hr      = MIN(moist_lg_d1hr   / char_moistfactor,1.0_wp)
      moist_10_100hr = MIN(RelFuelMoisture / char_moistfactor,1.0_wp)
    ENDIF

    ! Combustion factors
    ! combustion completeness of live grass
    IF (leaf_moisture <= 0.18_wp) THEN
      CombustCompleteness_lg = 1._wp * (1.0_wp - MINER_TOT)
    ELSEIF ((leaf_moisture > 0.18_wp) .AND. (leaf_moisture <= 0.73_wp)) THEN
      CombustCompleteness_lg = (1.11_wp -  0.62_wp * leaf_moisture) * (1._wp - MINER_TOT)
    ELSEIF ((leaf_moisture > 0.73_wp) .AND. (leaf_moisture <= 1._wp)) THEN
      CombustCompleteness_lg = (2.45_wp - 2.45_wp * leaf_moisture) * (1._wp - MINER_TOT)
    ENDIF
    ! combustion completeness of 1hr fuels
    IF (moist_1hr <= 0.18_wp) THEN
      CombustCompleteness_1hr = 1.0_wp * (1.0_wp - MINER_TOT)
    ELSEIF ((moist_1hr > 0.18_wp) .AND. (moist_1hr <= 0.73_wp)) THEN
      CombustCompleteness_1hr = (1.11_wp - 0.62_wp * moist_1hr) * (1.0_wp - MINER_TOT)
    ELSEIF ((moist_1hr > 0.73_wp) .AND. (moist_1hr <= 1._wp)) THEN
      CombustCompleteness_1hr = (2.45_wp - 2.45_wp * moist_1hr) * (1.0_wp - MINER_TOT)
    ENDIF
    ! 10hr fuel combustion completeness
    IF (moist_10_100hr <= 0.13_wp) THEN
      CombustCompleteness_10hr = 1._wp * (1._wp - MINER_TOT)
    ELSEIF ((moist_10_100hr > 0.13_wp) .AND. (moist_10_100hr <= 0.51_wp)) THEN
      CombustCompleteness_10hr = (1.09_wp - 0.72_wp * moist_10_100hr) * (1._wp - MINER_TOT)
    ELSEIF ((moist_10_100hr > 0.51_wp) .AND. (moist_10_100hr <= 1._wp)) THEN
      CombustCompleteness_10hr = (1.47_wp - 1.47_wp * moist_10_100hr) * (1._wp - MINER_TOT)
    ENDIF
    ! 100hr fuel combustion completeness
    CombustCompleteness_100hr = MIN(0.45_wp,(0.98_wp - 0.85_wp * moist_10_100hr)) * (1._wp - MINER_TOT)
    ! 1000hr fuel combustion completeness
    CombustCompleteness_1000hr = MIN(0.45_wp,(0.8_wp - 0.8_wp  * moist_10_100hr)) * (1._wp - MINER_TOT)

    ! Tiled fuel consumption
    fc_1hr    = CombustCompleteness_1hr    * fuel_1hr
    fc_10hr   = CombustCompleteness_10hr   * fuel_10hr
    fc_100hr  = CombustCompleteness_100hr  * fuel_100hr
    fc_1000hr = CombustCompleteness_1000hr * fuel_1000hr

    ! Aggregation
    fuel_consum1000perAreaBurned = fc_1hr + fc_10hr + fc_100hr + fc_1000hr
    fuel_consum =  burned_fract * (fc_1hr + fc_10hr + fc_100hr)

  END SUBROUTINE calc_fire_spitfire_fuel_consumption

  ! vegetation height (Shevlikova et al. 2009, GBC), h_max=24.19 m, hl=0.19 m2/kgC
  SUBROUTINE calc_vegetation_height(woody_pft,veg_fract_correction,c_green,c_reserve,c_woods,veg_height)
  USE mo_disturbance_constants, ONLY: hl, hmax
  USE mo_carbon_constants,      ONLY: molarMassC_kg

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_vegetation_height'
  LOGICAL,  INTENT(in)  :: woody_pft
  REAL(wp), INTENT(in)  :: veg_fract_correction(:)
  REAL(wp), INTENT(in)  :: c_green(:)
  REAL(wp), INTENT(in)  :: c_reserve(:)
  REAL(wp), INTENT(in)  :: c_woods(:)
  REAL(wp), INTENT(out) :: veg_height(:)

    IF (woody_pft) THEN
      veg_height(:) = hmax * (1.0_wp - EXP(-hl * veg_fract_correction(:) &
                                     * (c_green(:) + c_woods(:) + c_reserve(:)) * molarMassC_kg))
    ELSE
      veg_height(:) = 0._wp
    ENDIF

  END SUBROUTINE calc_vegetation_height

  ! Fire intensity according to eq. 15 in TH2010
  PURE ELEMENTAL SUBROUTINE calc_fire_spitfire_fire_intensity(burned_fract,fuel_consum,ROS,I_surface)
  USE mo_disturbance_constants, ONLY: HeatContentFuel

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fire_intensity'
  REAL(wp), INTENT(in)  :: burned_fract
  REAL(wp), INTENT(in)  :: fuel_consum
  REAL(wp), INTENT(in)  :: ROS
  REAL(wp), INTENT(out) :: I_surface

    I_surface = 0._wp
    IF ((burned_fract > EPSILON(1._wp)) .AND. (fuel_consum > EPSILON(1._wp))) THEN
      I_surface = HeatContentFuel * (fuel_consum / 1000._wp / burned_fract) * (ROS / 60._wp)
    ENDIF

  END SUBROUTINE calc_fire_spitfire_fire_intensity

  SUBROUTINE calc_fire_spitfire_fire_mortality(nc, woody_pft,                       &
                            crown_length_pft, bark_par1, bark_par2, RCK, mort_prob, &
                            flame_length_f, flaming2total,                          &
                            I_surface, fuel_consum,gamma,vegetation_height,         &
                            mortality_crownkill,mortality_camb,scorch_height,postfire_mortality)
  USE mo_disturbance_constants, ONLY: allom2, allom3

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_fire_mortality'
  INTEGER,            INTENT(in)  :: nc
  LOGICAL,            INTENT(in)  :: woody_pft
  REAL(wp),           INTENT(in)  :: crown_length_pft
  REAL(wp),           INTENT(in)  :: bark_par1
  REAL(wp),           INTENT(in)  :: bark_par2
  REAL(wp),           INTENT(in)  :: RCK
  REAL(wp),           INTENT(in)  :: mort_prob
  REAL(wp),           INTENT(in)  :: flame_length_f
  REAL(wp),           INTENT(in)  :: flaming2total            ! factor to calculate total burning time from flaming phase duration
  REAL(wp),           INTENT(in)  :: I_surface(:)             ! fire intensity
  REAL(wp),           INTENT(in)  :: fuel_consum(:)
  REAL(wp),           INTENT(in)  :: gamma(:)                 ! gamma:max reaction velocity*mineral dampening*moisture dampening
  REAL(wp),           INTENT(in)  :: vegetation_height(:)
  REAL(wp),           INTENT(out) :: mortality_crownkill(:)
  REAL(wp),           INTENT(out) :: mortality_camb(:)      ! fraction killed by cambial kill [0,1]
  REAL(wp),           INTENT(out) :: scorch_height(:)       ! [m]
  REAL(wp),           INTENT(out) :: postfire_mortality(:)

  ! local variables
  REAL(wp) :: FireResidenceTime(nc)                  ! residence time of fire
  REAL(wp) :: crownkill_class(5)
  REAL(wp) :: mortality_crownkill_class(5)
  REAL(wp) :: mortality_camb_class(5)
  REAL(wp) :: tau_c,bark_thickness,DBH,stem_length  ! tau_c: critical time of fire residence
  REAL(wp) :: crown_length
  REAL(wp) :: height
  INTEGER  :: i,clss

    ! initializations
    scorch_height      (:) = 0._wp
    FireResidenceTime  (:) = 0._wp
    mortality_crownkill(:) = 1._wp
    mortality_camb     (:) = 1._wp

    ! scorch height, FireResidenceTime residence time of fire (compare reference (2) below)
    IF (woody_pft) THEN
      DO i=1,nc
        IF ((vegetation_height(i) > 0._wp) .AND. (gamma(i) > 0._wp)) THEN
          scorch_height(i) = flame_length_f * (I_surface(i)**0.667_wp)
          ! FireResidenceTime: residence time of fire
          FireResidenceTime(i) = flaming2total * fuel_consum(i) / 1000._wp / gamma(i)
          ! crownkill
          ! Version using crownkill classes
          crownkill_class          (1:5) = 0._wp
          mortality_crownkill_class(1:5) = 0._wp
          mortality_camb_class     (1:5) = 0._wp
          DO clss=1,5
            height       = vegetation_height(i) * (REAL(clss,wp) * 0.25_wp + 0.25_wp)
            crown_length = crown_length_pft * height
            stem_length  = height - crown_length
            ! crown kill [0, 1]
            IF (scorch_height(i) < stem_length) THEN
              crownkill_class(clss) = 0._wp
            ELSEIF ((scorch_height(i) >= stem_length) .AND. (scorch_height(i) < height)) THEN
              crownkill_class(clss) = (scorch_height(i) - stem_length) / crown_length
            ELSE
              crownkill_class(clss) = 1._wp
            ENDIF
            mortality_crownkill_class(clss) = RCK * crownkill_class(clss)**mort_prob
            ! mortality cambial damage
            ! diameter at breast height
            DBH            = (height / allom2)**(1._wp / allom3) * 100._wp
            ! bark_thickness
            bark_thickness = bark_par1 * DBH + bark_par2
            tau_c          = 2.9_wp * bark_thickness**2._wp ! critical time for cambial damage
            ! eq. 19 in TH2010
            IF     (FireResidenceTime(i) / tau_c >= 2._wp) THEN
              mortality_camb_class(clss) = 1._wp
            ELSEIF (FireResidenceTime(i) / tau_c >  0.22_wp) THEN
              mortality_camb_class(clss) = (0.563_wp * FireResidenceTime(i) / tau_c ) - 0.125_wp
            ELSE
              mortality_camb_class(clss) = 0._wp
            ENDIF
          ENDDO ! class
        ENDIF
      ENDDO ! nc
      postfire_mortality(:) = mortality_crownkill(:) + mortality_camb(:) &
                            - mortality_crownkill(:) * mortality_camb(:)
    ELSE
      postfire_mortality(:) = 1._wp
    ENDIF

  END SUBROUTINE calc_fire_spitfire_fire_mortality

  PURE ELEMENTAL SUBROUTINE calc_fire_spitfire_pool_consumption(woody_pft,fract_green_active,burned_frac, &
                                                                CombustCompleteness_lg,                   &
                                                                CombustCompleteness_1hr,                  &
                                                                CombustCompleteness_10hr,                 &
                                                                CombustCompleteness_100hr,                &
                                                                CombustCompleteness_1000hr,               &
                                                                c_green,                                  &
                                                                c_acid_ag1,                               &
                                                                c_water_ag1,                              &
                                                                c_ethanol_ag1,                            &
                                                                c_nonsoluble_ag1,                         &
                                                                c_acid_ag2,                               &
                                                                c_water_ag2,                              &
                                                                c_ethanol_ag2,                            &
                                                                c_nonsoluble_ag2,                         &
                                                                fract_1hr_wood,                           &
                                                                fract_10hr_wood,                          &
                                                                fract_100hr_wood,                         &
                                                                fract_1000hr_wood,                        &
                                                                burned_litter_wood,                       &
                                                                burned_litter_green,                      &
                                                                burned_c_green,                           &
                                                                burned_acid_ag1,                          &
                                                                burned_water_ag1,                         &
                                                                burned_ethanol_ag1,                       &
                                                                burned_nonsoluble_ag1,                    &
                                                                burned_acid_ag2,                          &
                                                                burned_water_ag2,                         &
                                                                burned_ethanol_ag2,                       &
                                                                burned_nonsoluble_ag2 )
  USE mo_carbon_constants,      ONLY: fract_green_aboveGround
  USE mo_disturbance_constants, ONLY: fract_1hr_wood_new, fract_10hr_wood_new, fract_100hr_wood_new, fract_1000hr_wood_new

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fire_spitfire_pool_consumption'
  LOGICAL,  INTENT(in)    :: woody_pft
  REAL(wp), INTENT(in)    :: fract_green_active
  REAL(wp), INTENT(in)    :: CombustCompleteness_1hr,CombustCompleteness_10hr,CombustCompleteness_100hr
  REAL(wp), INTENT(in)    :: CombustCompleteness_1000hr,CombustCompleteness_lg
  REAL(wp), INTENT(in)    :: burned_frac
  REAL(wp), INTENT(in)    :: c_green
  REAL(wp), INTENT(in)    :: c_acid_ag1,       c_acid_ag2
  REAL(wp), INTENT(in)    :: c_water_ag1,      c_water_ag2
  REAL(wp), INTENT(in)    :: c_ethanol_ag1,    c_ethanol_ag2
  REAL(wp), INTENT(in)    :: c_nonsoluble_ag1, c_nonsoluble_ag2
  REAL(wp), INTENT(inout) :: fract_1hr_wood,fract_10hr_wood
  REAL(wp), INTENT(inout) :: fract_100hr_wood,fract_1000hr_wood
  REAL(wp), INTENT(out)   :: burned_acid_ag1, burned_acid_ag2, burned_water_ag1, burned_water_ag2
  REAL(wp), INTENT(out)   :: burned_ethanol_ag1,burned_ethanol_ag2,burned_nonsoluble_ag1,burned_nonsoluble_ag2
  REAL(wp), INTENT(out)   :: burned_litter_green, burned_litter_wood, burned_c_green

  ! local variables
  REAL(wp) :: burned_acid_leaves,        burned_water_leaves
  REAL(wp) :: burned_ethanol_leaves,     burned_nonsoluble_leaves
  REAL(wp) :: burned_acid_ag1_1hr,       burned_acid_ag2_1hr
  REAL(wp) :: burned_water_ag1_1hr,      burned_water_ag2_1hr
  REAL(wp) :: burned_ethanol_ag1_1hr,    burned_ethanol_ag2_1hr
  REAL(wp) :: burned_nonsoluble_ag1_1hr, burned_nonsoluble_ag2_1hr
  REAL(wp) :: burned_acid_ag2_10hr
  REAL(wp) :: burned_water_ag2_10hr
  REAL(wp) :: burned_ethanol_ag2_10hr
  REAL(wp) :: burned_nonsoluble_ag2_10hr
  REAL(wp) :: burned_acid_ag2_100hr
  REAL(wp) :: burned_water_ag2_100hr
  REAL(wp) :: burned_ethanol_ag2_100hr
  REAL(wp) :: burned_nonsoluble_ag2_100hr
  REAL(wp) :: burned_acid_ag2_1000hr
  REAL(wp) :: burned_water_ag2_1000hr
  REAL(wp) :: burned_ethanol_ag2_1000hr
  REAL(wp) :: burned_nonsoluble_ag2_1000hr
  REAL(wp) :: wood_fac_1hr     , wood_fac_10hr
  REAL(wp) :: wood_fac_100hr   , wood_fac_1000hr
  REAL(wp) :: litter_wood !, litter_green
  REAL(wp) :: burned_wood !, burned_green
!  REAL(wp) :: burned_leaves
  REAL(wp) :: burned_wood_1hr  , burned_wood_10hr
  REAL(wp) :: burned_wood_100hr, burned_wood_1000hr

    ! Initialization
    burned_c_green               = 0._wp
    burned_acid_ag1_1hr          = 0._wp
    burned_water_ag1_1hr         = 0._wp
    burned_ethanol_ag1_1hr       = 0._wp
    burned_nonsoluble_ag1_1hr    = 0._wp
    burned_acid_leaves           = 0._wp
    burned_water_leaves          = 0._wp
    burned_ethanol_leaves        = 0._wp
    burned_nonsoluble_leaves     = 0._wp
    burned_acid_ag2_1hr          = 0._wp
    burned_water_ag2_1hr         = 0._wp
    burned_ethanol_ag2_1hr       = 0._wp
    burned_nonsoluble_ag2_1hr    = 0._wp
    burned_acid_ag2_10hr         = 0._wp
    burned_water_ag2_10hr        = 0._wp
    burned_ethanol_ag2_10hr      = 0._wp
    burned_nonsoluble_ag2_10hr   = 0._wp
    burned_acid_ag2_100hr        = 0._wp
    burned_water_ag2_100hr       = 0._wp
    burned_ethanol_ag2_100hr     = 0._wp
    burned_nonsoluble_ag2_100hr  = 0._wp
    burned_acid_ag2_1000hr       = 0._wp
    burned_water_ag2_1000hr      = 0._wp
    burned_ethanol_ag2_1000hr    = 0._wp
    burned_nonsoluble_ag2_1000hr = 0._wp

    wood_fac_1hr    = CombustCompleteness_1hr    * fract_1hr_wood    * burned_frac
    wood_fac_10hr   = CombustCompleteness_10hr   * fract_10hr_wood   * burned_frac
    wood_fac_100hr  = CombustCompleteness_100hr  * fract_100hr_wood  * burned_frac
    wood_fac_1000hr = CombustCompleteness_1000hr * fract_1000hr_wood * burned_frac

    ! compute combusted litter in pools
    IF (woody_pft) THEN
      burned_acid_leaves           = c_acid_ag1       * fract_green_active * CombustCompleteness_1hr * burned_frac
      burned_water_leaves          = c_water_ag1      * fract_green_active * CombustCompleteness_1hr * burned_frac
      burned_ethanol_leaves        = c_ethanol_ag1    * fract_green_active * CombustCompleteness_1hr * burned_frac
      burned_nonsoluble_leaves     = c_nonsoluble_ag1 * fract_green_active * CombustCompleteness_1hr * burned_frac
      burned_acid_ag2_1hr          = c_acid_ag2       * wood_fac_1hr
      burned_water_ag2_1hr         = c_water_ag2      * wood_fac_1hr
      burned_ethanol_ag2_1hr       = c_ethanol_ag2    * wood_fac_1hr
      burned_nonsoluble_ag2_1hr    = c_nonsoluble_ag2 * wood_fac_1hr
      burned_acid_ag2_10hr         = c_acid_ag2       * wood_fac_10hr
      burned_water_ag2_10hr        = c_water_ag2      * wood_fac_10hr
      burned_ethanol_ag2_10hr      = c_ethanol_ag2    * wood_fac_10hr
      burned_nonsoluble_ag2_10hr   = c_nonsoluble_ag2 * wood_fac_10hr
      burned_acid_ag2_100hr        = c_acid_ag2       * wood_fac_100hr
      burned_water_ag2_100hr       = c_water_ag2      * wood_fac_100hr
      burned_ethanol_ag2_100hr     = c_ethanol_ag2    * wood_fac_100hr
      burned_nonsoluble_ag2_100hr  = c_nonsoluble_ag2 * wood_fac_100hr
      burned_acid_ag2_1000hr       = c_acid_ag2       * wood_fac_1000hr
      burned_water_ag2_1000hr      = c_water_ag2      * wood_fac_1000hr
      burned_ethanol_ag2_1000hr    = c_ethanol_ag2    * wood_fac_1000hr
      burned_nonsoluble_ag2_1000hr = c_nonsoluble_ag2 * wood_fac_1000hr
    ELSE
      ! combustion of live grass
      burned_c_green               = fract_green_aboveGround * c_green * CombustCompleteness_lg  * burned_frac
      burned_acid_ag1_1hr          = c_acid_ag1                        * CombustCompleteness_1hr * burned_frac
      burned_water_ag1_1hr         = c_water_ag1                       * CombustCompleteness_1hr * burned_frac
      burned_ethanol_ag1_1hr       = c_ethanol_ag1                     * CombustCompleteness_1hr * burned_frac
      burned_nonsoluble_ag1_1hr    = c_nonsoluble_ag1                  * CombustCompleteness_1hr * burned_frac
   ENDIF

    ! sum fuels for litter according to yasso classes
    ! ... burned green
    burned_acid_ag1       = burned_acid_ag1_1hr       + burned_acid_leaves
    burned_water_ag1      = burned_water_ag1_1hr      + burned_water_leaves
    burned_ethanol_ag1    = burned_ethanol_ag1_1hr    + burned_ethanol_leaves
    burned_nonsoluble_ag1 = burned_nonsoluble_ag1_1hr + burned_nonsoluble_leaves

    burned_litter_green   = burned_acid_ag1    + burned_water_ag1 + burned_ethanol_ag1 + burned_nonsoluble_ag1

    ! ... burned woody
    burned_acid_ag2       = burned_acid_ag2_1hr         + burned_acid_ag2_10hr         &
                          + burned_acid_ag2_100hr       + burned_acid_ag2_1000hr
    burned_water_ag2      = burned_water_ag2_1hr        + burned_water_ag2_10hr        &
                          + burned_water_ag2_100hr      + burned_water_ag2_1000hr
    burned_ethanol_ag2    = burned_ethanol_ag2_1hr      + burned_ethanol_ag2_10hr      &
                          + burned_ethanol_ag2_100hr    + burned_ethanol_ag2_1000hr
    burned_nonsoluble_ag2 = burned_nonsoluble_ag2_1hr   + burned_nonsoluble_ag2_10hr   &
                          + burned_nonsoluble_ag2_100hr + burned_nonsoluble_ag2_1000hr

    burned_litter_wood    = burned_acid_ag2  + burned_water_ag2  + burned_ethanol_ag2 + burned_nonsoluble_ag2

    ! ... and total
!    burned_leaves      = burned_acid_leaves + burned_water_leaves + burned_ethanol_leaves + burned_nonsoluble_leaves
!    litter_green       = c_acid_ag1 + c_water_ag1 + c_ethanol_ag1  + c_nonsoluble_ag1
!    burned_green       = burned_acid_ag1 + burned_water_ag1 + burned_ethanol_ag1  + burned_nonsoluble_ag1
    litter_wood        = c_acid_ag2             + c_water_ag2             + c_ethanol_ag2             + c_nonsoluble_ag2
    burned_wood        = burned_acid_ag2        + burned_water_ag2        + burned_ethanol_ag2        + burned_nonsoluble_ag2
    burned_wood_1hr    = burned_acid_ag2_1hr    + burned_water_ag2_1hr    + burned_ethanol_ag2_1hr    + burned_nonsoluble_ag2_1hr
    burned_wood_10hr   = burned_acid_ag2_10hr   + burned_water_ag2_10hr   + burned_ethanol_ag2_10hr   + burned_nonsoluble_ag2_10hr
    burned_wood_100hr  = burned_acid_ag2_100hr  + burned_water_ag2_100hr  + burned_ethanol_ag2_100hr  + burned_nonsoluble_ag2_100hr
    burned_wood_1000hr = burned_acid_ag2_1000hr + burned_water_ag2_1000hr + burned_ethanol_ag2_1000hr + burned_nonsoluble_ag2_1000hr

    IF (woody_pft .AND. ((litter_wood - burned_wood) > EPSILON(1._wp))) THEN
      fract_1hr_wood    = (litter_wood * fract_1hr_wood    - burned_wood_1hr   ) / (litter_wood - burned_wood)
      fract_10hr_wood   = (litter_wood * fract_10hr_wood   - burned_wood_10hr  ) / (litter_wood - burned_wood)
      fract_100hr_wood  = (litter_wood * fract_100hr_wood  - burned_wood_100hr ) / (litter_wood - burned_wood)
      fract_1000hr_wood = (litter_wood * fract_1000hr_wood - burned_wood_1000hr) / (litter_wood - burned_wood)
    ELSE
      fract_1hr_wood    = fract_1hr_wood_new
      fract_10hr_wood   = fract_10hr_wood_new
      fract_100hr_wood  = fract_100hr_wood_new
      fract_1000hr_wood = fract_1000hr_wood_new
    ENDIF

  END SUBROUTINE calc_fire_spitfire_pool_consumption

  SUBROUTINE fire_spitfire_relocate_carbon(nc, woody_pft, natural_pft,              &
    !with_nitrogen, &
                                           LeafLit_coef,             WoodLit_coef,  &
                                           veg_fract_correction, cf_burned,         &
                                           postfire_mortality, mortality_crownkill, &
                                           c_green, c_woods, c_reserve,             &
                                           c_acid_ag1,                              &
                                           c_water_ag1,                             &
                                           c_ethanol_ag1,                           &
                                           c_nonsoluble_ag1,                        &
                                           c_acid_ag2,                              &
                                           c_water_ag2,                             &
                                           c_ethanol_ag2,                           &
                                           c_nonsoluble_ag2,                        &
                                           c_acid_bg1,                              &
                                           c_water_bg1,                             &
                                           c_ethanol_bg1,                           &
                                           c_nonsoluble_bg1,                        &
                                           c_humus_1,                               &
                                           c_acid_bg2,                              &
                                           c_water_bg2,                             &
                                           c_ethanol_bg2,                           &
                                           c_nonsoluble_bg2,                        &
                                           c_humus_2,                               &
                                           fract_1hr_wood, fract_10hr_wood, fract_100hr_wood, fract_1000hr_wood, &
                                           burned_c_green, burned_c_reserve, burned_c_wood,                      &
                                           burned_c_litter_green_ag, burned_c_litter_wood_ag,                    &
                                           carbon_2_GreenLitterPools,carbon_2_WoodLitterPools,                   &
                                           carbon_2_atmos                                                        &
!                                           !Nitrogen Pools and fluxes
!                                           Npool_green, Npool_woods, Npool_mobile, SMINN_pool,      &
!                                           Npool_litter_green_ag,Npool_litter_green_bg,             &
!                                           Npool_litter_wood_ag,Npool_litter_wood_bg,               &
!                                           nitrogen_2_GreenLitterPools, nitrogen_2_WoodLitterPools, &
!                                           nitrogen_2_atmos, nitrogen_2_sminn &
                                           )
!  USE mo_cbal_parameters,  ONLY: cn_woods, cn_litter_wood
  USE mo_exception,             ONLY: message_text
  USE mo_carbon_constants,      ONLY: fract_wood_aboveGround, fract_green_aboveGround
  USE mo_disturbance_constants, ONLY: fract_1hr_wood_new, fract_10hr_wood_new, fract_100hr_wood_new, fract_1000hr_wood_new

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fire_spitfire_relocate_carbon'
  INTEGER,                  INTENT(in)    :: nc                 ! Vector length
  LOGICAL,                  INTENT(in)    :: woody_pft
  LOGICAL,                  INTENT(in)    :: natural_pft
  ! YASSO parameter
  REAL(wp),                 INTENT(in)    :: LeafLit_coef(:)
  REAL(wp),                 INTENT(in)    :: WoodLit_coef(:)
!  LOGICAL,                  INTENT(in)    :: with_nitrogen
  REAL(wp),                 INTENT(in)    :: veg_fract_correction(:)  ! Sparseness of vegetation)
  REAL(wp),                 INTENT(in)    :: cf_burned(:)     ! Frac of the veg area burned since the last call of this routine
  REAL(wp),                 INTENT(in)    :: postfire_mortality(:)
  REAL(wp),                 INTENT(in)    :: mortality_crownkill(:)
  ! vegetation C-pools
  REAL(wp),                 INTENT(inout) :: c_green(:)           ! green carbon pool [mol(C)/m2(canopy)]
  REAL(wp),                 INTENT(inout) :: c_woods(:)           ! wood carbon pool [mol(C)/m2(canopy)]
  REAL(wp),                 INTENT(inout) :: c_reserve(:)         ! reserve carbon pool [mol(C)/m2(canopy)]
  ! YASSO pools
  REAL(wp),                 INTENT(inout) :: c_acid_ag1      (:)
  REAL(wp),                 INTENT(inout) :: c_water_ag1     (:)
  REAL(wp),                 INTENT(inout) :: c_ethanol_ag1   (:)
  REAL(wp),                 INTENT(inout) :: c_nonsoluble_ag1(:)
  REAL(wp),                 INTENT(inout) :: c_acid_ag2      (:)
  REAL(wp),                 INTENT(inout) :: c_water_ag2     (:)
  REAL(wp),                 INTENT(inout) :: c_ethanol_ag2   (:)
  REAL(wp),                 INTENT(inout) :: c_nonsoluble_ag2(:)
  REAL(wp),                 INTENT(inout) :: c_acid_bg1      (:)
  REAL(wp),                 INTENT(inout) :: c_water_bg1     (:)
  REAL(wp),                 INTENT(inout) :: c_ethanol_bg1   (:)
  REAL(wp),                 INTENT(inout) :: c_nonsoluble_bg1(:)
  REAL(wp),                 INTENT(inout) :: c_humus_1       (:)
  REAL(wp),                 INTENT(inout) :: c_acid_bg2      (:)
  REAL(wp),                 INTENT(inout) :: c_water_bg2     (:)
  REAL(wp),                 INTENT(inout) :: c_ethanol_bg2   (:)
  REAL(wp),                 INTENT(inout) :: c_nonsoluble_bg2(:)
  REAL(wp),                 INTENT(inout) :: c_humus_2       (:)
  REAL(wp),                 INTENT(inout) :: burned_c_green(:), burned_c_reserve(:), burned_c_wood(:)
  REAL(wp),                 INTENT(inout) :: burned_c_litter_green_ag(:), burned_c_litter_wood_ag(:)
  REAL(wp),                 INTENT(inout) :: fract_1hr_wood(:), fract_10hr_wood(:), fract_100hr_wood(:), fract_1000hr_wood(:)
  ! C-fluxes
  REAL(wp),                 INTENT(inout) :: carbon_2_GreenLitterPools(:)  ! carbon to green litter pools[mol(C)/m2(vegetated area)]
  REAL(wp),                 INTENT(inout) :: carbon_2_WoodLitterPools(:)   ! carbon to wood  litter pools[mol(C)/m2(vegetated area)]
  REAL(wp),                 INTENT(out)   :: carbon_2_atmos(:)             ! C to atmosphere [mol(C)/m2(vegetated area)/timestep]
  ! Nitrogen
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_green(:,:)              !! Green nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_woods(:,:)              !! Wood nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_mobile(:,:)             !! Mobile nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: SMINN_pool(:,:)               !! Mineral nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_litter_green_ag(:,:)    !! Above ground green litter N pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_litter_green_bg(:,:)    !! Below ground green litter N pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_litter_wood_ag(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: Npool_litter_wood_bg(:,:)     !! Wood litter nitrogen pool [mol(N)/m^2(canopy)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: nitrogen_2_GreenLitterPools(:)!! Amount of N relocated by fire to the green litter
!                                                                           !!    pools [mol(N)/m^2(vegetated area)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: nitrogen_2_WoodLitterPools(:) !! Amount of nitrogen relocated by fire and wind damage
!                                                                           !!    to the wood litter pools [mol(N)/m^2(veg area)]
!  REAL(dp),OPTIONAL,        INTENT(out)   :: nitrogen_2_atmos(:)           !! Amount of N immediately released by fire and wind
!                                                                           !!    damage to the atmosphere [mol(N)/m^2(veg area)]
!  REAL(dp),OPTIONAL,        INTENT(inout) :: nitrogen_2_sminn(:)           !! Amount of N relocated by fire and wind damage to the
                                                                           !!    soil mineral N pool [mol(N)/m^2(vegetated area)]

  ! local variables
  REAL(wp) :: c_total_bf(nc)
  REAL(wp) :: c_total_af(nc)
!  REAL(wp) :: npools_total_bf(nc)
!  REAL(wp) :: npools_total_af(nc)
  REAL(wp) :: carbon_2_GreenLitterPools_tiled(nc)
  REAL(wp) :: carbon_2_WoodLitterPools_tiled(nc)
  REAL(wp) :: green_2_litter_green_ag(nc),green_2_litter_green_bg(nc)
  REAL(wp) :: reserve_2_litter_green_ag(nc),reserve_2_litter_green_bg(nc)
  REAL(wp) :: wood_2_litter_wood_ag(nc),wood_2_litter_wood_bg(nc)
  REAL(wp) :: wood_ck_2_litter_wood_ag, wood_camb_2_litter_wood_ag
!  REAL(wp) :: green_ck_2_litter_green_ag,
  REAL(wp) :: green_camb_2_litter_green_ag
!  REAL(wp) :: leaves_2_green_litter_pool,
  REAL(wp) :: f1hr_2_wood_litter_pool, f10hr_2_wood_litter_pool
  REAL(wp) :: f100hr_2_wood_litter_pool, f1000hr_2_wood_litter_pool
!  REAL(wp) :: leaves_old_pool
  REAL(wp) :: mortality_nocrown, cc_ck_wood, cc_ck_green, cc_green
  REAL(wp) :: f1hr_old_pool, f10hr_old_pool, f100hr_old_pool, f1000hr_old_pool
  REAL(wp) :: c_litter_wood
  REAL(wp) :: sum_fractions(nc), rem_frac(nc)
  REAL(wp) :: scale_fac(nc), sum_scale_fac
  INTEGER  :: i
!  ! Temporary nitrogen variables. All /m^2(canopy)
!  REAL(wp) :: N_2_litter_green_bg(nidx), N_2_litter_green_ag(nidx)
!  REAL(wp) :: N_2_litter_wood_bg(nidx),  N_2_litter_wood_ag (nidx)
!  REAL(wp) :: N_release(nidx), N_2_SMINN(nidx)
!  REAL(wp) :: N_litter_green_burned(nidx), N_litter_wood_burned(nidx)
!  REAL(wp) :: N_green_burned(nidx), N_mobile_burned(nidx), N_wood_burned(nidx)
!  REAL(wp) :: C_litter_green_ag(nidx), C_litter_wood_ag(nidx)

    ! initialisation
    carbon_2_GreenLitterPools(:) = 0._wp
    carbon_2_WoodLitterPools (:) = 0._wp
    green_2_litter_green_ag  (:) = 0._wp
    reserve_2_litter_green_ag(:) = 0._wp
    wood_2_litter_wood_ag    (:) = 0._wp
!    IF (with_nitrogen) THEN
!      nitrogen_2_GreenLitterPools(:)   = 0._dp
!      nitrogen_2_WoodLitterPools (:)   = 0._dp
!      nitrogen_2_atmos           (:)   = 0._dp
!      nitrogen_2_sminn           (:)   = 0._dp
!      N_2_litter_green_ag        (:) = 0._dp
!      N_2_litter_wood_ag         (:) = 0._dp
!      N_litter_green_burned      (:) = 0._dp
!      N_litter_wood_burned       (:) = 0._dp
!      N_green_burned             (:) = 0._dp
!      N_mobile_burned            (:) = 0._dp
!      N_wood_burned              (:) = 0._dp
!      IF (with_yasso) THEN
!        C_litter_wood_ag (:,:) = YCpool_acid_ag2   (:,:) + YCpool_water_ag2     (:,:) &
!                               + YCpool_ethanol_ag2(:,:) + YCpool_nonsoluble_ag2(:,:)
!        C_litter_green_ag(:,:) = YCpool_acid_ag1   (:,:) + YCpool_water_ag1     (:,:) &
!                               + YCpool_ethanol_ag1(:,:) + YCpool_nonsoluble_ag1(:,:)
!      ELSE
!        C_litter_wood_ag (:,:) = Cpool_litter_wood_ag (:,:)
!        C_litter_green_ag(:,:) = Cpool_litter_green_ag(:,:)
!      ENDIF
!    ENDIF

    ! Total carbon and nitrogen which must be preserved
    c_total_bf(:) = c_green(:) + c_woods(:) + c_reserve(:)    &
                  + c_acid_ag1      (:) + c_acid_ag2      (:) &
                  + c_water_ag1     (:) + c_water_ag2     (:) &
                  + c_ethanol_ag1   (:) + c_ethanol_ag2   (:) &
                  + c_nonsoluble_ag1(:) + c_nonsoluble_ag2(:) &
                  + c_acid_bg1      (:) + c_acid_bg2      (:) &
                  + c_water_bg1     (:) + c_water_bg2     (:) &
                  + c_ethanol_bg1   (:) + c_ethanol_bg2   (:) &
                  + c_nonsoluble_bg1(:) + c_nonsoluble_bg2(:) &
                  + c_humus_1       (:) + c_humus_2       (:)
!    IF (with_nitrogen) THEN
!      npools_total_bf(:,:) = Npool_green(:,:) + Npool_woods(:,:) + Npool_mobile (:,:) &
!                           + Npool_litter_green_ag(:,:) +   Npool_litter_wood_ag(:,:) &
!                           + Npool_litter_green_bg(:,:) +   Npool_litter_wood_bg(:,:) &
!                           + SMINN_pool(:,:)
!    ENDIF

    ! below ground carbon, same for woody and grass types, postfire_mortality=1 for grasses (initialisation)
    ! also the update of below ground litter pools is the same: litter below ground is not burned
    green_2_litter_green_bg  (:) = (1._wp - fract_green_aboveGround) * cf_burned(:) * postfire_mortality * c_green  (:)
    reserve_2_litter_green_bg(:) = (1._wp - fract_green_aboveGround) * cf_burned(:) * postfire_mortality * c_reserve(:)
    wood_2_litter_wood_bg    (:) = (1._wp - fract_wood_aboveGround ) * cf_burned(:) * postfire_mortality * c_woods  (:)

!    IF (with_nitrogen) THEN
!      ! Handling of below ground ex. living N now belonging to burned vegetation
!      N_2_litter_green_bg(:,:) = (1._dp-frac_green_aboveGround)*cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
!                               * (Npool_green(:,:) + Npool_mobile(:,:))
!      N_2_litter_wood_bg (:,:) = (1._dp-frac_wood_aboveGround) *cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
!                               *  Cpool_woods(:,:)/cn_litter_wood
!      ! Excess wood N go into SMINN. More terms are added later
!      N_2_SMINN          (:,:) = (1._dp-frac_wood_aboveGround) *cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:) &
!                               *  Cpool_woods(:,:)*(1._dp/cn_woods-1._dp/cn_litter_wood)
!    ENDIF
!
    ! combustion coefficient for living biomass of killed trees (see TH2010 appendix B2)
    ! cc for wood pool when crown kill occurs: above ground, 100% of twigs and small branches (1hr), 5% big branches (10hr)
    cc_ck_wood  = fract_wood_aboveGround  * (fract_1hr_wood_new + 0.05_wp * fract_10hr_wood_new)
    ! combustion coeffiecient for green pool (leaves and sap wood):
    ! all leaves, 1hr, and 5% of 10hr fuels above ground are burned
!    cc_ck_green = 1.0_dp
    cc_ck_green = fract_green_aboveGround

    DO i=1,nc
      IF (cf_burned(i) > 0.5_wp * EPSILON(1._wp)) THEN

        mortality_nocrown = postfire_mortality(i) - mortality_crownkill(i)

        IF (woody_pft .AND. natural_pft) THEN

          ! add biomass burned of woody species:crown_kill: all leaves, all 1, and 5% 10 hour woods
          burned_c_green  (i) = mortality_crownkill(i) * cf_burned(i) * c_green  (i) * cc_ck_green
          burned_c_wood   (i) = mortality_crownkill(i) * cf_burned(i) * c_woods  (i) * cc_ck_wood
          burned_c_reserve(i) = mortality_crownkill(i) * cf_burned(i) * c_reserve(i) * fract_green_aboveGround

          ! carbon affected by fire and transferred into the wood litter pool - diagnostic needed for C conversion test in JSBACH
          wood_ck_2_litter_wood_ag       = mortality_crownkill(i) * c_woods(i) * fract_wood_aboveGround * cf_burned(i) &
                                         * (0.95_wp * fract_10hr_wood_new + fract_100hr_wood_new + fract_1000hr_wood_new)
          wood_camb_2_litter_wood_ag     = mortality_nocrown * c_woods(i) * cf_burned(i) * fract_wood_aboveGround
          wood_2_litter_wood_ag(i)       = wood_ck_2_litter_wood_ag + wood_camb_2_litter_wood_ag

          ! carbon transferred into the green litter pool above ground
          green_camb_2_litter_green_ag = mortality_nocrown * c_green(i)   * fract_green_aboveGround * cf_burned(i)
          reserve_2_litter_green_ag(i) = mortality_nocrown * c_reserve(i) * fract_green_aboveGround * cf_burned(i)
          green_2_litter_green_ag  (i) = green_camb_2_litter_green_ag ! crownkill kills all leaves

          ! available wood litter
          c_litter_wood = c_acid_ag2(i) + c_water_ag2(i) + c_ethanol_ag2(i) + c_nonsoluble_ag2(i)

          ! new 1hr fuel to wood litter pool
          f1hr_2_wood_litter_pool     = fract_1hr_wood_new * wood_camb_2_litter_wood_ag
          f1hr_old_pool               = fract_1hr_wood(i) * c_litter_wood

          ! new 10hr fuel to wood litter pool
          f10hr_2_wood_litter_pool    = fract_10hr_wood_new * wood_camb_2_litter_wood_ag                         &
                                      + mortality_crownkill(i) * c_woods(i)                                      &
                                        * fract_wood_aboveGround  * cf_burned(i) * fract_10hr_wood_new * 0.95_wp
          f10hr_old_pool              = fract_10hr_wood(i) * c_litter_wood

          ! new 100hr fuel to wood litter pool
          f100hr_2_wood_litter_pool   = fract_100hr_wood_new * wood_camb_2_litter_wood_ag              &
                                      + mortality_crownkill(i) * c_woods(i)                            &
                                        * fract_wood_aboveGround * cf_burned(i) * fract_100hr_wood_new
          f100hr_old_pool             = fract_100hr_wood(i) * c_litter_wood

          ! new 1000hr fuel to wood litter pool
          f1000hr_2_wood_litter_pool  = fract_1000hr_wood_new * wood_camb_2_litter_wood_ag              &
                                      + mortality_crownkill(i) * c_woods(i)                             &
                                        * fract_wood_aboveGround * cf_burned(i) * fract_1000hr_wood_new
          f1000hr_old_pool            = fract_1000hr_wood(i) * c_litter_wood

          ! Update fuel fractions
          IF ((wood_2_litter_wood_ag(i) + c_litter_wood) > 0._wp) THEN
            fract_1hr_wood   (i) = (f1hr_2_wood_litter_pool    + f1hr_old_pool)    / (wood_2_litter_wood_ag(i) + c_litter_wood)
            fract_10hr_wood  (i) = (f10hr_2_wood_litter_pool   + f10hr_old_pool)   / (wood_2_litter_wood_ag(i) + c_litter_wood)
            fract_100hr_wood (i) = (f100hr_2_wood_litter_pool  + f100hr_old_pool)  / (wood_2_litter_wood_ag(i) + c_litter_wood)
            fract_1000hr_wood(i) = (f1000hr_2_wood_litter_pool + f1000hr_old_pool) / (wood_2_litter_wood_ag(i) + c_litter_wood)
          ELSE
            fract_1hr_wood   (i) = fract_1hr_wood_new
            fract_10hr_wood  (i) = fract_10hr_wood_new
            fract_100hr_wood (i) = fract_100hr_wood_new
            fract_1000hr_wood(i) = fract_1000hr_wood_new
          ENDIF

        END IF ! Natural woody pft

        ! combustion completeness of grasses, burned_c_green is computed in fuel_consumtion
        IF (.NOT. woody_pft) THEN
          green_2_litter_green_ag(i) = c_green(i) * fract_green_aboveGround * cf_burned(i) - burned_c_green(i)
          IF (burned_c_green(i) > EPSILON(1._wp)) THEN
            cc_green = burned_c_green(i) /  (c_green(i) * fract_green_aboveGround * cf_burned(i))
            reserve_2_litter_green_ag(i) = c_reserve(i) * fract_green_aboveGround * cf_burned(i) * (1._wp-cc_green)
            burned_c_reserve         (i) = c_reserve(i) * fract_green_aboveGround * cf_burned(i) *        cc_green
          ELSE
            reserve_2_litter_green_ag(i) = c_reserve(i) * fract_green_aboveGround * cf_burned(i)
            burned_c_reserve         (i) = 0._wp
          ENDIF
        ENDIF ! not woody and not crop
      ENDIF ! PFT burning
    ENDDO

    carbon_2_GreenLitterPools_tiled(:) = green_2_litter_green_ag  (:) + green_2_litter_green_bg  (:)  &
                                       + reserve_2_litter_green_ag(:) + reserve_2_litter_green_bg(:)
    carbon_2_WoodLitterPools_tiled (:) = wood_2_litter_wood_ag    (:) + wood_2_litter_wood_bg    (:)

    ! litter consumed and live biomass burned is transferred to atmosphere
    carbon_2_atmos(:) = burned_c_litter_wood_ag (:) &
                      + burned_c_litter_green_ag(:) &
                      + burned_c_green          (:) &
                      + burned_c_wood           (:) &
                      + burned_c_reserve        (:)

!    IF (with_nitrogen) THEN
!      ! Nitrogen to be released from burned material
!      N_litter_wood_burned (:,:) = state_FTH%burned_Cpool_litter_wood_ag (kidx0:kidx1,:)/MAX(1.e-30_dp,C_litter_wood_ag (:,:)) &
!                                                                                        *              Npool_litter_wood_ag (:,:)
!      N_litter_green_burned(:,:) = state_FTH%burned_Cpool_litter_green_ag(kidx0:kidx1,:)/MAX(1.e-30_dp,C_litter_green_ag(:,:)) &
!                                                                                        *              Npool_litter_green_ag(:,:)
!      N_green_burned       (:,:) = state_FTH%burned_Cpool_green(kidx0:kidx1,:) / MAX(1.e-30_dp,Cpool_green(:,:)) * Npool_green (:,:)
!      N_mobile_burned      (:,:) = state_FTH%burned_Cpool_green(kidx0:kidx1,:) / MAX(1.e-30_dp,Cpool_green(:,:)) * Npool_mobile(:,:)
!      N_wood_burned        (:,:) = state_FTH%burned_Cpool_wood (kidx0:kidx1,:) / MAX(1.e-30_dp,Cpool_woods(:,:)) * Npool_woods (:,:)
!      N_release(:,:) = N_litter_wood_burned(:,:) + N_litter_green_burned(:,:) &
!                     + N_green_burned      (:,:) + N_mobile_burned      (:,:) + N_wood_burned(:,:)
!      ! Nitrogen from living pools to above ground litter pools
!      N_2_litter_green_ag(:,:) = (Npool_green(:,:) + Npool_mobile(:,:)) * frac_green_aboveGround * cf_burned(:,:) &
!                               * state_FTH%postfire_mortality(kidx0:kidx1,:) - N_green_burned(:,:) - N_mobile_burned(:,:)
!      N_2_litter_wood_ag (:,:) = (Cpool_woods(:,:) * cf_burned(:,:) * state_FTH%postfire_mortality(kidx0:kidx1,:) &
!                                  * frac_wood_aboveGround - state_FTH%burned_Cpool_wood(kidx0:kidx1,:)) / cn_litter_wood
!      N_2_SMINN(:,:) = N_2_SMINN(:,:) + Npool_woods(:,:) * frac_wood_aboveGround * cf_burned(:,:) &
!                     * state_FTH%postfire_mortality(kidx0:kidx1,:) - N_wood_burned (:,:) - N_2_litter_wood_ag(:,:)
!    ENDIF

    ! Update pools
    c_reserve(:) = c_reserve(:) - (reserve_2_litter_green_bg(:) + reserve_2_litter_green_ag(:) + burned_c_reserve(:))
    c_woods  (:) = c_woods  (:) - (carbon_2_WoodLitterPools (:) + burned_c_wood(:))
    c_green  (:) = c_green  (:) - (green_2_litter_green_ag  (:) + green_2_litter_green_bg  (:) + burned_c_green(:))
    ! Aboveground green pools
    rem_frac(:)=(1._wp-burned_c_litter_green_ag(:)/MAX(1.e-30_wp,c_acid_ag1(:)+c_water_ag1(:)+c_ethanol_ag1(:)+c_nonsoluble_ag1(:)))
    c_acid_ag1      (:) = c_acid_ag1      (:)*rem_frac(:)+LeafLit_coef(1)*(green_2_litter_green_ag(:)+reserve_2_litter_green_ag(:))
    c_water_ag1     (:) = c_water_ag1     (:)*rem_frac(:)+LeafLit_coef(2)*(green_2_litter_green_ag(:)+reserve_2_litter_green_ag(:))
    c_ethanol_ag1   (:) = c_ethanol_ag1   (:)*rem_frac(:)+LeafLit_coef(3)*(green_2_litter_green_ag(:)+reserve_2_litter_green_ag(:))
    c_nonsoluble_ag1(:) = c_nonsoluble_ag1(:)*rem_frac(:)+LeafLit_coef(4)*(green_2_litter_green_ag(:)+reserve_2_litter_green_ag(:))
    ! Aboveground wood pools
    rem_frac(:)=(1._wp-burned_c_litter_wood_ag(:)/MAX(1.e-30_wp,c_acid_ag2(:)+c_water_ag2(:)+c_ethanol_ag2(:)+c_nonsoluble_ag2(:)))
    c_acid_ag2      (:) = c_acid_ag2      (:) * rem_frac(:) + WoodLit_coef(1) * wood_2_litter_wood_ag(:)
    c_water_ag2     (:) = c_water_ag2     (:) * rem_frac(:) + WoodLit_coef(2) * wood_2_litter_wood_ag(:)
    c_ethanol_ag2   (:) = c_ethanol_ag2   (:) * rem_frac(:) + WoodLit_coef(3) * wood_2_litter_wood_ag(:)
    c_nonsoluble_ag2(:) = c_nonsoluble_ag2(:) * rem_frac(:) + WoodLit_coef(4) * wood_2_litter_wood_ag(:)
    ! Belowground green pools
    c_acid_bg1      (:) = c_acid_bg1      (:) + LeafLit_coef(1) * (green_2_litter_green_bg(:) + reserve_2_litter_green_bg(:))
    c_water_bg1     (:) = c_water_bg1     (:) + LeafLit_coef(2) * (green_2_litter_green_bg(:) + reserve_2_litter_green_bg(:))
    c_ethanol_bg1   (:) = c_ethanol_bg1   (:) + LeafLit_coef(3) * (green_2_litter_green_bg(:) + reserve_2_litter_green_bg(:))
    c_nonsoluble_bg1(:) = c_nonsoluble_bg1(:) + LeafLit_coef(4) * (green_2_litter_green_bg(:) + reserve_2_litter_green_bg(:))
    c_humus_1       (:) = c_humus_1       (:) + LeafLit_coef(5) * (green_2_litter_green_bg(:) + reserve_2_litter_green_bg(:))
    ! Belowground wood pools
    c_acid_bg2      (:) = c_acid_bg2      (:) + WoodLit_coef(1) * wood_2_litter_wood_bg(:)
    c_water_bg2     (:) = c_water_bg2     (:) + WoodLit_coef(2) * wood_2_litter_wood_bg(:)
    c_ethanol_bg2   (:) = c_ethanol_bg2   (:) + WoodLit_coef(3) * wood_2_litter_wood_bg(:)
    c_nonsoluble_bg2(:) = c_nonsoluble_bg2(:) + WoodLit_coef(4) * wood_2_litter_wood_bg(:)
    c_humus_2       (:) = c_humus_2       (:) + WoodLit_coef(5) * wood_2_litter_wood_bg(:)

!    IF (with_nitrogen) THEN
!      Npool_green          (:,:) = Npool_green          (:,:) * (1._dp - cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:))
!      Npool_mobile         (:,:) = Npool_mobile         (:,:) * (1._dp - cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:))
!      Npool_woods          (:,:) = Npool_woods          (:,:) * (1._dp - cf_burned(:,:)*state_FTH%postfire_mortality(kidx0:kidx1,:))
!      SMINN_pool           (:,:) = SMINN_pool           (:,:) + N_2_SMINN(:,:) + (1._dp-fire_TH_options%frac_N_2_atm)*N_release(:,:)
!      Npool_litter_green_ag(:,:) = Npool_litter_green_ag(:,:) + N_2_litter_green_ag(:,:) - N_litter_green_burned(:,:)
!      Npool_litter_green_bg(:,:) = Npool_litter_green_bg(:,:) + N_2_litter_green_bg(:,:)
!      Npool_litter_wood_ag (:,:) = Npool_litter_wood_ag (:,:) + N_2_litter_wood_ag (:,:) - N_litter_wood_burned (:,:)
!      Npool_litter_wood_bg (:,:) = Npool_litter_wood_bg (:,:) + N_2_litter_wood_bg (:,:)
!    ENDIF

    ! bring to total tile area
    carbon_2_atmos           (:) = carbon_2_atmos           (:) * veg_fract_correction(:)
    carbon_2_GreenLitterPools(:) = carbon_2_GreenLitterPools(:) * veg_fract_correction(:)
    carbon_2_WoodLitterPools (:) = carbon_2_WoodLitterPools (:) * veg_fract_correction(:)

!    IF (with_nitrogen) THEN
!      nitrogen_2_atmos(:) = SUM( N_release(:,:) *        fire_TH_options%frac_N_2_atm &
!                                * veg_fract_correction(:,:) * cover_fract(:,:),DIM=2)
!      nitrogen_2_sminn(:) = SUM((N_release(:,:) * (1._dp-fire_TH_options%frac_N_2_atm) + N_2_SMINN(:,:)) &
!                                * veg_fract_correction(:,:) * cover_fract(:,:),DIM=2)
!      nitrogen_2_GreenLitterPools(:) = SUM((N_2_litter_green_ag(:,:) + N_2_litter_green_bg(:,:)) &
!                                           * veg_fract_correction(:,:) * cover_fract(:,:),DIM=2)
!      nitrogen_2_WoodLitterPools (:) = SUM((N_2_litter_wood_ag (:,:) + N_2_litter_wood_bg (:,:)) &
!                                           * veg_fract_correction(:,:) * cover_fract(:,:),DIM=2)
!    ENDIF

    ! normalize fractions
    sum_fractions    (:) = fract_1hr_wood   (:) + fract_10hr_wood(:) + fract_100hr_wood(:) + fract_1000hr_wood(:)
    fract_1hr_wood   (:) = fract_1hr_wood   (:) / sum_fractions  (:)
    fract_10hr_wood  (:) = fract_10hr_wood  (:) / sum_fractions  (:)
    fract_100hr_wood (:) = fract_100hr_wood (:) / sum_fractions  (:)
    fract_1000hr_wood(:) = fract_1000hr_wood(:) / sum_fractions  (:)

    ! Diagnose sum of carbon pools after relocate (without slow pool) for C mass balance check
    c_total_af(:) = c_green(:) + c_woods(:) + c_reserve(:)    &
                  + c_acid_ag1      (:) + c_acid_ag2      (:) &
                  + c_water_ag1     (:) + c_water_ag2     (:) &
                  + c_ethanol_ag1   (:) + c_ethanol_ag2   (:) &
                  + c_nonsoluble_ag1(:) + c_nonsoluble_ag2(:) &
                  + c_acid_bg1      (:) + c_acid_bg2      (:) &
                  + c_water_bg1     (:) + c_water_bg2     (:) &
                  + c_ethanol_bg1   (:) + c_ethanol_bg2   (:) &
                  + c_nonsoluble_bg1(:) + c_nonsoluble_bg2(:) &
                  + c_humus_1       (:) + c_humus_2       (:)
    c_total_af(:) = c_total_af(:) * veg_fract_correction(:)
    c_total_bf(:) = c_total_bf(:) * veg_fract_correction(:)

    WHERE(c_total_bf(:) > 1.0_wp)
      scale_fac(:) = (carbon_2_atmos(:) + c_total_af(:) - c_total_bf(:)) / c_total_bf(:)
    ELSEWHERE
      scale_fac(:) = 0._wp
    ENDWHERE
    sum_scale_fac = SUM(scale_fac)
    IF (ABS(sum_scale_fac) > 1.e-13_wp) THEN
      WRITE(message_text,*) 'Carbon conservation in relocation rather low (',sum_scale_fac,')'
      CALL message(TRIM(routine),TRIM(message_text))
    END IF

!    IF (with_nitrogen) THEN
!      npools_total_af(:,:) = Npool_green(:,:) + Npool_woods(:,:) + Npool_mobile (:,:) &
!                           + Npool_litter_green_ag(:,:) +   Npool_litter_wood_ag(:,:) &
!                           + Npool_litter_green_bg(:,:) +   Npool_litter_wood_bg(:,:) &
!                           + SMINN_pool(:,:)
!      npools_total_af(:,:) = npools_total_af(:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
!      npools_total_bf(:,:) = npools_total_bf(:,:) * veg_fract_correction(:,:) * cover_fract(:,:)
!      scale_fac(:,:) = 0._dp
!      DO itile=1,ntiles
!        WHERE(cover_fract(:,itile) > 0.5_dp*fract_small)
!          scale_fac(:,1) = scale_fac(:,1) + npools_total_af(:,itile) - npools_total_bf(:,itile)
!          scale_fac(:,2) = scale_fac(:,2) + npools_total_bf(:,itile)
!        ENDWHERE
!      ENDDO
!      WHERE(scale_fac(:,2) > 0._dp)
!        scale_fac(:,3) = (scale_fac(:,1) + nitrogen_2_atmos(:)) / scale_fac(:,2)
!      ENDWHERE
!      sum_scale_fac = SUM(scale_fac(:,3))
!      IF(abs(sum_scale_fac) > 1.e-13_dp) THEN
!        WRITE(message_text,*) 'Nitrogen conservation in relocation rather low (',sum_scale_fac,')'
!        CALL message('relocate_carbon_fire_thonicke',TRIM(message_text))
!      END IF
!
!    ENDIF

  END SUBROUTINE fire_spitfire_relocate_carbon

#endif
END MODULE mo_disturb_process
