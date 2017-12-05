!>
!! @brief Contains the routines for the disturb processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!!
MODULE mo_disturb_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_lctlib_class, ONLY: t_lctlib_element
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_parallel,     ONLY: Get_omp_thread
  dsl4jsb_Use_processes  DISTURB_
  dsl4jsb_Use_config(DISTURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: disturbed_frac, get_relative_humidity_air

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  ! Driver for calculation of damaged surface area (presently by fire and wind breaks)
  !
  ! History:
  ! 2011-04-29 - StW - Original code
  SUBROUTINE disturbed_frac(tile,                     & ! in
                            ! ics, ice, iblk,         & ! in
                            dist_type,                & ! in <- DIST_FIRE or DIST_WINDBREAK
                            pft_mask,                 & ! in <- If this pft should have disturbance then true
                            is_woody,                 & ! in <- would not be necessary for WINDBREAK
                            no_children,              & ! in
                            fract_st,                 & ! in <- would not be necessary for WINDBREAK
                            cover_frac_pot,           & ! in <- would not be necessary for DIST_FIRE 
                                                           !    <- in JSBACH3: name-switch to act_fpc
                            veg_frac_correction_st,   & ! in <- would not be necessary for WINDBREAK
                            q_rel_air,                & ! in <- would not be necessary for WINDBREAK
                            wind_10m,                 & ! in <- would not be necessary for WINDBREAK
                            !mem_process_ptr,         & ! in R: only to find FUNCTION Average_over_children
                            max_wind_10m_act,         & ! inout <- would not be necessary for WINDBREAK
                            prev_day_max_wind_10m,    & ! inout <- would not be necessary for DIST_FIRE
                            max_wind_10m,             & ! inout <- would not be necessary for DIST_FIRE
                            q_rel_air_climbuf,        & ! inout <- would not be necessary for WINDBREAK
                            q_rel_air_climbuf_yDay,   & ! out   <- would not be necessary for WINDBREAK
                            dist_frac,              & ! inout  <- burned_frac or damaged_frac
                            Ycpool_acid_ag1_st,     & ! in, optional, only for fire
                            YCpool_water_ag1,       & ! in, optional, only for fire
                            YCpool_ethanol_ag1,     & ! in, optional, only for fire
                            YCpool_nonsoluble_ag1,  & ! in, optional, only for fire
                            YCpool_acid_ag2,        & ! in, optional, only for fire
                            YCpool_water_ag2,       & ! in, optional, only for fire
                            YCpool_ethanol_ag2,     & ! in, optional, only for fire
                            YCpool_nonsoluble_ag2,  & ! in, optional, only for fire
                            fuel                    & ! inout, optional, only for fire
                            )

    USE mo_jsb_time,                ONLY: is_time_experiment_start, time_steps_per_day, is_newday
    USE mo_disturbance_constants,   ONLY: persist_rel_hum, persist_wind_10m, DIST_FIRE, DIST_WINDBREAK
    !USE mo_carbon_memory_class,     ONLY: t_carbon_memory
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_model_class,         ONLY: t_jsb_model

    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile ! R: see X1

    ! Arguments
    !INTEGER, INTENT(in)                     :: ics, ice, iblk
    INTEGER, INTENT(in)                     :: dist_type               ! Sub-disturbance type (one of the DIST_ constants)
    LOGICAL, INTENT(in)                     :: pft_mask(:)             ! True for pfts to which the disturb process should be applied
    LOGICAL, INTENT(in)                     :: is_woody(:)
    INTEGER, INTENT(in)                     :: no_children             ! Number of existing sub tiles
    REAL(wp), INTENT(in)                    :: fract_st(:,:)           ! Dimension: (ics:ice,no_children).
    REAL(wp), INTENT(in)                    :: cover_frac_pot(:,:)
    REAL(wp), INTENT(in)                    :: veg_frac_correction_st(:,:)

    REAL(wp), INTENT(in)                    :: q_rel_air(:)            ! R: In JS3 this is coming from: mo_soil.f90-2967:
    REAL(wp), INTENT(in)                    :: wind_10m(:)
    !TYPE(t_carbon_memory), POINTER, INTENT(in) :: mem_process_ptr      ! Points to tile%mem(CARBON_)%p% ; only to find
                                                                       ! FUNCTION Average_over_children
    REAL(wp), INTENT(inout)                 :: max_wind_10m_act(:)
    REAL(wp), INTENT(inout)                 :: prev_day_max_wind_10m(:)
    REAL(wp), INTENT(inout)                 :: max_wind_10m(:)
    REAL(wp), INTENT(inout)                 :: q_rel_air_climbuf(:)      ! Needed because we have to keep it in mind for
                                                                         ! timesteps within a day
    REAL(wp), INTENT(out)                   :: q_rel_air_climbuf_yDay(:) ! Smoothend relative humidity at new_day
                                                                         ! (needed as cbalone forcing)

    REAL(wp), INTENT(inout)                 :: dist_frac(:,:)          ! Fraction damaged by wind or fire, excluding survivors
    REAL(wp), OPTIONAL, INTENT(in)          :: Ycpool_acid_ag1_st(:,:) ! The follwoing Yasso C pools are all on canopy area.
                                                                       ! Dimension: (ics:ice,no_children).
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_water_ag1(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_ethanol_ag1(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_nonsoluble_ag1(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_acid_ag2(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_water_ag2(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_ethanol_ag2(:,:)
    REAL(wp), OPTIONAL, INTENT(in)          :: YCpool_nonsoluble_ag2(:,:)
    REAL(wp), OPTIONAL, INTENT(inout)       :: fuel(:)         ! Amount of fuel on this parent-tile on which the fire was based

    ! Local variables
    LOGICAL                                 :: lstart
    REAL(wp)                                :: persist ! Factor (as fraction) determining the persistance of q_rel_air_climbuf
                                                       ! relative to the time-step-actual q_rel_air
    dsl4jsb_Def_config(DISTURB_)
    CLASS(t_jsb_model), POINTER             :: model

    ! Necessary for dsl4jsb_Config below!
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(DISTURB_)


    ! ---------------------------
    ! Go

    dist_frac(:,:) = 0._wp
    lstart = is_time_experiment_start()

    ! R: This may be later in a climbuf process for dynveg and disturbance.
    ! Maximum daily wind speed and relative air humidity smoothed in time
    persist = persist_rel_hum ** (1._wp/REAL(time_steps_per_day(),wp))
    q_rel_air_climbuf(:) = q_rel_air_climbuf(:) * persist  +  MIN(q_rel_air(:),100._wp) * (1._wp - persist)

    IF (.NOT. is_newday() .OR. lstart) THEN
      WHERE (wind_10m(:) > max_wind_10m_act(:))
         max_wind_10m_act(:) = wind_10m(:)
      END WHERE
    ELSE
      prev_day_max_wind_10m(:) = max_wind_10m_act(:)
      max_wind_10m(:) = max_wind_10m(:) * persist_wind_10m + max_wind_10m_act(:) * (1._wp - persist_wind_10m)
      max_wind_10m_act(:) = wind_10m(:)
      q_rel_air_climbuf_yDay(:) = q_rel_air_climbuf(:)
    END IF

    SELECT CASE (dist_type)
      CASE (DIST_FIRE)
        SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

                CALL burned_frac_jsbach(tile, & !ics, ice, iblk,              & ! in
                                        pft_mask(:),                 & ! in
                                        is_woody(:),                 & ! in
                                        no_children,                 & ! in
                                        fract_st(:,:),               & ! in
                                        veg_frac_correction_st(:,:), & ! in
                                        q_rel_air_climbuf(:),        & ! in
                                        !mem_process_ptr,             & ! in
                                        dist_frac(:,:),              & ! inout <- bruned_frac in the called procedure
                                        fuel(:),                     & ! out
                                        Ycpool_acid_ag1_st(:,:),     & ! in
                                        YCpool_water_ag1(:,:),       & ! in
                                        YCpool_ethanol_ag1(:,:),     & ! in
                                        YCpool_nonsoluble_ag1(:,:),  & ! in
                                        YCpool_acid_ag2(:,:),        & ! in
                                        YCpool_water_ag2(:,:),       & ! in
                                        YCpool_ethanol_ag2(:,:),     & ! in
                                        YCpool_nonsoluble_ag2(:,:)   & ! in
                                        )

          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_frac','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_frac','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_frac
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown fire algorithm')
        END SELECT
        ! Correct burned fpc for survivors dependent on fuel_consumption algorithm when dynveg = TRUE
!         SELECT CASE (dsl4jsb_Config(DISTURB_)%fuel_algorithm)
!           CASE (0,1) ! No/jsbach fuel consumption. With and without survivors - data "copy" by initial pointer assignments
!           CASE (2)   ! Arora & Boer fuel consumption - survivors considered
!             CALL finish('disturbed_frac','Arora & Boer algorithm not implemented yet.')
!           CASE (3)!! Thonicke spitfire
!             CALL finish('disturbed_frac','Thornicke spit fire not implemented yet.')
!           CASE (4)   ! Read (GFED) fuel consumption - survivors considered 
!         END SELECT

      CASE (DIST_WINDBREAK)
        SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm
            CALL broken_woody_frac_jsbach(tile, &
                                          pft_mask(:),                & ! in
                                          no_children,                & ! in
                                          cover_frac_pot(:,:),        & ! in
                                          prev_day_max_wind_10m(:),   & ! in
                                          max_wind_10m(:),            & ! in
                                          dist_frac(:,:)              & ! inout <- damaged_frac in the called procedure
                                          )

          ! Insert other wind break algorithms here
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown windbreak algorithm')
        END SELECT
      CASE DEFAULT
        CALL finish('disturbed_fpc','Unknown disturbance type')
    END SELECT

  END SUBROUTINE disturbed_frac




  ! Calculate burned area. Woody and grass types call the subroutine with different parameters 
  SUBROUTINE burned_frac_jsbach(tile,  & !ics, ice, iblk,         & ! in
                                pft_mask,                 & ! in
                                is_woody,                 & ! in
                                no_children,              & ! in
                                fract_st,                 & ! in
                                veg_frac_correction_st,   & ! in
                                q_rel_air_climbuf,        & ! in
                                !mem_process_ptr,          & ! in
                                burned_frac_st,           & ! inout
                                fuel,                     & ! out
                                Ycpool_acid_ag1_st,       & ! in
                                YCpool_water_ag1_st,      & ! in
                                YCpool_ethanol_ag1_st,    & ! in
                                YCpool_nonsoluble_ag1_st, & ! in
                                YCpool_acid_ag2_st,       & ! in
                                YCpool_water_ag2_st,      & ! in
                                YCpool_ethanol_ag2_st,    & ! in
                                YCpool_nonsoluble_ag2_st  & ! in
                                )

     ! See X1
     USE mo_jsb_class,        ONLY: Get_model
     USE mo_jsb_memory_class, ONLY: t_jsb_memory
     USE mo_jsb_model_class,  ONLY: t_jsb_model

    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile ! R: see X1
    ! Input Arguments
    !INTEGER, INTENT(in)   :: ics, ice, iblk
    LOGICAL, INTENT(in)   :: pft_mask(:)   ! True for pfts to which the disturbance should be applied.
    LOGICAL, INTENT(in)   :: is_woody(:)
    INTEGER, INTENT(in)   :: no_children
    REAL(wp), INTENT(in)  :: fract_st(:,:) !  Dimension: (ics:ice,no_children).
    REAL(wp), INTENT(in)  :: veg_frac_correction_st(:,:) 
    REAL(wp), INTENT(in)  :: q_rel_air_climbuf(:)
    !TYPE(t_carbon_memory), POINTER, INTENT(in) :: mem_process_ptr

    REAL(wp), INTENT(in)  :: Ycpool_acid_ag1_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_water_ag1_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_ethanol_ag1_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_nonsoluble_ag1_st(:,:)
    REAL(wp), INTENT(in)  :: Ycpool_acid_ag2_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_water_ag2_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_ethanol_ag2_st(:,:)
    REAL(wp), INTENT(in)  :: YCpool_nonsoluble_ag2_st(:,:)

    ! Output Arguments
    REAL(wp), INTENT(INOUT) :: burned_frac_st(:,:) ! = dist_frac from calling procedure
    REAL(wp), INTENT(OUT)   :: fuel(:)             ! Amount of fuel on this parent-tile on which the fire was based

    ! Locals
    INTEGER  :: itile, nc
    INTEGER  :: no_omp_thread
    REAL(wp) :: delta_time_yr, &
                 fire_tau
                !fire_rel_hum_threshold, & ! maximal relative humidity for fire
                !fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                !fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                !fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                !fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                !fire_tau_grass,         & ! return period of fire for grass PFT [year] at 0% relative humidity

    REAL(wp), allocatable ::  fire_minimum(:),        &
                              fire_rel_hum_threshold_expanded(:)

    TYPE(t_jsb_model), POINTER :: model

    ! See X1
    dsl4jsb_Def_config(DISTURB_)
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(DISTURB_)

    no_omp_thread = Get_omp_thread()
    nc      = model%options(no_omp_thread)%nc

    ! ---------------------------
    ! Go

    delta_time_yr = 1._wp / 365._wp

    ! R: Note, in JS3 to calculate fuel, average_tiles is done in combination with the conversion to gridbox area.
    !    The average on the gridbox (including desert areas!) is then corrected so that the result is the average
    !     over the land tiles area, including vegetation, bare soil, glacier and lake!
    !
    !    In JS4 we should exclude glacier and lake in the average because this can only impair the parametrization
    !    of fire-fuel! In JS3 this was only done because we had no tiles/fractions of glacier and lake. Thus the
    !    results of fire on glacier and lake grid boxes was zero anyhow.
    !
    !    In JS4 the concept of vegetation tiles shall include uncovered "ground" areas (e.g. bare soil, stone desert).
    !    Assumed e.g. in radiation in calc_sky_view_fractions or in the variable veg_frac_correction.
    !    Therefore, we should assume desert tiles as separated from the veg tile (well separated) and not as well mixed,
    !    since the "bare ground" in reality that is well mixed with vegetation is already part of our concept of veg tiles!
    !    => in JS4 we should exclude desert tiles in the calculation of fuel. => fuel should be calculated on veg tile area!
    !
    !    In each case the thresholds for fire have to be adapted to JS4.
    !
    ! R: Normally we could average over the children with this routine:
    ! CALL mem_process_ptr%Average_over_children_real2d(ics, ice, iblk, no_children, mem_process_ptr,
    !      mem_process_ptr%Ycpool_acid_ag1_st, fract_st(:,:),    var_averaged(:) )
    !    However, in this case we do it better like that:
    fuel(:) = SUM(                                                                    &
                  (   Ycpool_acid_ag1_st   (:,:)   + YCpool_water_ag1_st     (:,:)    & ! C pools related to child 
                    + YCpool_ethanol_ag1_st(:,:)   + YCpool_nonsoluble_ag1_st(:,:)    & ! tile canopy area ...
                    + YCpool_acid_ag2_st   (:,:)   + YCpool_water_ag2_st     (:,:)    &
                    + YCpool_ethanol_ag2_st(:,:)   + YCpool_nonsoluble_ag2_st(:,:)  ) &
                  * veg_frac_correction_st(:,:)                                       & ! ... now related to whole child tile area...
                  * fract_st(:,:)                                                     & ! ... wheightet by the relative fracts of the
                                                                                        !     children to the parent tile...
                  , DIM=2 )                                                             ! ... and summed up over the children ...
                                                                                        ! ... results in the average C pool related
                                                                                        !     to the parent (=veg) tile area...
                                                                                        ! ... because the fract_st sum over DIM=2 is
                                                                                        ! always = 1...
    fire_rel_hum_threshold_expanded(1:nc) = dsl4jsb_Config(DISTURB_)%fire_rel_hum_threshold
    DO itile = 1,no_children
       ! See X1. This should be moved to the interface and give the results as arguments
       IF ( is_woody(no_children) ) THEN  ! Woody Type
         fire_minimum(1:nc) = dsl4jsb_Config(DISTURB_)%fire_minimum_woody
         fire_tau        = dsl4jsb_Config(DISTURB_)%fire_tau_woody
       ELSE                             ! all other types
         fire_minimum(1:nc) = dsl4jsb_Config(DISTURB_)%fire_minimum_grass
         fire_tau        = dsl4jsb_Config(DISTURB_)%fire_tau_grass
       ENDIF

       IF (pft_mask(itile)) THEN
          burned_frac_st(:,itile) = fire_minimum * delta_time_yr

          WHERE ( q_rel_air_climbuf(:) < dsl4jsb_Config(DISTURB_)%fire_rel_hum_threshold   &
                 .AND. fuel(:) > dsl4jsb_Config(DISTURB_)%fire_litter_threshold )
            burned_frac_st(:,itile) = ( fire_minimum(:)                                                     &
                                      + ( (fire_rel_hum_threshold_expanded(:) - q_rel_air_climbuf(:))       &
                                          / (fire_tau * dsl4jsb_Config(DISTURB_)%fire_rel_hum_threshold)  ) &
                                       ) * delta_time_yr
          END WHERE
        ENDIF
        !ENDDO
    ENDDO

  END SUBROUTINE burned_frac_jsbach

  ! Calculate windbreak for woody types
  SUBROUTINE broken_woody_frac_jsbach(tile, & ! in
                                      pft_mask,              & ! in <- = is woody()
                                      no_children,           & ! in
                                      cover_frac_pot_st,     & ! in
                                      prev_day_max_wind_10m, & ! in
                                      max_wind_10m,          & ! in
                                      damaged_frac_st        & ! inout
                                      )
    ! See X1
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile ! R: see X1

    ! Arguments
    LOGICAL , INTENT(in)    :: pft_mask(:)
    INTEGER , INTENT(in)    :: no_children
    REAL(wp), INTENT(in)    :: cover_frac_pot_st(:,:)
    REAL(wp), INTENT(in)    :: prev_day_max_wind_10m(:)
    REAL(wp), INTENT(in)    :: max_wind_10m(:)
    REAL(wp), INTENT(INOUT) :: damaged_frac_st(:,:)

    ! Locals
    INTEGER  :: itile, i, no_omp_thread
    REAL(wp) :: delta_time_yr

    ! See X1
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(DISTURB_)
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(DISTURB_)
    no_omp_thread = Get_omp_thread()

    ! ---------------------------
    ! Go

    delta_time_yr = 1._wp / 365._wp
    DO itile = 1,no_children
       IF (pft_mask(itile)) THEN
         ! where windbreak should occur
         DO i = 1,model%options(no_omp_thread)%nc
            IF ( cover_frac_pot_st(i,itile) > EPSILON(1._wp) &
                 .AND.  prev_day_max_wind_10m(i) > max_wind_10m(i) * dsl4jsb_Config(DISTURB_)%wnd_threshold ) THEN
               damaged_frac_st(i,itile) =  MIN(1._wp - EPSILON(1._wp) / cover_frac_pot_st(i,itile),            &
                                           delta_time_yr * dsl4jsb_Config(DISTURB_)%wnd_damage_scale           &
                                           * prev_day_max_wind_10m(i) ** 3._wp / max_wind_10m(i) )
            ENDIF
         ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE broken_woody_frac_jsbach



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
           CALL finish(TRIM(routine), 'air_qsat nearly 1.')
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
    USE mo_exception,            ONLY: message_text
    USE mo_echam_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2
    INTEGER,  INTENT(IN) :: nc
    REAL(wp), INTENT(IN) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sat_specific_humidity'

    INTEGER :: i
    REAL(wp) :: vtmpc1
    INTEGER, allocatable  :: it(:)
    REAL(wp), allocatable :: tluc(:)

    !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
    DO i = 1,nc
      IF ((temp(i) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)) .OR. (temp(i) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp))) THEN
        WRITE (message_text,*) &
          & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
          & '- temperature was: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
    ENDDO

    vtmpc1 = 0.6077686814143877_wp
    it = NINT(temp*1000._wp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
      qsat = tluc / (pressure - vtmpc1*tluc)
    ELSEWHERE
      qsat = HUGE(1._wp)
    END WHERE

  END SUBROUTINE sat_specific_humidity

#endif
END MODULE mo_disturb_process
