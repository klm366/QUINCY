!>
!! @brief module to determine the phenological state of a plant 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-04-11)
!!
MODULE mo_veg_phenology

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message_text, finish

  IMPLICIT NONE

  PUBLIC :: update_phenology

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> determine onset and end of the growing season for some PFT
  !! 
  !! PFT: evergreen, raingreen, cold-deciduous trees and perennial grasses \n
  !! 
  !! @todo OBS: currently uncalibrated, in particular rain green phenology \n
  !!  Known Issues: strange interannual variability of grass senescences \n
  !!  Dev: Have more cost-benefit in decision to start/stop growing season \n
  !!      what about seasonality beyond off/on behaviour? -> implies some collaboration between
  !!      turnover rates and phenology -> should the two be linked? 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_phenology(tile, phenology_type)

    USE mo_kind,               ONLY: eps4
    USE mo_base_constants,     ONLY: dt, one_day
    USE mo_base_memory,        ONLY: t_slm_tile
    USE mo_veg_memory
    USE mo_atmland_constants,  ONLY: Tzero
    USE mo_atmland_memory,     ONLY: atm2land
    USE mo_veg_parameters,     ONLY: ievergreen,iraingreen,isummergreen,iperennial, & 
                                 gdd_req_max,gdd_t_air_threshold,k_gdd_dormance,t_air_senescence, & 
                                 beta_soil_flush_grass,beta_soil_flush_tree,beta_soil_senescence_tree, & 
                                 min_leaf_age
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile      !< one tile
    INTEGER,INTENT(in) :: phenology_type                 !< type of phenology
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER               :: plant     ! plant variables per tile
    REAL(wp) :: hlp1                                     ! helper variable for tmp values

    !  bi-directional pointer
    plant => tile%plant


    !> 1.0 check growing season has started or ended
    !! 
    IF(.NOT.plant%growing_season) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(plant%moving_average%beta_soil_gs > beta_soil_flush_tree) plant%growing_season = .TRUE.
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max / exp ( k_gdd_dormance * plant%nd_dormance ) 
            IF(plant%gdd > hlp1 ) plant%growing_season = .TRUE.
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            IF(plant%moving_average%t_air_week > gdd_t_air_threshold .AND. &
               plant%moving_average%beta_soil_gs > beta_soil_flush_grass) plant%growing_season = .TRUE.
         CASE DEFAULT            ! abord if here
            WRITE(message_text,'(a,i4)') 'invalid phenology type ',phenology_type
            CALL finish("mo_veg_phenology:update_phenology",message_text) 
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(plant%moving_average%beta_soil_gs < beta_soil_senescence_tree .AND. plant%mean_leaf_age > min_leaf_age ) & 
               plant%growing_season = .FALSE.
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
                                 ! leaves have been out a given number of days
            IF(plant%moving_average%t_air_week < plant%moving_average%t_air_month .AND. &
                                plant%moving_average%t_air_week < t_air_senescence .AND. & 
                                plant%mean_leaf_age > min_leaf_age) plant%growing_season = .FALSE.
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((plant%moving_average%t_air_week < Tzero .OR. plant%moving_average%beta_soil_gs < eps4 .OR. &
                plant%moving_average%gpp_week < plant%moving_average%maint_respiration_week).AND. &
                plant%mean_leaf_age > min_leaf_age) &
                plant%growing_season = .FALSE.
         CASE DEFAULT            ! abord if here
            WRITE(message_text,'(a,i4)') 'invalid phenology type ',phenology_type
            CALL finish("mo_veg_phenology:update_phenology",message_text) 
       END SELECT
    ENDIF   

    !> 2.1 gdd counter
    !! 
    IF(atm2land%t_air > gdd_t_air_threshold) THEN
       plant%gdd = plant%gdd + ((atm2land%t_air - gdd_t_air_threshold) * dt/one_day)
    ENDIF
    IF(plant%moving_average%t_air_week < Tzero) plant%gdd = 0.0_wp

    !> 2.2 nd_dormance counter
    !! 
    IF (plant%growing_season) THEN 
      plant%gdd         = 0.0_wp
      plant%nd_dormance = 0.0_wp
    ELSE
      plant%nd_dormance = plant%nd_dormance + dt/one_day 
    ENDIF

  END SUBROUTINE update_phenology 

END MODULE mo_veg_phenology
