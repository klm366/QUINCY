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

  SUBROUTINE update_phenology(phenology_type,gpp_week,maint_respiration_week,mean_leaf_age, & 
                              t_air_week,t_air_month,beta_soil_gs,gdd,nd_dormance,growing_season) 
    ! 
    ! Routine to determine onset and end of the growing season for evergreen, raingreen
    !  cold-deciduous trees and perennial grasses
    ! OBS: currently uncalibrated, in particular rain green phenology
    ! Known Issues: strange interannual variability of grass senescences
    ! Dev: Have more cost-benefit in decision to start/stop growing season
    !      what about seasonality beyond off/on behaviour? -> implies some collaboration between
    !      turnover rates and phenology -> should the two be linked? 
    !
    USE mo_kind, ONLY: eps4
    USE mo_atmland_constants, ONLY: Tzero
    USE mo_veg_parameters, ONLY: ievergreen,iraingreen,isummergreen,iperennial, & 
                                 gdd_req_max,gdd_t_air_threshold,k_gdd_dormance,t_air_senescence, & 
                                 beta_soil_flush_grass,beta_soil_flush_tree,beta_soil_senescence_tree, & 
                                 min_leaf_age
    IMPLICIT NONE
    ! 0.1 Input
    INTEGER,INTENT(in) :: phenology_type                 ! type of phenology
    REAL(wp),INTENT(in) :: gpp_week                      ! weekly gross primary production (~µmol/m2/s)
    REAL(wp),INTENT(in) :: maint_respiration_week        ! weekly maintenance respiration (~µmol/m2/s)
    REAL(wp),INTENT(in) :: mean_leaf_age                 ! average age of foliage (days)
    REAL(wp),INTENT(in) :: t_air_week                    ! weekly avearge air temperature (K)
    REAL(wp),INTENT(in) :: t_air_month                   ! monthly average air temperature (K)
    REAL(wp),INTENT(in) :: beta_soil_gs                  ! weekly soil moisture stress (--)
    REAL(wp),INTENT(in) :: gdd                           ! growing degree days (days)
    REAL(wp),INTENT(in) :: nd_dormance                   ! number of dormancy days (days)
    ! 0.2 Inout
    LOGICAL, INTENT(inout) :: growing_season             ! whether or not a plant is in the growing season
    ! 0.4 Local
    REAL(wp) :: hlp1 

    IF(.NOT.growing_season) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(beta_soil_gs > beta_soil_flush_tree) growing_season = .TRUE.
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max / exp ( k_gdd_dormance * nd_dormance ) 
            IF(gdd > hlp1 ) growing_season = .TRUE.
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            IF(t_air_week > gdd_t_air_threshold .AND. beta_soil_gs > beta_soil_flush_grass) & 
               growing_season = .TRUE.
         CASE DEFAULT            ! abord if here
            WRITE(message_text,'(a,i)') 'invalid phenology type ',phenology_type
            CALL finish("mo_veg_phenology:update_phenology",message_text) 
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(beta_soil_gs < beta_soil_senescence_tree .AND. mean_leaf_age > min_leaf_age ) & 
               growing_season = .FALSE.
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
                                 ! leaves have been out a given number of days
            IF(t_air_week < t_air_month .AND. t_air_week < t_air_senescence .AND. mean_leaf_age > min_leaf_age) & 
               growing_season = .FALSE.
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((t_air_week < Tzero .OR. beta_soil_gs < eps4 .OR. gpp_week < maint_respiration_week) & 
                .AND. mean_leaf_age > min_leaf_age) & 
               growing_season = .FALSE.
         CASE DEFAULT            ! abord if here
            WRITE(message_text,'(a,i)') 'invalid phenology type ',phenology_type
            CALL finish("mo_veg_phenology:update_phenology",message_text) 
       END SELECT
    ENDIF   

  END SUBROUTINE update_phenology 

END MODULE mo_veg_phenology
