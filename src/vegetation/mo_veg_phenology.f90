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

  USE mo_kind,          ONLY: wp, eps4


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
  SUBROUTINE update_phenology(tile, lctlib)

    USE mo_base_constants,          ONLY: dt, one_day
    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_shared_constants,        ONLY: Tzero

    USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_atm2land_memory_class,   ONLY: t_atm2land_memory

    USE mo_veg_memory
    USE mo_veg_constants,      ONLY: ievergreen,iraingreen,isummergreen,iperennial, & 
                                 gdd_req_max,gdd_t_air_threshold,k_gdd_dormance,t_air_senescence, & 
                                 beta_soil_flush_grass,beta_soil_flush_tree,beta_soil_senescence_tree, & 
                                 min_leaf_age
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile      !< one tile
    TYPE(t_lctlib_element),INTENT(in)       :: lctlib    !< land-cover-type library - parameters for one pft
      ! 0.2 Local
    TYPE(t_slm_plant),POINTER               :: plant     ! plant variables per tile
    REAL(wp) :: hlp1                                     ! helper variable for tmp values
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND_)
    REAL(wp), POINTER                :: t_air
    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)            ! dsl4jsb_Get_memory(ATM2LAND_)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    t_air   => atm2land__mem%t_air%ptr


    ! bi-directional pointer
    plant => tile%plant


    !> 1.0 check growing season has started or ended
    !! 
    IF(.NOT.plant%growing_season) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (lctlib%phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(plant%beta_soil_gs_tphen_mavg > beta_soil_flush_tree) plant%growing_season = .TRUE.
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max / exp ( k_gdd_dormance * plant%nd_dormance ) 
            IF(plant%gdd > hlp1 ) plant%growing_season = .TRUE.
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            IF(plant%t_air_tphen_mavg > gdd_t_air_threshold .AND. &
               plant%beta_soil_gs_tphen_mavg > beta_soil_flush_grass) plant%growing_season = .TRUE.
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (lctlib%phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(plant%beta_soil_gs_tphen_mavg < beta_soil_senescence_tree .AND. plant%mean_leaf_age > min_leaf_age ) & 
               plant%growing_season = .FALSE.
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
                                 ! leaves have been out a given number of days
            IF(plant%t_air_week_mavg < plant%t_air_month_mavg .AND. &
                                plant%t_air_tphen_mavg < t_air_senescence .AND. & 
                                plant%mean_leaf_age > min_leaf_age) plant%growing_season = .FALSE.
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((plant%t_air_tphen_mavg < Tzero .OR. plant%beta_soil_gs_tphen_mavg < eps4 .OR. &
                plant%gpp_tlabile_mavg < plant%maint_respiration_tlabile_mavg).AND. &
                plant%mean_leaf_age > min_leaf_age) &
                plant%growing_season = .FALSE.
       END SELECT
    ENDIF   

    !> 2.1 gdd counter
    !! 
    IF(t_air > gdd_t_air_threshold) THEN
       plant%gdd = plant%gdd + ((t_air - gdd_t_air_threshold) * dt/one_day)
    ENDIF
    IF(plant%t_air_week_mavg < Tzero) plant%gdd = 0.0_wp

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
