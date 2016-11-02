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

#include "dsl_jsbach.inc"

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
  SUBROUTINE update_phenology(tile)

    USE mo_base_constants,          ONLY: dt, one_day
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_base_memory,             ONLY: t_slm_tile, t_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_shared_constants,        ONLY: Tzero
    USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_veg_constants,           ONLY: ievergreen,iraingreen,isummergreen,iperennial, & 
                                      gdd_req_max,gdd_t_air_threshold,k_gdd_dormance,t_air_senescence, & 
                                      beta_soil_flush_grass,beta_soil_flush_tree,beta_soil_senescence_tree, & 
                                      min_leaf_age

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)        :: tile         !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp) :: hlp1                                      ! helper variable for tmp values
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory),      POINTER :: veg               ! local pointer to veg__mem for SLM usability
    REAL(wp), POINTER                :: t_air
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)
!      TYPE IS(t_atm2land_memory)
!        atm2land__mem => mm
!    END SELECT
!    SELECT TYPE (mm => tile%mem(VEG_)%p)
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg     => veg__mem 
    t_air   => atm2land__mem%t_air%ptr

    model   => Get_model(tile%owner_model_id)
    lctlib  => model%lctlib(tile%lcts%lib_id)


    !> 1.0 check growing season has started or ended
    !! 
    IF(.NOT.veg%growing_season%ptr) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (lctlib%phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(veg%beta_soil_gs_tphen_mavg%ptr > beta_soil_flush_tree) veg%growing_season%ptr = .TRUE.
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max / exp ( k_gdd_dormance * veg%nd_dormance%ptr ) 
            IF(veg%gdd%ptr > hlp1) veg%growing_season%ptr = .TRUE.
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            IF(veg%t_air_tphen_mavg%ptr > gdd_t_air_threshold .AND. &
               veg%beta_soil_gs_tphen_mavg%ptr > beta_soil_flush_grass) veg%growing_season%ptr = .TRUE.
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (lctlib%phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(veg%beta_soil_gs_tphen_mavg%ptr < beta_soil_senescence_tree .AND. veg%mean_leaf_age%ptr > min_leaf_age ) & 
               veg%growing_season%ptr = .FALSE.
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
                                 ! leaves have been out a given number of days
            IF(veg%t_air_week_mavg%ptr < veg%t_air_month_mavg%ptr .AND. &
                                veg%t_air_tphen_mavg%ptr < t_air_senescence .AND. & 
                                veg%mean_leaf_age%ptr > min_leaf_age) veg%growing_season%ptr = .FALSE.
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((veg%t_air_tphen_mavg%ptr < Tzero .OR. veg%beta_soil_gs_tphen_mavg%ptr < eps4 .OR. &
                veg%gpp_tlabile_mavg%ptr < veg%maint_respiration_tlabile_mavg%ptr).AND. &
                veg%mean_leaf_age%ptr > min_leaf_age) &
                veg%growing_season%ptr = .FALSE.
       END SELECT
    ENDIF   

    !> 2.1 gdd counter
    !! 
    IF(t_air > gdd_t_air_threshold) THEN
       veg%gdd%ptr = veg%gdd%ptr + ((t_air - gdd_t_air_threshold) * dt/one_day)
    ENDIF
    IF(veg%t_air_week_mavg%ptr < Tzero) veg%gdd%ptr = 0.0_wp

    !> 2.2 nd_dormance counter
    !! 
    IF (veg%growing_season%ptr) THEN 
      veg%gdd%ptr         = 0.0_wp
      veg%nd_dormance%ptr = 0.0_wp
    ELSE
      veg%nd_dormance%ptr = veg%nd_dormance%ptr + dt/one_day 
    ENDIF

  END SUBROUTINE update_phenology 

END MODULE mo_veg_phenology
