!>
!! @brief module to calculate vegetation dynamics 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-04-14)
!!
MODULE mo_veg_dynamics

#include "dsl_jsbach.inc"

  USE mo_kind,                 ONLY: wp, eps4


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_dynamics

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculate establishment and mortality of a population (currently static mortality rate only) 
  !!
  !! a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !!
  !! flux\%litterfall & flux\%growth is set to zero each timestep in the plant model driver 
  !-----------------------------------------------------------------------------------------------------
  !PURE ELEMENTAL SUBROUTINE update_veg_dynamics(tile, lctlib) -> currently not working *slmdev
  SUBROUTINE update_veg_dynamics(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants
    USE mo_shared_constants,      ONLY: Tzero
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_veg_constants
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)     :: tile           !< current tile 
    ! 0.2 Local
    TYPE(t_model),          POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                                  :: mortality_rate
    REAL(wp)                                  :: establishment_rate
    REAL(wp)                                  :: boc_cbalance
    REAL(wp)                                  :: tc
    REAL(wp)                                  :: rfr_ratio
    REAL(wp)                                  :: flim_light,flim_moist,flim_temp
    INTEGER :: icanop
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg => veg__mem 

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)


    !> 1.0 background mortality 
    !! select rate and convert from 1/year to 1/timestep
    SELECT CASE (lctlib%growthform)
      CASE (itree)
         mortality_rate = ((background_mort_rate_tree /  one_year) / one_day) * dt
      CASE (igrass)
         mortality_rate = ((background_mort_rate_grass / one_year) / one_day) * dt
    ENDSELECT 

    IF(veg_config%flag_tile_veg_dynamics) THEN 

       !> 2.0 self-thinning related mortality
       !! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost
       !! in mol / m2 LAI / year
       !! transformed to a mortality estimate using a Weibull function
       boc_cbalance = veg%an_boc_tvegdyn_mavg%ptr * one_day * one_year / 1.e6 & 
                    - ( 1._wp + fresp_growth ) * lctlib%tau_leaf/lctlib%sla
       IF(boc_cbalance < -eps4 )THEN
          mortality_rate = mortality_rate & 
                           + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance ) ** k_mort_light)) & 
                           / one_year / one_day * dt
       ENDIF

       !> 3.0 establishment rate
       !!
       !> 3.1 light limitation
       !! depending on the red to far-red ratio of the light at the bottom of the canopy
       tc = 0.0_wp
       DO icanop=1,ncanop
         tc = tc + veg%nl_tot_cl%ptr(icanop,1) * veg%fnl_chl_cl%ptr(icanop,1) * & 
                   veg%lai_cl%ptr(icanop,1)
       ENDDO
       rfr_ratio = toc_rfr_ratio * exp( k_r2fr_chl * tc) 
       flim_light = 1._wp-exp(-(lctlib%lambda_est_light*rfr_ratio)**lctlib%k_est_light)

       !> 3.2 temperature and soil moisture limitation 
       !! OBS: should be top soil moisture!
       IF(veg%t_air_week_mavg%ptr - Tzero < 0.0_wp.OR..NOT.veg%growing_season%ptr) THEN 
          flim_temp=0.0_wp
       ELSE
          flim_temp = 1.0_wp - exp(-(lambda_est_temp * (veg%t_air_week_mavg%ptr - Tzero))**k_est_temp)
       ENDIF
       flim_moist = 1.0_wp - exp(-(lambda_est_moist * tile%soil%theta)**k_est_moist) 
   
       !> 3.3 actual establishment rate and associated matter flux from seed bed to reserve pool
       !!
       establishment_rate = flim_light * flim_temp * flim_moist * & 
                            lctlib%tau_seed_est / one_year / one_day * dt
       veg%establishment  = veg%pool%seed_bed * establishment_rate

       !> 4.0 diagnose change in the number of individuals (#/m2)
       !!
       IF(lctlib%growthform == itree) THEN
          veg%delta_dens_ind%ptr = veg%establishment%carbon%ptr / lctlib%seed_size - & 
                                   mortality_rate * veg%dens_ind%ptr 
       END IF

    ENDIF

    !> 5.0 calculate litter fall from mortality
    !!
    !! for all tissues but seed_bed
!! previous code  *slmdev
!    CALL mortality%init(mortality_rate)   ! set mortality to mortality_rate
!    CALL mortality%seed_bed%init(0.0_wp)  ! set seed_bed mortality to zero as its unaffected by actual plant death
!    litter_flux_mortality = mortality * tile%plant%pool
!    tile%plant%litterfall = tile%plant%litterfall + litter_flux_mortality
!! new code replacing above lines and save using: TYPE(t_slm_veg_pool) :: mortality, litter_flux_mortality
    veg%litterfall%leaf         = veg%litterfall%leaf         + veg%pool%leaf * mortality_rate
    veg%litterfall%fine_root    = veg%litterfall%fine_root    + veg%pool%fine_root * mortality_rate
    veg%litterfall%coarse_root  = veg%litterfall%coarse_root  + veg%pool%coarse_root * mortality_rate
    veg%litterfall%sap_wood     = veg%litterfall%sap_wood     + veg%pool%sap_wood * mortality_rate
    veg%litterfall%heart_wood   = veg%litterfall%heart_wood   + veg%pool%heart_wood * mortality_rate
    veg%litterfall%labile       = veg%litterfall%labile       + veg%pool%labile * mortality_rate
    veg%litterfall%reserve      = veg%litterfall%reserve      + veg%pool%reserve * mortality_rate
    veg%litterfall%fruit        = veg%litterfall%fruit        + veg%pool%fruit * mortality_rate

  END SUBROUTINE update_veg_dynamics

END MODULE mo_veg_dynamics
