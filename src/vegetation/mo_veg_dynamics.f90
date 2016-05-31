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

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message_text, finish 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_veg_dynamics

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculate establishment and mortality of a population (currently static mortality rate only) 
  !!
  !! a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !!
  !! flux%litterfall & flux%growth is set to zero each timestep in the plant model driver 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_dynamics(plant, veg_param)

    USE mo_base_constants
    USE mo_shared_memory,    ONLY: t_element, element_zero
    USE mo_veg_memory
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_param),INTENT(in)          :: veg_param            !< list of vegetation parameters
    TYPE(t_slm_plant),INTENT(inout)           :: plant                !< current plant  
    ! 0.2 Local
    TYPE(t_element)                           :: mortality_spread_element
    TYPE(t_slm_veg_pool)                      :: mortality_spread_pool
    TYPE(t_slm_veg_pool)                      :: litter_from_mortality
    REAL(wp)                                  :: mortality

    ! select rate and convert from 1/year to 1/timestep
    SELECT CASE (veg_param%growthform(plant%pftid))
      CASE (itree)
         mortality = ((background_mort_rate_tree /  one_year) / one_day) * dt
      CASE (igrass)
         mortality = ((background_mort_rate_grass / one_year) / one_day) * dt
      CASE DEFAULT
          WRITE(message_text,'(a,i4)') 'invalid growthform',veg_param%growthform(plant%pftid)
          CALL finish("mo_veg_dynamics:update_iveg_dynamics",message_text) 
    ENDSELECT 
    ! 
    ! apply rate to all elements of a pool
    mortality_spread_element%carbon     = mortality 
    mortality_spread_element%nitrogen   = mortality
    mortality_spread_element%phosphorus = mortality
    mortality_spread_element%carbon13   = mortality
    mortality_spread_element%carbon14   = mortality
    mortality_spread_element%nitrogen15 = mortality
    ! 
    ! apply rate to all pools
    mortality_spread_pool%leaf        = mortality_spread_element
    mortality_spread_pool%fine_root   = mortality_spread_element
    mortality_spread_pool%coarse_root = mortality_spread_element
    mortality_spread_pool%sap_wood    = mortality_spread_element
    mortality_spread_pool%heart_wood  = mortality_spread_element
    mortality_spread_pool%fruit       = mortality_spread_element
    mortality_spread_pool%labile      = mortality_spread_element
    mortality_spread_pool%reserve     = mortality_spread_element
    mortality_spread_pool%seed_bed    = element_zero 

    ! calculate litter fall from mortality
    litter_from_mortality = mortality_spread_pool * plant%pool
    plant%flux%litterfall = plant%flux%litterfall + litter_from_mortality
  END SUBROUTINE update_veg_dynamics

END MODULE mo_veg_dynamics
