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
  USE mo_veg_config_class, ONLY: veg_config

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
  !PURE ELEMENTAL SUBROUTINE update_veg_dynamics(tile, veg_param)
  SUBROUTINE update_veg_dynamics(tile, veg_param)

    USE mo_base_memory, ONLY: t_slm_tile
    USE mo_base_constants
    USE mo_atmland_constants, ONLY: Tzero
    USE mo_shared_memory,    ONLY: init_element
    USE mo_veg_memory
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_param),INTENT(in)          :: veg_param            !< list of vegetation parameters
    TYPE(t_slm_tile),INTENT(inout)            :: tile                 !< current tile 
    ! 0.2 Local
    TYPE(t_slm_veg_pool)                      :: mortality
    TYPE(t_slm_veg_pool)                      :: litter_flux_mortality
    REAL(wp)                                  :: mortality_rate
    REAL(wp)                                  :: establishment_rate
    REAL(wp)                                  :: boc_cbalance
    REAL(wp)                                  :: tc
    REAL(wp)                                  :: rfr_ratio
    REAL(wp)                                  :: flim_light,flim_moist,flim_temp
    INTEGER :: icanop

    !> 1.0 background mortality 
    !! select rate and convert from 1/year to 1/timestep
    SELECT CASE (veg_param%growthform(tile%plant%pftid))
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
       boc_cbalance = tile%plant%moving_average%boc_an * one_day * one_year / 1.e6 & 
                    - ( 1._wp + fresp_growth ) * veg_param%tau_leaf(tile%plant%pftid)/veg_param%sla(tile%plant%pftid)
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
         tc = tc + tile%plant%canopy_layer(icanop)%nl_tot*tile%plant%canopy_layer(icanop)%fnl_chl
       ENDDO
       rfr_ratio = toc_rfr_ratio * exp( k_r2fr_chl * tc) 
       flim_light = 1._wp-exp(-(veg_param%lambda_est_light(tile%plant%pftid)*rfr_ratio)**veg_param%k_est_light(tile%plant%pftid))
   
       !> 3.2 temperature and soil moisture limitation 
       !! OBS: should be top soil moisture!
       flim_temp = 1._wp-exp(-(lambda_est_temp * (tile%plant%moving_average%t_air_week - Tzero))**k_est_temp)
       flim_moist = 1.0_wp - exp(-(lambda_est_moist * tile%soil%theta)**k_est_moist) 
   
       !> 3.3 actual establishment rate and associated matter flux from seed bed to reserve pool
       !!
       establishment_rate = flim_light * flim_temp * flim_moist * & 
                            veg_param%tau_seed_est(tile%plant%pftid) / one_year / one_day * dt
       tile%plant%flux%establishment = establishment_rate * tile%plant%pool%seed_bed 

       !> 4.0 diagnose change in the number of individuals (#/m2)
       !!
       IF(veg_param%growthform(tile%plant%pftid)==itree) &
       tile%plant%flux%delta_dens_ind = & 
              tile%plant%flux%establishment%carbon/veg_param%seed_size(tile%plant%pftid) & 
              - mortality_rate * tile%plant%dens_ind 

    ENDIF

    !> 5.0 calculate litter fall from mortality
    !!
    CALL init_veg_pool(mortality,mortality_rate)
!@todo: veg_dynamics: it would probably make things easier if the seed bed was part of the plant, but not the pool structure?
    CALL init_element(mortality%seed_bed,0.0_wp)
    litter_flux_mortality = mortality * tile%plant%pool
    tile%plant%flux%litterfall = tile%plant%flux%litterfall + litter_flux_mortality
   
  END SUBROUTINE update_veg_dynamics

END MODULE mo_veg_dynamics
