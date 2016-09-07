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
  !PURE ELEMENTAL SUBROUTINE update_veg_dynamics(tile, lctlib)
  SUBROUTINE update_veg_dynamics(tile, lctlib)

    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_base_constants
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_shared_constants,      ONLY: Tzero
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_veg_constants
    USE mo_veg_memory

    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_veg_memory_class,      ONLY: t_veg_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)         :: lctlib         !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout)            :: tile           !< current tile 
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
! --------------------------------
! @slmdev @ slmtemp see ticket #79
    ! 0.3 Memory allocation
    TYPE(t_veg_memory), POINTER         :: veg__mem   ! dsl4jsb_Def_memory(VEG_)
    TYPE(t_veg_memory), POINTER         :: plant1     ! local pointer to veg__mem for SLM usability
    REAL(wp), POINTER :: height
    REAL(wp), POINTER :: diameter

    SELECT TYPE (mm => tile%mem(VEG_)%p)       ! dsl4jsb_Get_memory(VEG_)
      TYPE IS(t_veg_memory)
        veg__mem => mm
    END SELECT
    plant1 => veg__mem          ! @slmdev  -> plant, vegetation, vegm, ... ?


  ! using the memory of the veg process exemplary
    plant1%height%ptr    = plant1%height%ptr + 0.001_wp
    plant1%diameter%ptr  = plant1%diameter%ptr + 0.001_wp

    plant1%pool%leaf%carbon%ptr    = plant1%pool%leaf%carbon%ptr    + 0.001_wp
    plant1%pool%leaf%nitrogen%ptr  = plant1%pool%leaf%nitrogen%ptr  + 0.001_wp
    plant1%growth%leaf%carbon%ptr     = plant1%growth%leaf%carbon%ptr     + 0.001_wp

    CALL test_routine(plant1%litterfall, &                ! adds 0.001_wp to litterfall%labile%carbon%ptr
                      plant1%litterfall%leaf%carbon%ptr)  ! adds 0.001_wp to litterfall%leaf%carbon%ptr

    ! just for testing purposes   @slmtemp
!    write(999,*) plant1%height%ptr, plant1%diameter%ptr, &
!                 plant1%pool%leaf%carbon%ptr, plant1%pool%leaf%nitrogen%ptr, &
!                 plant1%growth%leaf%carbon%ptr, &
!                 plant1%litterfall%labile%carbon%ptr, plant1%litterfall%leaf%carbon%ptr, ' |-| '
! @slmdev @ slmtemp see ticket #79
! --------------------------------



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
       boc_cbalance = tile%plant%boc_an_mavg * one_day * one_year / 1.e6 & 
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
         tc = tc + tile%plant%nl_tot_cl(icanop) * tile%plant%fnl_chl_cl(icanop)
       ENDDO
       rfr_ratio = toc_rfr_ratio * exp( k_r2fr_chl * tc) 
       flim_light = 1._wp-exp(-(lctlib%lambda_est_light*rfr_ratio)**lctlib%k_est_light)
   
       !> 3.2 temperature and soil moisture limitation 
       !! OBS: should be top soil moisture!
       flim_temp = 1._wp-exp(-(lambda_est_temp * (tile%plant%t_air_week_mavg - Tzero))**k_est_temp)
       flim_moist = 1.0_wp - exp(-(lambda_est_moist * tile%soil%theta)**k_est_moist) 
   
       !> 3.3 actual establishment rate and associated matter flux from seed bed to reserve pool
       !!
       establishment_rate = flim_light * flim_temp * flim_moist * & 
                            lctlib%tau_seed_est / one_year / one_day * dt
       tile%plant%establishment = establishment_rate * tile%plant%pool%seed_bed 

       !> 4.0 diagnose change in the number of individuals (#/m2)
       !!
       IF(lctlib%growthform==itree) &
       tile%plant%delta_dens_ind = & 
              tile%plant%establishment%carbon/lctlib%seed_size & 
              - mortality_rate * tile%plant%dens_ind 

    ENDIF

    !> 5.0 calculate litter fall from mortality
    !!
    CALL mortality%init(mortality_rate)   ! set mortality to mortality_rate
!@todo: veg_dynamics: it would probably make things easier if the seed bed was part of the plant, but not the pool structure?
    CALL mortality%seed_bed%init(0.0_wp)  ! set seed_bed to zero
    litter_flux_mortality = mortality * tile%plant%pool
    tile%plant%litterfall = tile%plant%litterfall + litter_flux_mortality
   
  END SUBROUTINE update_veg_dynamics



!! <<< temporary test code - non-functional to slm - see ticket #79 - Jan Engel 20160901 >>> !!
  subroutine test_routine(litterfall, test_var)

    use mo_veg_memory_class,  only: t_slm_struct_veg_pool

    type(t_slm_struct_veg_pool), INTENT(inout) :: litterfall
    REAL(wp), INTENT(inout) :: test_var

    test_var = test_var + 0.001_wp
    litterfall%labile%carbon%ptr =litterfall%labile%carbon%ptr + 0.001_wp

  end subroutine test_routine
!! <<< temporary test code - non-functional to slm - see ticket #79 - Jan Engel 20160901 >>> !!



END MODULE mo_veg_dynamics
