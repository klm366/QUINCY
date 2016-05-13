!>
!! @brief Contains the interface for plant growth and turnover processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_interface

  USE mo_kind

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_vegetation 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> update vegetation pools
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegetation(tile, veg_param)

    USE mo_base_memory,        ONLY: t_slm_soil, t_slm_tile
    USE mo_atmland_constants,  ONLY: Tzero
    USE mo_atmland_memory
    USE mo_land_constants
    USE mo_land_memory,        ONLY: add_elements, subtract_elements
    USE mo_veg_memory
    USE mo_veg_parameters,     ONLY: t_slm_veg_param, gdd_t_air_threshold
    USE mo_veg_photosynthesis, ONLY: update_canopy_fluxes
    USE mo_veg_phenology,      ONLY: update_phenology
    USE mo_veg_turnover
    USE mo_veg_growth
    USE mo_veg_dynamics,       ONLY: update_veg_dynamics
    USE mo_veg_respiration,    ONLY: calc_maintenance_respiration

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_veg_param),INTENT(in)        :: veg_param     !< vegetation type parameters
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER         :: soil           ! soil variables per tile
    TYPE(t_slm_plant),POINTER        :: plant          ! plant variables per tile
    ! 0.2 Local
    REAL(wp) :: hlp1

    !  bi-directional pointer
    plant => tile%plant
    soil  => tile%soil

    !> 1.0 update vegetation parameter switches
    !! 
    !! @todo we might find a better (more general) solution for this
    CALL plant%moving_average%update_veg_switches(atm2land%t_air, gdd_t_air_threshold, Tzero)

    !> 1.1 update phenological state
    !! 
    !! 
    CALL update_phenology(veg_param%phenology_type(plant%pftid), &
                          plant%moving_average%gpp_week, &
                          plant%moving_average%maint_respiration_week, &
                          plant%mean_leaf_age, &
                          plant%moving_average%t_air_week, &
                          plant%moving_average%t_air_month, &
                          plant%moving_average%beta_soil_gs, &
                          plant%moving_average%gdd, &
                          plant%moving_average%nd_dormance, &
                          plant%growing_season)  

    !> 1.2 photosynthesis and leaf respiration 
    !! 
    !! @todo may move elsewhere
    CALL update_canopy_fluxes(tile, veg_param)

    !> 2.0 Calculate growth
    !! 
    !! 
    CALL  calc_growth(atm2land%t_air,soil%theta,plant,veg_param) 

    !> 2.5 Calculate maintenance respiration
    !! 
    !! 
    CALL  calc_maintenance_respiration(atm2land%t_air,atm2land%t_air,plant%pool%sap_wood%nitrogen,plant%pool%fine_root%nitrogen, &
            plant%pool%coarse_root%nitrogen,plant%pool%labile%carbon,plant%pool%labile%carbon13, & 
            plant%pool%labile%carbon14,plant%flux%maint_respiration_leaf,&
            plant%flux%maint_respiration,plant%flux%maint_respiration_c13,plant%flux%maint_respiration_c14)

    !> 3.0 Calculate turnover from tissue aging and senescence
    !! 
    !! 
    CALL calc_turnover(plant,veg_param)
    CALL update_veg_dynamics(plant,veg_param)

    !> 4.1 Add carbon from photosynthesis to the labile pool
    !! 
    !! 
    plant%flux%npp = plant%flux%gpp - plant%flux%growth_respiration - plant%flux%maint_respiration 
    plant%flux%npp_c13 = plant%flux%gpp_c13 - plant%flux%growth_respiration_c13 - plant%flux%maint_respiration_c13 
    plant%flux%npp_c14 = plant%flux%gpp_c14 - plant%flux%growth_respiration_c14 - plant%flux%maint_respiration_c14 
    plant%flux%growth%labile%carbon = plant%flux%growth%labile%carbon + plant%flux%npp * dt / 1.e6
    plant%flux%growth%labile%carbon13 = plant%flux%growth%labile%carbon13 + & 
       plant%flux%npp_c13 * dt / 1.e6
    plant%flux%growth%labile%carbon14 = plant%flux%growth%labile%carbon14 + & 
       plant%flux%npp_c14 * dt / 1.e6

    !> 4.2 Add nutrients from uptake calculations
    !! 
    !! 
    plant%flux%growth%labile%nitrogen = plant%flux%growth%labile%nitrogen + &
       plant%flux%uptake_n * dt / 1.e6
    plant%flux%growth%labile%nitrogen15 = plant%flux%growth%labile%nitrogen15 + & 
       plant%flux%uptake_n15 * dt / 1.e6
    plant%flux%growth%labile%phosphorus = plant%flux%growth%labile%phosphorus + &
       plant%flux%uptake_p * dt / 1.e6

    !> 5.1 update labile pool and reserve 
    !! 
    !! 
    plant%pool%labile = add_elements(plant%pool%labile,plant%flux%reserve_use)
    plant%pool%reserve = subtract_elements(plant%pool%reserve,plant%flux%reserve_use)

    !> 5.2 update fluxes from retranslocation of leaves
    !! 
    !! 
    plant%pool%leaf%nitrogen = plant%pool%leaf%nitrogen - plant%flux%leaf_resorption_n
    plant%pool%labile%nitrogen = plant%pool%labile%nitrogen + plant%flux%leaf_resorption_n
    plant%pool%leaf%phosphorus = plant%pool%leaf%phosphorus - plant%flux%leaf_resorption_p
    plant%pool%labile%phosphorus = plant%pool%labile%phosphorus + plant%flux%leaf_resorption_p
    plant%pool%leaf%nitrogen15 = plant%pool%leaf%nitrogen15 - plant%flux%leaf_resorption_n15
    plant%pool%labile%nitrogen15 = plant%pool%labile%nitrogen15 + plant%flux%leaf_resorption_n15

    !> 6.0 Update plant biogeochemical pools with growth and litterfall
    !! 
    !! 
    plant%pool = subtract_veg_pools(plant%pool,plant%flux%litterfall)
    plant%pool = add_veg_pools(plant%pool,plant%flux%growth)

    !> 7.0 Update plant diagnostics
    !! 
    !! 
    plant%lai = plant%pool%leaf%carbon * veg_param%sla(plant%pftid)
    CALL calc_allometry(plant%pool%sap_wood%carbon,plant%pool%heart_wood%carbon, & 
                          veg_param%wood_density(plant%pftid),veg_param%sla(plant%pftid), &
                          veg_param%allom_k1(plant%pftid),veg_param%allom_k2(plant%pftid), & 
                          veg_param%k_latosa(plant%pftid),veg_param%k_rtos(plant%pftid), &
                          veg_param%growthform(plant%pftid), & 
                          plant%diameter,plant%height,hlp1,hlp1,hlp1,hlp1)
    IF(plant%pool%leaf%carbon .GT. eps8) THEN 
       plant%mean_leaf_age = & 
          (plant%mean_leaf_age + dt/one_day) * (1._wp - plant%flux%growth%leaf%carbon/plant%pool%leaf%carbon)
    ELSE 
       plant%mean_leaf_age = 0.0_wp
    ENDIF

  END SUBROUTINE update_vegetation

END MODULE mo_veg_interface
