!>
!! @brief Contains the routines for the Task update_canopy_fluxes (photosynthesis related processes)
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719): svn move vegetation/mo_veg_photosynthesis.f90 assimilation/mo_assimi_update_canopy_fluxes.f90
!!
MODULE mo_veg_update_veg_pools

  USE mo_kind
  USE mo_jsb_lctlib_class,     ONLY: t_lctlib_element


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_pools

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update vegetation pools
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_pools(tile, lctlib)

    USE mo_base_memory,         ONLY: t_slm_tile
    USE mo_base_config_class,   ONLY: base_config
    USE mo_shared_constants        ! only lambda_C14 ?
    USE mo_veg_memory              ! only t_slm_plant ?
    USE mo_veg_canopy,          ONLY: calc_canopy_layers
    USE mo_veg_growth,          ONLY: calc_allometry

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)        :: lctlib        !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER        :: plant          ! plant variables per tile
    REAL(wp)                         :: hlp1           ! helper variable

    !  bi-directional pointer
    plant => tile%plant

    !> 1.0 for plant and land model, update vegetation pools with fluxes
    !! 
    !!
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('plant','land')
   
       !> 1.1 Add carbon from photosynthesis to the labile pool
       !! 
       !! 
       plant%npp = plant%gpp - plant%growth_respiration - plant%maint_respiration - & 
                   plant%n_processing_respiration
       plant%npp_c13 = plant%gpp_c13 - plant%growth_respiration_c13 - plant%maint_respiration_c13 - & 
                   plant%n_processing_respiration_c13
       plant%npp_c14 = plant%gpp_c14 - plant%growth_respiration_c14 - plant%maint_respiration_c14 - &
                   plant%n_processing_respiration_c14
       plant%growth%labile%carbon = plant%growth%labile%carbon + plant%npp * dt / 1000000.0_wp
       plant%growth%labile%carbon13 = plant%growth%labile%carbon13 + & 
          plant%npp_c13 * dt / 1000000.0_wp
       plant%growth%labile%carbon14 = plant%growth%labile%carbon14 + & 
          plant%npp_c14 * dt / 1000000.0_wp
   
       !> 1.2 Add nutrients from uptake calculations
       !! 
       !! 
       plant%growth%labile%nitrogen = plant%growth%labile%nitrogen + &
          ( plant%uptake_n + plant%n_fixation ) * dt / 1000000.0_wp
       plant%growth%labile%nitrogen15 = plant%growth%labile%nitrogen15 + & 
          ( plant%uptake_n15 + plant%n_fixation_n15) * dt / 1000000.0_wp
       plant%growth%labile%phosphorus = plant%growth%labile%phosphorus + &
          plant%uptake_p * dt / 1000000.0_wp
   
       !> 1.3 update labile pool and reserve 
       !! 
       plant%pool%labile  = plant%pool%labile  + plant%reserve_use
       plant%pool%reserve = plant%pool%reserve - plant%reserve_use
   
       !> 1.4 update fluxes from continuous nutrient turnover and resorption 
       !!
       plant%pool%labile%nitrogen = plant%pool%labile%nitrogen + plant%recycling_leaf_n & 
            + plant%recycling_fine_root_n + plant%recycling_sap_wood_n
       plant%pool%leaf%nitrogen = plant%pool%leaf%nitrogen - plant%recycling_leaf_n
       plant%pool%fine_root%nitrogen = plant%pool%fine_root%nitrogen - plant%recycling_fine_root_n
       plant%pool%sap_wood%nitrogen = plant%pool%sap_wood%nitrogen - plant%recycling_sap_wood_n
   
       plant%pool%labile%phosphorus = plant%pool%labile%phosphorus + plant%recycling_leaf_p & 
           + plant%recycling_fine_root_p + plant%recycling_sap_wood_p
       plant%pool%leaf%phosphorus = plant%pool%leaf%phosphorus - plant%recycling_leaf_p
       plant%pool%fine_root%phosphorus = plant%pool%fine_root%phosphorus - plant%recycling_fine_root_p
       plant%pool%sap_wood%phosphorus = plant%pool%sap_wood%phosphorus - plant%recycling_sap_wood_p
   
       plant%pool%labile%nitrogen15 = plant%pool%labile%nitrogen15 + plant%recycling_leaf_n15 & 
           + plant%recycling_fine_root_n15 + plant%recycling_sap_wood_n15
       plant%pool%leaf%nitrogen15 = plant%pool%leaf%nitrogen15 - plant%recycling_leaf_n15
       plant%pool%fine_root%nitrogen15 = plant%pool%fine_root%nitrogen15 - plant%recycling_fine_root_n15
       plant%pool%sap_wood%nitrogen15 = plant%pool%sap_wood%nitrogen15 - plant%recycling_sap_wood_n15
   
       !> 1.5 Update plant biogeochemical pools with growth, litterfall and establishment
       !! 
       plant%pool = plant%pool - plant%litterfall
       plant%pool = plant%pool + plant%growth
       plant%pool%labile  = plant%pool%labile  + plant%establishment
       plant%pool%seed_bed = plant%pool%seed_bed - plant%establishment
   
       !> 1.6 Radioactive decay of C14 in living plant tissue
       !!
       plant%pool%leaf%carbon14         = (1._wp - lambda_C14 * dt ) * plant%pool%leaf%carbon14
       plant%pool%fine_root%carbon14    = (1._wp - lambda_C14 * dt ) * plant%pool%fine_root%carbon14
       plant%pool%sap_wood%carbon14     = (1._wp - lambda_C14 * dt ) * plant%pool%sap_wood%carbon14
       plant%pool%heart_wood%carbon14   = (1._wp - lambda_C14 * dt ) * plant%pool%heart_wood%carbon14
       plant%pool%coarse_root%carbon14  = (1._wp - lambda_C14 * dt ) * plant%pool%coarse_root%carbon14
       plant%pool%fruit%carbon14        = (1._wp - lambda_C14 * dt ) * plant%pool%fruit%carbon14
       plant%pool%labile%carbon14       = (1._wp - lambda_C14 * dt ) * plant%pool%labile%carbon14
       plant%pool%reserve%carbon14      = (1._wp - lambda_C14 * dt ) * plant%pool%reserve%carbon14
       plant%pool%seed_bed%carbon14     = (1._wp - lambda_C14 * dt ) * plant%pool%seed_bed%carbon14
   
       !> 2.0 Update plant diagnostics given the above fluxes
       !!
       !! 
       plant%lai = plant%pool%leaf%carbon * lctlib%sla
       plant%dens_ind = plant%dens_ind + plant%delta_dens_ind
!@todo: check for consistency with Silvia's branch
       CALL calc_allometry(plant%pool%sap_wood%carbon,plant%pool%heart_wood%carbon,plant%dens_ind,plant%pool%labile%carbon,&
            plant%pool%labile%nitrogen, plant%pool%labile%phosphorus, lctlib%wood_density, & 
            lctlib%sla,lctlib%allom_k1,lctlib%allom_k2, & 
            lctlib%k_latosa,lctlib%k_rtos,lctlib%growthform, & 
            plant%growth_cn_plant_mavg,plant%growth_cp_plant_mavg,plant%unit_npp_mavg,&
            plant%unit_uptake_n_talloc_mavg, &
            plant%unit_uptake_p_talloc_mavg, plant%allometry_npp_mavg,hlp1, plant%diameter,plant%height,hlp1,hlp1, &
            hlp1,hlp1,hlp1, hlp1)  
   
       IF(plant%pool%leaf%carbon .GT. eps8) THEN 
           plant%mean_leaf_age = & 
             (plant%mean_leaf_age + dt/one_day) * (1._wp - plant%growth%leaf%carbon/plant%pool%leaf%carbon)
         ELSE 
           plant%mean_leaf_age = 0.0_wp
       END IF

    END SELECT

    !> 3.0 Given current vegetation pools, update canopy layers 
    !!
    !!
    CALL calc_canopy_layers(lctlib,plant%pool%leaf%nitrogen,plant%lai, &
                             plant%par_sunlit_mavg_cl,&
                             plant%par_shaded_mavg_cl, & 
                             plant%fleaf_sunlit_mavg_cl,&
                             plant%t_air_mavg, &
                             plant%press_srf_mavg, &
                             plant%co2_mixing_ratio_mavg, &
                             plant%ga_mavg, &
                             plant%beta_air_mavg, &
                             plant%beta_soil_ps_mavg, &
                             plant%beta_soil_gs_mavg,& 
                             plant%t_jmax_opt_mavg, &
                             plant%lai_cl, &
                             plant%cumm_lai_cl, &
                             plant%nl_tot_cl, &
                             plant%fnl_rub_cl, &
                             plant%fnl_et_cl, &
                             plant%fnl_pepc_cl, &
                             plant%fnl_chl_cl, &
                             plant%fnl_oth_cl)

  END SUBROUTINE update_veg_pools


END MODULE mo_veg_update_veg_pools
