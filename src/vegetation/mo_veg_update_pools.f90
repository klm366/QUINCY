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
MODULE mo_veg_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_jsb_lctlib_class,    ONLY: t_lctlib_element


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_update_pools'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update vegetation pools
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_pools(tile)

    USE mo_jsb_class,           ONLY: Get_model
    USE mo_base_memory,         ONLY: t_slm_tile, t_model
    USE mo_base_constants,      ONLY: dtime, one_day
    USE mo_base_config_class,   ONLY: base_config
    USE mo_jsb_process_class,   ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_shared_constants,    ONLY: lambda_C14, rhoh2o, grav
    USE mo_veg_canopy,          ONLY: calc_canopy_layers
    USE mo_veg_growth,          ONLY: calc_diameter_from_woody_biomass,calc_height_from_diameter
    USE mo_veg_constants,       ONLY: itree

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout) :: tile         !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib           !< land-cover-type library - parameter across pft's
    REAL(wp)                              :: hlp1,hlp2        ! helper variable
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)                       ! TYPE(t_assimi_memory),      POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                          !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                         !     TYPE(t_soil_memory),     POINTER :: soil__mem
    TYPE(t_assimi_memory),   POINTER      :: assimi   ! local pointer to assimi__mem for SLM usability
    TYPE(t_veg_memory),      POINTER      :: veg      ! local pointer to veg__mem for SLM usability
    TYPE(t_rad_memory),      POINTER      :: rad      ! local pointer to rad__mem for SLM usability
    TYPE(t_soil_memory),     POINTER      :: soil     ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)

    assimi    => assimi__mem 
    veg       => veg__mem 
    rad       => rad__mem
    soil      => soil__mem 

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !> 1.0 for plant and land model, update vegetation pools with fluxes
    !! 
    !!
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('plant','land')
   
       !> 1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !! 
       !! 
       veg%npp%ptr      = assimi%gross_assimilation%ptr     - & 
           veg%growth_respiration%ptr - veg%maint_respiration%ptr - veg%n_processing_respiration%ptr
       veg%npp_c13%ptr  = assimi%gross_assimilation_C13%ptr - &
           veg%growth_respiration_c13%ptr - veg%maint_respiration_c13%ptr - veg%n_processing_respiration_c13%ptr
       veg%npp_c14%ptr  = assimi%gross_assimilation_C14%ptr - & 
           veg%growth_respiration_c14%ptr - veg%maint_respiration_c14%ptr - veg%n_processing_respiration_c14%ptr
       veg%growth%labile%carbon%ptr   = veg%growth%labile%carbon%ptr    + veg%npp%ptr     * dtime / 1000000.0_wp
       veg%growth%labile%carbon13%ptr = veg%growth%labile%carbon13%ptr  + veg%npp_c13%ptr * dtime / 1000000.0_wp
       veg%growth%labile%carbon14%ptr = veg%growth%labile%carbon14%ptr  + veg%npp_c14%ptr * dtime / 1000000.0_wp
   
       !> 1.2 Add nutrients from uptake calculations
       !! 
       !! 
       veg%growth%labile%nitrogen%ptr   = veg%growth%labile%nitrogen%ptr + & 
          ( veg%uptake_nh4%ptr + veg%uptake_no3%ptr + veg%n_fixation%ptr) * dtime / 1000000.0_wp
       veg%growth%labile%nitrogen15%ptr = veg%growth%labile%nitrogen15%ptr + & 
          ( veg%uptake_nh4_n15%ptr + veg%uptake_no3_n15%ptr + veg%n_fixation_n15%ptr) * dtime / 1000000.0_wp
       veg%growth%labile%phosphorus%ptr = veg%growth%labile%phosphorus%ptr + & 
          veg%uptake_p%ptr * dtime / 1000000.0_wp
   
       !> 1.3 update labile pool and reserve 
       !! 
       veg%pool%labile  = veg%pool%labile  + veg%reserve_use
       veg%pool%reserve = veg%pool%reserve - veg%reserve_use
   
       !> 1.4 update fluxes from continuous nutrient turnover and resorption 
       !!
       veg%pool%labile%nitrogen%ptr       = veg%pool%labile%nitrogen%ptr      + veg%recycling_leaf_n%ptr + & 
                                                    veg%recycling_fine_root_n%ptr + veg%recycling_heart_wood_n%ptr
       veg%pool%leaf%nitrogen%ptr         = veg%pool%leaf%nitrogen%ptr        - veg%recycling_leaf_n%ptr
       veg%pool%fine_root%nitrogen%ptr    = veg%pool%fine_root%nitrogen%ptr   - veg%recycling_fine_root_n%ptr
       veg%pool%heart_wood%nitrogen%ptr   = veg%pool%heart_wood%nitrogen%ptr  - veg%recycling_heart_wood_n%ptr
   
       veg%pool%labile%phosphorus%ptr     = veg%pool%labile%phosphorus%ptr    + veg%recycling_leaf_p%ptr + & 
                                                    veg%recycling_fine_root_p%ptr + veg%recycling_heart_wood_p%ptr
       veg%pool%leaf%phosphorus%ptr       = veg%pool%leaf%phosphorus%ptr      - veg%recycling_leaf_p%ptr
       veg%pool%fine_root%phosphorus%ptr  = veg%pool%fine_root%phosphorus%ptr - veg%recycling_fine_root_p%ptr
       veg%pool%heart_wood%phosphorus%ptr = veg%pool%heart_wood%phosphorus%ptr - veg%recycling_heart_wood_p%ptr
   
       veg%pool%labile%nitrogen15%ptr     = veg%pool%labile%nitrogen15%ptr    + veg%recycling_leaf_n15%ptr + & 
                                                    veg%recycling_fine_root_n15%ptr + veg%recycling_heart_wood_n15%ptr
       veg%pool%leaf%nitrogen15%ptr       = veg%pool%leaf%nitrogen15%ptr      - veg%recycling_leaf_n15%ptr
       veg%pool%fine_root%nitrogen15%ptr  = veg%pool%fine_root%nitrogen15%ptr - veg%recycling_fine_root_n15%ptr
       veg%pool%heart_wood%nitrogen15%ptr = veg%pool%heart_wood%nitrogen15%ptr - veg%recycling_heart_wood_n15%ptr
   
       !> 1.5 Update plant biogeochemical pools with growth, litterfall and establishment
       !! 
       veg%pool = veg%pool - veg%litterfall
       veg%pool = veg%pool + veg%growth
       veg%pool%labile    = veg%pool%labile   + veg%establishment
       veg%pool%seed_bed  = veg%pool%seed_bed - veg%establishment
   
       !> 1.6 Radioactive decay of C14 in living plant tissue
       !!
       veg%pool%leaf%carbon14%ptr         = veg%pool%leaf%carbon14%ptr        * (1._wp - lambda_C14 * dtime )
       veg%pool%fine_root%carbon14%ptr    = veg%pool%fine_root%carbon14%ptr   * (1._wp - lambda_C14 * dtime )
       veg%pool%sap_wood%carbon14%ptr     = veg%pool%sap_wood%carbon14%ptr    * (1._wp - lambda_C14 * dtime )
       veg%pool%heart_wood%carbon14%ptr   = veg%pool%heart_wood%carbon14%ptr  * (1._wp - lambda_C14 * dtime )
       veg%pool%coarse_root%carbon14%ptr  = veg%pool%coarse_root%carbon14%ptr * (1._wp - lambda_C14 * dtime )
       veg%pool%fruit%carbon14%ptr        = veg%pool%fruit%carbon14%ptr       * (1._wp - lambda_C14 * dtime )
       veg%pool%labile%carbon14%ptr       = veg%pool%labile%carbon14%ptr      * (1._wp - lambda_C14 * dtime )
       veg%pool%reserve%carbon14%ptr      = veg%pool%reserve%carbon14%ptr     * (1._wp - lambda_C14 * dtime )
       veg%pool%seed_bed%carbon14%ptr     = veg%pool%seed_bed%carbon14%ptr    * (1._wp - lambda_C14 * dtime )
   
       !> 2.0 Update plant diagnostics given the above fluxes
       !!
       !! 
       veg%lai%ptr      = veg%pool%leaf%carbon%ptr * lctlib%sla
       veg%dens_ind%ptr = veg%dens_ind%ptr + veg%delta_dens_ind%ptr
       IF(lctlib%growthform == itree)THEN
          veg%diameter%ptr  = calc_diameter_from_woody_biomass(lctlib                         , &
                                                               veg%pool%sap_wood%carbon%ptr   , &
                                                               veg%pool%heart_wood%carbon%ptr , &
                                                               veg%dens_ind%ptr)
          veg%height%ptr    = calc_height_from_diameter(lctlib, veg%diameter%ptr) 
       ENDIF 

       IF(veg%pool%leaf%carbon%ptr > eps8) THEN 
           veg%mean_leaf_age%ptr = & 
             (veg%mean_leaf_age%ptr + dtime/one_day) * (1._wp - veg%growth%leaf%carbon%ptr / veg%pool%leaf%carbon%ptr)
         ELSE 
           veg%mean_leaf_age%ptr = 0.0_wp
       END IF

       veg%dphi%ptr = soil%w_soil_root_pot%ptr(1,1) - lctlib%phi_leaf_min - grav * rhoh2o * veg%height%ptr * 1.e-6_wp

    END SELECT

    !> 3.0 Given current vegetation pools, update canopy layers 
    !!
    !!
    CALL calc_canopy_layers( lctlib, &
                             veg%pool%leaf%nitrogen%ptr, &
                             veg%lai%ptr, &
                             rad%ppfd_sunlit_tfrac_mavg_cl%ptr(:,1),&
                             rad%ppfd_shaded_tfrac_mavg_cl%ptr(:,1), & 
                             veg%fleaf_sunlit_tfrac_mavg_cl%ptr(:,1),&
                             veg%t_air_tfrac_mavg%ptr, &
                             veg%t_air_tacclim_mavg%ptr, &
                             veg%press_srf_tfrac_mavg%ptr, &
                             veg%co2_mixing_ratio_tfrac_mavg%ptr, &
                             veg%ga_tfrac_mavg%ptr, &
                             assimi%beta_air_tfrac_mavg%ptr, &
                             assimi%beta_soil_ps_tfrac_mavg%ptr, &
                             veg%beta_sinklim_ps_tfrac_mavg%ptr, &
                             assimi%beta_soil_gs_tfrac_mavg%ptr, & 
                             veg%t_jmax_opt_mavg%ptr, &
                             veg%lai_cl%ptr(:,1), &
                             veg%cumm_lai_cl%ptr(:,1), &
                             veg%leaf_nitrogen_cl%ptr(:,1), &
                             veg%fn_rub_cl%ptr(:,1), &
                             veg%fn_et_cl%ptr(:,1), &
                             veg%fn_pepc_cl%ptr(:,1), &
                             veg%fn_chl_cl%ptr(:,1), &
                             veg%fn_oth_cl%ptr(:,1))

    !> 4.0 analysis of mass conservation: sum elements across pools into veg\%total_biomass
    !!
    !!
    veg%total_biomass =   veg%pool%leaf        + &     !< elements of leaf
                          veg%pool%fine_root   + &     !< elements of fine_root
                          veg%pool%coarse_root + &     !< elements of coarse_root
                          veg%pool%sap_wood    + &     !< elements of sap_wood
                          veg%pool%heart_wood  + &     !< elements of heart_wood
                          veg%pool%labile      + &     !< elements of labile
                          veg%pool%reserve     + &     !< elements of reserve
                          veg%pool%fruit       + &     !< elements of fruit
                          veg%pool%seed_bed            !< elements of seed_bed

  END SUBROUTINE update_veg_pools

END MODULE mo_veg_update_pools
