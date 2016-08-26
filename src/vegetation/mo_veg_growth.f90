!>
!! @brief Contains the routines for the plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_growth

  USE mo_kind
  USE mo_exception,             ONLY: message_text, finish
  USE mo_base_constants
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_veg_memory
  USE mo_veg_constants
  USE mo_veg_config_class,      ONLY: veg_config


  PRIVATE

  PUBLIC :: calc_growth, calc_allometry

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> calculate plant growth (flux\%growth)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_growth (tile, lctlib)

    USE mo_base_memory,             ONLY: t_slm_tile
    USE mo_veg_memory
    USE mo_shared_memory,           ONLY: t_element
    USE mo_jsb_process_class,       ONLY: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, HD, RAD, ASSIMI, PHENO
    USE mo_atm2land_memory_class,   ONLY: t_atm2land_memory

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile      !< one tile
    TYPE(t_lctlib_element),INTENT(in)       :: lctlib    !< land-cover-type library - parameters for one pft
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER               :: plant     ! plant variables per tile
    TYPE(t_slm_veg_pool) :: falloc                                 ! fractional allocation to pools
    REAL(wp) :: growth_req_n,growth_req_p                          ! N and P reqiurement of growth (mol N / P per mol C)
    REAL(wp) :: target_lai                                         ! LAI implied from sap wood area 
    REAL(wp) :: tau_labile_pot                                     ! labile pool turnover rate assuming only temperature and moisture
                                                                   ! constraints (no stoichiometry)
    ! target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    REAL(wp) :: target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit
    REAL(wp) :: diameter,height,dheight,leaf2sapwood_mass_ratio,leaf2root_mass_ratio
    REAl(wp) :: k1_opt, k2_opt
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND)
    REAL(wp), POINTER                :: t_air
    SELECT TYPE (mm => tile%mem(ATM2LAND)%p)            ! dsl4jsb_Get_memory(ATM2LAND)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    t_air   => atm2land__mem%t_air%ptr

    ! bi-directional pointer
    plant => tile%plant


    !> 1.0 calculate potential growth rate relative to labile pool size
    !! 
    !! 
    tau_labile_pot = meristem_activity(t_air,tile%soil%theta,plant%growing_season)
  
    !> 2.0 update target C:N:P for each tissue type
    !! 
    !! 
    CALL update_stoichiometry(plant,lctlib%growthform, & 
                              target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, & 
                              target_cn_fruit,target_np_fine_root,target_np_coarse_root, & 
                              target_np_sap_wood,target_np_fruit)

    !> 3.0 calculate current allometry 
    !! 
    !! 
    CALL calc_allometry(plant%pool%sap_wood%carbon,plant%pool%heart_wood%carbon,plant%dens_ind,plant%moving_average%allometry_labile%carbon,&
         plant%moving_average%allometry_labile%nitrogen, plant%moving_average%allometry_labile%phosphorus,&
         lctlib%wood_density, lctlib%sla,lctlib%allom_k1,lctlib%allom_k2, & 
         lctlib%k_latosa,lctlib%k_rtos,lctlib%growthform, & 
         plant%moving_average%growth_cn_plant,plant%moving_average%growth_cp_plant, plant%moving_average%unit_npp,&
         plant%moving_average%unit_uptake_n,plant%moving_average%unit_uptake_p, &
         plant%moving_average%allometry_npp,tau_labile_pot,diameter,height,dheight,leaf2sapwood_mass_ratio, & 
         leaf2root_mass_ratio,target_lai,k1_opt, k2_opt)  

    !> 4.0 calculate net flux between reserve and labile pool 
    !! 
    !! 
    plant%flux%reserve_use = calc_reserve_useage(tau_labile_pot,plant%growing_season, &
                              plant%moving_average%gpp_week, &
                              plant%moving_average%maint_respiration_week, &
                              plant%moving_average%growth_req_n_week, &
                              plant%moving_average%growth_req_p_week, &
                              lctlib%sla, &
                              lctlib%tau_leaf,lctlib%tau_fine_root, & 
                              plant%target_cn_leaf, target_cn_fine_root, & 
                              plant%target_np_leaf, target_np_fine_root, plant%moving_average%leaf2root_mass_ratio, &
                              plant%lai,target_lai, &
                              plant%pool%labile,plant%pool%reserve)
    !CALL calc_nutrient_recycling(plant, target_cn_fine_root, target_np_fine_root)

    !> 5.0 calculate allocation coefficients for structural plant pools (leaves, wood, fine and coarse root, fruits)
    !! 
    !! 
    CALL calc_allocation_fraction(tau_labile_pot,plant%pool%leaf%carbon,plant%pool%sap_wood%carbon, &
         plant%pool%coarse_root%carbon,plant%pool%fine_root%carbon,plant%pool%labile%carbon,&
         lctlib%growthform, &
         height,dheight,leaf2sapwood_mass_ratio,leaf2root_mass_ratio, & 
         lctlib%k_crtos,k1_opt,k2_opt,lctlib%k_max_falloc_fruit,& 
         plant%target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
         target_cn_fruit,plant%target_np_leaf,target_np_fine_root,target_np_coarse_root, &
         target_np_sap_wood,target_np_fruit,plant%flux%reserve_use%carbon, & 
         falloc,plant%growth_req_n,plant%growth_req_p)

    !> 6.0 partition labile carbon into growth, growth_respiration, and maintenance 
    !! 
    CALL calc_partitioning(tau_labile_pot, falloc, plant, t_air)

  END SUBROUTINE calc_growth 


  !-----------------------------------------------------------------------------------------------------
  !> Calculate partitioning of labile pool outflow into maintenance, growth respiration and growth 
  !! Follows the concept of Amthor: Growth = Yg * (P - Rm), where\n
  !!   Yg = growth efficiency of biomass production (1-fresp_growth)\n
  !!   P  = current turnover of labile pool (Amthor's production)\n
  !!   Rm = maintenance respiration
  !!
  !! In the current implementation, maintenance is given preference over growth, but 
  !! nevertheless, maintenance is limited to the available C   
  !!
  !! Input: plant
  !! 
  !! Output: growth of tissue types, maintenance respiration and growth respiration 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_partitioning(tau_labile_pot, falloc, plant, t_air) 

    USE mo_shared_memory,    ONLY: t_element
    USE mo_veg_respiration,  ONLY: calc_maintenance_respiration

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in)                     :: tau_labile_pot       !< potential labile pool turnover (1/day) 
    REAL(wp),INTENT(in)                     :: t_air 
    TYPE(t_slm_veg_pool), INTENT(in)        :: falloc               !< allocation coefficients 
    TYPE(t_slm_plant),INTENT(inout),TARGET  :: plant                !< plant structure
    ! 0.2 Local
    REAL(wp)                  :: tau_labile                ! actual turnover time of labile pool (1/day)
    TYPE(t_element),POINTER   :: growth
    TYPE(t_element),POINTER   :: labile_pool

    labile_pool => plant%pool%labile
    growth => plant%flux%growth%labile

    !> 1.0 Determine C use for maintenance 
    !! 
    !! based on tissue N and Temperature, constrained by default turnover of the labile pool
    !! to avoid that maitenance respiration consumes all C in times of low GPP 

    ! potential rate of maintenance respiration (µmol C / m2 / s)
    plant%flux%maint_respiration_pot = & 
          calc_maintenance_respiration(t_air,t_air,plant%pool%sap_wood%nitrogen, & 
                                       plant%pool%fine_root%nitrogen,plant%pool%coarse_root%nitrogen) + & 
          plant%flux%maint_respiration_leaf
    ! actual rate of maintentance respiration (µmol C / m2 / s), constrained to 
    ! support leaf respiration, but not more than the current labile pool turnover
    plant%flux%maint_respiration = MAX(plant%flux%maint_respiration_leaf, & 
                                   MIN(tau_labile_default / one_day * labile_pool%carbon * 1.e6, &
                                   plant%flux%maint_respiration_pot))

    !> 2.0 Determine C, N, P use for growth
    !!
    !> 2.1 potential C growth rate (mol C / m2 / timestep) 
    !! is determined from the labile turnover minus the carbon used for maintenance, i
    !! accounting for the fraction of growth respiration required to achieve this growth, to ensure that
    !! growth does not empty the labile pool
    growth%carbon = (labile_pool%carbon * tau_labile_pot * dt/one_day -  plant%flux%maint_respiration * dt / 1.e6 ) * & 
                    1._wp /  (1._wp + fresp_growth) 
    IF(growth%carbon < eps8) growth%carbon = 0.0_wp

    !> 2.2 Elemental limitation of growth by N & P (mol [N|P] / m2 / timestep)
    !! is determined by the minimum turnover of the labile N and P pools and the stoichiometric growth requirement
    !! given the potential growth rate of C
    growth%nitrogen   = MIN(plant%growth_req_n * growth%carbon, & 
                            labile_pool%nitrogen * tau_labile_pot * dt/one_day) 
    growth%phosphorus = MIN(plant%growth_req_p * growth%nitrogen, & 
                            labile_pool%phosphorus * tau_labile_pot * dt/one_day) 
    ! adjust N growth rate to P growth rate if P is limiting 
    IF(growth%nitrogen > growth%phosphorus / plant%growth_req_p) & 
        growth%nitrogen = growth%phosphorus / plant%growth_req_p
    ! adjust C growth rate to N growth rate if N is limiting 
    IF(growth%carbon > growth%nitrogen / plant%growth_req_n) & 
        growth%carbon = growth%nitrogen / plant%growth_req_n
 
    ! Diagnose growth respiration (µmol C / m2 / s)
    plant%flux%growth_respiration = fresp_growth * growth%carbon * 1.e6 / dt

    !> 3.0 Calculate isotopic signature of respiration fluxes, and growth
    !!
    IF(labile_pool%carbon > eps8) THEN
       ! growth (in mol C / m2 / timestep)
       growth%carbon13 = growth%carbon * labile_pool%carbon13/labile_pool%carbon
       growth%carbon14 = growth%carbon * labile_pool%carbon14/labile_pool%carbon
       ! respiration (in µmol C / m2 / s)
       plant%flux%maint_respiration_c13 = plant%flux%maint_respiration * labile_pool%carbon13/labile_pool%carbon
       plant%flux%maint_respiration_c14 = plant%flux%maint_respiration * labile_pool%carbon14/labile_pool%carbon
       plant%flux%growth_respiration_c13 = plant%flux%growth_respiration * labile_pool%carbon13/labile_pool%carbon
       plant%flux%growth_respiration_c14 = plant%flux%growth_respiration * labile_pool%carbon14/labile_pool%carbon
    ENDIF
    IF(labile_pool%nitrogen > eps8) THEN
       ! growth (in mol  / m2 / timestep)
       growth%nitrogen15 = growth%nitrogen * labile_pool%nitrogen15/labile_pool%nitrogen
    ENDIF

    !> 4.0 partition growth to tissue types
    !! 
    plant%flux%growth%leaf        = falloc%leaf        * growth
    plant%flux%growth%fine_root   = falloc%fine_root   * growth 
    plant%flux%growth%coarse_root = falloc%coarse_root * growth
    plant%flux%growth%sap_wood    = falloc%sap_wood    * growth
    plant%flux%growth%fruit       = falloc%fruit       * growth

    !> 4.1 invert flux of labile to be out of the labile pool
    !! 
    plant%flux%growth%labile%carbon     = -1.0_wp * growth%carbon 
    plant%flux%growth%labile%nitrogen   = -1.0_wp * growth%nitrogen
    plant%flux%growth%labile%phosphorus = -1.0_wp * growth%phosphorus
    plant%flux%growth%labile%carbon13   = -1.0_wp * growth%carbon13
    plant%flux%growth%labile%carbon14   = -1.0_wp * growth%carbon14
    plant%flux%growth%labile%nitrogen15 = -1.0_wp * growth%nitrogen15

  END SUBROUTINE calc_partitioning 


  !-----------------------------------------------------------------------------------------------------
  !> Calculates allocation fractions for each tissue type and element
  !! 
  !! Input: target stoichiometry of the tissue types, plant allometry, and current pool sizes
  !! 
  !! Output:\n 
  !!         (1) Allocation coefficients for the distribution of biomass for growth;\n
  !!         (2) total N and P required to sustain a unit of growth (in mol); \n
  !!
  !! @todo: calc_allocation_fraction: check calculation of b_inc. Currently a fixed value for stability,
  !!        Needs to be adjusted to better match prescriebed allometry and to account for falloc_fruit
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_allocation_fraction(tau_labile_pot, &
                                      leaf_carbon, &
                                      sap_wood_carbon, &
                                      coarse_root_carbon, &
                                      root_carbon, &
                                      labile_carbon, &
                                      growthform, &
                                      height, &
                                      dheight, &
                                      leaf2sapwood_mass_ratio, &
                                      leaf2root_mass_ratio, &
                                      coarse_root2sapwood_ratio, & 
                                      k1_opt, &
                                      k2_opt, &
                                      k_max_falloc_fruit, &
                                      target_cn_leaf, &
                                      target_cn_fine_root, &
                                      target_cn_coarse_root, & 
                                      target_cn_sap_wood, &
                                      target_cn_fruit, &
                                      target_np_leaf, &
                                      target_np_fine_root, & 
                                      target_np_coarse_root, &
                                      target_np_sap_wood, &
                                      target_np_fruit, &
                                      reserve_use_carbon, &
                                      falloc, &
                                      growth_req_n, & 
                                      growth_req_p)

    USE mo_veg_memory,     ONLY: t_slm_veg_pool
    USE mo_veg_constants,  ONLY: itree, igrass, lambda_fruit_alloc, weibull_k_fruit_alloc, min_falloc_fruit, offset_fruit_alloc

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: tau_labile_pot             !< potential labile pool turnover rate (1/day) without stoichiometric limit
    REAL(wp), INTENT(in) :: leaf_carbon                !< carbon in leaves (mol/individual)
    REAL(wp), INTENT(in) :: sap_wood_carbon            !< carbon in sap wood (mol/individual)
    REAL(wp), INTENT(in) :: coarse_root_carbon         !< carbon in heart wood (mol/individual)
    REAL(wp), INTENT(in) :: root_carbon                !< carbon in roots (mol/individual) 
    REAL(wp), INTENT(in) :: labile_carbon              !< carbon in labile pool (mol/individual) 
    INTEGER, INTENT(in)  :: growthform                 !< plant growth form 
    REAL(wp), INTENT(in) :: height                     !< plant height  
    REAL(wp), INTENT(in) :: dheight                    !< derivative of height per unit sapwood mass change 
    REAL(wp), INTENT(in) :: leaf2sapwood_mass_ratio    !< leaf to shoot mass ratio
    REAL(wp), INTENT(in) :: leaf2root_mass_ratio       !< leaf to root mass ratio
    REAL(wp), INTENT(in) :: coarse_root2sapwood_ratio  !< sapwood to coarse root mass ratio 
    REAL(wp), INTENT(in) :: k1_opt                     !< coefficient for optimal leaf:root ratio
    REAL(wp), INTENT(in) :: k2_opt                     !< coefficient for optimal leaf:root ratio
    REAL(wp), INTENT(in) :: k_max_falloc_fruit         !< maximum fraction allocation to fruits 
    REAL(wp), INTENT(in) :: target_cn_leaf             !< target C:N ratio for tissue type 
    REAL(wp), INTENT(in) :: target_cn_fine_root        !< target C:N ratio for tissue type 
    REAL(wp), INTENT(in) :: target_cn_coarse_root      !< target C:N ratio for tissue type 
    REAL(wp), INTENT(in) :: target_cn_sap_wood         !< target C:N ratio for tissue type 
    REAL(wp), INTENT(in) :: target_cn_fruit            !< target C:N ratio for tissue type 
    REAL(wp), INTENT(in) :: target_np_leaf             !< target N:P ratio for tissue type 
    REAL(wp), INTENT(in) :: target_np_fine_root        !< target N:P ratio for tissue type 
    REAL(wp), INTENT(in) :: target_np_coarse_root      !< target N:P ratio for tissue type 
    REAL(wp), INTENT(in) :: target_np_sap_wood         !< target N:P ratio for tissue type 
    REAL(wp), INTENT(in) :: target_np_fruit            !< target N:P ratio for tissue type 
    REAL(wp), INTENT(in) :: reserve_use_carbon         !< current reserve usage (mol C / m2 / timestep)

    TYPE(t_slm_veg_pool), INTENT(out) :: falloc        !< allocation coefficients 
    REAL(wp), INTENT(out)             :: growth_req_n  !< N requirement for unit growth in C and N
    REAL(wp), INTENT(out)             :: growth_req_p  !< P requirement for unit growth in C and N
    ! 0.2 Local
    REAL(wp) :: b_inc,                   &             ! total and compartment-wise carbon increment (mol/dt)
                all_inc,                 &
                leaf_carbon_inc,         &
                root_carbon_inc,         &
                sap_wood_carbon_inc,     &
                coarse_root_carbon_inc 
    REAL(wp) :: leaf2wood_mass_ratio,    &
                wood_carbon
    REAL(wp) :: hlp1

    !> 1.0 Allocation factors for C (this should be calculated based on various principles of plant growth and allometry)
    !! 
    !! 

    ! SZ: uncomment this to revert to fixed constant allocation
    ! falloc%leaf%carbon = 0.3_wp
    ! falloc%fine_root%carbon = 0.3_wp
    ! falloc%coarse_root%carbon = 0.1_wp
    ! falloc%sap_wood%carbon = 0.2_wp
    ! falloc%fruit%carbon = 0.1_wp 

    ! 
    ! seed for allocation coefficients is the potential growth rate per time step.
    ! given by the average turnover of the labile pool 
    b_inc = tau_labile_pot * labile_carbon * dt/one_day

    !> 1.1 Calculate allocation fraction to reproduction
    !!  Fruit allocation depends on whether or not the plant is able to build reserves
    !!  as soon as the carbon balance exceeds a PFT specific offset, the maximal fruit growth rate is achieved
    hlp1                = MAX(0.0_wp, reserve_use_carbon / dt * 1000000.0_wp + offset_fruit_alloc)
    falloc%fruit%carbon = min_falloc_fruit + ( k_max_falloc_fruit - min_falloc_fruit ) * & 
                          exp(-(lambda_fruit_alloc * hlp1)**weibull_k_fruit_alloc)

    !> 1.2 Calculate allocation fractions to functional organs
    !!
    !! 
    SELECT CASE(growthform)
      !> 1.2 trees
      !!
      CASE (itree)
       ! carbon increment
       IF(b_inc > eps8) THEN

          !! needed for stability - to be checked BLARPP! 
          !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> check stability issues with b_inc
          b_inc = 0.2_wp * dt / one_day 

          leaf2wood_mass_ratio = leaf2sapwood_mass_ratio / (1._wp + coarse_root2sapwood_ratio)
          wood_carbon = sap_wood_carbon + coarse_root_carbon

          !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> 
          !! rewrite "IF (veg_config\%flag_optimal_allocation)" such that instead of kopt1 and 2 leaf2root_mass ratio can be given
          IF (veg_config%flag_optimal_allocation) THEN

             !> 1.2a Calculate allocation fractions assuming optimal relationship between leaf and root allocation
             !! 
             leaf_carbon_inc = MAX((b_inc - k1_opt + root_carbon + sap_wood_carbon - dheight/leaf2sapwood_mass_ratio * &
                  (b_inc-k1_opt + root_carbon) - height/leaf2sapwood_mass_ratio*leaf_carbon) / &
                  (1.0_wp + k2_opt + dheight/leaf2sapwood_mass_ratio * &
                  (b_inc - (1.0_wp + k2_opt) * leaf_carbon - k1_opt + root_carbon) + &
                  height/leaf2sapwood_mass_ratio), 0.0_wp)
             sap_wood_carbon_inc = MAX((( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * height - sap_wood_carbon) / & 
                  (1._wp - ( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * dheight ), 0.0_wp)        
             IF (sap_wood_carbon_inc > 0.0_wp) THEN
                coarse_root_carbon_inc = MAX((sap_wood_carbon + sap_wood_carbon_inc) * coarse_root2sapwood_ratio & 
                     - coarse_root_carbon, 0.0_wp)
                root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
             ELSE
                leaf_carbon_inc = MAX((b_inc - k1_opt + root_carbon)/(1+k2_opt), 0.0_wp)
                root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
                sap_wood_carbon_inc = 0.0_wp
                coarse_root_carbon_inc = 0.0_wp
             ENDIF
          ELSE 

 	         !> 1.2b solving the allometric rations between leaf and root, leaf and sapwood and sapwood and coarse roots 
             !! simulataneously with the height-to-diameter relationship 
             leaf_carbon_inc = MAX(( b_inc + leaf_carbon * dheight / leaf2wood_mass_ratio * & 
                            ( leaf_carbon / leaf2root_mass_ratio - b_inc - root_carbon ) + wood_carbon + root_carbon - & 
                            leaf_carbon * (height / leaf2wood_mass_ratio + 1.0_wp / leaf2root_mass_ratio )) / & 
                            (1.0_wp + ( b_inc - (1.0_wp + 2.0_wp/leaf2root_mass_ratio) * leaf_carbon + root_carbon ) & 
                            * dheight / leaf2wood_mass_ratio + height / leaf2wood_mass_ratio + 1._wp / leaf2root_mass_ratio ), & 
                            0.0_wp) 
   
             sap_wood_carbon_inc = MAX((( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * height - sap_wood_carbon) / & 
                                (1._wp - ( leaf_carbon + leaf_carbon_inc ) / leaf2sapwood_mass_ratio * dheight ), 0.0_wp)          
             IF(sap_wood_carbon_inc > 0.0_wp) THEN
                coarse_root_carbon_inc = MAX((sap_wood_carbon + sap_wood_carbon_inc) * coarse_root2sapwood_ratio & 
                                           - coarse_root_carbon, 0.0_wp)
                root_carbon_inc = MAX(( leaf_carbon + leaf_carbon_inc ) / leaf2root_mass_ratio - root_carbon, 0.0_wp)
             ELSE ! special case of no wood growth (phenology etc.): redistribute available carbon 
                leaf_carbon_inc = MAX(( b_inc - leaf_carbon / leaf2root_mass_ratio + root_carbon ) / & 
                                    (1.0_wp + 1._wp/ leaf2root_mass_ratio),0.0_wp)
                root_carbon_inc = MAX(0.0_wp,MIN(b_inc - leaf_carbon_inc, & 
                                       ( leaf_carbon_inc + leaf_carbon ) / leaf2root_mass_ratio - root_carbon))
                sap_wood_carbon_inc = 0.0_wp
                coarse_root_carbon_inc = 0.0_wp
             ENDIF
          ENDIF
        ENDIF
      !> 1.2 grass
      !!
      CASE(igrass)
         IF (veg_config%flag_optimal_allocation) THEN
            leaf_carbon_inc = MAX((b_inc - k1_opt - leaf_carbon/leaf2sapwood_mass_ratio + root_carbon + &
                 sap_wood_carbon) / (1 + k2_opt + 1/leaf2sapwood_mass_ratio), 0.0_wp)
            root_carbon_inc = MAX(k1_opt + k2_opt * leaf_carbon_inc - root_carbon, 0.0_wp)
            sap_wood_carbon_inc = ( leaf_carbon_inc + leaf_carbon ) / leaf2sapwood_mass_ratio - sap_wood_carbon
            coarse_root_carbon_inc = 0.0_wp
         ELSE
            leaf_carbon_inc = MAX( ( b_inc - leaf_carbon / leaf2root_mass_ratio + sap_wood_carbon + root_carbon ) / & 
                 (1.0_wp + 1._wp / leaf2sapwood_mass_ratio + 1._wp / leaf2root_mass_ratio ) , 0.0_wp) 
            root_carbon_inc = MAX(MIN( b_inc, ( leaf_carbon + leaf_carbon_inc ) / leaf2root_mass_ratio - root_carbon), 0.0_wp)
            sap_wood_carbon_inc = ( leaf_carbon_inc + leaf_carbon ) / leaf2sapwood_mass_ratio - sap_wood_carbon
            coarse_root_carbon_inc = 0.0_wp
         ENDIF
      !> 1.2 default - write and finish
      !!
      CASE DEFAULT
       !! Write statement
       !! 
       WRITE(message_text,'(a,i4)') 'invalid growthform type ',growthform
       CALL finish("mo_veg_growth",message_text)
    END SELECT

    !> 1.3 Convert biomass increments to allocation fractions
    !! taking account of the previously calculated fraction of growth that goes to fruits
    !!
    IF(b_inc > eps8) THEN
       all_inc = leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc+coarse_root_carbon_inc
       IF(all_inc > eps8) THEN 
          falloc%leaf%carbon        = (1._wp - falloc%fruit%carbon ) * leaf_carbon_inc / all_inc
          falloc%fine_root%carbon   = (1._wp - falloc%fruit%carbon ) * root_carbon_inc / all_inc
          falloc%coarse_root%carbon = (1._wp - falloc%fruit%carbon ) * coarse_root_carbon_inc / all_inc
          falloc%sap_wood%carbon    = (1._wp - falloc%fruit%carbon ) * sap_wood_carbon_inc/ all_inc
       ENDIF
    ENDIF

    IF(b_inc > eps8) THEN

       !> 2.0 Allocation factors for N, derived from the target C:N ratios and the C partitioning coefficients
       !!  
 
       !> 2.1 N required to sustain a unit of C growth
       !!
       growth_req_n = falloc%leaf%carbon / target_cn_leaf + & 
                      falloc%fine_root%carbon / target_cn_fine_root + &
                      falloc%coarse_root%carbon / target_cn_coarse_root + &
                      falloc%sap_wood%carbon / target_cn_sap_wood + & 
                      falloc%fruit%carbon / target_cn_fruit 
   
       !> 2.2 Adjustment of N allocation factors
       !!
       IF(growth_req_n > eps8)THEN
          falloc%leaf%nitrogen        = falloc%leaf%carbon / target_cn_leaf / growth_req_n
          falloc%fine_root%nitrogen   = falloc%fine_root%carbon / target_cn_fine_root / growth_req_n 
          falloc%coarse_root%nitrogen = falloc%coarse_root%carbon / target_cn_coarse_root / growth_req_n
          falloc%sap_wood%nitrogen    = falloc%sap_wood%carbon / target_cn_sap_wood / growth_req_n
          falloc%fruit%nitrogen       = falloc%fruit%carbon / target_cn_fruit / growth_req_n
       ENDIF

       !> 3.0 Allocation factors for P, derived from the target N:P ratios and the N partitioning coefficients
       !! 
       !! 
       ! P required to sustain a unit of N growth
       growth_req_p = falloc%leaf%nitrogen / target_np_leaf + & 
                      falloc%fine_root%nitrogen / target_np_fine_root + &
                      falloc%coarse_root%nitrogen / target_np_coarse_root + &
                      falloc%sap_wood%nitrogen / target_np_sap_wood + & 
                      falloc%fruit%nitrogen / target_np_fruit 
   
       ! Adjustment of P allocation factors
       IF(growth_req_p > eps8)THEN
          falloc%leaf%phosphorus = falloc%leaf%nitrogen / target_np_leaf / growth_req_p
          falloc%fine_root%phosphorus = falloc%fine_root%nitrogen / target_np_fine_root / growth_req_p 
          falloc%coarse_root%phosphorus = falloc%coarse_root%nitrogen / target_np_coarse_root / growth_req_p
          falloc%sap_wood%phosphorus = falloc%sap_wood%nitrogen / target_np_sap_wood / growth_req_p
          falloc%fruit%phosphorus = falloc%fruit%nitrogen / target_np_fruit / growth_req_p
       ENDIF
   
       !> 4.0 Allocation factors for isotopes follow the main isotope 
       !! 
       !! 
       falloc%leaf%carbon13 = falloc%leaf%carbon
       falloc%fine_root%carbon13 = falloc%fine_root%carbon
       falloc%coarse_root%carbon13 = falloc%coarse_root%carbon
       falloc%sap_wood%carbon13 = falloc%sap_wood%carbon
       falloc%fruit%carbon13 = falloc%fruit%carbon
   
       falloc%leaf%carbon14 = falloc%leaf%carbon
       falloc%fine_root%carbon14 = falloc%fine_root%carbon
       falloc%coarse_root%carbon14 = falloc%coarse_root%carbon
       falloc%sap_wood%carbon14 = falloc%sap_wood%carbon
       falloc%fruit%carbon14 = falloc%fruit%carbon
   
       falloc%leaf%nitrogen15 = falloc%leaf%nitrogen
       falloc%fine_root%nitrogen15 = falloc%fine_root%nitrogen
       falloc%coarse_root%nitrogen15 = falloc%coarse_root%nitrogen
       falloc%sap_wood%nitrogen15 = falloc%sap_wood%nitrogen
       falloc%fruit%nitrogen15 = falloc%fruit%nitrogen
    ENDIF

  END SUBROUTINE calc_allocation_fraction

  ! ------------------------------------------------------------------------------------------
  !> Routine to calculate the target stoichiometry of the newly growth tissue
  !! 
  !! Includes optimal CN ratio calculation \n
  !! Output: target stoichiometries
  ! --------------------------------------------------------------------------------------------
  SUBROUTINE update_stoichiometry(plant,growthform,target_cn_fine_root,target_cn_coarse_root, & 
                              target_cn_sap_wood,target_cn_fruit,target_np_fine_root, & 
                              target_np_coarse_root,target_np_sap_wood,target_np_fruit)

    USE mo_assimi_constants,     ONLY: delta_n_leaf

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)              :: growthform             !< plant growth form
    TYPE(t_slm_plant), INTENT(inout) :: plant                  !< plant structure
    REAL(wp), INTENT(out)            :: target_cn_fine_root    !< fine root target C:N
    REAL(wp), INTENT(out)            :: target_cn_coarse_root  !< coarse root target C:N
    REAL(wp), INTENT(out)            :: target_cn_sap_wood     !< sap wood target C:N
    REAL(wp), INTENT(out)            :: target_cn_fruit        !< fruit target C:N
    REAL(wp), INTENT(out)            :: target_np_fine_root    !< fine root target N:P
    REAL(wp), INTENT(out)            :: target_np_coarse_root  !< coarse root target N:P
    REAL(wp), INTENT(out)            :: target_np_sap_wood     !< sap wood target N:P
    REAL(wp), INTENT(out)            :: target_np_fruit        !< fruit target N:P


    !> 1.0 calculate optimal allocation fractions
    !! 
    IF (veg_config%flag_optimal_Nleaf.AND.plant%growing_season) THEN
       ! update leaf CN target. leaf_cn_direction has values of -1, 0 or 1
       plant%target_cn_leaf = plant%target_cn_leaf * (1.0_wp - plant%leaf_cn_direction * delta_n_leaf)
    ENDIF

    !> 2.0 determine C:N of other plant tissues
    !!
    target_cn_fine_root = plant%target_cn_leaf / root2leaf_cn
    SELECT CASE(growthform)
      CASE (itree) ! trees
         target_cn_coarse_root = plant%target_cn_leaf / wood2leaf_cn
         target_cn_sap_wood = plant%target_cn_leaf / wood2leaf_cn 
      CASE (igrass) ! grasses 
         target_cn_coarse_root = plant%target_cn_leaf / root2leaf_cn
         target_cn_sap_wood = plant%target_cn_leaf / root2leaf_cn 
    END SELECT
    target_cn_fruit = plant%target_cn_leaf / root2leaf_cn

    !> 3.0 determine N:P of other plant tissues
    !!
    target_np_fine_root = plant%target_np_leaf / root2leaf_np
    SELECT CASE(growthform)
      CASE (itree) ! trees
         target_np_coarse_root = plant%target_np_leaf / wood2leaf_np
         target_np_sap_wood = plant%target_np_leaf / wood2leaf_np
      CASE (igrass) ! grasses 
         target_np_coarse_root = plant%target_np_leaf / root2leaf_np
         target_np_sap_wood = plant%target_np_leaf / root2leaf_np
    END SELECT
    target_np_fruit = plant%target_np_leaf / root2leaf_np

  END SUBROUTINE update_stoichiometry

  !-----------------------------------------------------------------------------------------------------
  !> calculates the activity of the meristem for growth 
  !!
  !! As this is proportional to the size of the labile pool, there is no need for a direct temperature 
  !! function of the turnover, as this will result from the dynamics of the labile pool. \nThis function
  !!  therefore only takes account of the impact of low temperatures (<5°C) and low moisture 
  !! (<20\% saturation). \nAdditionally, phenological dormancy is accounted for.
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION meristem_activity(t_air,theta,growing_season) RESULT(tau_labile_pot)

    USE mo_shared_constants,  ONLY: Tzero
    USE mo_veg_constants,     ONLY: lambda_tau_labile_temp, weibull_k_tau_labile_temp, & 
                                    lambda_tau_labile_theta,weibull_k_tau_labile_theta,tau_labile_default
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in)     :: t_air                  !< air temperature (K) 
    REAL(wp),INTENT(in)     :: theta                  !< soil moisture scalar (%) 
    LOGICAL,INTENT(in)      :: growing_season         !< is the plant growing? 
    REAL(wp)                :: tau_labile_pot         !< current potential labile pool turnover (1/day)
    ! 0.2 Local
    REAL(wp)                :: ftemp                  ! temperature limit on cell construction rate
    REAL(wp)                :: fmoist                 ! moisture limit on cell construction rate
  
    IF(t_air > Tzero) THEN 
        ftemp = 1.0_wp - exp(-(lambda_tau_labile_temp * (t_air - Tzero))**weibull_k_tau_labile_temp)
      ELSE 
        ftemp = 0.0_wp
    ENDIF
    fmoist = 1.0_wp - exp(-(lambda_tau_labile_theta * theta)**weibull_k_tau_labile_theta) 
    IF(growing_season)THEN
        tau_labile_pot = ftemp * fmoist * tau_labile_default
      ELSE 
        tau_labile_pot = 0.0_wp
    ENDIF

  END FUNCTION meristem_activity 

  !-----------------------------------------------------------------------------------------------------
  !> calculates the current allometry of a tree 
  !!
  !! diagnoses diameter and height from woody biomass by assuming that the trunk is a cylinder following
  !!  a given allometric relationship: H = k_allom1 * D ** k_allom2
  !!
  !! calculates leaf to sapwood and root mass ratios
  !!
  !! calculates maximum LAI given current root and sapwood mass
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_allometry(sap_wood_carbon,heart_wood_carbon,dens_ind, labile_carbon, labile_nitrogen, labile_phosphorus,&
       wood_density,sla,allom_k1, allom_k2,k_latosa,k_rtos,growthform,cn_plant, cp_plant,  dnpp, dN_up, dP_up, npp, tau_labile_pot,&
       diameter,height,dheight,leaf2sapwood_mass_ratio,leaf2root_mass_ratio,target_lai, k1_opt, k2_opt)
    
    USE mo_base_constants, ONLY: pi
    USE mo_veg_constants,  ONLY: min_height, itree, igrass, sm2lm_grass, target_lai_max
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in)  :: sap_wood_carbon         !< woody biomass of an individual tree (mol C)
    REAL(wp),INTENT(in)  :: heart_wood_carbon       !< woody biomass of an individual tree (mol C)
    REAL(wp),INTENT(in)  :: dens_ind                !< density of populations (#indivuduals /m2) 
    REAL(wp),INTENT(in)  :: wood_density            !< wood density (mol C / m3)
    REAl(wp),INTENT(in)  :: labile_carbon           !< labile carbon 
    REAL(wp),INTENT(in)  :: labile_nitrogen         !< labile nitrogen
    REAL(wp),INTENT(in)  :: labile_phosphorus       !< labile phosphorus
    REAL(wp),INTENT(in)  :: sla                     !< specific leaf area (mol C / m2)
    REAL(wp),INTENT(in)  :: allom_k1                !< coeffcient in H-D relationship
    REAL(wp),INTENT(in)  :: allom_k2                !< exponent in H-D relationship
    REAL(wp),INTENT(in)  :: k_latosa                !< leaf area to sapwood area ratio 
    REAL(wp),INTENT(in)  :: k_rtos                  !< root to sapwood mass ratio
    INTEGER,INTENT(in)   :: growthform              !< plant growth form 
    REAl(wp), INTENT(in) :: cn_plant                !< C:N of whole plant growth increment
    REAl(wp), INTENT(in) :: cp_plant                !< C:P of whole plant growth increment
    REAl(wp), INTENT(in) :: dnpp                    !< unit increase in NPP for unti increase in leaf mass
    REAl(wp), INTENT(in) :: dN_up                   !< unit N uptake per unit root mass
    REAl(wp), INTENT(in) :: dP_up                   !< unit P uptake per unit root mass
    REAl(wp), INTENT(in) :: npp                     !< net primary productivity
    REAL(wp), INTENT(in) :: tau_labile_pot          !< labile pool turnover rate
    REAL(wp),INTENT(out) :: diameter                !< diameter (m)
    REAL(wp),INTENT(out) :: height                  !< height (m)
    REAL(wp),INTENT(out) :: dheight                 !< change in height (m) resulting from a small change in woody biomass
    REAL(wp),INTENT(out) :: leaf2sapwood_mass_ratio !< leaf to sapwood mass ratio 
    REAL(wp),INTENT(out) :: leaf2root_mass_ratio    !< leaf to root mass ratio
    REAL(wp),INTENT(out) :: target_lai              !< lai implied by current available supporting structure 
    REAL(Wp),INTENT(out) :: k1_opt, &               !< coefficient in optimized leaf:root ratio
                            k2_opt                  !< coefficient in optimized leaf:root ratio
    ! 0.2 Local
    REAL(wp)             :: diameter2               ! diameter (m) resulting from a infitisemal small increase in woody biomas 

    !> 1.0 case growthform
    !! 
    SELECT CASE (growthform)
      !> 1.1 trees
      !! 
      CASE (itree)
        ! Diameter and Height calculated assuming that 
        !   all woody biomass forms a cylinder: diameter = SQRT(4 * wood volume / height / pi)
        !   height = allom_k1 * diameter ^ allom_k2
        IF(dens_ind > eps8) THEN
           diameter = ( ( sap_wood_carbon + heart_wood_carbon ) / dens_ind / & 
                        (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
           height = MAX(min_height,allom_k1 * (diameter ** allom_k2))

           ! Calculate change of height given a small change in woody biomass
           diameter2 = ( ( sap_wood_carbon + heart_wood_carbon + eps4 ) / dens_ind / & 
                         (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
           dheight = MAX(min_height,allom_k1 * (diameter2 ** allom_k2)) - height
        ELSE 
           diameter = 0.0_wp
           height = 0.0_wp
           dheight = 0.0_wp
        ENDIF
 
        ! scaling factor for leaf to wood mass, derived from pipe-model (Shinozaki et al. 1964, J J Ecol)
        ! leaf area = k_latosa * sap_wood area
        ! => Cl * sla = k * Cs / wood_density / height
        ! => Cl = k_latosa / sla / wood_density / height * Cs = leaf2sapwood_mass_ratio * Cs / height
        leaf2sapwood_mass_ratio = k_latosa / sla / wood_density

        ! full leaf area implied by current sapwood area
        target_lai = sap_wood_carbon * k_latosa / wood_density / height

        IF (veg_config%flag_optimal_allocation) THEN
           !> Calculate allocation fractions assuming optimal relationship between leaf and root allocation
           !! potential carbon growth rate = potential nittrogen growth rate * plant cn ratio
           !! NPP = gN * R * cn_plant
           !! Or similar relationships for P
           IF ( (npp+(labile_carbon - labile_nitrogen*cn_plant)*1.e6*tau_labile_pot/one_day)/(dN_up * cn_plant) >= &
                 (npp+(labile_carbon - labile_phosphorus*cp_plant)*1.e6*tau_labile_pot/one_day)/(dP_up * cp_plant)) THEN
              ! if N is more limiting than P
              k1_opt = (npp+(labile_carbon - labile_nitrogen*cn_plant)*1.e6*tau_labile_pot/one_day)/(dN_up * cn_plant)
              k2_opt = dnpp/(dN_up * cn_plant)
           ELSE
              k1_opt = (npp+(labile_carbon - labile_phosphorus*cp_plant)*1.e6*tau_labile_pot/one_day)/(dP_up * cp_plant)
              k2_opt = dnpp/(dP_up * cp_plant)
           ENDIF

        ELSE             
           ! scaling factor for root to leaf mass
           ! ad 1: derived from hydraulic homeostatis of root and wood hydraulic resistance (Magnani et al. 2000, PCE)
           ! Note: currently still lacks a moisture response!
           ! sapwood mass / height =  k_rtos * fine_root_mass
           ! Cr = Cs / k_rtos / height
           ! let Cs = Cl / leaf2sapwood_mass_ratio * height
           ! => Cl = Cr * leaf2sapwood_mass_ratio * k_rtos = Cl = Cr * leaf2root_mass_ratio 
           ! ad 2: nitrogen/phosphorus sensitivity still needed
           leaf2root_mass_ratio = k_rtos * leaf2sapwood_mass_ratio
        ENDIF 
      !> 1.1 grass
      !! 
      CASE (igrass)
         ! leaf to stem mass ratio
         leaf2sapwood_mass_ratio = 1._wp / sm2lm_grass
         ! leaf:root ratio
         IF (veg_config%flag_optimal_allocation) THEN
            !> Calculate allocation fractions assuming optimal relationship between leaf and root allocation
            !! potential carbon growth rate = potential nittrogen growth rate * plant cn ratio
            !! NPP = gN * R * cn_plant
            !! Or similar relationships for P
           IF ( (npp+(labile_carbon - labile_nitrogen*cn_plant)*1.e6*tau_labile_pot/one_day)/(dN_up * cn_plant) >= &
                 (npp+(labile_carbon - labile_phosphorus*cp_plant)*1.e6*tau_labile_pot/one_day)/(dP_up * cp_plant)) THEN
              ! if N is more limiting than P
              k1_opt = (npp+(labile_carbon - labile_nitrogen*cn_plant)*1.e6*tau_labile_pot/one_day)/(dN_up * cn_plant)
              k2_opt = dnpp/(dN_up * cn_plant)
           ELSE
              k1_opt = (npp+(labile_carbon - labile_phosphorus*cp_plant)*1.e6*tau_labile_pot/one_day)/(dP_up * cp_plant)
              k2_opt = dnpp/(dP_up * cp_plant)
           ENDIF

         ELSE      
            leaf2root_mass_ratio = sm2lm_grass * k_rtos 
         ENDIF
         ! full leaf area - should depend on carbon balance of lowest leaf - works in principle but needs some filtering
         target_lai = 6.0_wp
    ENDSELECT

  END SUBROUTINE calc_allometry 

  !----------------------------------------------------------------------------------------------------------
  !> Calculate flux between reserve and labile pools for each element 
  !! 
  !! Input:\n
  !!    (1) weekly GPP and maintenance respiration, and N,P growth requirements to calculate
  !!        target labile pool size;\n
  !!    (2) target LAI, leaf-to-root mass ratio and target C:N:P of foliage and fine roots to calculate
  !!        target reserve pool size;\n
  !!    (3) current LAI, labile and reserve pools, growing season status and potential turnover rate
  !!          of the labile pool;\n    
  !!  
  !! Output: net flux between reserve and labile pools 
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_reserve_useage(tau_labile_pot, growing_season, gpp_week, rm_week, & 
                               growth_req_n_week, growth_req_p_week, sla, tau_leaf, tau_fine_root, & 
                               target_cn_leaf, target_cn_fine_root, target_np_leaf, target_np_fine_root, & 
                               leaf2root_mass_ratio, lai, target_lai, labile_pool, reserve_pool) & 
                               RESULT (reserve_use)

    USE mo_shared_memory,    ONLY: t_element

    IMPLICIT NONE
    ! 0.1 InOut
    LOGICAL,INTENT(in)          :: growing_season            !< does the plant intend to grow?
    REAL(wp),INTENT(in)         :: tau_labile_pot            !< potential labile pool turnover (1/day) 
    REAL(wp),INTENT(in)         :: gpp_week                  !< last week's GPP (µmol CO2/m2/s) 
    REAL(wp),INTENT(in)         :: rm_week                   !< last week's Rm (µmol CO2/m2/s)
    REAL(wp),INTENT(in)         :: growth_req_n_week         !< last week's requirement for unit growth in C and N
    REAL(wp),INTENT(in)         :: growth_req_p_week         !< last week's rP requirement for unit growth in C and N
    REAL(wp),INTENT(in)         :: sla                       !< specific leaf area
    REAL(wp),INTENT(in)         :: tau_leaf                  !< default turnover time of leaves (1/year) 
    REAL(wp),INTENT(in)         :: tau_fine_root             !< default turnover time of fine roots (1/year)
    REAL(wp),INTENT(in)         :: target_cn_leaf            !< target C:N and N:P rations for tissue type 
    REAL(wp),INTENT(in)         :: target_cn_fine_root       !< target C:N and N:P rations for tissue type 
    REAL(wp),INTENT(in)         :: target_np_leaf            !< target C:N and N:P rations for tissue type 
    REAL(wp),INTENT(in)         :: target_np_fine_root       !< target C:N and N:P rations for tissue type 
    REAL(wp),INTENT(in)         :: leaf2root_mass_ratio      !< leaf to root mass ratio
    REAL(wp),INTENT(in)         :: lai                       !< current LAI 
    REAL(wp),INTENT(in)         :: target_lai                !< target LAI 
    TYPE(t_element),INTENT(in)  :: labile_pool               !< the plant's labile pool size
    TYPE(t_element),INTENT(in)  :: reserve_pool              !< the plant's reserve pool size
    TYPE(t_element)             :: reserve_use               !< current resource use rate (mol/time step)
                                                             !!  (+ve: net flux to labile / -ve: net flux to reserve)
    ! 0.2 Local
    REAL(wp)             :: tau_labile                       ! actual turnover time of labile pool (1/day)
    REAL(wp)             :: phi_phenology                    ! relative pull from reserve due to phenology
    REAL(wp)             :: phi_maintenance                  ! relative pull from reserve due to maintenance resp. 
    REAL(wp)             :: phi_storage                      ! relative pull from labile due to low reserve
    REAL(wp)             :: target_labile_pool_carbon        ! desirable size (mol C/m2) of the labile pool
    REAL(wp)             :: target_labile_pool_nitrogen      ! desirable size (mol N/m2) of the labile pool
    REAL(wp)             :: target_labile_pool_phosphorus    ! desirable size (mol P/m2) of the labile pool
    REAL(wp)             :: target_reserve_pool_carbon       ! desirable size (mol C/m2) of the reserve pool
    REAL(wp)             :: target_reserve_pool_nitrogen     ! desirable size (mol N/m2) of the reserve pool
    REAL(wp)             :: target_reserve_pool_phosphorus   ! desirable size (mol P/m2) of the reserve pool
    REAL(wp)             :: hlp1                             ! helper variable

    !> 1.0 target size of the labile and reserve pools
    !!

    !> 1.1 labile pool
    !! 
    !! This is taken as the average idaily input + maintenance requirements of the last week 
    !!  scaled to the default turnover time of the labile pool
    !! 
    !! for N and P this corresponds to the current growth C:N:P  
    target_labile_pool_carbon = MAX( gpp_week , rm_week ) * one_day / 1.e6 / tau_labile_default
    target_labile_pool_nitrogen = labile_pool%carbon * growth_req_n_week
    target_labile_pool_phosphorus = labile_pool%nitrogen * growth_req_p_week

    !> 1.2 reserve pool
    !!
    !! For the reserve pool, this is the C:N:P needed to regrow one year worth of folige and roots,
    !! accounting for losses to respiration and the nutrients currently required to grow these
    target_reserve_pool_carbon = target_lai / sla * & 
           ( MIN(1.0_wp,tau_leaf) + tau_fine_root / leaf2root_mass_ratio ) * (1._wp + fresp_growth)
    target_reserve_pool_nitrogen = target_lai / sla * & 
           ( MIN(1.0_wp,tau_leaf) / target_cn_leaf + tau_fine_root / leaf2root_mass_ratio / target_cn_fine_root) 
    target_reserve_pool_phosphorus = target_lai / sla * & 
           ( MIN(1.0_wp,tau_leaf) / target_cn_leaf / target_np_leaf & 
             + tau_fine_root / leaf2root_mass_ratio / target_cn_fine_root / target_np_fine_root)  

    !> 2.0 carbon reserve flow
    !! 
    !! calculate the pull from reserve to labile (phi_phenology,phi_maintenance) or
    !! 
    !! push from labile to reserve (phi_storage) depending on the difference of the LAI/labile/reserve from their targets
    !! 
    !! phi_storage is adjusted by phi_phenology to avoid storage building to interfere with flushing (push phi_storage towards 
    !!  reserve at times of low LAI or labile pool to give priority to growth and maintentance)
    phi_phenology = 0._wp
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    hlp1 = MIN(MAX(min_lai,target_lai),target_lai_max)
    IF(lai < hlp1 .AND. growing_season) & 
       phi_phenology = exp(-(lambda_phiphen * lai / hlp1 )**weibull_k_phiphen) & 
                       - exp(-(lambda_phiphen)**weibull_k_phiphen) 

    hlp1 = MAX(0.0_wp,labile_pool%carbon)
    IF(target_labile_pool_carbon > eps8 .AND. hlp1 < target_labile_pool_carbon) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_carbon)**weibull_k_phimaint) & 
                       - exp(-(lambda_phimaint)**weibull_k_phimaint)

    hlp1 = MAX(0.0_wp,reserve_pool%carbon)
    IF(target_reserve_pool_carbon > eps8 .AND. hlp1 < target_reserve_pool_carbon) & 
       phi_storage = exp(-(lambda_phistore * hlp1 / target_reserve_pool_carbon)**weibull_k_phistore) & 
                       - exp(-(lambda_phistore)**weibull_k_phistore)

    phi_storage = phi_storage * ( 1.0_wp - phi_phenology) 
    IF(phi_maintenance > 0.75_wp) phi_storage = phi_storage * 4.0_wp * ( 1.0_wp - phi_maintenance) 

    !> 2.2 convert phi to turnover times (1/timestep)
    !! 
    !! taking account that flushing can only happen if the meristem is active
    phi_phenology = phi_phenology * tau_labile_pot * dt/one_day   
    phi_maintenance = phi_maintenance * tau_labile_default * dt/one_day
    phi_storage = phi_storage * tau_labile_default * dt/one_day

    reserve_use%carbon   = ((phi_phenology + phi_maintenance) * reserve_pool%carbon - & 
                           phi_storage * labile_pool%carbon) 
    reserve_use%carbon13 = ((phi_phenology + phi_maintenance) * reserve_pool%carbon13 - & 
                             phi_storage * labile_pool%carbon13) 
    reserve_use%carbon14 = ((phi_phenology + phi_maintenance) * reserve_pool%carbon14 - & 
                             phi_storage * labile_pool%carbon14) 

    !> 3.0 net transfer of N from labile to reserve
    !! 
    !! 
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    hlp1 = MAX(0.0_wp,labile_pool%nitrogen)
    IF(target_labile_pool_nitrogen > eps8 .AND. hlp1 < target_labile_pool_nitrogen) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_nitrogen)**weibull_k_phimaint) & 
                       - exp(-(lambda_phimaint)**weibull_k_phimaint)
    IF(target_reserve_pool_nitrogen > eps8 .AND. reserve_pool%nitrogen < target_reserve_pool_nitrogen) & 
       phi_storage = exp(-(lambda_phistore * reserve_pool%nitrogen / target_reserve_pool_nitrogen)**weibull_k_phistore) & 
                       - exp(-(lambda_phistore)**weibull_k_phistore)
       
    reserve_use%nitrogen = (phi_maintenance * reserve_pool%nitrogen - & 
                           phi_storage * labile_pool%nitrogen) * tau_labile_default * dt/one_day
    reserve_use%nitrogen15 = (phi_maintenance * reserve_pool%nitrogen15 - & 
                           phi_storage * labile_pool%nitrogen15) * tau_labile_default * dt/one_day

    !> 4.0 net transfer of P from labile to reserve
    !! 
    !! 
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    hlp1 = MAX(0.0_wp,labile_pool%phosphorus)
    IF(target_labile_pool_phosphorus > eps8 .AND. hlp1 < target_labile_pool_phosphorus) & 
       phi_maintenance = exp(-(lambda_phimaint * hlp1 / target_labile_pool_phosphorus)**weibull_k_phimaint) & 
                       - exp(-(lambda_phimaint)**weibull_k_phimaint)
    IF(target_reserve_pool_phosphorus > eps8 .AND. reserve_pool%phosphorus < target_reserve_pool_phosphorus) & 
       phi_storage = exp(-(lambda_phistore * reserve_pool%phosphorus / target_reserve_pool_phosphorus)**weibull_k_phistore) & 
                       - exp(-(lambda_phistore)**weibull_k_phistore)

    reserve_use%phosphorus = (phi_maintenance * reserve_pool%phosphorus - & 
                           phi_storage * labile_pool%phosphorus) * tau_labile_default * dt/one_day

  END FUNCTION calc_reserve_useage

  !----------------------------------------------------------------------------------------------------------
  !> Calculate nutrient flux from enzyme turnover and other biochemical turnover in leaves and roots
  !! at the time scale of tau_canopy_acclimation
  !!
  !! Calculations only during growing season \n
  !! Input: plant status and leaf and root C:N:P target ratios \n
  !! Output: net flux between labile and leaves/roots (+ve: flux to labile, -ve: flux from labile)
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_nutrient_recycling(plant, target_cn_fine_root, target_np_fine_root)

    USE mo_veg_constants, ONLY: tau_canopy_acclimation, tau_labile_default

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout) :: plant                    !< plant object
    REAL(wp), INTENT(in)             :: target_cn_fine_root      !< target C:N of fine_roots
    REAL(wp), INTENT(in)             :: target_np_fine_root      !< target N:P of fine_roots 
    ! 0.2 Local
    REAL(wp)                         :: labile_dn15              ! labile pools 15N to N molar mixing ratio
    REAL(wp)                         :: scal                     ! helper variable
    REAL(wp)                         :: recycling_leaf,recycling_fine_root  ! more helpers

    IF(plant%growing_season)THEN 
       !> 1.0 replace overturning N from leaves and roots by new N from the labile pool using the target C:N
       !!
       recycling_leaf = ( plant%pool%leaf%nitrogen - & 
                  plant%pool%leaf%carbon / plant%target_cn_leaf ) * tau_canopy_acclimation * dt / one_day
       recycling_fine_root = ( plant%pool%fine_root%nitrogen - &
                  plant%pool%fine_root%carbon / target_cn_fine_root ) * tau_canopy_acclimation * dt / one_day
       !> 1.1 if N is addded to leaves or fine root (negative flux), make sure there is enough labile N to support this
       !! 
       IF( -1.0_wp * (recycling_leaf + recycling_fine_root) > eps8) THEN
          scal = -1.0_wp * ( tau_labile_default * dt / one_day * plant%pool%labile%nitrogen ) / & 
                           ( recycling_leaf + recycling_fine_root ) 
          IF(scal < 1.0_wp) THEN
             recycling_leaf = recycling_leaf * scal
             recycling_fine_root = recycling_fine_root * scal
          ENDIF
       ENDIF
       plant%flux%recycling_leaf_n      = plant%flux%recycling_leaf_n      + recycling_leaf
       plant%flux%recycling_fine_root_n = plant%flux%recycling_fine_root_n + recycling_fine_root
   
       !> 2.0 replace overturning P from leaves and roots by new N from the labile pool using the target N:P
       !!
       recycling_leaf = ( plant%pool%leaf%phosphorus - & 
                  plant%pool%leaf%nitrogen / plant%target_np_leaf ) * tau_canopy_acclimation * dt / one_day
       recycling_fine_root = ( plant%pool%fine_root%phosphorus - &
                  plant%pool%fine_root%nitrogen / target_np_fine_root ) * tau_canopy_acclimation * dt / one_day
       !> 2.1 if P is addded to leaves or fine root, make sure there is enough labile P to support this
       !! 
       IF( -1.0_wp * ( recycling_leaf + recycling_fine_root ) > eps8) THEN
          scal = -1.0_wp * ( tau_labile_default * dt / one_day * plant%pool%labile%phosphorus ) / & 
                           ( recycling_leaf + recycling_fine_root ) 
          IF(scal < 1.0_wp) THEN
             recycling_leaf = recycling_leaf * scal
             recycling_fine_root = recycling_fine_root * scal
          ENDIF
       ENDIF
       plant%flux%recycling_leaf_p        = plant%flux%recycling_leaf_p        + recycling_leaf
       plant%flux%recycling_fine_root_p   = plant%flux%recycling_fine_root_p   + recycling_fine_root
   
       !> 3.0 replace overturning 15N from leaves and roots by new 15N from the labile pool using the target C:N
       !! and the labile pool's 15N molar mixing ratio
       !!
       IF(plant%pool%labile%nitrogen > eps8) THEN
          labile_dn15 = plant%pool%labile%nitrogen15/plant%pool%labile%nitrogen
       ELSE 
          labile_dn15 = 0.0_wp
       ENDIF
       recycling_leaf = ( labile_dn15 * plant%pool%leaf%carbon / plant%target_cn_leaf - & 
                  labile_dn15 * plant%pool%leaf%carbon ) * tau_canopy_acclimation * dt / one_day
       recycling_fine_root = ( plant%pool%fine_root%nitrogen15 - & 
                  labile_dn15 * plant%pool%fine_root%carbon / target_cn_fine_root ) * tau_canopy_acclimation * dt / one_day
       !> 3.1 if N is addded to leaves or fine root, make sure there is enough labile N to support this
       !! 
       IF( -1.0_wp * ( recycling_leaf + recycling_fine_root ) > eps8) THEN
          scal = -1.0_wp * ( tau_labile_default * dt / one_day * plant%pool%labile%nitrogen15 ) / & 
                           ( recycling_leaf + recycling_fine_root ) 
          IF(scal < 1.0_wp) THEN
             recycling_leaf = recycling_leaf * scal
             recycling_fine_root = recycling_fine_root * scal
          ENDIF
       ENDIF
       plant%flux%recycling_leaf_n15      = plant%flux%recycling_leaf_n15      + recycling_leaf
       plant%flux%recycling_fine_root_n15 = plant%flux%recycling_fine_root_n15 + recycling_fine_root
    ENDIF

  END SUBROUTINE calc_nutrient_recycling

END MODULE mo_veg_growth
