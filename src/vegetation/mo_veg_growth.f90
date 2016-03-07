!!
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
  USE mo_base_memory

  PRIVATE

  PUBLIC :: calc_growth

CONTAINS

  SUBROUTINE calc_growth(plant,veg_param,growth,growth_respiration,growth_respiration_c13,growth_respiration_c14)
    !
    ! Routine to calculate plant growth 
    ! 
    USE mo_veg_parameters
    USE mo_base_memory, ONLY: ElementMul

    IMPLICIT NONE

    TYPE(PlantClass), INTENT(in) :: plant
    TYPE(VegetationParameterClass),INTENT(in) :: veg_param 
    TYPE(VegPool),INTENT(out) :: growth
    REAL(wp),INTENT(out) :: growth_respiration,growth_respiration_c13,growth_respiration_c14 
    TYPE(VegPool) :: falloc
    !
    ! local variables
    !
    ! target C:N and N:P ratios for various tissue types (mol/mol)
    REAL(wp) :: target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit
    ! N and P reqiurement of growth (mol N / P per mol C)
    REAL(wp) :: nreq_growth,preq_growth

    CALL update_stoichiometry(target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, & 
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, & 
                              target_np_sap_wood,target_np_fruit)
    CALL calc_allocation_fraction(plant%pool%leaf%carbon,plant%pool%sap_wood%carbon,plant%pool%heart_wood%carbon, & 
                              plant%pool%fine_root%carbon,veg_param%sla(plant%pftid),veg_param%k_latosa(plant%pftid), &
                              veg_param%k_rtos(plant%pftid),veg_param%wood_density(plant%pftid), & 
                              veg_param%allom_k1(plant%pftid),veg_param%allom_k2(plant%pftid), & 
                              veg_param%growthform(plant%pftid), & 
                              target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit,falloc,nreq_growth,preq_growth)
    CALL meristem_activity(nreq_growth,preq_growth,plant%pool%labile,growth%labile)

    growth_respiration = fresp_growth * growth%labile%carbon * 1.e6 / dt
    growth_respiration_c13 = fresp_growth * growth%labile%carbon13 * 1.e6 / dt
    growth_respiration_c14 = fresp_growth * growth%labile%carbon14 * 1.e6 / dt

    growth%leaf = ElementMul ( falloc%leaf, growth%labile )
    growth%fine_root = ElementMul ( falloc%fine_root, growth%labile ) 
    growth%coarse_root = ElementMul ( falloc%coarse_root, growth%labile ) 
    growth%sap_wood = ElementMul ( falloc%sap_wood, growth%labile ) 
    growth%fruit = ElementMul ( falloc%fruit, growth%labile ) 

    ! invert flux of labile to be out of the labile pool
    growth%labile%carbon = -1.0_wp * growth%labile%carbon 
    growth%labile%nitrogen = -1.0_wp * growth%labile%nitrogen
    growth%labile%phosphorus = -1.0_wp * growth%labile%phosphorus
    growth%labile%carbon13 = -1.0_wp * growth%labile%carbon13
    growth%labile%carbon14 = -1.0_wp * growth%labile%carbon14
    growth%labile%nitrogen15 = -1.0_wp * growth%labile%nitrogen15

  END SUBROUTINE calc_growth 

  SUBROUTINE meristem_activity(nreq_growth,preq_growth,labile_pool,growth)
    ! 
    ! Routine to calculate meristem activity and thus labile pool turnover
    ! OPS: Dummy only
    ! Input: mol N and P required per mol C growth, labile pool size
    ! Output: growth rate
    ! 
    USE mo_land_constants
    USE mo_base_memory
    USE mo_veg_parameters
    IMPLICIT NONE
    ! 0. Input
    REAL(wp), INTENT(in)     :: nreq_growth,preq_growth     ! N and P requirement for unit growth in C and N
    TYPE(Element)            :: labile_pool                 ! the plant's labile pool size
    ! 0.2 Output
    TYPE(Element)            :: growth                      ! current growth rate (mol/time step)
    ! 0.3 Local
    REAL(wp)                 :: tau_labile                  ! actual growth (1/day)
    !
    !  1.0 Initial guess of growth rate, should be influenced by moisture and temperature
    !
    tau_labile = tau_labile_default * dt/one_day

    ! 2.0 Elemental limitation of growth (this should be replaced at a later stage - here to conserve mass)
    !
    ! 2.1 N
    growth%carbon = tau_labile * labile_pool%carbon
    growth%nitrogen = tau_labile * labile_pool%nitrogen
    IF(nreq_growth * growth%carbon > eps8)THEN
        IF(growth%nitrogen < nreq_growth * growth%carbon) & 
               tau_labile = tau_labile * growth%nitrogen / ( nreq_growth * growth%carbon ) 
    ELSE 
        tau_labile = zero
    ENDIF 

    ! 2.2 P
    !
    growth%nitrogen = tau_labile * labile_pool%nitrogen
    growth%phosphorus = tau_labile * labile_pool%phosphorus
    IF(preq_growth * growth%phosphorus > eps8)THEN
        IF(growth%phosphorus < preq_growth * growth%nitrogen ) & 
              tau_labile = tau_labile * growth%phosphorus / ( preq_growth * growth%nitrogen ) 
    ELSE 
        tau_labile = zero
    ENDIF 

    ! 3.0 Calculate final growth rate
    !
    growth%carbon = tau_labile * labile_pool%carbon
    growth%nitrogen = MIN(tau_labile * labile_pool%nitrogen, nreq_growth * growth%carbon)
    growth%phosphorus = MIN(tau_labile * labile_pool%phosphorus, preq_growth * growth%nitrogen)
    growth%carbon13 = tau_labile * labile_pool%carbon13
    growth%carbon14 = tau_labile * labile_pool%carbon14
    growth%nitrogen15 = tau_labile * labile_pool%nitrogen15

  END SUBROUTINE meristem_activity

  SUBROUTINE calc_allocation_fraction(leaf_carbon,sap_wood_carbon,heart_wood_carbon,root_carbon,sla, &
                              k_latosa,k_rtos,wood_density,allom_k1,allom_k2,growthform, &
                              target_cn_leaf,target_cn_fine_root,target_cn_coarse_root, & 
                              target_cn_sap_wood,target_cn_fruit,target_np_leaf,target_np_fine_root, & 
                              target_np_coarse_root,target_np_sap_wood,target_np_fruit,falloc,nreq_growth,& 
                              preq_growth)
    !
    ! Routine to calculate allocation fractions for each element
    ! Input: target stoichiometry of the tissue types
    ! Output: 1) Partitioning coefficients for the distribution of biomass for growth
    !         2) total N and P required to sustain a unit of growth (in mol)
    ! 
    USE mo_land_constants
    USE mo_base_memory
    USE mo_base_constants, ONLY: pi
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: leaf_carbon           ! carbon in leaves (mol/individual)
    REAL(wp), INTENT(in) :: sap_wood_carbon       ! carbon in sap wood (mol/individual)
    REAL(wp), INTENT(in) :: heart_wood_carbon     ! carbon in heart wood (mol/individual)
    REAL(wp), INTENT(in) :: root_carbon           ! carbon in roots (mol/individual) 
    REAL(wp), INTENT(in) :: sla                   ! specific leaf area
    REAL(wp), INTENT(in) :: k_latosa              ! leaf area to sapwood area ratio
    REAL(wp), INTENT(in) :: k_rtos                ! root tradeoff parameter
    REAL(wp), INTENT(in) :: wood_density          ! specific wood density 
    REAL(wp), INTENT(in) :: allom_k1              ! parameter 1 in height diameter relationship
    REAL(wp), INTENT(in) :: allom_k2              ! parameter 2 in height diameter relationship
    INTEGER, INTENT(in) :: growthform             ! plant growth form 
    REAL(wp), INTENT(in) :: target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit        ! target C:N and N:P rations for tissue types
    ! 0.2 In- and Output
    TYPE(VegPool), INTENT(inout) :: falloc                              ! allocation coefficients 
    ! 0.3 Output
    REAL(wp), INTENT(out) :: nreq_growth, preq_growth                   ! mol N/P needed per mol growth of C
    ! 0.4 Local
    REAL(wp) height,dheight,diameter,ddiameter                          ! height, diameter and their derivative
    REAL(wp) b_inc,leaf_carbon_inc,root_carbon_inc,sap_wood_carbon_inc  ! total and compartment-wise carbon increment (mol/dt)
    REAL(wp) ltos,ltor                                                  ! leaf to sapwood and leaf to root mass ratio

    !
    ! 1.0 Allocation factors for C (this should be calculated based on various principles of plant growth and allometry)
    ! 
    ! SZ: uncomment this to revert to fixed constant allocation
    ! falloc%leaf%carbon = 0.3_wp
    ! falloc%fine_root%carbon = 0.3_wp
    ! falloc%coarse_root%carbon = 0.1_wp
    ! falloc%sap_wood%carbon = 0.2_wp
    ! falloc%fruit%carbon = 0.1_wp 

    IF(growthform==itree)THEN
       ! 
       ! seed for allocation coefficients is the average growth rate per time step.
       ! for now assume a standard growth rate of 2.4 g / m2 / day
       b_inc = 0.2_wp / one_day * dt 
       !
       ! Diagnose current height and diameter from wood mass assuming a zylindric shape and 
       ! height = k1 * D ^k2
       !  
       diameter = ((sap_wood_carbon + heart_wood_carbon) / & 
                  (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
       height = MAX(3.0_wp,allom_k1 * (diameter ** allom_k2))
       ! derivative of height against biomass increase for approximation of height growth effect on partitioning
       ddiameter = ((sap_wood_carbon + heart_wood_carbon + b_inc) / &
                  (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
       dheight = allom_k1 * ddiameter ** allom_k2 - height
       ! scaling factor for leaf to wood mass, derived from pipe-model (Shinozaki et al. 1964, J J Ecol)
       ! leaf area = k_latosa * sap_wood area
       ! => Cl * sla = k * Cs / wood_density / height
       ! => Cl = k_latosa / sla / wood_density / height * Cs = ltos * Cs / height
       ltos = k_latosa / sla / wood_density
       ! scaling factor for root to leaf mass
       ! ad 1: derived from hydraulic homeostatis of root and wood hydraulic resistance (Magnani et al. 2000, PCE)
       ! Note: currently still lacks a moisture response!
       ! sapwood mass / height =  k_rtos * fine_root_mass
       ! Cr = Cs / k_rtos / height
       ! let Cs = Cl / ltos * height
       ! => Cl = Cr * ltos * k_rtos = Cl = Cr * ltor 
       ! ad 2: nitrogen/phosphorus sensitivity still needed
       ltor = k_rtos * ltos

       ! solving the above three equations simulatenously
       leaf_carbon_inc = MAX(( b_inc + leaf_carbon * dheight / ltos * & 
                         ( leaf_carbon / ltor - b_inc - root_carbon ) + sap_wood_carbon + root_carbon - & 
                         leaf_carbon * (height / ltos + 1.0_wp / ltor )) / & 
                         (1.0_wp + ( b_inc - (1.0_wp + 2.0_wp/ltor) * leaf_carbon + root_carbon ) & 
                         * dheight / ltos + height / ltos + 1._wp / ltor ), 0.0_wp) 

       sap_wood_carbon_inc = MAX((( leaf_carbon + leaf_carbon_inc ) / ltos * height - sap_wood_carbon) / & 
                             (1._wp - ( leaf_carbon + leaf_carbon_inc ) / ltos * dheight ), 0.0_wp)          
       IF(sap_wood_carbon_inc > 0.0_wp) THEN
          root_carbon_inc = MAX(( leaf_carbon + leaf_carbon_inc ) / ltor - root_carbon, 0.0_wp)
       ELSE ! special case of no wood growth (phenology etc.): redistribute available carbon 
          leaf_carbon_inc = MAX(( b_inc - leaf_carbon / ltor + root_carbon ) / (1.0_wp + 1._wp/ ltor),0.0_wp)
          root_carbon_inc = MAX(MIN(b_inc - leaf_carbon_inc,( leaf_carbon_inc + leaf_carbon ) / ltor - root_carbon),&
                            0.0_wp)
       ENDIF

       ! convert to allocation fractions 
       falloc%leaf%carbon = leaf_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
       falloc%fine_root%carbon = root_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
       falloc%coarse_root%carbon = 0.0_wp
       falloc%sap_wood%carbon = sap_wood_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
       falloc%fruit%carbon = 0.0_wp 

    ENDIF
    !
    ! 2.0 Allocation factors for N, derived from the target C:N ratios and the C partitioning coefficients
    !  
    ! N required to sustain a unit of C growth
    nreq_growth = falloc%leaf%carbon / target_cn_leaf + & 
                  falloc%fine_root%carbon / target_cn_fine_root + &
                  falloc%coarse_root%carbon / target_cn_coarse_root + &
                  falloc%sap_wood%carbon / target_cn_sap_wood + &
                  falloc%fruit%carbon / target_cn_fruit

    ! Adjustment of N allocation factors
    IF(nreq_growth > eps8)THEN
       falloc%leaf%nitrogen = falloc%leaf%carbon / target_cn_leaf / nreq_growth
       falloc%fine_root%nitrogen = falloc%fine_root%carbon / target_cn_fine_root / nreq_growth 
       falloc%coarse_root%nitrogen = falloc%coarse_root%carbon / target_cn_coarse_root / nreq_growth
       falloc%sap_wood%nitrogen = falloc%sap_wood%carbon / target_cn_sap_wood / nreq_growth
       falloc%fruit%nitrogen = falloc%fruit%carbon / target_cn_fruit / nreq_growth
    ENDIF
    !
    ! 3.0 Allocation factors for P, derived from the target N:P ratios and the N partitioning coefficients
    !  
    ! P required to sustain a unit of N growth
    preq_growth = falloc%leaf%nitrogen / target_np_leaf + & 
                  falloc%fine_root%nitrogen / target_np_fine_root + &
                  falloc%coarse_root%nitrogen / target_np_coarse_root + &
                  falloc%sap_wood%nitrogen / target_np_sap_wood + &
                  falloc%fruit%nitrogen / target_np_fruit

    ! Adjustment of P allocation factors
    IF(preq_growth > eps8)THEN
       falloc%leaf%phosphorus = falloc%leaf%nitrogen / target_np_leaf / preq_growth
       falloc%fine_root%phosphorus = falloc%fine_root%nitrogen / target_np_fine_root / preq_growth 
       falloc%coarse_root%phosphorus = falloc%coarse_root%nitrogen / target_np_coarse_root / preq_growth
       falloc%sap_wood%phosphorus = falloc%sap_wood%nitrogen / target_np_sap_wood / preq_growth
       falloc%fruit%phosphorus = falloc%fruit%nitrogen / target_np_fruit / preq_growth
    ENDIF

    !
    ! 4.0 Allocation factors for isotopes follow the main isotope 
    ! 
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

  END SUBROUTINE calc_allocation_fraction

  SUBROUTINE update_stoichiometry(target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, & 
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, & 
                              target_np_sap_wood,target_np_fruit)
    ! 
    ! Routine to calculate the target stoichiometry of the newly growth tissue
    ! for now: prescribed C:N:P - should be part of the optimality approach
    ! Output: target stoichiometries
    ! 
    USE mo_land_constants
    IMPLICIT NONE
    ! 0.2: Output
    REAL(wp), INTENT(out) ::  target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit
                              ! target C:N and N:P rations (Ops! this is in mol:mol)

    ! for now prescribed and based on OCN
    target_cn_leaf = 25.0_wp * molar_mass_N14 / molar_mass_C12  ! CN adjusted for g/g to mol/mol
    target_cn_fine_root = target_cn_leaf / 0.85_wp 
    target_cn_coarse_root = target_cn_leaf / 0.125_wp 
    target_cn_sap_wood = target_cn_leaf / 0.125_wp
    target_cn_fruit = target_cn_leaf / 0.85_wp 

    ! Mock allocation fractions for P, assuming a univeral aim of N:P = 14 ("optimal") 
    target_np_leaf = 14.0_wp * molar_mass_P31 / molar_mass_N14  ! NP adjusted for g/g to mol/mol
    target_np_fine_root = target_np_leaf / 0.85_wp
    target_np_coarse_root = target_np_leaf / 0.125_wp
    target_np_sap_wood = target_np_leaf / 0.125_wp
    target_np_fruit = target_np_leaf / 0.85_wp

  END SUBROUTINE update_stoichiometry

END MODULE mo_veg_growth
