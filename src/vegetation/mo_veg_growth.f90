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

  SUBROUTINE calc_growth(t_air,theta,plant,veg_param,growth,growth_respiration,growth_respiration_c13,growth_respiration_c14,reserve_use)
    !
    ! Routine to calculate plant growth 
    ! 
    USE mo_base_memory, ONLY: ElementMul
    USE mo_veg_parameters

    IMPLICIT NONE

    ! 0.1 Input
    REAL(wp), INTENT(in) :: t_air                              ! ambient air temperature (K)
    REAL(wp), INTENT(in) :: theta                              ! soil moisture (%)
    TYPE(PlantClass), INTENT(in) :: plant                      ! plant structure
    TYPE(VegetationParameterClass),INTENT(in) :: veg_param     ! vegetation type parameters
    ! 0.3 Output
    TYPE(VegPool),INTENT(out) :: growth                        ! growth rate (mol m-2 timestep-1)
    REAL(wp),INTENT(out) :: growth_respiration,growth_respiration_c13,growth_respiration_c14
                                                               ! growth respiration rate 
                                                               ! (µmol m-2 s-1) 
    TYPE(Element),INTENT(out) :: reserve_use                   ! net reserve use (mol m-2 timestep-1) 
    ! 0.4 Local
    TYPE(VegPool) :: falloc                                    ! fractional allocation to pools
    REAL(wp) :: nreq_growth,preq_growth                        ! N and P reqiurement of growth (mol N / P per mol C)
    REAL(wp) :: target_lai                                     ! LAI implied from sap wood area 
    REAL(wp) :: target_reserve_pool                            ! target reserve pool size (mol m-2)
    REAL(wp) :: tau_labile_pot                                 ! labile pool turnover rate assuming only temperature and moisture
                                                               ! constraints (no stoichiometry)
    ! target C:N and N:P ratios for various tissue types (mol/mol)
    REAL(wp) :: target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit
    ! 
    ! 1.0 calculate potential growth rate relative to labile pool size
    !
    tau_labile_pot = meristem_activity(t_air,theta)
    ! 
    ! 2.0 update target C:N:P for each tissue type
    !
    CALL update_stoichiometry(veg_param%cn_leaf(plant%pftid),veg_param%np_leaf(plant%pftid), & 
                              target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, & 
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, & 
                              target_np_sap_wood,target_np_fruit)
    !
    ! 3.0 calculate allocation coefficients for structural plant pools (leaves, wood, root)
    ! 
    CALL calc_allocation_fraction(tau_labile_pot,plant%pool%leaf%carbon,plant%pool%sap_wood%carbon, & 
                              plant%pool%heart_wood%carbon,plant%pool%fine_root%carbon,plant%pool%labile%carbon,& 
                              veg_param%sla(plant%pftid),veg_param%k_latosa(plant%pftid), &
                              veg_param%k_rtos(plant%pftid),veg_param%wood_density(plant%pftid), & 
                              veg_param%allom_k1(plant%pftid),veg_param%allom_k2(plant%pftid), & 
                              veg_param%tau_leaf(plant%pftid),veg_param%tau_fine_root(plant%pftid), & 
                              veg_param%growthform(plant%pftid), & 
                              target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit,falloc,nreq_growth,preq_growth, &
                              target_lai,target_reserve_pool)
    ! 
    ! 4.0 partition labile carbon into growth, reserve and fruit growth (not yet implemented).
    ! OBS: Should maintenance respiration be here as well??
    ! 
    CALL calc_partitioning(tau_labile_pot,plant%accumulated%gpp_week,plant%accumulated%maint_respiration_week, & 
                              plant%lai,target_lai,target_reserve_pool,nreq_growth,preq_growth, & 
                              plant%pool%labile,plant%pool%reserve,growth%labile,reserve_use)
    !
    ! 5.0 diagnose growth respiration and partition growth to tissue types, convert to µmol m-2 s-1
    ! 
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

  SUBROUTINE calc_partitioning(tau_labile_pot,gpp_week,rm_week,lai,target_lai,target_reserve_pool, & 
                               nreq_growth,preq_growth,labile_pool,reserve_pool,growth,reserve_use)
    ! 
    ! Routine to calculate meristem activity and thus labile pool turnover
    ! OPS: Dummy only
    ! Input: mol N and P required per mol C growth, labile pool size
    ! Output: growth rate
    ! 
    USE mo_base_memory
    USE mo_land_constants
    USE mo_veg_parameters
    IMPLICIT NONE
    ! 0. Input
    REAL(wp),INTENT(in)       :: tau_labile_pot            ! potential labile pool turnover (1/day) 
    REAL(wp),INTENT(in)       :: gpp_week                  ! last weeks GPP (µmol CO2/m2/s) 
    REAL(wp),INTENT(in)       :: rm_week                   ! last weeks Rm (µmol CO2/m2/s)
    REAL(wp),INTENT(in)       :: lai                       ! current LAI 
    REAL(wp),INTENT(in)       :: target_lai                ! target LAI 
    REAL(wp),INTENT(in)       :: target_reserve_pool       ! desirable size (mol C/m2) of the reserve pool 
    REAL(wp),INTENT(in)       :: nreq_growth,preq_growth   ! N and P requirement for unit growth in C and N
    TYPE(Element),INTENT(in)  :: labile_pool               ! the plant's labile pool size
    TYPE(Element),INTENT(in)  :: reserve_pool              ! the plant's reserve pool size
    ! 0.2 Output
    TYPE(Element),INTENT(out) :: growth                    ! current growth rate (mol/time step)
    TYPE(Element),INTENT(out) :: reserve_use               ! current resource use rate (mol/time step) (+ve: net flux to labile
                                                           !                                            -ve: net flux to reserve
    ! 0.3 Local
    REAL(wp)                  :: tau_labile                ! actual turnover time of labile pool (1/day)
    REAL(wp)                  :: phi_phenology,phi_maintenance,phi_storage 
                                                           ! relative need for C for phenology, maintenance and storage 
    REAL(wp)                  :: target_labile_pool        ! desirable size (mol/m2) of the labile pool
    !
    ! 1. desired size of the labile pool, taken as the average input + maintenance requirements of the last week
    !    scaled to the default turnover time of the labile pool 
    !
    target_labile_pool= ( gpp_week + rm_week ) * one_day / 1.e6 / tau_labile_default
    !
    ! 2. mobilisation of reserve: calculate the pull from reserve to labile (phi_phenology,phi_maintenance) or
    !    push from labile to reserve (phi_storage) depending on the difference of the LAI/labile/reserve from
    !    their targets. phi_storage is adjusted by phi_phenology to avoid storage building to interfere with flushing 
    !
    phi_phenology = 0._wp
    phi_maintenance = 0._wp
    phi_storage = 0._wp
    IF(target_lai > eps8.AND.lai < target_lai ) & 
       phi_phenology = exp(-(lambda_phiphen * lai / target_lai )**weibull_k_phiphen) & 
                       - exp(-(lambda_phiphen)**weibull_k_phiphen) 
    IF(target_labile_pool > eps8 .AND. labile_pool%carbon < target_labile_pool) & 
       phi_maintenance = exp(-(lambda_phimaint * labile_pool%carbon / target_labile_pool)**weibull_k_phimaint) & 
                       - exp(-(lambda_phimaint)**weibull_k_phimaint)
    IF(target_reserve_pool > eps8 .AND. reserve_pool%carbon < target_reserve_pool) & 
       phi_storage = exp(-(lambda_phistore * reserve_pool%carbon / target_reserve_pool)**weibull_k_phistore) & 
                       - exp(-(lambda_phistore)**weibull_k_phistore)
    phi_storage = phi_storage * ( 1.0_wp - phi_phenology) 
    !
    ! 2.1 convert phi to turnover times (1/timestep), taking account that flushing can only happen if the meristem is active
    ! 
    phi_phenology = phi_phenology * tau_labile_pot * dt/one_day   
    phi_maintenance = phi_maintenance * tau_labile_default * dt/one_day
    phi_storage = phi_storage * tau_labile_default * dt/one_day
    !
    ! 2.2 net transfer of carbon from labile to reserve
    ! 
    reserve_use%carbon = ((phi_phenology + phi_maintenance) * reserve_pool%carbon - & 
                           phi_storage * labile_pool%carbon) 
    reserve_use%carbon13 = ((phi_phenology + phi_maintenance) * reserve_pool%carbon13 - & 
                             phi_storage * labile_pool%carbon13) 
    reserve_use%carbon14 = ((phi_phenology + phi_maintenance) * reserve_pool%carbon14 - & 
                             phi_storage * labile_pool%carbon14) 
    !
    ! 2.3 net transfer of N and P from labile to reserve (still missing)
    !
    ! 3.0 Elemental limitation of growth (this should be replaced at a later stage - here to conserve mass)
    !
    ! 3.1 N
    tau_labile=tau_labile_pot * dt / one_day
    growth%carbon = tau_labile * labile_pool%carbon
    growth%nitrogen = tau_labile * labile_pool%nitrogen
    IF(nreq_growth * growth%carbon > eps8)THEN
        IF(growth%nitrogen < nreq_growth * growth%carbon) & 
               tau_labile = tau_labile * growth%nitrogen / ( nreq_growth * growth%carbon ) 
    ELSE 
        tau_labile = zero
    ENDIF 

    ! 3.2 P
    !
    growth%nitrogen = tau_labile * labile_pool%nitrogen
    growth%phosphorus = tau_labile * labile_pool%phosphorus
    IF(preq_growth * growth%phosphorus > eps8)THEN
        IF(growth%phosphorus < preq_growth * growth%nitrogen ) & 
              tau_labile = tau_labile * growth%phosphorus / ( preq_growth * growth%nitrogen ) 
    ELSE 
        tau_labile = zero
    ENDIF 

    ! 4.0 Calculate final growth rate
    !
    growth%carbon = tau_labile * labile_pool%carbon
    growth%nitrogen = MIN(tau_labile * labile_pool%nitrogen, nreq_growth * growth%carbon)
    growth%phosphorus = MIN(tau_labile * labile_pool%phosphorus, preq_growth * growth%nitrogen)
    growth%carbon13 = tau_labile * labile_pool%carbon13
    growth%carbon14 = tau_labile * labile_pool%carbon14
    growth%nitrogen15 = tau_labile * labile_pool%nitrogen15

  END SUBROUTINE calc_partitioning 

  SUBROUTINE calc_allocation_fraction(tau_labile_pot,leaf_carbon,sap_wood_carbon,heart_wood_carbon,root_carbon,labile_carbon, &
                              sla,k_latosa,k_rtos,wood_density,allom_k1,allom_k2,tau_leaf,tau_fine_root,growthform, &
                              target_cn_leaf,target_cn_fine_root,target_cn_coarse_root, & 
                              target_cn_sap_wood,target_cn_fruit,target_np_leaf,target_np_fine_root, & 
                              target_np_coarse_root,target_np_sap_wood,target_np_fruit,falloc,nreq_growth,& 
                              preq_growth,target_lai,target_reserve_pool)
    !
    ! Routine to calculate allocation fractions for each element
    ! Input: target stoichiometry of the tissue types
    ! Output: 1) Allocation coefficients for the distribution of biomass for growth
    !         2) total N and P required to sustain a unit of growth (in mol)
    !         3) implied target LAI and reserve pool size
    ! 
    USE mo_land_constants
    USE mo_base_memory
    USE mo_base_constants, ONLY: pi
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: tau_labile_pot        ! potential labile pool turnover rate (1/day) without stoichiometric limit
    REAL(wp), INTENT(in) :: leaf_carbon           ! carbon in leaves (mol/individual)
    REAL(wp), INTENT(in) :: sap_wood_carbon       ! carbon in sap wood (mol/individual)
    REAL(wp), INTENT(in) :: heart_wood_carbon     ! carbon in heart wood (mol/individual)
    REAL(wp), INTENT(in) :: root_carbon           ! carbon in roots (mol/individual) 
    REAL(wp), INTENT(in) :: labile_carbon         ! carbon in labile pool (mol/individual) 
    REAL(wp), INTENT(in) :: sla                   ! specific leaf area
    REAL(wp), INTENT(in) :: k_latosa              ! leaf area to sapwood area ratio
    REAL(wp), INTENT(in) :: k_rtos                ! root tradeoff parameter
    REAL(wp), INTENT(in) :: wood_density          ! specific wood density 
    REAL(wp), INTENT(in) :: allom_k1              ! parameter 1 in height diameter relationship
    REAL(wp), INTENT(in) :: allom_k2              ! parameter 2 in height diameter relationship
    REAL(wp), INTENT(in) :: tau_leaf              ! default turnover time of leaves (1/year) 
    REAL(wp), INTENT(in) :: tau_fine_root         ! default turnover time of fine roots (1/year)
    INTEGER, INTENT(in) :: growthform             ! plant growth form 
    REAL(wp), INTENT(in) :: target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit        ! target C:N and N:P rations for tissue types
    ! 0.2 In- and Output
    TYPE(VegPool), INTENT(inout) :: falloc                              ! allocation coefficients 
    ! 0.3 Output
    REAL(wp), INTENT(out) :: nreq_growth, preq_growth                   ! mol N/P needed per mol growth of C
    REAL(wp), INTENT(out) :: target_lai                                 ! implied LAI from sapwood area
    REAL(wp), INTENT(out) :: target_reserve_pool                        ! desirable size of the reserve pool 
                                                                        ! (replace one year's growth of leaf and fine root)
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

    ! 
    ! seed for allocation coefficients is the potential growth rate per time step.
    ! given by the average turnover of the labile pool 
    b_inc = tau_labile_pot * labile_carbon * dt/one_day
    !
    IF(growthform==itree)THEN
       !
       ! Diagnose current height and diameter from wood mass assuming a zylindric shape and 
       ! height = k1 * D ^k2
       !  
       diameter = ((sap_wood_carbon + heart_wood_carbon) / & 
                  (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
       height = MAX(3.0_wp,allom_k1 * (diameter ** allom_k2))
       ! scaling factor for leaf to wood mass, derived from pipe-model (Shinozaki et al. 1964, J J Ecol)
       ! leaf area = k_latosa * sap_wood area
       ! => Cl * sla = k * Cs / wood_density / height
       ! => Cl = k_latosa / sla / wood_density / height * Cs = ltos * Cs / height
       ltos = k_latosa / sla / wood_density
       ! full leaf area implied by current sapwood area
       target_lai = (sap_wood_carbon + heart_wood_carbon) * k_latosa / wood_density / height
       ! scaling factor for root to leaf mass
       ! ad 1: derived from hydraulic homeostatis of root and wood hydraulic resistance (Magnani et al. 2000, PCE)
       ! Note: currently still lacks a moisture response!
       ! sapwood mass / height =  k_rtos * fine_root_mass
       ! Cr = Cs / k_rtos / height
       ! let Cs = Cl / ltos * height
       ! => Cl = Cr * ltos * k_rtos = Cl = Cr * ltor 
       ! ad 2: nitrogen/phosphorus sensitivity still needed
       ltor = k_rtos * ltos
       ! implied target of the reserve storage to regrow one year worth of folige and roots, accounting for 
       ! losses to respiration to grow these.
       target_reserve_pool = target_lai / sla * ( tau_leaf + tau_fine_root / ltor ) * (1._wp + fresp_growth)
       IF(b_inc > eps8) THEN
          b_inc = 0.2_wp * dt / one_day ! needed for stability - to be checked BLARPP! 
          ! derivative of height against biomass increase for approximation of height growth effect on partitioning
          ddiameter = ((sap_wood_carbon + heart_wood_carbon + b_inc) / &
                     (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))
          dheight = allom_k1 * ddiameter ** allom_k2 - height
   
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
          IF(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc>eps8) THEN 
             falloc%leaf%carbon = leaf_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
             falloc%fine_root%carbon = root_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
             falloc%coarse_root%carbon = 0.0_wp
             falloc%sap_wood%carbon = sap_wood_carbon_inc/(leaf_carbon_inc+sap_wood_carbon_inc+root_carbon_inc)
          ENDIF
       ENDIF
    ENDIF

    IF(b_inc > eps8) THEN
       !
       ! 2.0 Allocation factors for N, derived from the target C:N ratios and the C partitioning coefficients
       !  
       ! N required to sustain a unit of C growth
       nreq_growth = falloc%leaf%carbon / target_cn_leaf + & 
                     falloc%fine_root%carbon / target_cn_fine_root + &
                     falloc%coarse_root%carbon / target_cn_coarse_root + &
                     falloc%sap_wood%carbon / target_cn_sap_wood 
   
       ! Adjustment of N allocation factors
       IF(nreq_growth > eps8)THEN
          falloc%leaf%nitrogen = falloc%leaf%carbon / target_cn_leaf / nreq_growth
          falloc%fine_root%nitrogen = falloc%fine_root%carbon / target_cn_fine_root / nreq_growth 
          falloc%coarse_root%nitrogen = falloc%coarse_root%carbon / target_cn_coarse_root / nreq_growth
          falloc%sap_wood%nitrogen = falloc%sap_wood%carbon / target_cn_sap_wood / nreq_growth
       ENDIF
       !
       ! 3.0 Allocation factors for P, derived from the target N:P ratios and the N partitioning coefficients
       !  
       ! P required to sustain a unit of N growth
       preq_growth = falloc%leaf%nitrogen / target_np_leaf + & 
                     falloc%fine_root%nitrogen / target_np_fine_root + &
                     falloc%coarse_root%nitrogen / target_np_coarse_root + &
                     falloc%sap_wood%nitrogen / target_np_sap_wood 
   
       ! Adjustment of P allocation factors
       IF(preq_growth > eps8)THEN
          falloc%leaf%phosphorus = falloc%leaf%nitrogen / target_np_leaf / preq_growth
          falloc%fine_root%phosphorus = falloc%fine_root%nitrogen / target_np_fine_root / preq_growth 
          falloc%coarse_root%phosphorus = falloc%coarse_root%nitrogen / target_np_coarse_root / preq_growth
          falloc%sap_wood%phosphorus = falloc%sap_wood%nitrogen / target_np_sap_wood / preq_growth
       ENDIF
   
       !
       ! 4.0 Allocation factors for isotopes follow the main isotope 
       ! 
       falloc%leaf%carbon13 = falloc%leaf%carbon
       falloc%fine_root%carbon13 = falloc%fine_root%carbon
       falloc%coarse_root%carbon13 = falloc%coarse_root%carbon
       falloc%sap_wood%carbon13 = falloc%sap_wood%carbon
   
       falloc%leaf%carbon14 = falloc%leaf%carbon
       falloc%fine_root%carbon14 = falloc%fine_root%carbon
       falloc%coarse_root%carbon14 = falloc%coarse_root%carbon
       falloc%sap_wood%carbon14 = falloc%sap_wood%carbon
   
       falloc%leaf%nitrogen15 = falloc%leaf%nitrogen
       falloc%fine_root%nitrogen15 = falloc%fine_root%nitrogen
       falloc%coarse_root%nitrogen15 = falloc%coarse_root%nitrogen
       falloc%sap_wood%nitrogen15 = falloc%sap_wood%nitrogen

    ENDIF

  END SUBROUTINE calc_allocation_fraction

  SUBROUTINE update_stoichiometry(cn_leaf,np_leaf,target_cn_leaf,target_cn_fine_root,target_cn_coarse_root, & 
                              target_cn_sap_wood,target_cn_fruit,target_np_leaf,target_np_fine_root, & 
                              target_np_coarse_root,target_np_sap_wood,target_np_fruit)
    ! 
    ! Routine to calculate the target stoichiometry of the newly growth tissue
    ! for now: prescribed C:N:P - should be part of the optimality approach
    ! Output: target stoichiometries
    ! 
    USE mo_veg_parameters, ONLY: root2leaf_cn,root2leaf_np,wood2leaf_cn,wood2leaf_np
    IMPLICIT NONE
    ! 0.1: Input
    REAL(wp), INTENT(in) :: cn_leaf
    REAL(wp), INTENT(in) :: np_leaf
    ! 0.2: Output
    REAL(wp), INTENT(out) ::  target_cn_leaf,target_cn_fine_root,target_cn_coarse_root,target_cn_sap_wood, &
                              target_cn_fruit,target_np_leaf,target_np_fine_root,target_np_coarse_root, &
                              target_np_sap_wood,target_np_fruit
                              ! target C:N and N:P rations (Ops! this is in mol:mol)

    ! for now prescribed and based on OCN
    target_cn_leaf = cn_leaf
    target_cn_fine_root = target_cn_leaf / root2leaf_cn
    target_cn_coarse_root = target_cn_leaf / wood2leaf_cn
    target_cn_sap_wood = target_cn_leaf / wood2leaf_cn 
    target_cn_fruit = target_cn_leaf / root2leaf_cn

    ! Mock allocation fractions for P, assuming a univeral aim of N:P = 14 ("optimal") 
    target_np_leaf = np_leaf 
    target_np_fine_root = target_np_leaf / root2leaf_np
    target_np_coarse_root = target_np_leaf / wood2leaf_np
    target_np_sap_wood = target_np_leaf / wood2leaf_np
    target_np_fruit = target_np_leaf / root2leaf_np

  END SUBROUTINE update_stoichiometry

  PURE ELEMENTAL FUNCTION meristem_activity(t_air,theta) RESULT(tau_labile_pot)
    ! 
    ! This function calculates the activity of the meristem for growth. As this is proportional to the size
    ! of the labile pool, there is no need for a direct temperature function of the turnover, as this will 
    ! result from the dynamic of the labile pool. This function therefore only takes account of the impact
    ! of low temperatures (<5°C) and low moisture (<20% saturation). Additionally, phenological dormancy
    ! is accounted for (OBS: not yet implemented!)
    !
    USE mo_atmland_constants, ONLY: Tzero
    USE mo_veg_parameters, ONLY: lambda_tau_labile_temp,weibull_k_tau_labile_temp,& 
                                 lambda_tau_labile_theta,weibull_k_tau_labile_theta,tau_labile_default
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp),INTENT(in)     :: t_air                  ! air temperature in K 
    REAL(wp),INTENT(in)     :: theta                  ! soil moisture scalar (%) 
    ! 0.3 Output
    REAL(wp)                :: tau_labile_pot             ! current potential labile pool turnover (1/day)
    ! 0.4 Local
    REAL(wp)                :: ftemp                  ! temperature limit on cell construction rate
    REAL(wp)                :: fmoist                 ! moisture limit on cell construction rate
  
    IF(t_air > Tzero) THEN 
       ftemp = 1.0_wp - exp(-(lambda_tau_labile_temp * (t_air - Tzero))**weibull_k_tau_labile_temp)
    ELSE 
       ftemp = 0.0_wp
    ENDIF
    fmoist = 1.0_wp - exp(-(lambda_tau_labile_theta * theta)**weibull_k_tau_labile_theta) 
    tau_labile_pot = ftemp * fmoist * tau_labile_default

  END FUNCTION meristem_activity 

END MODULE mo_veg_growth
