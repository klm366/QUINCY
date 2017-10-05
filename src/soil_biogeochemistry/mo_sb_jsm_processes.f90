!>
!! @brief Contains the science routines for the JSM version in SLM 
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170124): adopted from JSM_interface and JSM_core
!!
MODULE mo_sb_jsm_processes

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps4, eps8
  USE mo_exception,       ONLY: message, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options
  USE mo_base_constants,  ONLY: dtime


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_transfer_litterfall_pools,calc_litter_turnover, calc_depolymerisation, & 
            calc_microbial_turnover,calc_microbial_growth, &
            calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax, &
            calc_om_sorption, &
            calc_inorganic_uptake, calc_po4_uptake_ECA, calc_flux_po4, calc_adsorption_po4

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_processes'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the transfer of vegetation litter fall to JSM litter pools
  !! still requires some thought as to the stoichiometry of the litter and the fraction of soluable litter 
  !! provides the litterfall interface between vegetation and soil_biogeochemistry
  !!
  !! Input: litterfall from vegetation, root distribution, soil layer depth
  !!
  !! Output: increase in soluable, polymeric, and woody litter
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_transfer_litterfall_pools(options,soil_depth_sl,root_fraction_sl,veg_litterfall, &
                                  formation_soluable_litter,formation_polymeric_litter,formation_woody_litter)
                                 
    USE mo_sb_constants
    dsl4slm_Use_memory_a_veg
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                    !< model options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: soil_depth_sl              !< depth of soil layer [m]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: root_fraction_sl           !< fraction of roots per layer 
    TYPE(t_veg_pool_a),                               INTENT(in)    :: veg_litterfall             !< litter input from vegetation [mol/m2/time step] 
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_soluable_litter
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_polymeric_litter
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_woody_litter
    ! ---------------------------
    ! 0.2 Local 
    INTEGER                                                 :: isoil
    REAL(wp), DIMENSION(options%nc)                         :: hlp1
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb)        :: hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_transfer_litterfall_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ... 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 pure soluable litter -> labile + reserve?
    !!
    formation_soluable_litter%carbon(:,1) = formation_soluable_litter%carbon(:,1) + &
        (veg_litterfall%labile%carbon + veg_litterfall%reserve%carbon) / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen(:,1) = formation_soluable_litter%nitrogen(:,1) + &
        (veg_litterfall%labile%nitrogen + veg_litterfall%reserve%nitrogen) / soil_depth_sl(:,1)
    formation_soluable_litter%phosphorus(:,1) = formation_soluable_litter%phosphorus(:,1) + &
        (veg_litterfall%labile%phosphorus + veg_litterfall%reserve%phosphorus) / soil_depth_sl(:,1)
    formation_soluable_litter%carbon13(:,1) = formation_soluable_litter%carbon13(:,1) + &
        (veg_litterfall%labile%carbon13 + veg_litterfall%reserve%carbon13) / soil_depth_sl(:,1)
    formation_soluable_litter%carbon14(:,1) = formation_soluable_litter%carbon14(:,1) + &
        (veg_litterfall%labile%carbon14 + veg_litterfall%reserve%carbon14) / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen15(:,1) = formation_soluable_litter%nitrogen15(:,1) + &
        (veg_litterfall%labile%nitrogen15 + veg_litterfall%reserve%nitrogen15) / soil_depth_sl(:,1)

    !> 2.0 mixed soluable and polymeric stuff
    !! 
    !> 2.1 leaves, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    hlp1 = f_soluable_leaf 
    formation_soluable_litter%carbon(:,1) = formation_soluable_litter%carbon(:,1) + &
        hlp1 * veg_litterfall%leaf%carbon / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon(:,1) = formation_polymeric_litter%carbon(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%carbon / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen(:,1) = formation_soluable_litter%nitrogen(:,1) + &
        hlp1 * veg_litterfall%leaf%nitrogen / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen(:,1) = formation_polymeric_litter%nitrogen(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%nitrogen / soil_depth_sl(:,1)
    formation_soluable_litter%phosphorus(:,1) = formation_soluable_litter%phosphorus(:,1) + &
        hlp1 * veg_litterfall%leaf%phosphorus / soil_depth_sl(:,1)
    formation_polymeric_litter%phosphorus(:,1) = formation_polymeric_litter%phosphorus(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%phosphorus / soil_depth_sl(:,1)
    formation_soluable_litter%carbon13(:,1) = formation_soluable_litter%carbon13(:,1) + &
        hlp1 * veg_litterfall%leaf%carbon13 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon13(:,1) = formation_polymeric_litter%carbon13(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%carbon13 / soil_depth_sl(:,1)
    formation_soluable_litter%carbon14(:,1) = formation_soluable_litter%carbon14(:,1) + &
        hlp1 * veg_litterfall%leaf%carbon14 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon14(:,1) = formation_polymeric_litter%carbon14(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%carbon14 / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen15(:,1) = formation_soluable_litter%nitrogen15(:,1) + &
        hlp1 * veg_litterfall%leaf%nitrogen15 / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen15(:,1) = formation_polymeric_litter%nitrogen15(:,1) + &
        (1._wp - hlp1) * veg_litterfall%leaf%nitrogen15 / soil_depth_sl(:,1)

    !> 2.2 fruit, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    !!
    hlp1 = f_soluable_fruit 
    formation_soluable_litter%carbon(:,1) = formation_soluable_litter%carbon(:,1) + &
        hlp1 * veg_litterfall%fruit%carbon / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon(:,1) = formation_polymeric_litter%carbon(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%carbon / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen(:,1) = formation_soluable_litter%nitrogen(:,1) + &
        hlp1 * veg_litterfall%fruit%nitrogen / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen(:,1) = formation_polymeric_litter%nitrogen(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%nitrogen / soil_depth_sl(:,1)
    formation_soluable_litter%phosphorus(:,1) = formation_soluable_litter%phosphorus(:,1) + &
        hlp1 * veg_litterfall%fruit%phosphorus / soil_depth_sl(:,1)
    formation_polymeric_litter%phosphorus(:,1) = formation_polymeric_litter%phosphorus(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%phosphorus / soil_depth_sl(:,1)
    formation_soluable_litter%carbon13(:,1) = formation_soluable_litter%carbon13(:,1) + &
        hlp1 * veg_litterfall%fruit%carbon13 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon13(:,1) = formation_polymeric_litter%carbon13(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%carbon13 / soil_depth_sl(:,1)
    formation_soluable_litter%carbon14(:,1) = formation_soluable_litter%carbon14(:,1) + &
        hlp1 * veg_litterfall%fruit%carbon14 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon14(:,1) = formation_polymeric_litter%carbon14(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%carbon14 / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen15(:,1) = formation_soluable_litter%nitrogen15(:,1) + &
        hlp1 * veg_litterfall%fruit%nitrogen15 / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen15(:,1) = formation_polymeric_litter%nitrogen15(:,1) + &
        (1._wp - hlp1) * veg_litterfall%fruit%nitrogen15 / soil_depth_sl(:,1)

    !> 2.3 seed_bed, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    !!
    hlp1 = f_soluable_seed_bed 
    formation_soluable_litter%carbon(:,1) = formation_soluable_litter%carbon(:,1) + &
        hlp1 * veg_litterfall%seed_bed%carbon / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon(:,1) = formation_polymeric_litter%carbon(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%carbon / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen(:,1) = formation_soluable_litter%nitrogen(:,1) + &
        hlp1 * veg_litterfall%seed_bed%nitrogen / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen(:,1) = formation_polymeric_litter%nitrogen(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%nitrogen / soil_depth_sl(:,1)
    formation_soluable_litter%phosphorus(:,1) = formation_soluable_litter%phosphorus(:,1) + &
        hlp1 * veg_litterfall%seed_bed%phosphorus / soil_depth_sl(:,1)
    formation_polymeric_litter%phosphorus(:,1) = formation_polymeric_litter%phosphorus(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%phosphorus / soil_depth_sl(:,1)
    formation_soluable_litter%carbon13(:,1) = formation_soluable_litter%carbon13(:,1) + &
        hlp1 * veg_litterfall%seed_bed%carbon13 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon13(:,1) = formation_polymeric_litter%carbon13(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%carbon13 / soil_depth_sl(:,1)
    formation_soluable_litter%carbon14(:,1) = formation_soluable_litter%carbon14(:,1) + &
        hlp1 * veg_litterfall%seed_bed%carbon14 / soil_depth_sl(:,1)
    formation_polymeric_litter%carbon14(:,1) = formation_polymeric_litter%carbon14(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%carbon14 / soil_depth_sl(:,1)
    formation_soluable_litter%nitrogen15(:,1) = formation_soluable_litter%nitrogen15(:,1) + &
        hlp1 * veg_litterfall%seed_bed%nitrogen15 / soil_depth_sl(:,1)
    formation_polymeric_litter%nitrogen15(:,1) = formation_polymeric_litter%nitrogen15(:,1) + &
        (1._wp - hlp1) * veg_litterfall%seed_bed%nitrogen15 / soil_depth_sl(:,1)

    !> 3.0 roots
    !!
    DO isoil = 1,options%nsoil_sb
       !> 3.1 fine roots
       !! still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
       hlp1 = f_soluable_fine_root 
       WHERE(soil_depth_sl(:,isoil) > eps8)
          formation_soluable_litter%carbon(:,isoil) = formation_soluable_litter%carbon(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon / soil_depth_sl(:,isoil)
          formation_polymeric_litter%carbon(:,isoil) = formation_polymeric_litter%carbon(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon / soil_depth_sl(:,isoil)
          formation_soluable_litter%nitrogen(:,isoil) = formation_soluable_litter%nitrogen(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%nitrogen / soil_depth_sl(:,isoil)
          formation_polymeric_litter%nitrogen(:,isoil) = formation_polymeric_litter%nitrogen(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%nitrogen / soil_depth_sl(:,isoil)
          formation_soluable_litter%phosphorus(:,isoil) = formation_soluable_litter%phosphorus(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%phosphorus / soil_depth_sl(:,isoil)
          formation_polymeric_litter%phosphorus(:,isoil) = formation_polymeric_litter%phosphorus(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%phosphorus / soil_depth_sl(:,isoil)
          formation_soluable_litter%carbon13(:,isoil) = formation_soluable_litter%carbon13(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon13 / soil_depth_sl(:,isoil)
          formation_polymeric_litter%carbon13(:,isoil) = formation_polymeric_litter%carbon13(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon13 / soil_depth_sl(:,isoil)
          formation_soluable_litter%carbon14(:,isoil) = formation_soluable_litter%carbon14(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon14 / soil_depth_sl(:,isoil)
          formation_polymeric_litter%carbon14(:,isoil) = formation_polymeric_litter%carbon14(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%carbon14 / soil_depth_sl(:,isoil)
          formation_soluable_litter%nitrogen15(:,isoil) = formation_soluable_litter%nitrogen15(:,isoil) + &
              hlp1 * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%nitrogen15 / soil_depth_sl(:,isoil)
          formation_polymeric_litter%nitrogen15(:,isoil) = formation_polymeric_litter%nitrogen15(:,isoil) + &
              (1._wp - hlp1) * root_fraction_sl(:,isoil) * veg_litterfall%fine_root%nitrogen15 / soil_depth_sl(:,isoil)

          !> 3.2 woody litter below ground 
          !!
          formation_woody_litter%carbon(:,isoil) = formation_woody_litter%carbon(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%carbon / soil_depth_sl(:,isoil)
          formation_woody_litter%nitrogen(:,isoil) = formation_woody_litter%nitrogen(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%nitrogen / soil_depth_sl(:,isoil)
          formation_woody_litter%phosphorus(:,isoil) = formation_woody_litter%phosphorus(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%phosphorus / soil_depth_sl(:,isoil)
          formation_woody_litter%carbon13(:,isoil) = formation_woody_litter%carbon13(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%carbon13 / soil_depth_sl(:,isoil)
          formation_woody_litter%carbon14(:,isoil) = formation_woody_litter%carbon14(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%carbon14 / soil_depth_sl(:,isoil)
          formation_woody_litter%nitrogen15(:,isoil) = formation_woody_litter%nitrogen15(:,isoil) + & 
             root_fraction_sl(:,isoil) * veg_litterfall%coarse_root%nitrogen15 / soil_depth_sl(:,isoil)
       ENDWHERE
    ENDDO

    !> 4.0 woody litter above ground 
    !!
    formation_woody_litter%carbon(:,1) = formation_woody_litter%carbon(:,1) + & 
        (veg_litterfall%sap_wood%carbon + veg_litterfall%heart_wood%carbon) / soil_depth_sl(:,1)
    formation_woody_litter%nitrogen(:,1) = formation_woody_litter%nitrogen(:,1) + & 
        (veg_litterfall%sap_wood%nitrogen + veg_litterfall%heart_wood%nitrogen) / soil_depth_sl(:,1)
    formation_woody_litter%phosphorus(:,1) = formation_woody_litter%phosphorus(:,1) + & 
        (veg_litterfall%sap_wood%phosphorus + veg_litterfall%heart_wood%phosphorus) / soil_depth_sl(:,1)
    formation_woody_litter%carbon13(:,1) = formation_woody_litter%carbon13(:,1) + & 
        (veg_litterfall%sap_wood%carbon13 + veg_litterfall%heart_wood%carbon13) / soil_depth_sl(:,1)
    formation_woody_litter%carbon14(:,1) = formation_woody_litter%carbon14(:,1) + & 
        (veg_litterfall%sap_wood%carbon14 + veg_litterfall%heart_wood%carbon14) / soil_depth_sl(:,1)
    formation_woody_litter%nitrogen15(:,1) = formation_woody_litter%nitrogen15(:,1) + & 
        (veg_litterfall%sap_wood%nitrogen15 + veg_litterfall%heart_wood%nitrogen15) / soil_depth_sl(:,1)

  END SUBROUTINE calc_transfer_litterfall_pools 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  ! 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of litter pools according to first order decay 
  !! provides a link between the soluable and woody litter pools introduced in SLM, and the
  !! DOM and polymeric pools of JSM
  !!
  !! Input: temperature and moisture dependence of litter decomposition, litter pool size, and its
  !!        turnover time at standard conditions 
  !! Output: transfer of litter to its sink pool (soluable -> DOM; woody -> polymeric litter) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_turnover(options               , &
                                  tau_litter_pool       , &
                                  rtm_depolymerisation  , &
                                  rmm_depolymerisation  , & 
                                  pool_litter           , &
                                  loss_litter           , &
                                  formation_sink        )

    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options               !< model options
    REAL(wp),                                         INTENT(in)    :: tau_litter_pool       !< turnover time of the litter pool [1/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier [-]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_litter           !< litter pool to decompose [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_litter           !< reduction in litter pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_sink        !< formation of SOM pool [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1 , &
                                                        decay_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ..


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calulate decay per timestep
    !!
    decay_rate(:,:) = 1._wp/tau_litter_pool * dtime

    !> 2.0 calculate source and sink terms
    !! assuming to change in composition or mass of the material to decompose
    !! 
    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon(:,:)
    loss_litter%carbon(:,:)         = loss_litter%carbon(:,:) + hlp1(:,:)
    formation_sink%carbon(:,:)      = formation_sink%carbon(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen(:,:)
    loss_litter%nitrogen(:,:)       = loss_litter%nitrogen(:,:) + hlp1(:,:)
    formation_sink%nitrogen(:,:)    = formation_sink%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%phosphorus(:,:)
    loss_litter%phosphorus(:,:)     = loss_litter%phosphorus(:,:) + hlp1(:,:)
    formation_sink%phosphorus(:,:)  = formation_sink%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon13(:,:)
    loss_litter%carbon13(:,:)       = loss_litter%carbon13(:,:) + hlp1(:,:)
    formation_sink%carbon13(:,:)    = formation_sink%carbon13(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%carbon14(:,:)
    loss_litter%carbon14(:,:)       = loss_litter%carbon14(:,:) + hlp1(:,:)
    formation_sink%carbon14(:,:)    = formation_sink%carbon14(:,:) + hlp1(:,:) 

    hlp1(:,:)                       = decay_rate(:,:) * pool_litter%nitrogen15(:,:)
    loss_litter%nitrogen15(:,:)     = loss_litter%nitrogen15(:,:) + hlp1(:,:)
    formation_sink%nitrogen15(:,:)  = formation_sink%nitrogen15(:,:) + hlp1(:,:) 

  END SUBROUTINE calc_litter_turnover

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate sinks of the bio-available inorganic nutrient pools
  !!  according to Zhu et al. (2017)
  !!
  !! Input: solute concentration, microbial biomass, fine root biomass, rate modifiers 
  !!
  !! Output: microbial and plant inorganic uptake rates
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_inorganic_uptake(solute, microbial_carbon, fine_root_carbon, &
                                             rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, & 
                                             vmax_mic_uptake, km_mic_uptake, vmax_plant_uptake, km_plant_uptake, & 
                                             mic_uptake, plant_uptake)

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)            :: solute                  !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp), INTENT(in)            :: microbial_carbon        !< microbial carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: fine_root_carbon        !< fine_root density (mol / m3)
    REAL(wp), INTENT(in)            :: rtm_mic_uptake          !< temperature rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rmm_mic_uptake          !< moisture rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rtm_hsc                 !< temperature rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rmm_hsc                 !< moisture rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: vmax_mic_uptake         !< maximium uptake rate of microbial biomass (µmol S / mol C / s)
    REAL(wp), INTENT(in)            :: km_mic_uptake           !< half-saturation constant of microbial nutrient uptake (mol S / m3)
    REAL(wp), INTENT(in)            :: vmax_plant_uptake       !< maximium uptake rate of microbial biomass (µmol S / mol C / s)
    REAL(wp), INTENT(in)            :: km_plant_uptake         !< half-saturation constant of microbial nutrient uptake (mol S / m3)
    REAL(wp), INTENT(out)           :: mic_uptake              !< microbial uptake flux (µmol / m3 / s)
    REAL(wp), INTENT(out)           :: plant_uptake            !< plant uptake flux (µmol / m3 / s)
    ! 0.2 Local
    REAL(wp)              :: hlp1,hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_uptake'

    !> 1.0 uptake
    !! co-limited rates of N uptake for microbes and plants 
    hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake
    hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake
    mic_uptake     = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake * microbial_carbon * solute / &
                         (hlp1 + solute + hlp1/hlp2 * fine_root_carbon + microbial_carbon)
    plant_uptake   = rtm_mic_uptake * rmm_mic_uptake * vmax_plant_uptake * fine_root_carbon * solute / &
                         (hlp2 + solute + hlp2/hlp1 * microbial_carbon + fine_root_carbon)

    !> 1.1 limit to possible solution
    !!
    hlp1 = (mic_uptake + plant_uptake) 
    IF(hlp1 > (solute * 1000000._wp/dtime) .AND. hlp1 > eps8)THEN
       mic_uptake   = mic_uptake   / hlp1 * MAX(solute-eps8,0.0_wp) * 1000000.0_wp
       plant_uptake = plant_uptake / hlp1 * MAX(solute-eps8,0.0_wp) * 1000000.0_wp
    ENDIF

  END SUBROUTINE calc_inorganic_uptake

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from sb_jsm_main
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate the flux rates of P processes
  !! weathering and occlusion according to Wang et al. 2010; 
  !! desorption (slow exchange) rate based on the ANIMO model (Groenendijk et al. 2005)
  !! biocehmical mineralization is calculated using MM-kinetics assuming Mic C is the limiting factor
  !! Input: rate modifiers and source pools
  !!
  !! Output: biochemical P mineralization rate and loss of source pools, gross P mineralization, 
  !!         weathering, occlusion and (slow) desorption
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_flux_po4(microbial_carbon, dom_carbon, dom_phosphorus, &
                                            rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &
                                            polymeric_litter_po4, dom_assoc_po4, &
                                            residue_po4, residue_assoc_po4, &
                                            loss_polymeric_litter_po4, loss_dom_assoc_po4, &
                                            loss_residue_po4, loss_residue_assoc_po4, assoc_slow, &
                                            biochem_min_po4, gross_mineralisation_po4, weathering_po4, occlusion_po4, slow_sorption_po4)


    USE mo_sb_constants

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)            :: microbial_carbon          !< microbial carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: dom_carbon                !< dom carbon density (mol / m3)
    REAL(wp), INTENT(in)            :: dom_phosphorus            !< dom phosphorus density (mol / m3)
    REAL(wp), INTENT(in)            :: rtm_mic_uptake            !< temperature rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rmm_mic_uptake            !< moisture rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)            :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
    REAL(wp), INTENT(in)            :: polymeric_litter_po4      !< PO4 in the polymeric litter pool [mol/m3]
    REAL(wp), INTENT(in)            :: dom_assoc_po4             !< PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), INTENT(in)            :: residue_po4               !< PO4 in the microbial necromass pool [mol/m3]
    REAL(wp), INTENT(in)            :: residue_assoc_po4         !< PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), INTENT(in)            :: assoc_slow                !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), INTENT(inout)         :: loss_polymeric_litter_po4 !< loss of PO4 in the polymeric litter pool [mol/m3]
    REAL(wp), INTENT(inout)         :: loss_dom_assoc_po4        !< loss of PO4 in the mineral associated dom pool [mol/m3]
    REAL(wp), INTENT(inout)         :: loss_residue_po4          !< loss of PO4 in the microbial necromass pool [mol/m3]
    REAL(wp), INTENT(inout)         :: loss_residue_assoc_po4    !< loss of PO4 in the mineral associated microbial necromass pool [mol/m3]
    REAL(wp), INTENT(out)           :: biochem_min_po4           !< phosphorus biochemical mineralization rate, [µmol P /m3 /s]
    REAL(wp), INTENT(out)           :: gross_mineralisation_po4  !< phosphorus gross(biological) mineralization rate, [µmol P /m3 /s]
    REAL(wp), INTENT(out)           :: weathering_po4            !< phosphorus weathering rate [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: occlusion_po4             !< phosphorus occlusion rate, assoc_slow --> occluded [µmol P / m3 / s]
    REAL(wp), INTENT(out)           :: slow_sorption_po4             !< phosphorus desorption rate, assoc_fast --> assoc_slow [µmol P / m3 / s]
    ! 0.2 Local
    REAL(wp)                        :: hlp1, hlp2, uptake_c, f_uptake, uptake_p
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_flux_po4'
    ! ---------------------------


    !> 1.0 Calculates the biochemical P mineralization, and the loss of the contributing pools.
    !! P is assumed to clevage out from poly litter, dom associated, residue and residue associated 
    hlp1 = (polymeric_litter_po4 + dom_assoc_po4 + residue_po4 + residue_assoc_po4)
    biochem_min_po4 = rtm_mic_uptake * rmm_mic_uptake * vmax_biochem_po4 * hlp1 * microbial_carbon / (km_mic_biochem_po4 + microbial_carbon)
    IF (hlp1 > eps8) THEN
        hlp2 = polymeric_litter_po4 / hlp1
        loss_polymeric_litter_po4 = biochem_min_po4 * hlp2 * dtime
        hlp2 = dom_assoc_po4 / hlp1
        loss_dom_assoc_po4 = biochem_min_po4 * hlp2 * dtime
        hlp2 = residue_po4 / hlp1
        loss_residue_po4 = biochem_min_po4 * hlp2 * dtime
        hlp2 = residue_assoc_po4 / hlp1
        loss_residue_assoc_po4 = biochem_min_po4 * hlp2 * dtime
    ELSE
        loss_polymeric_litter_po4 = 0.0_wp
        loss_dom_assoc_po4        = 0.0_wp
        loss_residue_po4          = 0.0_wp
        loss_residue_assoc_po4    = 0.0_wp
    ENDIF
    biochem_min_po4 = biochem_min_po4 * 1.e6_wp

    !> 2.0 calculate the gross P mineralization, based on the potential depolymerization of DOM
    !!
    hlp1 = dom_carbon / (km_mic_uptake * rtm_hsc * rmm_hsc + dom_carbon )
    uptake_c = vmax_mic_uptake * rtm_mic_uptake * rmm_mic_uptake * hlp1 * microbial_carbon
    IF(dom_carbon > eps8) THEN
       f_uptake = uptake_c / dom_carbon
    ELSE
       f_uptake = 0.0_wp
    ENDIF
    uptake_p = f_uptake * dom_phosphorus
    gross_mineralisation_po4 = (1._wp - microbial_pue) * uptake_p

    !> 3.0 Pi fluxes between different mineral P pools
    !! weathering, occlusion and slow exchange between sorbed and labile P pools
    weathering_po4 = rtm_hsc * rmm_hsc * k_weath_mineral * 1.e6_wp
    occlusion_po4 = assoc_slow * k_occlude_po4 * 1.e6_wp
    slow_sorption_po4 = k_desorpt_po4 * (eq_sorp_mineral_po4 - assoc_slow) * 1.e6_wp

  END SUBROUTINE calc_flux_po4

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate sinks of the bio-available inorganic po4 pools
  !! ECA approach based on Zhu et al. (2016) and Zhu et al. (2017)
  !!
  !! Input: Pi solute concentration, microbial biomass, fine root biomass, soil Pi pools, rate modifiers
  !!
  !! Output: microbial and plant inorganic uptake rates
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_po4_uptake_ECA(  microbial_carbon, fine_root_carbon, solute, &
                                             rtm_mic_uptake, rmm_mic_uptake, rtm_hsc, rmm_hsc, &
                                             qmax_po4, assoc_fast, mic_uptake, plant_uptake, flag_sb_ECA_upt_po4)

    USE mo_sb_constants

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)    :: microbial_carbon          !< microbial carbon density (mol / m3)
    REAL(wp), INTENT(in)    :: fine_root_carbon          !< fine_root density (mol / m3)
    REAL(wp), INTENT(in)    :: rtm_mic_uptake            !< temperature rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)    :: rmm_mic_uptake            !< moisture rate modifier for microbial (and plant) uptake
    REAL(wp), INTENT(in)    :: rtm_hsc                   !< temperature rate modifier for half_saturation constants
    REAL(wp), INTENT(in)    :: rmm_hsc                   !< moisture rate modifier for half_saturation constants
    REAL(wp), INTENT(in)    :: solute                    !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp), INTENT(in)    :: qmax_po4                  !< sorption capacity of po4 to mineral soil [mol P / m3]
    REAL(wp), INTENT(in)    :: assoc_fast                !< fast minerally associated PO4 pool [mol/m3]
    LOGICAL,  INTENT(in)    :: flag_sb_ECA_upt_po4       !< uptake option for P, TRUE: ECA; FALSE: non-ECA
    REAL(wp), INTENT(out)   :: mic_uptake                !< microbial uptake flux (µmol / m3 / s)
    REAL(wp), INTENT(out)   :: plant_uptake              !< plant uptake flux (µmol / m3 / s)


    ! 0.2 Local
    REAL(wp)                                             :: hlp1, hlp2, avail_sorption_sites
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_po4_uptake_ECA'
    ! ---------------------------


    !> define uptake scheme to use, 1. ECA; 2. non-ECA
    IF (flag_sb_ECA_upt_po4) THEN
    !> 1.0 uptake
    !! co-limited rates of P uptake for microbes and plants

      !> 1.1 ECA approach
      !! 
        avail_sorption_sites = qmax_po4 - assoc_fast 
        hlp1 = km_mic_uptake_po4 + solute + f_mic_carrier_po4 * microbial_carbon + &
               km_mic_uptake_po4 / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon + &
               km_mic_uptake_po4 / km_adsorpt_po4 * avail_sorption_sites
        hlp2 = km_plant_uptake_po4 + solute + f_plant_carrier_po4 * fine_root_carbon + &
               km_plant_uptake_po4 / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon + &
               km_plant_uptake_po4 / km_adsorpt_po4 * avail_sorption_sites

        mic_uptake     = rtm_mic_uptake * rmm_mic_uptake * k_mic_uptake_po4 * solute * f_mic_carrier_po4 * &
                          microbial_carbon / hlp1 * 1.e6_wp
        plant_uptake   = rtm_mic_uptake * rmm_mic_uptake * k_plant_uptake_po4 * solute * f_plant_carrier_po4 * &
                          fine_root_carbon / hlp2 * 1.e6_wp

    ELSE
      !> 1.2 non-ECA approach
      !! 
        hlp1 = rtm_hsc * rmm_hsc * km_mic_uptake_po4
        hlp2 = rtm_hsc * rmm_hsc * km_plant_uptake_po4
        mic_uptake     = rtm_mic_uptake * rmm_mic_uptake * vmax_mic_uptake_po4 * microbial_carbon * solute / &
                         (hlp1 + solute + hlp1/hlp2 * fine_root_carbon + microbial_carbon)
        plant_uptake   = rtm_mic_uptake * rmm_mic_uptake * vmax_plant_uptake_po4 * fine_root_carbon * solute / &
                         (hlp2 + solute + hlp2/hlp1 * microbial_carbon + fine_root_carbon)

    ENDIF
  END SUBROUTINE calc_po4_uptake_ECA


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate the adsorption rate of po4 based on
  !!      1)the ECA equation for adsorption (Zhu et al. 2016); 
  !!        the potential adsorption is defined by Langmuir kinetics, Van der Zee et al. 1989
  !!      2)the Langmuir kinetics, Van der Zee et al. 1989, where
  !!        the Langmuir equilibrium is either maintained or not
  !!      3)the Langmuir equilibrium assumption according to Yang et al. 2014
  !!       
  !! Input: uptake, immobilization, desorption, occlusion, weathering rates; po4 pool sizes
  !!
  !! Output: adsorption rate
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_adsorption_po4(options, solute, assoc_fast, assoc_slow, occluded, qmax_po4, &
                                        weathering_po4, occlusion_po4, slow_sorption_po4, percolation_sl, bulk_dens_corr_sl, &
                                        gross_mineralisation_po4, biochem_min_po4, microbial_carbon, fine_root_carbon, &
                                        mic_uptake, plant_uptake, fast_adsorption_po4, &
                                        sb_adsorp_scheme, flag_sb_CLE_po4, flag_sb_LK_po4)

    USE mo_sb_constants

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options), INTENT(in)                             :: options           !< model options
    CHARACTER(len=*), INTENT(in)                                     :: sb_adsorp_scheme  !< eca_full; 2: eca_part (only uptake, not for adsorption); 3: eca_none
    LOGICAL,          INTENT(in)                                     :: flag_sb_CLE_po4   !< T: Constant Langmuir Equilibrium; F: non-Equilibrium
    LOGICAL,          INTENT(in)                                     :: flag_sb_LK_po4    !< T: Langmuir kinetics for adsorption; F: adsorption calculated based on other fluxes
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: solute            !< concentration of the solvent to be taken up (mol / m3)
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: assoc_fast        !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: assoc_slow        !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: occluded          !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: qmax_po4          !< sorption capacity of po4 to mineral soil [mol P / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: percolation_sl    !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: bulk_dens_corr_sl !< bulk density of the mineral soil [kg / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: microbial_carbon  !< microbial carbon density (mol / m3)
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: fine_root_carbon  !< fine_root density (mol / m3)
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: gross_mineralisation_po4     !< phosphorus biochemical mineralization rate, [µmol P /m3 /s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: biochem_min_po4   !< phosphorus biochemical mineralization rate, [µmol P /m3 /s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: weathering_po4    !< phosphorus weathering rate [µmol P / m3 / s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: occlusion_po4     !< phosphorus occlusion rate, assoc_slow --> occluded [µmol P / m3 / s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)     :: slow_sorption_po4 !< phosphorus desorption rate, assoc_fast --> assoc_slow [µmol P / m3 / s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)    :: fast_adsorption_po4 !< phosphorus adsorption rate, solute --> assoc_fast [µmol P / m3 / s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout)  :: mic_uptake        !< microbial uptake flux (µmol / m3 / s)
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout)  :: plant_uptake      !< plant uptake flux (µmol / m3 / s)



    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb)                :: hlp1,hlp2,hlp3,hlp4,hlp5,hlp6,hlp9,scal
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb)                :: perco_po4, desorp_f_po4, adsorpt_f_po4, &
                                                                       lab_po4, avail_sorption_sites
    INTEGER                                                         :: isoil

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_adsorption_po4'
    ! ---------------------------

    !> 0.3 percolation flux and available sorption sites for po4
    perco_po4(:,1) = - percolation_sl(:,1) * solute(:,1)/1000_wp
    DO isoil = 2 , options%nsoil_sb
        perco_po4(:,isoil) = percolation_sl(:,isoil-1) /1000_wp * solute(:,isoil-1) - percolation_sl(:,isoil) / &
                              1000_wp * solute(:,isoil)
    ENDDO
    avail_sorption_sites = qmax_po4 - assoc_fast


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 ECA option for adsorption rate
    SELECT CASE (sb_adsorp_scheme)
      CASE("eca_full")
    !> 2.0 option 1: full ECA approach
    !!
        IF (flag_sb_CLE_po4) THEN
        !> 2.1 By default, full ECA only applies to Langmuir kinetics; option with Constant Equilibrium
        !! potential adsorption is defined as the fast adsorption rate of Langmuir kinetics in Van der Zee et al. 1989
        !! constant Langmuir equilibrium is maintained by FUNCTION calc_CLE_po4
            hlp4 = km_adsorpt_po4 * bulk_dens_corr_sl
            hlp5 = avail_sorption_sites
            hlp2 = hlp4 + solute + hlp5 + &
                   hlp4 / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon + &
                   hlp4 / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon
            hlp3 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            adsorpt_f_po4 = hlp3 * solute / hlp2
            desorp_f_po4 = assoc_fast * k_desorp_f_po4

            hlp6 = ((weathering_po4 + biochem_min_po4 + gross_mineralisation_po4)/1.e6_wp + desorp_f_po4 + perco_po4(:,:)) * dtime

            hlp9 = ((mic_uptake + plant_uptake)/1.e6_wp + adsorpt_f_po4) * dtime

        !! Langmuir equilibrium considering possible limiting sources 
            scal = 1.0_wp
            lab_po4 = 0.0_wp
            WHERE ((solute + hlp6 - hlp9) < eps8)
                hlp1 = assoc_fast + (- desorp_f_po4 - slow_sorption_po4/1.e6_wp)*dtime
                scal = calc_CLE_po4("insufficient_p_input",solute, hlp4, hlp5, hlp1, hlp6, hlp9, adsorpt_f_po4*dtime)
                adsorpt_f_po4 = adsorpt_f_po4 * scal
                fast_adsorption_po4 = adsorpt_f_po4 - desorp_f_po4
                mic_uptake = mic_uptake * scal
                plant_uptake = plant_uptake * scal
            ELSEWHERE
                hlp1 = assoc_fast + (- desorp_f_po4 - slow_sorption_po4/1.e6_wp)*dtime
                scal = calc_CLE_po4("sufficient_p_input",solute, hlp4, hlp5, hlp1, hlp6, hlp9, adsorpt_f_po4*dtime)
                adsorpt_f_po4 = adsorpt_f_po4 * scal
                fast_adsorption_po4 = adsorpt_f_po4 - desorp_f_po4
                mic_uptake = mic_uptake * scal
                plant_uptake = plant_uptake * scal
                scal = -scal -1._wp
            ENDWHERE
        ELSE
        !> 2.2 By default, full ECA only applies to Langmuir kinetics; option with non-Equilibrium
        !! Langmuir equilibrium is not constantly maintained
            hlp4 = km_adsorpt_po4 * bulk_dens_corr_sl
            hlp5 = avail_sorption_sites
            hlp2 = hlp4 + solute + hlp5 + &
                   hlp4 / km_plant_uptake_po4 * f_plant_carrier_po4 * fine_root_carbon + &
                   hlp4 / km_mic_uptake_po4 * f_mic_carrier_po4 * microbial_carbon
            hlp3 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            adsorpt_f_po4 = hlp3 * solute / hlp2 
            desorp_f_po4 = assoc_fast * k_desorp_f_po4

        !> 2.2.1 limit to possible P sources
        !!
            hlp1 = (mic_uptake + plant_uptake)/1.e6_wp + adsorpt_f_po4
            hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4)/1.e6_wp + desorp_f_po4 + perco_po4(:,:)

            WHERE((hlp1 + solute / dtime) > hlp6 .AND. solute > eps8)
               mic_uptake    = mic_uptake    / hlp1 * hlp6
               plant_uptake  = plant_uptake  / hlp1 * hlp6
               adsorpt_f_po4 = adsorpt_f_po4 / hlp1 * hlp6
            ENDWHERE
        ENDIF




      CASE("eca_part")
      !> 3.0 option 2: partial-ECA, adsorption is not determined by ECA
      !!
        IF (flag_sb_CLE_po4 .AND. flag_sb_LK_po4) THEN
        !> 3.1 option with Langmuir kinetics and constant equilibrium
        !! @to be finished...
            hlp5 = avail_sorption_sites
            adsorpt_f_po4 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            desorp_f_po4 = assoc_fast * k_desorp_f_po4


        ELSEIF (flag_sb_CLE_po4 .AND. (.NOT. flag_sb_LK_po4)) THEN
        !> 3.2 option with constant equilibrium, i.e. the same as Wang et al. 2007 and so forth
        !! 
            hlp4 = km_adsorpt_po4 * bulk_dens_corr_sl
            hlp2 = qmax_po4 * hlp4 / ((solute + hlp4)**2.0_wp + qmax_po4 * hlp4)
            hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4 - &
                    (mic_uptake + plant_uptake + slow_sorption_po4))/1.e6_wp + perco_po4(:,:)
            WHERE (hlp6 * dtime + solute + assoc_fast < eps8 .AND. solute > eps8)
                hlp1 = mic_uptake + plant_uptake
                hlp3 = (hlp6 + hlp1) * dtime + solute + assoc_fast
                mic_uptake = mic_uptake / hlp1 * hlp3
                plant_uptake = plant_uptake / hlp1 * hlp3
                hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4 - &
                        (mic_uptake + plant_uptake + slow_sorption_po4))/1.e6_wp + perco_po4(:,:)
            ENDWHERE
            fast_adsorption_po4 = hlp2 * hlp6 + slow_sorption_po4

        ELSEIF (flag_sb_LK_po4 .AND. (.NOT. flag_sb_CLE_po4)) THEN
        !> 3.3 option with Langmuir kinetics but not constant equilibrium
        !!
            hlp5 = avail_sorption_sites
            adsorpt_f_po4 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            desorp_f_po4 = assoc_fast * k_desorp_f_po4

        !> 3.3.1 limit to possible solution (P sources)
            hlp1 = (mic_uptake + plant_uptake)/1.e6_wp + adsorpt_f_po4
            hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4)/1.e6_wp + desorp_f_po4 + perco_po4(:,:)

            WHERE((hlp1 + solute / dtime) > hlp6 .AND. solute > eps8)
               mic_uptake    = mic_uptake    / hlp1 * hlp6
               plant_uptake  = plant_uptake  / hlp1 * hlp6
               adsorpt_f_po4 = adsorpt_f_po4 / hlp1 * hlp6
            ENDWHERE
        ENDIF

      CASE("eca_none")
      !> 4.0 option 3: non-ECA, this only applies to upt_opt=2 in SUB calc_po4_uptake_ECA
      !! 
        IF (flag_sb_CLE_po4 .AND. flag_sb_LK_po4) THEN
        !> 4.1 option with Langmuir kinetics and constant equilibrium
        !! @to be finished...
            hlp5 = avail_sorption_sites
            adsorpt_f_po4 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            desorp_f_po4 = assoc_fast * k_desorp_f_po4


        ELSEIF (flag_sb_CLE_po4 .AND. (.NOT. flag_sb_LK_po4)) THEN
        !> 4.2 option with constant equilibrium, i.e. the same as Wang et al. 2007 and so forth
        !! 
            hlp4 = km_adsorpt_po4 * bulk_dens_corr_sl
            hlp2 = qmax_po4 * hlp4 / ((solute + hlp4)**2.0_wp + qmax_po4 * hlp4)
            hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4 - &
                    (mic_uptake + plant_uptake + slow_sorption_po4))/1.e6_wp + perco_po4(:,:)
            WHERE (hlp6 * dtime + solute + assoc_fast < eps8 .AND. solute > eps8)
                hlp1 = mic_uptake + plant_uptake
                hlp3 = (hlp6 + hlp1) * dtime + solute + assoc_fast
                mic_uptake = mic_uptake / hlp1 * hlp3
                plant_uptake = plant_uptake / hlp1 * hlp3
                hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4 - &
                        (mic_uptake + plant_uptake + slow_sorption_po4))/1.e6_wp + perco_po4(:,:)
            ENDWHERE
            fast_adsorption_po4 = hlp2 * hlp6 + slow_sorption_po4

        ELSEIF (flag_sb_LK_po4 .AND. (.NOT. flag_sb_CLE_po4)) THEN
        !> 4.3 option with Langmuir kinetics but not constant equilibrium
        !!
            hlp5 = avail_sorption_sites
            adsorpt_f_po4 = (k_adsorpt_f_po4/bulk_dens_corr_sl) * solute * hlp5
            desorp_f_po4 = assoc_fast * k_desorp_f_po4

        !> 4.3.1 limit to possible P sources
        !!
            hlp1 = (mic_uptake + plant_uptake)/1.e6_wp + adsorpt_f_po4
            hlp6 = (weathering_po4 + biochem_min_po4 + gross_mineralisation_po4)/1.e6_wp + desorp_f_po4 + perco_po4(:,:)

            WHERE((hlp1 + solute / dtime) > hlp6 .AND. solute > eps8)
               mic_uptake    = mic_uptake    / hlp1 * hlp6
               plant_uptake  = plant_uptake  / hlp1 * hlp6
               adsorpt_f_po4 = adsorpt_f_po4 / hlp1 * hlp6
            ENDWHERE
        ENDIF


    END SELECT

    !! convert the units
    fast_adsorption_po4   = fast_adsorption_po4 * 1.e6_wp

  END SUBROUTINE calc_adsorption_po4
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of the microbial pool
  !! corresponds to the term pi * Cb in Ahrens et al. 2015, eq. 7
  !!
  !! Input: microbial biomass
  !!
  !! Output: transfer of microbial biomass to its sink pools (soluable -> DOM; polymeric -> residue)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_turnover(options,pool_microbial,loss_microbial,formation_dom,formation_residue)

    USE mo_sb_constants,          ONLY: k_mort_microbial, min_mic_biomass, f_soluable_residue
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),       INTENT(in)    :: options            !< model options
    TYPE(t_sb_element_a),           INTENT(in)    :: pool_microbial     !< microbial biomass [mol/m3]
    TYPE(t_sb_element_a),           INTENT(inout) :: loss_microbial     !< microbial death [mol/m3/timestep]
    TYPE(t_sb_element_a),           INTENT(inout) :: formation_dom      !< DOM formation from microbial death [mol/m3/timestep]
    TYPE(t_sb_element_a),           INTENT(inout) :: formation_residue  !< residue formation from microbial death [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1 , &
                                                        mortality_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate microbial mortality rate, limited to minimum microbial biomass density (-> dormant microbes)
    !!
    mortality_rate(:,:) = k_mort_microbial * dtime
    hlp1(:,:)           = mortality_rate(:,:) * pool_microbial%carbon(:,:) 
    WHERE(pool_microbial%carbon(:,:) - hlp1(:,:) < min_mic_biomass)
      mortality_rate(:,:) = 1._wp - (min_mic_biomass / pool_microbial%carbon(:,:))
    ENDWHERE

    !> 2.0 calculate matter fluxes associated with microbial death
    !! a fraction turns to DOM, the remainder becomes residue
    !! the assumption is made that the composition of this material is not changed
    !!@todo: define C:N:P stoichiometries of DOM and residues from microbial death 
    !!
    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon(:,:)
    loss_microbial%carbon(:,:)        = loss_microbial%carbon(:,:) + hlp1(:,:)
    formation_dom%carbon(:,:)         = formation_dom%carbon(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon(:,:)     = formation_residue%carbon(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen(:,:)
    loss_microbial%nitrogen(:,:)      = loss_microbial%nitrogen(:,:) + hlp1(:,:)
    formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%nitrogen(:,:)   = formation_residue%nitrogen(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%phosphorus(:,:)
    loss_microbial%phosphorus(:,:)    = loss_microbial%phosphorus(:,:) + hlp1(:,:)
    formation_dom%phosphorus(:,:)     = formation_dom%phosphorus(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%phosphorus(:,:) = formation_residue%phosphorus(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon13(:,:)
    loss_microbial%carbon13(:,:)      = loss_microbial%carbon13(:,:) + hlp1(:,:)
    formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon13(:,:)   = formation_residue%carbon13(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%carbon14(:,:)
    loss_microbial%carbon14(:,:)      = loss_microbial%carbon14(:,:) + hlp1(:,:)
    formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%carbon14(:,:)   = formation_residue%carbon14(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

    hlp1(:,:)                         = mortality_rate(:,:) * pool_microbial%nitrogen15(:,:)
    loss_microbial%nitrogen15(:,:)    = loss_microbial%nitrogen15(:,:) + hlp1(:,:)
    formation_dom%nitrogen15(:,:)     = formation_dom%nitrogen15(:,:) + f_soluable_residue * hlp1(:,:) 
    formation_residue%nitrogen15(:,:) = formation_residue%nitrogen15(:,:) + (1._wp - f_soluable_residue) * hlp1(:,:) 

  END SUBROUTINE calc_microbial_turnover


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate growth of the microbial pool 
  !! corresponds to Ahrens et al. 2015, eq. 2 
  !! extended to account for N and P limitation of microbial growth given a prescribed 
  !! stoichometry
  !!
  !! Input: DOM and microbial pools, microbial nutrient uptake rate, temperature and moisture rate modifiers 
  !!
  !! Output: DOM consumption, microbial growth, hetertrophic respiration and net mineralisation
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_growth(options,pool_dom,pool_microbial, & 
                                   rtm_mic_uptake,rmm_mic_uptake,rtm_hsc,rmm_hsc, &
                                   mic_uptake_nh4,mic_uptake_nh4_n15, &
                                   mic_uptake_no3,mic_uptake_no3_n15,mic_uptake_po4, &
                                   loss_dom,formation_microbial, &
                                   het_respiration,het_respiration_c13,het_respiration_c14, &
                                   net_mineralisation_nh4,net_mineralisation_nh4_n15, &
                                   net_mineralisation_no3,net_mineralisation_no3_n15,net_mineralisation_po4)

    USE mo_shared_isotope, ONLY: calc_fractionation
    USE mo_sb_constants,   ONLY: vmax_mic_uptake,km_mic_uptake,microbial_cue,microbial_nue,microbial_pue, & 
                                 microbial_cn,microbial_np,eta_ammonification
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options             !< model options
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom            !< DOM [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_microbial      !< microbial biomass [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_mic_uptake      !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_mic_uptake      !< moisture rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_hsc             !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_hsc             !< moisture rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mic_uptake_nh4      !< microbial uptake rate [µmol N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mic_uptake_nh4_n15  !< microbial uptake rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mic_uptake_no3      !< microbial uptake rate [µmol N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mic_uptake_no3_n15  !< microbial uptake rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mic_uptake_po4      !< microbial uptake rate [µmol P/m3/s]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom            !< microbial death [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_microbial !< DOM formation from microbial death [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: het_respiration     !< respiration rate [µmol C/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: het_respiration_c13 !< respiration rate [µmol 13C/m3/s] 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: het_respiration_c14 !< respiration rate [µmol 14C/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: net_mineralisation_nh4     !< net mineralisation rate [µmol N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: net_mineralisation_nh4_n15 !< net mineralisation rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: net_mineralisation_no3     !< net mineralisation rate [µmol N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: net_mineralisation_no3_n15 !< net mineralisation rate [µmol 15N/m3/s]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: net_mineralisation_po4     !< net mineralisation rate [µmol P/m3/s]
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1             , &
                                                        f_uptake         , &
                                                        uptake_dom  
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: growth_c,&
                                                        growth_n,&
                                                        growth_p,&
                                                        growth,&
                                                        uptake_c,&
                                                        uptake_n,&
                                                        uptake_p 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: scal
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: microbial_cue_eff,&
                                                        microbial_nue_eff,&
                                                        microbial_pue_eff 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: gross_mineralisation_nh4,&
                                                        gross_mineralisation_nh4_n15,&
                                                        gross_mineralisation_po4 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_growth'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...
    

    !> 1.0 potential microbial DOM consumption (µmol / m3 / s) 
    !! 
    hlp1          = pool_dom%carbon(:,:) / (km_mic_uptake * rtm_hsc * rmm_hsc + pool_dom%carbon(:,:) )
    uptake_c(:,:) = vmax_mic_uptake * rtm_mic_uptake(:,:) * rmm_mic_uptake(:,:) * hlp1(:,:) * pool_microbial%carbon(:,:) 
    WHERE(pool_dom%carbon(:,:) > eps8)
       f_uptake(:,:) = uptake_c(:,:) / pool_dom%carbon(:,:)
    ELSEWHERE
       f_uptake(:,:) = 0.0_wp
    ENDWHERE
    uptake_n(:,:) = f_uptake(:,:) * pool_dom%nitrogen(:,:) 
    uptake_p(:,:) = f_uptake(:,:) * pool_dom%phosphorus(:,:)

    !> 1.1 microbial growth rates in µmol C / m3 / s, given by C, N and P availability, and the rate-limiting growth rate
    !! 
    growth_c(:,:) = microbial_cue * uptake_c(:,:) 
    growth_n(:,:) = (microbial_nue * uptake_n(:,:) + mic_uptake_nh4(:,:) + mic_uptake_no3(:,:)) * microbial_cn
    growth_p(:,:) = (microbial_pue * uptake_p(:,:) + mic_uptake_po4(:,:)) * microbial_np * microbial_cn
    WHERE(growth_c(:,:) <= growth_n(:,:) .AND. growth_c(:,:) <= growth_p(:,:))
       ! carbon uptake is limiting
       growth(:,:) = growth_c(:,:)
    ELSEWHERE(growth_n(:,:) < growth_c(:,:) .AND. growth_n(:,:) <= growth_p(:,:))
       ! nitrogen uptake is limiting
       growth(:,:) = growth_n(:,:)
    ELSEWHERE(growth_p(:,:) < growth_c(:,:) .AND. growth_p(:,:) < growth_n(:,:))
       ! phosphorus uptake is limiting
       growth(:,:) = growth_p(:,:)
    ENDWHERE

    !> 2.0 Adjust carbon-use efficiency such that the growth is stoichiometrically balanced
    !!
    microbial_cue_eff(:,:) = microbial_cue
    WHERE(growth(:,:) < growth_c(:,:) .AND. uptake_c(:,:) > eps8)
       ! excess carbon is going to be respired
       microbial_cue_eff(:,:) = growth / uptake_c(:,:)
    ELSEWHERE(growth(:,:) < growth_c(:,:))
       microbial_cue_eff(:,:) = 0.0_wp
    ENDWHERE
 
    !> 2.1 Adjust nitrogen-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
    !!     and calculate gross mineralisation rate (and its isotopic signature)
    microbial_nue_eff(:,:) = microbial_nue
    WHERE(growth(:,:) < growth_n(:,:)) 
       ! inorganic N immobilisation flux in agreement with C growth
       hlp1(:,:) = MAX(growth / microbial_cn - microbial_nue * uptake_n(:,:),0.0_wp)
       WHERE(mic_uptake_nh4(:,:) + mic_uptake_no3(:,:) > eps8)
          scal(:,:) = hlp1(:,:) / (mic_uptake_nh4(:,:) + mic_uptake_no3(:,:)) 
       ELSEWHERE
          scal(:,:) = 0.0_wp
       ENDWHERE
       mic_uptake_nh4(:,:)     = mic_uptake_nh4(:,:)     * scal(:,:) 
       mic_uptake_nh4_n15(:,:) = mic_uptake_nh4_n15(:,:) * scal(:,:)
       mic_uptake_no3(:,:)     = mic_uptake_no3(:,:)     * scal(:,:) 
       mic_uptake_no3_n15(:,:) = mic_uptake_no3_n15(:,:) * scal(:,:)
       WHERE(uptake_n(:,:) > eps8)
          microbial_nue_eff(:,:) = ( growth / microbial_cn - mic_uptake_nh4(:,:) - mic_uptake_no3(:,:) ) / uptake_n(:,:)
       ENDWHERE
    ENDWHERE
    gross_mineralisation_nh4(:,:)     = (1._wp - microbial_nue_eff(:,:)) * uptake_n(:,:) 
    gross_mineralisation_nh4_n15(:,:) = gross_mineralisation_nh4(:,:) * calc_fractionation(pool_dom%nitrogen(:,:), & 
                                                                                           pool_dom%nitrogen15(:,:),& 
                                                                                           eta_ammonification)

    !> 2.2 Adjust phosphorus-use efficiency and immobilisation flux such that the growth is stoichiometrically balanced
    !!     and calculate gross mineralisation rate 
    microbial_pue_eff(:,:) = microbial_pue
    WHERE(growth(:,:) < growth_p(:,:))  
       ! inorganic P immobilisation flux, in agreement with C and N uptake
       mic_uptake_po4(:,:) = MAX(growth / microbial_cn / microbial_np - microbial_pue * uptake_p(:,:),0.0)
       ! In the case that the source C:N:P is too tight for microbes, cause gross mineralisation
       WHERE(uptake_p(:,:) > eps8)
          microbial_pue_eff(:,:) = ( growth / microbial_cn / microbial_np - mic_uptake_po4(:,:) ) / uptake_p(:,:)
       ENDWHERE
    ENDWHERE
    gross_mineralisation_po4(:,:) = (1._wp - microbial_pue_eff(:,:)) * uptake_p(:,:)

    !> 3.0 update fluxes (note that loss and formation of organic matter is in mol/m3/timestep, whereas fluxes are in µmol/m3/s)
    !! carbon
    uptake_dom(:,:)                 = f_uptake(:,:) * pool_dom%carbon(:,:)
    loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + uptake_dom(:,:) * dtime / 1000000._wp
    formation_microbial%carbon(:,:) = formation_microbial%carbon(:,:) + & 
                                        microbial_cue_eff(:,:) * uptake_dom(:,:) * dtime / 1000000._wp 
    het_respiration(:,:)            = het_respiration(:,:) + & 
                                        (1._wp - microbial_cue_eff(:,:)) * uptake_dom(:,:) 

    !! nitrogen
    net_mineralisation_nh4(:,:)       = gross_mineralisation_nh4(:,:) - mic_uptake_nh4(:,:)
    net_mineralisation_no3(:,:)       = - mic_uptake_no3(:,:)
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%nitrogen(:,:)
    loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%nitrogen(:,:) = formation_microbial%nitrogen(:,:) + & 
                                          (microbial_nue_eff(:,:) * uptake_dom(:,:) + & 
                                          mic_uptake_nh4(:,:) + mic_uptake_no3(:,:)) * dtime / 1000000.0_wp 
   
    !! phosphorus
    net_mineralisation_po4(:,:)         = gross_mineralisation_po4(:,:) - mic_uptake_po4(:,:)
    uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%phosphorus(:,:)
    loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%phosphorus(:,:) = formation_microbial%phosphorus(:,:) + & 
                                            (microbial_pue_eff(:,:) * uptake_dom(:,:) + & 
                                            mic_uptake_po4(:,:)) * dtime / 1000000.0_wp 
   
    !! carbon13
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon13(:,:)
    loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + uptake_dom(:,:) * dtime / 1000000._wp 
    formation_microbial%carbon13(:,:) = formation_microbial%carbon13(:,:) + & 
                                          microbial_cue_eff(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp 
    het_respiration_c13(:,:)          = het_respiration_c13(:,:) + (1._wp - microbial_cue_eff(:,:) ) * uptake_dom(:,:)

    !! carbon14
    uptake_dom(:,:)                   = f_uptake(:,:) * pool_dom%carbon14(:,:)
    loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%carbon14(:,:) = formation_microbial%carbon14(:,:) + & 
                                          microbial_cue_eff(:,:) * uptake_dom(:,:) * dtime / 1000000.0_wp
    het_respiration_c14(:,:)          = het_respiration_c14(:,:) + (1._wp - microbial_cue_eff(:,:) ) * uptake_dom(:,:)

    !! nitrogen15
    net_mineralisation_nh4_n15(:,:)     = gross_mineralisation_nh4_n15(:,:) - mic_uptake_nh4_n15(:,:)
    net_mineralisation_no3_n15(:,:)     = - mic_uptake_no3_n15(:,:)
    uptake_dom(:,:)                     = f_uptake(:,:) * pool_dom%nitrogen15(:,:)
    loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + uptake_dom(:,:) * dtime / 1000000.0_wp 
    formation_microbial%nitrogen15(:,:) = formation_microbial%nitrogen15(:,:) + & 
                                            (microbial_nue * uptake_dom(:,:) + & 
                                            mic_uptake_nh4_n15(:,:) + mic_uptake_no3_n15(:,:)) * dtime / 1000000.0_wp

  END SUBROUTINE calc_microbial_growth


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate depolymerisation rate of polymeric litter and microbial residue 
  !! corresponds Ahrens et al. 2015, eq. 1 
  !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  !!
  !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for
  !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass
  !!
  !! Output: transfer of residue/polymeric pool to DOM
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_depolymerisation(options,vmax_depolymerisation,rtm_depolymerisation,rmm_depolymerisation, &
                                   rtm_hsc,rmm_hsc,microbial_carbon,pool_polymeric,loss_polymeric,formation_dom)

    dsl4slm_Use_memory_a_sb

    USE mo_sb_constants,          ONLY: km_depolymerisation

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options               !< model options
    REAL(wp),                                         INTENT(in)    :: vmax_depolymerisation !< maximum depolymerisation rate
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier for depolymerisation
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier for depolymerisation
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_hsc               !< temperature rate modifier for half-saturation constants
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_hsc               !< moisture rate modifier for half-saturation constants
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: microbial_carbon      !< microbial biomass carbon [mol C/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_polymeric        !< polymeric/residue pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_polymeric        !< loss of polymeric/residue pool [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom         !< formation of DOM [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1                   , &
                                                        hlp2                   , &
                                                        depolymerisation_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...
    
    

    !> 1.0 determine potential depolymerisation rate [mol / m3 / timestep]
    !!
    hlp1(:,:) = 1._wp   
    hlp1(:,:) = microbial_carbon(:,:) / (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:)) 
    hlp2(:,:) = vmax_depolymerisation * rtm_depolymerisation * rmm_depolymerisation * hlp1(:,:) * dtime * & 
                 pool_polymeric%carbon(:,:) 
    
    !> 2.0 convert to turnover rate [1 / timestep] 
    !!
    WHERE(pool_polymeric%carbon(:,:) > eps8)
        depolymerisation_rate(:,:) = hlp2(:,:) / pool_polymeric%carbon(:,:)
    ELSEWHERE
        depolymerisation_rate(:,:) = 0.0_wp
    ENDWHERE

    !> 3.0 apply to all elements in this pool (should be a simple operator overload, but does currently not work safely)
    !!
    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon(:,:)
    loss_polymeric%carbon(:,:)      = loss_polymeric%carbon(:,:) + hlp2(:,:) 
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen(:,:)
    loss_polymeric%nitrogen(:,:)    = loss_polymeric%nitrogen(:,:) + hlp2(:,:) 
    formation_dom%nitrogen(:,:)     = formation_dom%nitrogen(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%phosphorus(:,:)
    loss_polymeric%phosphorus(:,:)  = loss_polymeric%phosphorus(:,:) + hlp2(:,:) 
    formation_dom%phosphorus(:,:)   = formation_dom%phosphorus(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon13(:,:)
    loss_polymeric%carbon13(:,:)    = loss_polymeric%carbon13(:,:) + hlp2(:,:) 
    formation_dom%carbon13(:,:)     = formation_dom%carbon13(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%carbon14(:,:)
    loss_polymeric%carbon14(:,:)    = loss_polymeric%carbon14(:,:) + hlp2(:,:) 
    formation_dom%carbon14(:,:)     = formation_dom%carbon14(:,:) + hlp2(:,:) 

    hlp2(:,:)                       = depolymerisation_rate(:,:) * pool_polymeric%nitrogen15(:,:)
    loss_polymeric%nitrogen15(:,:)  = loss_polymeric%nitrogen15(:,:) + hlp2(:,:) 
    formation_dom%nitrogen15(:,:)   = formation_dom%nitrogen15(:,:) + hlp2(:,:) 

  END SUBROUTINE calc_depolymerisation


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the bulk density of the soil corrected by its organic component 
  !! derived from Federer et al. 1993, Can J. For. Res. 23 1026-1032 
  !!
  !! Input: carbon content by layer 
  !!
  !! Output: corrected bulk density of the soil 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl) RESULT(bulk_dens_corr_sl)

    USE mo_sb_constants,          ONLY: rho_bulk_org 
    USE mo_shared_constants,      ONLY: molar_mass_C 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_soil_carbon_sl     !< soil carbon contributing to density [mol C / m3] 
    REAL(wp), INTENT(in) :: bulk_dens_sl             !< bulk density of the mineral soil [kg / m3]
    REAL(wp)             :: bulk_dens_corr_sl       !< combined bulk density of the mineral and organic soil [kg / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: fract_org
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_density_correction'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

 
    !> 1.0 calculate organic mass fraction per soil layer 
    !!
    IF(bulk_soil_carbon_sl * molar_mass_C + bulk_dens_sl > eps8)THEN
       fract_org = bulk_soil_carbon_sl * molar_mass_C/1000._wp / & 
                   (bulk_soil_carbon_sl * molar_mass_C/1000._wp + bulk_dens_sl)
    ELSE 
       fract_org = 0.0_wp
    ENDIF

    !> 2.0 actual bulk density of the soil 
    !!
    bulk_dens_corr_sl = bulk_dens_sl * rho_bulk_org / & 
                          (fract_org * bulk_dens_sl + (1._wp - fract_org ) * rho_bulk_org)

  END FUNCTION calc_bulk_density_correction 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the maximum sorption capacity for organic material of a soil layer
  !! derived from COMISSION, the implict assumption here is that no association happens on organic material
  !!
  !! Input: bulk density of the soil (mineral and organic combined) 
  !!
  !! Output: sorption capacity for organic material 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax(bulk_dens_corr_sl, maxcapacity) RESULT(qmax)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_dens_corr_sl     !< combined bulk density of the mineral and organic soil [kg / m3]
    REAL(wp), INTENT(in) :: maxcapacity           !< maximium sorption capacity per kg soil
    REAL(wp)             :: qmax                  !< maximum sorption capacity of the soil [kg C / m3] or [mol P / m3]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax'

    qmax = bulk_dens_corr_sl * maxcapacity

  END FUNCTION calc_qmax


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates bulk soil carbon used for the calculation of soil bulk density and sorption capacity 
  !!
  !! Input: soil pools 
  !!
  !! Output: total_bulk_soil_carbon 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_bulk_soil_carbon(options, sb_pool) RESULT (bulk_soil_carbon_sl)

    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    TYPE(t_sb_pool_a),                     INTENT(in)    :: sb_pool                 !< pools [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb)     :: bulk_soil_carbon_sl
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: nsoil_sb, isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_soil_carbon'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Get local variables from options argument
    nsoil_sb  = options%nsoil_sb

    DO isoil = 1,nsoil_sb
       bulk_soil_carbon_sl(:,isoil) = &
            sb_pool%microbial%carbon(:,isoil) + &
            sb_pool%residue%carbon(:,isoil) + &
            sb_pool%residue_assoc%carbon(:,isoil) + &
            sb_pool%dom%carbon(:,isoil) + &
            sb_pool%dom_assoc%carbon(:,isoil) + &
            sb_pool%soluable_litter%carbon(:,isoil) + &
            sb_pool%polymeric_litter%carbon(:,isoil)
    !   IF(isoil > 1 ) THEN       
         bulk_soil_carbon_sl(:,isoil) = &
            bulk_soil_carbon_sl(:,isoil) + &
            sb_pool%woody_litter%carbon(:,isoil) 
    !   ENDIF
    ENDDO

  END FUNCTION calc_bulk_soil_carbon


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates po4 adsorption rate based on the Langmuir equilibrium condition
  !!
  !! Input: potential po4 fluxes
  !!
  !! Output: actual fast_adsorption_po4
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_CLE_po4(solution_selection, solute, Km, Smax, opt1, Pin, Pup, opt3) RESULT (opt2)
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(*),INTENT(in)        :: solution_selection  !< for the selection of solution
    REAL(wp), INTENT(in)           :: solute        !< solute po4 concentration (mol / m3)
    REAL(wp), INTENT(in)           :: Km            !< half-saturation constant in Langmuir (mol / m3)
    REAL(wp), INTENT(in)           :: Smax          !< max po4 sorption capacity in Langmuir (mol / m3)
    REAL(wp), INTENT(in)           :: opt1          !< assoc_fast po4 (mol / m3)
    REAL(wp), INTENT(in), OPTIONAL :: opt3          !< potential po4 adsorption or desorption2 (mol / m3)
    REAL(wp), INTENT(in), OPTIONAL :: Pin           !< po4 total input of the timestep (mol / m3)
    REAL(wp), INTENT(in), OPTIONAL :: Pup           !< po4 total potential uptake of the timestep (mol / m3)
    REAL(wp)                       :: opt2          !< scaling factor for plant, microbe and adsorption, 0<opt<1
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: hlp1,a1,a2,a3,a4,a5,k1,k2   !< helper variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_CLE_po4'

    SELECT CASE (solution_selection)
    !> CASE 1 and 2 reserved for future
    !! 
      CASE ("insufficient_p_input")
    !> Quadratic solution for Langmuir Equilibrium: Smax*Solute/(Km+solute)=assoc_fast
    !! real solution between 0 and 1 is the scaling factor for P demand, otherwise no downscaling
    !! for case: solute + Pin < Pup
        a3 = Smax * (solute + Pin)
        a2 = opt1 * Pup
        a1 = opt1 * (Km + solute + Pin)
        a4 = opt3 * (Km + solute + Pin) + Smax * Pup
        a5 = Pup * opt3
        k1 = (a4-a2+sqrt((a4-a2)*(a4-a2)-4._wp*a5*(a3-a1)))/(2._wp*a5)
        k2 = (a4-a2-sqrt((a4-a2)*(a4-a2)-4._wp*a5*(a3-a1)))/(2._wp*a5)
        IF ( k1<(solute + Pin)/Pup .AND. k1>eps8) THEN
            opt2 = k1
        ELSEIF (k2<(solute + Pin)/Pup .AND. k2>eps8) THEN
            opt2 = k2
        ELSE
            opt2 = 0.0_wp
        ENDIF
      CASE ("sufficient_p_input")
    !> Quadratic solution for Langmuir Equilibrium: Smax*Solute/(Km+solute)=assoc_fast
    !! 
        a3 = Smax * (solute + Pin)
        a2 = opt1 * Pup
        a1 = opt1 * (Km + solute + Pin)
        a4 = opt3 * (Km + solute + Pin) + Smax * Pup
        a5 = Pup * opt3
        k1 = (a4-a2+sqrt((a4-a2)*(a4-a2)-4._wp*a5*(a3-a1)))/(2._wp*a5)
        k2 = (a4-a2-sqrt((a4-a2)*(a4-a2)-4._wp*a5*(a3-a1)))/(2._wp*a5)
        IF ( k1<1._wp .AND. k1>eps8) THEN
            opt2 = k1
        ELSEIF (k2<1._wp .AND. k2>eps8) THEN
            opt2 = k2
        ELSEIF (k1>1._wp .OR. k2>1._wp) THEN
            opt2 = 1._wp
        ELSE
            opt2 = 0.0_wp
        ENDIF
      END SELECT

  END FUNCTION calc_CLE_po4


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates sorption and desorption of organic material
  !! corresponds to the sorption and desorption terms in Ahrens et al. 2015, eq. 5 & 6
  !!
  !! Input: sorption temperature and moisture modifiers, sorption pools
  !!
  !! Output: sorption and desorption rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_om_sorption(options,qmax_org,rtm_sorption,rmm_sorption,rtm_desorption, &
                              rmm_desorption,pool_dom,pool_dom_assoc,pool_residue, & 
                              pool_residue_assoc,formation_dom,formation_dom_assoc,formation_residue, &
                              formation_residue_assoc,loss_dom,loss_dom_assoc,loss_residue,loss_residue_assoc)

    USE mo_sb_constants,          ONLY: k_adsorpt_dom,k_desorpt_dom,k_adsorpt_det,k_desorpt_det 
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: qmax_org 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_sorption            !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_sorption            !< moisture rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_desorption          !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rmm_desorption          !< moisture rate modifier [-]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom                !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue            !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_dom_assoc          !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue_assoc      !< DOM pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom           !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_dom_assoc     !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue       !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: formation_residue_assoc !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom                !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_dom_assoc          !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue            !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: loss_residue_assoc      !< rate of DOM transport [mol/m3/timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1                   , &
                                                        avail_sorption_sites   , &
                                                        sorption_rate          , &
                                                        desorption_rate        
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_om_sorption'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...


    !> 1.0 Determine available sorption sites [mol C m-3]
    !!
    avail_sorption_sites(:,:) = qmax_org(:,:) - (pool_dom_assoc%carbon(:,:) + pool_residue_assoc%carbon(:,:))

    !> 2.0 mineral association of DOM
    sorption_rate(:,:)   = k_adsorpt_dom * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) 
    desorption_rate(:,:) = k_desorpt_dom * rtm_desorption(:,:) * rmm_desorption(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon(:,:)
    formation_dom_assoc%carbon(:,:) = formation_dom_assoc%carbon(:,:) + hlp1(:,:) 
    loss_dom%carbon(:,:)            = loss_dom%carbon(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon(:,:)
    formation_dom%carbon(:,:)       = formation_dom%carbon(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon(:,:)      = loss_dom_assoc%carbon(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen(:,:)
    formation_dom_assoc%nitrogen(:,:) = formation_dom_assoc%nitrogen(:,:) + hlp1(:,:) 
    loss_dom%nitrogen(:,:)            = loss_dom%nitrogen(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen(:,:)
    formation_dom%nitrogen(:,:)       = formation_dom%nitrogen(:,:) + hlp1(:,:) 
    loss_dom_assoc%nitrogen(:,:)      = loss_dom_assoc%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%phosphorus(:,:)
    formation_dom_assoc%phosphorus(:,:) = formation_dom_assoc%phosphorus(:,:) + hlp1(:,:) 
    loss_dom%phosphorus(:,:)            = loss_dom%phosphorus(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%phosphorus(:,:)
    formation_dom%phosphorus(:,:)       = formation_dom%phosphorus(:,:) + hlp1(:,:) 
    loss_dom_assoc%phosphorus(:,:)      = loss_dom_assoc%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon13(:,:)
    formation_dom_assoc%carbon13(:,:) = formation_dom_assoc%carbon13(:,:) + hlp1(:,:) 
    loss_dom%carbon13(:,:)            = loss_dom%carbon13(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon13(:,:)
    formation_dom%carbon13(:,:)       = formation_dom%carbon13(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon13(:,:)      = loss_dom_assoc%carbon13(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%carbon14(:,:)
    formation_dom_assoc%carbon14(:,:) = formation_dom_assoc%carbon14(:,:) + hlp1(:,:) 
    loss_dom%carbon14(:,:)            = loss_dom%carbon14(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%carbon14(:,:)
    formation_dom%carbon14(:,:)       = formation_dom%carbon14(:,:) + hlp1(:,:) 
    loss_dom_assoc%carbon14(:,:)      = loss_dom_assoc%carbon14(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_dom%nitrogen15(:,:)
    formation_dom_assoc%nitrogen15(:,:) = formation_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 
    loss_dom%nitrogen15(:,:)            = loss_dom%nitrogen15(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_dom_assoc%nitrogen15(:,:)
    formation_dom%nitrogen15(:,:)       = formation_dom%nitrogen15(:,:) + hlp1(:,:) 
    loss_dom_assoc%nitrogen15(:,:)      = loss_dom_assoc%nitrogen15(:,:) + hlp1(:,:) 

    !> 3.0 mineral association of microbial residue 
    sorption_rate(:,:)   = k_adsorpt_det * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) 
    desorption_rate(:,:) = k_desorpt_det * rtm_desorption(:,:) * rmm_desorption(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon(:,:)
    formation_residue_assoc%carbon(:,:) = formation_residue_assoc%carbon(:,:) + hlp1(:,:) 
    loss_residue%carbon(:,:)            = loss_residue%carbon(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon(:,:)
    formation_residue%carbon(:,:)       = formation_residue%carbon(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon(:,:)      = loss_residue_assoc%carbon(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen(:,:)
    formation_residue_assoc%nitrogen(:,:) = formation_residue_assoc%nitrogen(:,:) + hlp1(:,:) 
    loss_residue%nitrogen(:,:)            = loss_residue%nitrogen(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen(:,:)
    formation_residue%nitrogen(:,:)       = formation_residue%nitrogen(:,:) + hlp1(:,:) 
    loss_residue_assoc%nitrogen(:,:)      = loss_residue_assoc%nitrogen(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%phosphorus(:,:)
    formation_residue_assoc%phosphorus(:,:) = formation_residue_assoc%phosphorus(:,:) + hlp1(:,:) 
    loss_residue%phosphorus(:,:)            = loss_residue%phosphorus(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%phosphorus(:,:)
    formation_residue%phosphorus(:,:)       = formation_residue%phosphorus(:,:) + hlp1(:,:) 
    loss_residue_assoc%phosphorus(:,:)      = loss_residue_assoc%phosphorus(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon13(:,:)
    formation_residue_assoc%carbon13(:,:) = formation_residue_assoc%carbon13(:,:) + hlp1(:,:) 
    loss_residue%carbon13(:,:)            = loss_residue%carbon13(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon13(:,:)
    formation_residue%carbon13(:,:)       = formation_residue%carbon13(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon13(:,:)      = loss_residue_assoc%carbon13(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%carbon14(:,:)
    formation_residue_assoc%carbon14(:,:) = formation_residue_assoc%carbon14(:,:) + hlp1(:,:) 
    loss_residue%carbon14(:,:)            = loss_residue%carbon14(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%carbon14(:,:)
    formation_residue%carbon14(:,:)       = formation_residue%carbon14(:,:) + hlp1(:,:) 
    loss_residue_assoc%carbon14(:,:)      = loss_residue_assoc%carbon14(:,:) + hlp1(:,:) 

    hlp1(:,:) = sorption_rate(:,:) * pool_residue%nitrogen15(:,:)
    formation_residue_assoc%nitrogen15(:,:) = formation_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 
    loss_residue%nitrogen15(:,:)            = loss_residue%nitrogen15(:,:) + hlp1(:,:) 
    hlp1(:,:) = desorption_rate(:,:) * pool_residue_assoc%nitrogen15(:,:)
    formation_residue%nitrogen15(:,:)       = formation_residue%nitrogen15(:,:) + hlp1(:,:) 
    loss_residue_assoc%nitrogen15(:,:)      = loss_residue_assoc%nitrogen15(:,:) + hlp1(:,:) 

  END SUBROUTINE calc_om_sorption 

END MODULE mo_sb_jsm_processes
