!>
!! @brief Contains the routines for the Task update_sb_pools
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_sb_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_pools


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_update_pools'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_base_config_class,     ONLY: base_config
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_sb_constants
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_soil_memory), POINTER            :: soil       ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)

    soil => soil__mem
    sb   => sb__mem

    !> 1.0 Calculate biophysical rate modifiers
    !!
    !! for the full and plant only model, update soil modifiers. In the soil model, these are read from a file 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

         sb%rtm_depolymerisation%ptr = calc_temperature_rate_modifier(soil%t_soil%ptr,ea_depolymerisation) 
         sb%rtm_mic_uptake%ptr       = calc_temperature_rate_modifier(soil%t_soil%ptr,ea_mic_uptake)
         sb%rtm_sorption%ptr         = calc_temperature_rate_modifier(soil%t_soil%ptr,ea_sorption) 
         sb%rtm_desorption%ptr       = calc_temperature_rate_modifier(soil%t_soil%ptr,ea_desorption) 
         sb%rtm_hsc%ptr              = calc_temperature_rate_modifier(soil%t_soil%ptr,ea_hsc) 
     
         sb%rmm_depolymerisation%ptr = calc_oxygen_rate_modifier(soil%w_soil_sl%ptr,soil%w_soil_fc_sl%ptr,kmr_depolymerisation) 
         sb%rmm_mic_uptake%ptr       = calc_oxygen_rate_modifier(soil%w_soil_sl%ptr,soil%w_soil_fc_sl%ptr,kmr_mic_uptake)
         sb%rmm_sorption%ptr         = calc_convert_bulk_to_liquid(soil%w_soil_sl%ptr,kmr_sorption) 
         sb%rmm_desorption%ptr       = calc_convert_bulk_to_liquid(soil%w_soil_sl%ptr,kmr_desorption) 
         sb%rmm_hsc%ptr              = calc_substrate_rate_modifier(soil%w_soil_sl%ptr,soil%w_soil_fc_sl%ptr,kmr_hsc) 
     
    ENDSELECT
 
    !> 2.0 Calculate soil pools 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

         CALL update_sb_pools_simple(tile)

    ENDSELECT

  END SUBROUTINE update_sb_pools

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task
  !
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools_simple(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime,one_day,one_year
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_sb_constants

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET    :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools_simple'
    INTEGER                               :: isoil
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    TYPE(t_soil_memory), POINTER            :: soil       ! local pointer to soil__mem for SLM usability
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)

    veg  => veg__mem
    soil => soil__mem
    sb   => sb__mem

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !> 1 Calculate soil litter pool decomposition.
    !!
    sb%het_respiration%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon%ptr 
    sb%net_mineralisation_n%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%nitrogen%ptr 
    sb%net_mineralisation_p%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%phosphorus%ptr 
    sb%het_respiration_c13%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon13%ptr 
    sb%het_respiration_c14%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon14%ptr 
    sb%net_mineralisation_n15%ptr = sb%rtm_depolymerisation%ptr/tau_litter/one_day/one_year * sb%pool%structural_litter%nitrogen15%ptr 

    ! > 1.1 subtract losses
    !!
    sb%pool%structural_litter%carbon%ptr = sb%pool%structural_litter%carbon%ptr - sb%het_respiration%ptr * dtime
    sb%pool%structural_litter%nitrogen%ptr = sb%pool%structural_litter%nitrogen%ptr - sb%net_mineralisation_n%ptr * dtime
    sb%pool%structural_litter%phosphorus%ptr = sb%pool%structural_litter%phosphorus%ptr - sb%net_mineralisation_p%ptr * dtime
    sb%pool%structural_litter%carbon13%ptr = sb%pool%structural_litter%carbon13%ptr - sb%het_respiration_c13%ptr * dtime
    sb%pool%structural_litter%carbon14%ptr = sb%pool%structural_litter%carbon14%ptr - sb%het_respiration_c14%ptr * dtime
    sb%pool%structural_litter%nitrogen15%ptr = sb%pool%structural_litter%nitrogen15%ptr - sb%net_mineralisation_n15%ptr * dtime

    !> 2 Update soil litter pool with litter fall
    !!
    sb%pool%structural_litter%carbon%ptr(:,1,:) = sb%pool%structural_litter%carbon%ptr(:,1,:) + & 
       (veg%litterfall%leaf%carbon%ptr + veg%litterfall%sap_wood%carbon%ptr + & 
        veg%litterfall%heart_wood%carbon%ptr + veg%litterfall%fruit%carbon%ptr + & 
        veg%litterfall%labile%carbon%ptr + veg%litterfall%reserve%carbon%ptr + & 
        veg%litterfall%seed_bed%carbon%ptr ) / soil%soil_depth_l%ptr(:,1,:)
    sb%pool%structural_litter%nitrogen%ptr(:,1,:) = sb%pool%structural_litter%nitrogen%ptr(:,1,:) +  & 
       (veg%litterfall%leaf%nitrogen%ptr + veg%litterfall%sap_wood%nitrogen%ptr + & 
        veg%litterfall%heart_wood%nitrogen%ptr + veg%litterfall%fruit%nitrogen%ptr + & 
        veg%litterfall%labile%nitrogen%ptr + veg%litterfall%reserve%nitrogen%ptr + & 
        veg%litterfall%seed_bed%nitrogen%ptr ) / soil%soil_depth_l%ptr(:,1,:)
    sb%pool%structural_litter%phosphorus%ptr(:,1,:) = sb%pool%structural_litter%phosphorus%ptr(:,1,:) + & 
       (veg%litterfall%leaf%phosphorus%ptr + veg%litterfall%sap_wood%phosphorus%ptr + & 
        veg%litterfall%heart_wood%phosphorus%ptr + veg%litterfall%fruit%phosphorus%ptr + & 
        veg%litterfall%labile%phosphorus%ptr + veg%litterfall%reserve%phosphorus%ptr + & 
        veg%litterfall%seed_bed%phosphorus%ptr ) / soil%soil_depth_l%ptr(:,1,:)
    sb%pool%structural_litter%carbon13%ptr(:,1,:) = sb%pool%structural_litter%carbon13%ptr(:,1,:) + & 
       (veg%litterfall%leaf%carbon13%ptr + veg%litterfall%sap_wood%carbon13%ptr + & 
        veg%litterfall%heart_wood%carbon13%ptr + veg%litterfall%fruit%carbon13%ptr + & 
        veg%litterfall%labile%carbon13%ptr + veg%litterfall%reserve%carbon13%ptr + & 
        veg%litterfall%seed_bed%carbon13%ptr ) / soil%soil_depth_l%ptr(:,1,:)
    sb%pool%structural_litter%carbon14%ptr(:,1,:) = sb%pool%structural_litter%carbon14%ptr(:,1,:) + & 
       (veg%litterfall%leaf%carbon14%ptr + veg%litterfall%sap_wood%carbon14%ptr + & 
        veg%litterfall%heart_wood%carbon14%ptr + veg%litterfall%fruit%carbon14%ptr + & 
        veg%litterfall%labile%carbon14%ptr + veg%litterfall%reserve%carbon14%ptr + & 
        veg%litterfall%seed_bed%carbon14%ptr ) / soil%soil_depth_l%ptr(:,1,:)
    sb%pool%structural_litter%nitrogen15%ptr(:,1,:) = sb%pool%structural_litter%nitrogen15%ptr(:,1,:) + &
       (veg%litterfall%leaf%nitrogen15%ptr + veg%litterfall%sap_wood%nitrogen15%ptr + & 
        veg%litterfall%heart_wood%nitrogen15%ptr + veg%litterfall%fruit%nitrogen15%ptr + & 
        veg%litterfall%labile%nitrogen15%ptr + veg%litterfall%reserve%nitrogen15%ptr + & 
        veg%litterfall%seed_bed%nitrogen15%ptr ) / soil%soil_depth_l%ptr(:,1,:)

    DO isoil = 1,1

       sb%pool%structural_litter%carbon%ptr(:,isoil,:) = sb%pool%structural_litter%carbon%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%carbon%ptr + veg%litterfall%coarse_root%carbon%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 
       sb%pool%structural_litter%nitrogen%ptr(:,isoil,:) = sb%pool%structural_litter%nitrogen%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%nitrogen%ptr + veg%litterfall%coarse_root%nitrogen%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 
       sb%pool%structural_litter%phosphorus%ptr(:,isoil,:) = sb%pool%structural_litter%phosphorus%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%phosphorus%ptr + veg%litterfall%coarse_root%phosphorus%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 
       sb%pool%structural_litter%carbon13%ptr(:,isoil,:) = sb%pool%structural_litter%carbon13%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%carbon13%ptr + veg%litterfall%coarse_root%carbon13%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 
       sb%pool%structural_litter%carbon14%ptr(:,isoil,:) = sb%pool%structural_litter%carbon14%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%carbon14%ptr + veg%litterfall%coarse_root%carbon14%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 
       sb%pool%structural_litter%nitrogen15%ptr(:,isoil,:) = sb%pool%structural_litter%nitrogen15%ptr(:,isoil,:) + & 
          ( veg%litterfall%fine_root%nitrogen15%ptr + veg%litterfall%coarse_root%nitrogen15%ptr ) / & ! * & 
          !soil%root_fraction_l%ptr(:,isoil,:) / & 
            soil%soil_depth_l%ptr(:,isoil,:) 

    ENDDO

    !> 3 Update mineral pool
    !!
    IF(veg_config%n_uptake_scheme /= 'prescribed') THEN
       sb%volatilisation_n%ptr = sb%net_mineralisation_n%ptr * floss_nmin * 1.e6_wp 
       sb%nitrification%ptr = sb%net_mineralisation_n%ptr * (1._wp - floss_nmin) * fnit_nmin * 1.e6_wp 
       sb%nh4_solute%ptr = sb%nh4_solute%ptr + sb%net_mineralisation_n%ptr * (1._wp - floss_nmin) * (1._wp - fnit_nmin) * dtime
       sb%no3_solute%ptr = sb%no3_solute%ptr + sb%net_mineralisation_n%ptr * (1._wp - floss_nmin) * fnit_nmin * dtime
       sb%nh4_solute%ptr(:,1,:) = sb%nh4_solute%ptr(:,1,:) - veg%uptake_nh4%ptr / soil%soil_depth_l%ptr(:,1,:) * dtime/1.e6_wp
       sb%no3_solute%ptr(:,1,:) = sb%no3_solute%ptr(:,1,:) - veg%uptake_no3%ptr / soil%soil_depth_l%ptr(:,1,:) * dtime/1.e6_wp
    ENDIF

    IF(veg_config%p_uptake_scheme /= 'prescribed') THEN
       sb%po4_solute%ptr = sb%po4_solute%ptr + sb%net_mineralisation_p%ptr * dtime
       sb%po4_solute%ptr(:,1,:) = sb%po4_solute%ptr(:,1,:) - veg%uptake_p%ptr / soil%soil_depth_l%ptr(:,1,:) * dtime/1.e6_wp
    ENDIF

    END SUBROUTINE update_sb_pools_simple


!--------------------------------------------------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates temperature response factor for various soil biogeochemical reactions
  !!  according to the Arrhenius equation as defined in
  !!  Wang, G., Post, W.M., Mayes, M.A., Frerichs, J.T., Sindhu, J., 2012.
  !!  Parameter estimation for models of ligninolytic and cellulolytic enzyme kinetics.
  !!  Soil Biology and Biochemistry 48, 28-38.
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil temperature, activiation energy 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_temperature_rate_modifier(t_soil,Ea) RESULT(rate_modifier)

    USE mo_shared_constants,    ONLY: r_gas
    USE mo_sb_constants,        ONLY: t_ref_tresponse 

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil        !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    IF(Ea>eps8)THEN
       rate_modifier = exp(-Ea/r_gas*(1._wp/t_soil - 1._wp/t_ref_tresponse))
    ELSE
       rate_modifier = 1._wp
    ENDIF   

  END FUNCTION calc_temperature_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Oxygen limitation of decomposition and uptake due to soil moisture response function 
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    USE mo_sb_constants,        ONLY: k1_afps 

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl     !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl  !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param         !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_oxygen_rate_modifier'

    IF(param>eps8) THEN
       IF(w_soil_fc_sl > eps8) THEN
          air_filled_pore_space = (w_soil_fc_sl - w_soil_sl ) / w_soil_fc_sl
          rate_modifier = (air_filled_pore_space**k1_afps) / param + (air_filled_pore_space**k1_afps)
       ELSE
          rate_modifier = 1._wp
       ENDIF
    ELSE 
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_oxygen_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Converts bulk density to liquid concentration 
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: conversion modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_convert_bulk_to_liquid(w_soil_sl,param) RESULT(factor)

    USE mo_kind,  ONLY: eps4

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: param               !< 
    REAL(wp)              :: factor              !< conversion factor [unitless] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_convert_bulk_to_liquid'

    IF(param>eps8) THEN
       IF(w_soil_sl > eps4) THEN
         factor = 1._wp / w_soil_sl 
       ELSE
          factor = 1._wp / eps4 
       ENDIF
    ELSE 
       factor = 1._wp
    ENDIF

  END FUNCTION calc_convert_bulk_to_liquid

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Substrate limitation of decomposition and uptake due to soil moisture response function 
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    USE mo_kind,  ONLY: eps4

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param               !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_substrate_rate_modifier'
    
    IF(param>eps8.AND.w_soil_fc_sl > eps8) THEN
       air_filled_pore_space = (w_soil_sl - w_soil_fc_sl ) / w_soil_fc_sl
       IF(air_filled_pore_space > eps4) THEN
          rate_modifier = 1._wp / air_filled_pore_space
       ELSE 
          rate_modifier = 1._wp / eps4
       ENDIF
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_substrate_rate_modifier

END MODULE mo_sb_update_pools
