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
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's

    CALL update_sb_pools_simple(tile)

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
    !! Needs temperature and moisture correction
    sb%het_respiration%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon%ptr 
    sb%net_mineralisation_n%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%nitrogen%ptr 
    sb%net_mineralisation_p%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%phosphorus%ptr 
    sb%het_respiration_c13%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon13%ptr 
    sb%het_respiration_c14%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%carbon14%ptr 
    sb%net_mineralisation_n15%ptr = 1.0_wp/tau_litter/one_day/one_year * sb%pool%structural_litter%nitrogen15%ptr 

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
    IF(veg_config%n_uptake_scheme=='prescribed') THEN
       sb%nh4_solute%ptr(:,1,:) = 0.2_wp
    ELSE
       sb%nh4_solute%ptr = sb%nh4_solute%ptr + sb%net_mineralisation_n%ptr
       sb%nh4_solute%ptr(:,1,:) = sb%nh4_solute%ptr(:,1,:) - veg%uptake_n%ptr / soil%soil_depth_l%ptr(:,1,:) * dtime/1.e6_wp
    ENDIF

    !sb%po4_solute%ptr = sb%po4_solute%ptr + sb%net_mineralisation_p%ptr
    !sb%po4_solute%ptr = sb%po4_solute%ptr - veg%uptake_p%ptr * dtime/1.e6_wp

    END SUBROUTINE update_sb_pools_simple

END MODULE mo_sb_update_pools
