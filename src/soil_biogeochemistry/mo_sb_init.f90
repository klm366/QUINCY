!>
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sb_init

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish,message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_config_class,     ONLY: base_config  
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_shared_isotope,        ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, calc_mixing_ratio_N15N14
    USE mo_shared_constants,      ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_sb_jsm_processes
    USE mo_sb_jsm_transport
    USE mo_veg_constants,         ONLY: itree
    USE mo_sb_constants

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(SOIL_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)
    dsl4slm_Use_memory_a_sb
    dsl4jsb_Use_memory(SOIL_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout)      :: tile       !< data structure for dynamic vegetation and soil variables
    TYPE(t_jsb_task_options), INTENT(in)         :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: isoil
    INTEGER                                      :: iblk, ics, ice, nc
    INTEGER                                      :: nsoil_sb
    REAL(wp)                                     :: hlp1,hlp2,hlp3,hlp4,hlp5,hlp6
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SOIL_)
    ! Declare pool structures in memory
#define Def_sb_pool
#include "slm_pools_a_incl.f90"
#undef Def_sb_pool
!! *slmtemp    
    TYPE(t_soil_memory),    POINTER              :: soil       ! local pointer to soil__mem for SLM usability
    TYPE(t_sb_memory),    POINTER              :: sb 
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(SOIL_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    ! ...
    ! Get process memory pool structures
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool

!! *slmtemp
    soil    => soil__mem 
    sb    => sb__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))

       ! Initial values for Soenke's four pool soil biogeochemical model
       CASE("simple_0d", "simple_1d")

          hlp1 = 100._wp ! g C m-2 metabolic litter
          hlp2 = 100._wp ! g C m-2 structural litter
          hlp3 = 1000._wp ! g C m-2 woody litter
          hlp4 = 4500._wp ! g C m-2 fast SOM
          hlp5 = 4500._wp ! g C m-2 slow SOM
          hlp6 = 0.6_wp   ! fraction of production not roots 

          IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
             sb_pool%soluable_litter%carbon(:,1) = hlp1 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             IF(lctlib%growthform==itree) &
                sb_pool%woody_litter%carbon(:,1) = hlp3 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb_pool%microbial%carbon(:,1) = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb_pool%residue%carbon(:,1) = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)

             sb%nh4_solute%ptr(ics:ice,1,iblk) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%no3_solute%ptr(ics:ice,1,iblk) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%po4_solute%ptr(ics:ice,1,iblk) = 0.01_wp / dsl4jsb_Config(SOIL_)%dz_water(1) 
          ELSE 
             sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             IF(lctlib%growthform==itree) &
                sb_pool%woody_litter%carbon(:,1)  = hlp3 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb_pool%microbial%carbon(:,1)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             sb_pool%residue%carbon(:,1)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * & 
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             DO isoil=2,options%nsoil_sb
                sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb_pool%microbial%carbon(:,isoil)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb_pool%residue%carbon(:,isoil)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * & 
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk)

                sb%nh4_solute%ptr(ics:ice,isoil,iblk) = 0.1_wp
                sb%no3_solute%ptr(ics:ice,isoil,iblk) = 0.1_wp 
                sb%po4_solute%ptr(ics:ice,isoil,iblk) = 0.01_wp 
             ENDDO
          ENDIF
 
          ! metabolic litter all other elements
          sb_pool%soluable_litter%nitrogen(:,:)   = sb_pool%soluable_litter%carbon(:,:) / & 
                                                                    (25._wp * molar_mass_N / molar_mass_C)
          sb_pool%soluable_litter%phosphorus(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%soluable_litter%carbon13(:,:)   = sb_pool%soluable_litter%carbon(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb_pool%soluable_litter%carbon14(:,:)   = sb_pool%soluable_litter%carbon(:,:) * & 
                                                                    calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb_pool%soluable_litter%nitrogen15(:,:) = sb_pool%soluable_litter%nitrogen(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool%polymeric_litter%nitrogen(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / & 
                                                                    (100._wp * molar_mass_N / molar_mass_C)
          sb_pool%polymeric_litter%phosphorus(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%polymeric_litter%carbon13(:,:)   = sb_pool%polymeric_litter%carbon(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb_pool%polymeric_litter%carbon14(:,:)   = sb_pool%polymeric_litter%carbon(:,:) * & 
                                                                    calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb_pool%polymeric_litter%nitrogen15(:,:) = sb_pool%polymeric_litter%nitrogen(:,:) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool%woody_litter%nitrogen(:,:)   = sb_pool%woody_litter%carbon(:,:) / & 
                                                                (260._wp * molar_mass_N / molar_mass_C)
          sb_pool%woody_litter%phosphorus(:,:) = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                                (14._wp * molar_mass_P / molar_mass_N) 
          sb_pool%woody_litter%carbon13(:,:)   = sb_pool%woody_litter%carbon(:,:) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb_pool%woody_litter%carbon14(:,:)   = sb_pool%woody_litter%carbon(:,:) * & 
                                                                calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb_pool%woody_litter%nitrogen15(:,:) = sb_pool%woody_litter%nitrogen(:,:) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (=microbial) all other elements
          sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / microbial_cn
          sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / microbial_np 
          sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb_pool%residue%nitrogen(:,:)   = sb_pool%residue%carbon(:,:) / microbial_cn
          sb_pool%residue%phosphorus(:,:) = sb_pool%residue%nitrogen(:,:) / microbial_np 
          sb_pool%residue%carbon13(:,:)   = sb_pool%residue%carbon(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb_pool%residue%carbon14(:,:)   = sb_pool%residue%carbon(:,:) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb_pool%residue%nitrogen15(:,:) = sb_pool%residue%nitrogen(:,:) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          sb%no3_n15_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

       CASE('jsm')
 
         ! set microbial biomass to minimum possible density
         sb_pool%microbial%carbon(:,:)     = min_mic_biomass
         sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / microbial_cn 
         sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / microbial_np
         sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
         sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
         sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:)  / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
     
         ! some initial values for testing
         sb%nh4_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%no3_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%no3_n15_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%po4_solute%ptr(ics:ice,:,iblk) = 0.01_wp 

         ! inital Qmax based on current soil carbon content
         sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk) = calc_bulk_soil_carbon(options, sb%pool)
         sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk) = calc_bulk_density_correction(sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk), & 
                                                                         soil%bulk_dens_sl%ptr(ics:ice,:,iblk))
         sb%qmax_org%ptr(ics:ice,:,iblk) = calc_qmax(sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk),qmax_org_mineral_default)
         sb%qmax_po4%ptr(ics:ice,:,iblk) = calc_qmax(sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk),qmax_po4_mineral)
 
         ! initialise diffusion by bioturbation
         sb%k_bioturb%ptr(ics:ice,:,iblk) = calc_bioturbation_rate(options,sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk)) 

		 ! initialise mineral soil P pools
         sb%po4_assoc_fast%ptr(ics:ice,:,iblk) = sb%po4_solute%ptr(ics:ice,:,iblk) * sb%qmax_po4%ptr(ics:ice,:,iblk) &
                                                 / (sb%po4_solute%ptr(ics:ice,:,iblk) + km_adsorpt_po4 * sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk))
         sb%po4_assoc_slow%ptr(ics:ice,:,iblk) = 70.0_wp / molar_mass_P
         sb%po4_occluded%ptr(ics:ice,:,iblk)   = 150.0_wp / molar_mass_P		 
		 
       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT

! Pass back values to process memory & deallocate pool structures 
#define Dealloc_sb_pool
#include "slm_pools_a_incl.f90"
#undef Dealloc_sb_pool

  END SUBROUTINE sb_init

END MODULE mo_sb_init
