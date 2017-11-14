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
    USE mo_sb_memory_class,       ONLY: sb_pool
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
    ! Declare pointers to variables in memory
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_solute
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: k_bioturb
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: ph_sl
    dsl4jsb_Real3D_onChunk :: lithology_class_sl
    dsl4jsb_Real3D_onChunk :: po4_primary
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
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     bulk_dens_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded) 
    dsl4jsb_Get_var3D_onChunk(SB_,       ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lithology_class_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)
    ! Get process memory pool structures
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 all soil models
    !! 
    ph_sl(:,:)              = dsl4jsb_Config(SB_)%soil_ph
    lithology_class_sl(:,:) = dsl4jsb_Config(SB_)%lithology_class


    !> 2.0 differentiate soil model init
    !! 
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

             nh4_solute(:,1) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             no3_solute(:,1) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             po4_solute(:,1) = 0.01_wp / dsl4jsb_Config(SOIL_)%dz_water(1) 
             
             ! initialise mineral soil P pools
             po4_primary(:,1) = 100.0_wp / molar_mass_P	
          ELSE 
             sb_pool%soluable_litter%carbon(:,1)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             sb_pool%polymeric_litter%carbon(:,1) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             IF(lctlib%growthform==itree) &
                sb_pool%woody_litter%carbon(:,1)  = hlp3 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb_pool%microbial%carbon(:,1)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             sb_pool%residue%carbon(:,1)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * & 
                                                                       ( hlp6 + (1.0_wp - hlp6) * root_fraction_sl(:,1)) 
             DO isoil=2,options%nsoil_sb
                sb_pool%soluable_litter%carbon(:,isoil)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%polymeric_litter%carbon(:,isoil) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%microbial%carbon(:,isoil)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil) 
                sb_pool%residue%carbon(:,isoil)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * & 
                                                                          (1.0_wp - hlp6) * root_fraction_sl(:,isoil)

                nh4_solute(:,isoil) = 0.1_wp
                no3_solute(:,isoil) = 0.1_wp 
                po4_solute(:,isoil) = 0.01_wp 
             ENDDO
             
             ! initialise mineral soil P pools
             po4_primary(:,1) = 100.0_wp / molar_mass_P	
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

          nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

       CASE('jsm')
 
         ! set microbial biomass to minimum possible density
         sb_pool%microbial%carbon(:,:)     = min_mic_biomass
         sb_pool%microbial%nitrogen(:,:)   = sb_pool%microbial%carbon(:,:) / microbial_cn 
         sb_pool%microbial%phosphorus(:,:) = sb_pool%microbial%nitrogen(:,:) / microbial_np
         sb_pool%microbial%carbon13(:,:)   = sb_pool%microbial%carbon(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
         sb_pool%microbial%carbon14(:,:)   = sb_pool%microbial%carbon(:,:) * calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
         sb_pool%microbial%nitrogen15(:,:) = sb_pool%microbial%nitrogen(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
     
         ! some initial values for testing
         nh4_solute(:,:) = 1.0_wp 
         nh4_n15_solute(:,:) = nh4_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         no3_solute(:,:) = 1.0_wp 
         no3_n15_solute(:,:) = no3_solute(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         po4_solute(:,:) = 0.01_wp 

         ! inital Qmax based on current soil carbon content
         bulk_soil_carbon_sl(:,:) = calc_bulk_soil_carbon(options, sb_pool)
         bulk_dens_corr_sl(:,:) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:), bulk_dens_sl(:,:))
         qmax_org(:,:) = calc_qmax(bulk_dens_corr_sl(:,:),qmax_org_mineral_default)
         qmax_po4(:,:) = calc_qmax(bulk_dens_corr_sl(:,:),qmax_po4_mineral)
 
         ! initialise diffusion by bioturbation
         k_bioturb(:,:) = calc_bioturbation_rate(options,bulk_dens_corr_sl(:,:)) 

         ! initialise mineral soil P pools
         po4_assoc_fast(:,:) = po4_solute(:,:) * qmax_po4(:,:) / (po4_solute(:,:) + km_adsorpt_po4 * bulk_dens_corr_sl(:,:))
         po4_assoc_slow(:,:) = 70.0_wp / molar_mass_P
         po4_occluded(:,:)   = 150.0_wp / molar_mass_P		 
         po4_primary(:,1)    = 100.0_wp / molar_mass_P			 
       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT

! Pass back values to process memory
#define Passback_sb_pool
#include "slm_pools_a_incl.f90"
#undef Passback_sb_pool

  END SUBROUTINE sb_init

END MODULE mo_sb_init
