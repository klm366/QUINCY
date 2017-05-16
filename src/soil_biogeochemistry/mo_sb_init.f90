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
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_sb_config_class,       ONLY: sb_config
    USE mo_sb_jsm_processes
    USE mo_sb_jsm_transport
    USE mo_sb_constants

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)
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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SOIL_)
!! *slmtemp    
    TYPE(t_sb_memory),    POINTER                :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_soil_memory),    POINTER              :: soil       ! local pointer to soil__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------


!! *slmtemp
    sb      => sb__mem 
    soil    => soil__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

 
    SELECT CASE(TRIM(sb_config%sb_model_scheme))

       CASE("simple_0d", "simple_1d")

          ! Initial values for Soenke's four pool soil biogeochemical model
          sb%pool%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = 500._wp / molar_mass_C / soil_config%dz_water(1)
          sb%pool%woody_litter%carbon%ptr(ics:ice,1,iblk) = 1000._wp / molar_mass_C / soil_config%dz_water(1)
          sb%pool%microbial%carbon%ptr(ics:ice,1,iblk) = 1000._wp / molar_mass_C / soil_config%dz_water(1)
          sb%pool%residue%carbon%ptr(ics:ice,1,iblk) = 5000._wp / molar_mass_C / soil_config%dz_water(1)
 
          ! non-woody litter all other elements
          sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (50._wp * molar_mass_N / molar_mass_C)
          sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) * & 
                                                                    calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                (220._wp * molar_mass_N / molar_mass_C)
          sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                (14._wp * molar_mass_P / molar_mass_N) 
          sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) * & 
                                                                calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (=microbial) all other elements
          sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / microbial_cn
          sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / microbial_np 
          sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) / microbial_cn
          sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) / microbial_np 
          sb%pool%residue%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb%nh4_solute%ptr(ics:ice,1,iblk) = 0.1_wp / soil_config%dz_water(1)
          sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          sb%po4_solute%ptr(ics:ice,1,iblk) = 0.01_wp / soil_config%dz_water(1) 

       CASE('jsm')
 
         ! set microbial biomass to minimum possible density
         sb%pool%microbial%carbon%ptr(ics:ice,:,iblk)     = min_mic_biomass
         sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / microbial_cn 
         sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / microbial_np
         sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
         sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) * calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
         sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)  / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
     
         ! some initial values for testing
         sb%nh4_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%no3_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%no3_n15_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%po4_solute%ptr(ics:ice,:,iblk) = 1.0_wp 

         ! inital Qmax based on current soil carbon content
         sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk) = calc_bulk_soil_carbon(options,soil_config%nsoil_water, sb%pool)
         sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk) = calc_bulk_density_correction(sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk), & 
                                                                         soil%bulk_dens_sl%ptr(ics:ice,:,iblk))
         sb%qmax_org%ptr(ics:ice,:,iblk) = calc_qmax_org(sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk))
 
         ! initialise diffusion by bioturbation
         sb%k_bioturb%ptr(ics:ice,:,iblk) = calc_bioturbation_rate(options,soil_config%nsoil_water,sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk)) 

       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT

  END SUBROUTINE sb_init

END MODULE mo_sb_init
