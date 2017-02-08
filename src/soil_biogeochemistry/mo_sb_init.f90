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
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_config_class,     ONLY: base_config  
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_sb_config_class,       ONLY: sb_config
    USE mo_sb_jsm_processes
    USE mo_sb_constants

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SOIL_)


  IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout)      :: tile       !< data structure for dynamic vegetation and soil variables
    TYPE(t_jsb_task_options), INTENT(in)         :: options      !< model options
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: isoil
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SOIL_)
    TYPE(t_sb_memory),    POINTER                :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_soil_memory),    POINTER              :: soil       ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)

    sb      => sb__mem 
    soil    => soil__mem 
    model   => Get_model(tile%owner_model_id)
    lctlib  => model%lctlib(tile%lcts%lib_id)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
 
    SELECT CASE(TRIM(base_config%slm_model_name))
!! BLARPP: still treating plant and land the same here. Dummy version only. Later land and soil need to go together
       CASE('land','plant')  
 
         ! Initial values for Siliva's simple one-bucket soil biogeochemical scheme
         ! Initialise surface litter pool to 0.5 mol / m2 N and P, and provide some initial NH4 and PO4
         sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = 8.45_wp / soil_config%dz_water(1)
         sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(:,1,iblk) 
         sb%nh4_solute%ptr(ics:ice,1,iblk) = 0.1_wp / soil_config%dz_water(1)
         sb%po4_solute%ptr(ics:ice,1,iblk) = 0.01_wp / soil_config%dz_water(1) 
 
       CASE('soil')
 
         ! initialise diffusion by bioturbation
         sb%k_bioturb%ptr(ics:ice,:,iblk) = k_diff_org
         ! set microbial biomass to minimum possible sensity (does this require adjustment with layer thickness, probably so...?
         sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) = min_mic_biomass
     
         ! some initial values for testing
         sb%pool%dom%carbon%ptr(ics:ice,1,iblk) = 1.0_wp
         sb%pool%microbial%carbon%ptr(ics:ice,1,iblk) = .1_wp
         sb%nh4_assoc%ptr(ics:ice,1,iblk) = 1.0_wp 

         ! inital Qmax based on current soil carbon content
         sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk) = calc_bulk_soil_carbon(options,soil_config%nsoil_water, sb%pool)
         sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk) = calc_bulk_density_correction(sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk), & 
                                                                         soil%bulk_dens_l%ptr(ics:ice,:,iblk))
         sb%qmax_org%ptr(ics:ice,:,iblk) = calc_qmax_org(sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk))
 
       CASE DEFAULT
         CALL finish(routine,'Invalid SLM model name. SB_ must only be used with plant, land or soil')
 
    END SELECT

  END SUBROUTINE sb_init

END MODULE mo_sb_init
